#!/usr/bin/env bash
# collect_manifest.sh — collect a filesystem manifest for mirror reconciliation
#
# Usage:
#   bash collect_manifest.sh [OPTIONS] /fs/mount/path /logical/root > manifest.jsonl
#
# Arguments:
#   /fs/mount/path   Actual filesystem path where the tree is mounted on this server.
#   /logical/root    Canonical path that this tree represents, shared across all
#                    servers being compared. Must be an absolute path. All record
#                    paths in the manifest will be relative to this logical root,
#                    making manifests from different servers directly comparable
#                    regardless of where each is mounted.
#
# Example:
#   Server A (mounted at /mnt/vol_alpha, represents /projects):
#     bash collect_manifest.sh /mnt/vol_alpha /projects > manifest_a.jsonl
#
#   Server B (mounted at /data/storage, represents /projects):
#     bash collect_manifest.sh /data/storage /projects > manifest_b.jsonl
#
#   Both manifests will have paths like "src/main.py", "docs/readme.txt", etc.,
#   anchored to /projects, making them directly comparable.
#
# Options:
#   -x PATH   Exclude path (relative to logical root, with or without leading /);
#             may be repeated
#   -d N      Maximum symlink chain depth before flagging as a loop (default: 10)
#   -h        Show this help
#
# Output: one JSON record per line (JSONL).
#
# The first line is always a header record (type "header"):
#   {
#     "type":         "header",
#     "logical_root": "/projects",
#     "fs_root":      "/mnt/vol_alpha",
#     "hostname":     "serverA",
#     "collected_at": 1700000000,
#     "excludes":     ["tmp", "cache"]
#   }
#
# Subsequent lines are one record per filesystem entry, with fields:
#   path          path relative to logical root (no leading slash)
#   dir           directory component (empty string for root-level entries)
#   name          filename component
#   type          "file" | "symlink" | "directory"
#   size          bytes for files; null for symlinks and directories
#   mtime         modification time as unix timestamp; null for symlinks
#   md5           hex content hash for files; null for symlinks and directories
#   target        symlink target string as written; null for files/dirs
#   target_exists true/false whether symlink resolves; null for files/dirs
#   target_type   "file"|"directory"|"broken"|"loop"|"unknown"; null for files/dirs
#   link_depth    hops to resolve (-1 if loop); null for files/dirs
#
# Platform notes:
#   Requires: bash 4+, find, stat, readlink, md5sum (Linux) or md5 (macOS), python3.
#   Tested on Ubuntu 22.04 (GNU coreutils) and macOS with BSD stat.
#   python3 is used solely for JSON serialisation — bash cannot safely escape
#   arbitrary filenames (spaces, quotes, newlines, backslashes) in JSON output.

set -euo pipefail

# ── Argument parsing ──────────────────────────────────────────────────────────

usage() {
    grep '^#' "$0" | sed 's/^# \?//'
    exit 0
}

MAX_LINK_DEPTH=10
EXCLUDES=()

while getopts ":x:d:h" opt; do
    case $opt in
        x) EXCLUDES+=("$OPTARG") ;;
        d) MAX_LINK_DEPTH="$OPTARG" ;;
        h) usage ;;
        :) echo "ERROR: option -$OPTARG requires an argument" >&2; exit 1 ;;
        \?) echo "ERROR: unknown option -$OPTARG" >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

if [[ $# -lt 2 ]]; then
    echo "ERROR: two positional arguments required" >&2
    echo "Usage: $0 [OPTIONS] /fs/mount/path /logical/root" >&2
    exit 1
fi

FS_ROOT="${1%/}"       # actual mount point on this server
LOGICAL_ROOT="${2%/}"  # canonical name shared across all servers

if [[ ! -d "$FS_ROOT" ]]; then
    echo "ERROR: '$FS_ROOT' is not a directory or does not exist" >&2
    exit 1
fi

if [[ "$LOGICAL_ROOT" != /* ]]; then
    echo "ERROR: logical_root '$LOGICAL_ROOT' must be an absolute path (start with /)" >&2
    exit 1
fi

# ── Platform detection ────────────────────────────────────────────────────────

if stat --version &>/dev/null 2>&1; then
    STAT_PLATFORM="gnu"
else
    STAT_PLATFORM="bsd"
fi

if command -v md5sum &>/dev/null; then
    MD5_CMD="md5sum"
elif command -v md5 &>/dev/null; then
    MD5_CMD="md5 -q"
else
    echo "ERROR: neither md5sum nor md5 found in PATH" >&2
    exit 1
fi

if ! command -v python3 &>/dev/null; then
    echo "ERROR: python3 is required for safe JSON serialisation" >&2
    exit 1
fi

# ── Stat helper ───────────────────────────────────────────────────────────────
# Stats the entry itself (not its target) — lstat() semantics.
# Prints two space-separated values: size mtime

get_size_mtime() {
    if [[ "$STAT_PLATFORM" == "gnu" ]]; then
        stat -c '%s %Y' "$1"
    else
        stat -f '%z %m' "$1"
    fi
}

# ── Symlink chain resolver ────────────────────────────────────────────────────
# Follows a symlink chain manually up to MAX_LINK_DEPTH hops, without using
# realpath (which can hang on certain loop configurations).
#
# Outputs exactly two lines to stdout:
#   line 1: depth as integer (-1 means loop)
#   line 2: target_type — "file" | "directory" | "broken" | "loop" | "unknown"
#
# NOTE: (( n++ )) exits with status 1 when n is 0 (the expression evaluates to
# the pre-increment value, which is falsy). Under set -e this would kill the
# subshell silently. All arithmetic uses (( expr )) || true to prevent this.

resolve_link() {
    local current="$1"
    local depth=0
    local target absdir

    while [[ $depth -le $MAX_LINK_DEPTH ]]; do

        if [[ ! -L "$current" ]]; then
            # Reached the final target — classify it
            if   [[ -f "$current" ]]; then echo "$depth"; echo "file";      return
            elif [[ -d "$current" ]]; then echo "$depth"; echo "directory"; return
            elif [[ -e "$current" ]]; then echo "$depth"; echo "unknown";   return
            else                           echo "$depth"; echo "broken";    return
            fi
        fi

        target=$(readlink "$current")    # one hop only, never recursive

        # Resolve relative targets against the directory containing the link
        if [[ "$target" != /* ]]; then
            absdir=$(cd "$(dirname "$current")" 2>/dev/null && pwd -P) || {
                echo "$depth"; echo "unknown"; return
            }
            target="${absdir}/${target}"
        fi

        # Normalise path in pure bash (no realpath — it follows links)
        while [[ "$target" == *"//"*  ]]; do target="${target//\/\//\/}";  done
        while [[ "$target" == *"/./"* ]]; do target="${target//\/.\//\/}"; done

        (( depth++ )) || true    # || true: (( 0 )) is falsy, would trigger set -e
        current="$target"
    done

    echo "-1"; echo "loop"
}

# ── JSON emitters ─────────────────────────────────────────────────────────────
# All variable data is passed as positional arguments to python3 so that
# special characters in filenames (newlines, quotes, backslashes, etc.) are
# handled by Python's json.dumps rather than by fragile bash string escaping.

emit_header() {
    local collected_at
    collected_at=$(date +%s)

    python3 - "$FS_ROOT" "$LOGICAL_ROOT" "$(hostname)" "$collected_at" \
               "${EXCLUDES[@]+"${EXCLUDES[@]}"}" <<'PYEOF'
import json, sys

args         = sys.argv[1:]
fs_root      = args[0]
logical_root = args[1]
hostname     = args[2]
collected_at = int(args[3])
excludes     = args[4:]      # empty list if no -x flags were given

print(json.dumps({
    "type":         "header",
    "logical_root": logical_root,
    "fs_root":      fs_root,
    "hostname":     hostname,
    "collected_at": collected_at,
    "excludes":     excludes,
}))
PYEOF
}

emit_record() {
    # Arguments (all strings; "null" is the sentinel for JSON null):
    #  1  rel_path    2  dir_part    3  name_part   4  ftype
    #  5  size        6  mtime       7  md5
    #  8  target      9  tgt_exists  10 tgt_type    11 link_depth
    python3 - "$@" <<'PYEOF'
import json, sys

def maybe_int(v):
    try:    return int(v)
    except: return None

def maybe_str(v):
    return None if v == "null" else v

def maybe_bool(v):
    if v == "true":  return True
    if v == "false": return False
    return None

_, path, dir_, name, ftype, size, mtime, md5, \
    target, tgt_exists, tgt_type, link_depth = sys.argv

print(json.dumps({
    "path":          path,
    "dir":           dir_,
    "name":          name,
    "type":          ftype,
    "size":          maybe_int(size),
    "mtime":         maybe_int(mtime),
    "md5":           maybe_str(md5),
    "target":        maybe_str(target),
    "target_exists": maybe_bool(tgt_exists),
    "target_type":   maybe_str(tgt_type),
    "link_depth":    maybe_int(link_depth),
}))
PYEOF
}

# ── Exclusion check ───────────────────────────────────────────────────────────

is_excluded() {
    local rel="$1"   # relative path, no leading slash
    local excl
    for excl in "${EXCLUDES[@]+"${EXCLUDES[@]}"}"; do
        excl="${excl#/}"    # tolerate leading slash in -x arguments
        if [[ "$rel" == "$excl" || "$rel" == "$excl/"* ]]; then
            return 0
        fi
    done
    return 1
}

# ── Main traversal ────────────────────────────────────────────────────────────
# find output is written to a temp file and read back with a redirect rather
# than a pipe. A pipeline runs the while loop in a subshell, which prevents
# mapfile from working correctly with process substitution in all bash versions.

TMPFILE=$(mktemp /tmp/manifest_find.XXXXXX)
trap 'rm -f "$TMPFILE"' EXIT

main() {
    emit_header

    find "$FS_ROOT" -mindepth 1 \( -type f -o -type l -o -type d \) -print0 \
        > "$TMPFILE"

    while IFS= read -r -d '' entry; do

        # Path relative to the filesystem root (no leading slash)
        rel="${entry#"$FS_ROOT"/}"
        dir_part=$(dirname "$rel")
        name_part=$(basename "$rel")

        # dirname returns "." for entries directly under the root
        [[ "$dir_part" == "." ]] && dir_part=""

        is_excluded "$rel" && continue

        # ── Symlink ──────────────────────────────────────────────────────────
        if [[ -L "$entry" ]]; then
            link_target=$(readlink "$entry")

            # resolve_link emits exactly two lines: depth, then target_type
            mapfile -t resolved < <(resolve_link "$entry")
            link_depth="${resolved[0]}"
            tgt_type="${resolved[1]}"

            if [[ "$tgt_type" == "broken" || "$tgt_type" == "loop" ]]; then
                tgt_exists="false"
            else
                tgt_exists="true"
            fi

            emit_record "$rel" "$dir_part" "$name_part" \
                "symlink" "null" "null" "null" \
                "$link_target" "$tgt_exists" "$tgt_type" "$link_depth"

        # ── Regular file ─────────────────────────────────────────────────────
        elif [[ -f "$entry" ]]; then
            read -r size mtime < <(get_size_mtime "$entry")
            md5=$($MD5_CMD "$entry" | awk '{print $1}')

            emit_record "$rel" "$dir_part" "$name_part" \
                "file" "$size" "$mtime" "$md5" \
                "null" "null" "null" "null"

        # ── Directory ────────────────────────────────────────────────────────
        elif [[ -d "$entry" ]]; then
            read -r _dsize mtime < <(get_size_mtime "$entry")

            emit_record "$rel" "$dir_part" "$name_part" \
                "directory" "null" "$mtime" "null" \
                "null" "null" "null" "null"
        fi

    done < "$TMPFILE"
}

main
