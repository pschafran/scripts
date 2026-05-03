#!/usr/bin/env python3
"""
Outsource session documentation writing to a local Ollama model.

Gathers git context from the working directory, reads the lab notebook entry
template, and asks the model to write:
  1. A lab notebook entry at /media/data/projects/{project}/notebook/
  2. A GitHub session note at ~/Notes/claude-sessions/YYYY-MM-DD_{project}.md

After the model writes both files, updates the project metadata and TOC via
lab_notebook_manager.py.

Usage:
    python ollama_write_docs.py --project repro_coexp --working-dir /media/data/projects/repro_coexp \
        --description "Analyzed co-expression communities and plotted results"

    # Read description from a file or stdin
    python ollama_write_docs.py --project repro_coexp --working-dir /path/to/project \
        --description - < session_notes.txt
"""

from ollama import chat
import argparse
import importlib.util
import json
import os
import sys
import subprocess
import uuid
from datetime import datetime
from pathlib import Path

NOTEBOOK_MANAGER = "/media/data/notebook/lab_notebook_manager.py"
NOTES_DIR = Path.home() / "Notes"
ENTRY_TEMPLATE = "/media/data/notebook/.entry_template.md"

# Import only the metadata/TOC functions — not main() — to avoid creating a
# duplicate skeleton entry file
_spec = importlib.util.spec_from_file_location("lab_notebook_manager", NOTEBOOK_MANAGER)
_nbm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_nbm)

parser = argparse.ArgumentParser()
parser.add_argument("--project", required=True, help="Project name (e.g. repro_coexp)")
parser.add_argument("--working-dir", required=True, help="Project working directory for git context")
parser.add_argument("--description", required=True,
                    help="Description of session work. Use '-' to read from stdin.")
parser.add_argument("--model", default="qwen3.6-128k", help="Ollama model to use")
parser.add_argument("--claude-model", default="claude-sonnet-4-6",
                    help="Claude model name to record as author in the entry")
parser.add_argument("--context-file", help="Path to the full Claude session transcript for Ollama to read")
parser.add_argument("--push", action="store_true", help="Commit and push Notes to GitHub after writing")
args = parser.parse_args()

# --- Gather inputs ---

if args.description == "-":
    description = sys.stdin.read().strip()
else:
    description = args.description

today = datetime.now().strftime("%Y-%m-%d")
session_id = str(uuid.uuid4())
session_id_short = session_id[:8]
working_dir = os.path.abspath(args.working_dir)

# Git context from the working directory
def run(cmd, cwd=None):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    return result.stdout.strip() if result.returncode == 0 else "(unavailable)"

git_log = run("git log --oneline -20", cwd=working_dir)
git_diff_stat = run("git diff --stat HEAD~1 2>/dev/null || git diff --stat HEAD", cwd=working_dir)
git_status = run("git status --short", cwd=working_dir)

# Templates
entry_template = Path(ENTRY_TEMPLATE).read_text() if Path(ENTRY_TEMPLATE).exists() else "(template not found)"

# Output paths
notebook_dir = Path(f"/media/data/projects/{args.project}/notebook")
entry_filename = f"entry_{today}_{session_id_short}.md"
entry_path = str(notebook_dir / entry_filename)
session_note_path = str(NOTES_DIR / "claude-sessions" / f"{today}_{args.project}.md")

# --- Transcript stripping ---

def _tool_summary(name: str, inp: dict) -> str:
    key_args = {
        'Read':       'file_path',
        'Edit':       'file_path',
        'Write':      'file_path',
        'Bash':       'command',
        'WebFetch':   'url',
        'Agent':      'description',
        'Glob':       'pattern',
        'Grep':       'pattern',
        'ToolSearch': 'query',
    }
    arg = inp.get(key_args.get(name, '')) or (next(iter(inp.values()), '') if inp else '')
    arg = str(arg)
    if len(arg) > 120:
        arg = arg[:117] + '...'
    return '[%s: %s]' % (name, arg)


def strip_transcript(path: str) -> str:
    """Parse a Claude Code JSONL session transcript and return a clean,
    readable conversation log. Keeps user messages, assistant text, and
    one-line tool call summaries. Strips metadata, thinking blocks,
    tool results, file attachments, and all bookkeeping entry types."""
    turns = []
    with open(path, 'r', errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue

            t = obj.get('type')
            ts = obj.get('timestamp', '')
            time_str = ts[11:16] if len(ts) >= 16 else ''  # HH:MM from ISO timestamp

            if t == 'user':
                content = obj.get('message', {}).get('content', '')
                if isinstance(content, list):
                    # List content mixes text and tool_result blocks; keep only text
                    parts = [c.get('text', '') for c in content if c.get('type') == 'text']
                    text = '\n'.join(p for p in parts if p).strip()
                else:
                    text = str(content).strip()
                if text:
                    turns.append('[%s] USER: %s' % (time_str, text))

            elif t == 'assistant':
                content = obj.get('message', {}).get('content', [])
                if not isinstance(content, list):
                    continue
                lines = []
                for block in content:
                    bt = block.get('type')
                    if bt == 'text':
                        text = block.get('text', '').strip()
                        if text:
                            lines.append(text)
                    elif bt == 'tool_use':
                        lines.append('  ' + _tool_summary(block.get('name', '?'), block.get('input', {})))
                    # thinking blocks: skip
                if lines:
                    turns.append('[%s] ASSISTANT:\n%s' % (time_str, '\n'.join(lines)))

            # All other types (permission-mode, file-history-snapshot, ai-title,
            # last-prompt, system, attachment): skip

    return '\n\n'.join(turns)


# --- Tools ---

def write_file(path: str, content: str) -> str:
    """Write content to a file, creating parent directories as needed.

    Args:
        path: Absolute path to write.
        content: Full text content of the file.

    Returns:
        str: Confirmation message or error.
    """
    try:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(content)
        return "Written: %s (%d bytes)" % (path, len(content))
    except Exception as e:
        return "Error writing %s: %s" % (path, e)

# --- Prompt ---

prompt = """You are writing session documentation for a plant genomics research lab.

Today's date: {today}
Project: {project}
Session ID: {session_id}
Claude model that did the work: {claude_model}
Working directory: {working_dir}

## Primary session notes

{description}
{context_file_section}

## Git context (recent commits)

{git_log}

## Files changed (git diff --stat)

{git_diff_stat}

## Uncommitted changes (git status)

{git_status}

---

## Your task

Write TWO files using the write_file tool.

### File 1 — Lab notebook entry

Path: {entry_path}

Follow this template exactly, filling in every section with specific, accurate
content based on the session description and git context above. Do not leave
placeholder text. Write in first person ("I worked on..."). Be specific about
file paths, commands run, and outputs generated.

TEMPLATE:
{entry_template}

### File 2 — GitHub session note

Path: {session_note_path}

Use this format:

# Claude Code Session - {today}

**Project:** {project}
**Location:** `{working_dir}`
**Duration:** [estimate based on description if not given]

## Summary
[2-3 sentences describing what was accomplished]

## Work Completed

### Files Modified
- `path/to/file` — [description of changes]

### Commands Executed
```bash
[key commands that were run]
```

## Key Decisions & Rationale
- **Decision:** [what was decided]
  - **Rationale:** [why]

## Technical Details
[Important technical information]

## Next Steps
- [ ] [specific action item]

## Tags
`#[topic]` `#[language]`

---

Write File 1 first, then File 2. Do not output any other text.
""".format(
    today=today,
    project=args.project,
    session_id=session_id,
    claude_model=args.claude_model,
    working_dir=working_dir,
    description=description,
    context_file_section=(
        "\n## Full session transcript\n\n%s\n" % strip_transcript(args.context_file)
    ) if args.context_file else "",
    git_log=git_log,
    git_diff_stat=git_diff_stat,
    git_status=git_status,
    entry_path=entry_path,
    entry_template=entry_template,
    session_note_path=session_note_path,
)

# --- Agentic loop ---

messages = [{'role': 'user', 'content': prompt}]

print("Generating documentation with %s..." % args.model)

tool_dispatch = {'write_file': write_file}

while True:
    response = chat(model=args.model, messages=messages, tools=[write_file])
    messages.append(response.message)

    if not response.message.tool_calls:
        break

    for tool_call in response.message.tool_calls:
        fn_name = tool_call.function.name
        fn_args = tool_call.function.arguments
        print("  %s -> %s" % (fn_name, fn_args.get('path', '')))
        if fn_name in tool_dispatch:
            result = tool_dispatch[fn_name](**fn_args)
        else:
            result = "Unknown tool: %s" % fn_name
        print("  %s" % result[:120])
        messages.append({'role': 'tool', 'content': result})

# --- Post-write: update metadata and TOC ---

if Path(entry_path).exists():
    print("\nUpdating project metadata and TOC...")
    _nbm.update_project_metadata(args.project, Path(entry_path), args.claude_model, description[:120])
    _nbm.generate_toc(args.project)
else:
    print("\nWarning: entry file was not written — skipping metadata update.")

# --- Optionally commit and push Notes ---

if args.push:
    print("\nCommitting Notes to GitHub...")
    cmds = [
        "git -C ~/Notes add claude-sessions/",
        "git -C ~/Notes commit -m 'Claude session: %s - %s'" % (args.project, today),
        "git -C ~/Notes push origin master",
    ]
    for cmd in cmds:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        print(result.stdout.strip() or result.stderr.strip())

print("\nDone.")
print("  Lab notebook entry: %s" % entry_path)
print("  Session note:       %s" % session_note_path)
