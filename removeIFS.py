#!/usr/bin/env python3
"""
removeIFS.py — Remove sequences containing in-frame stop codons from protein FASTA files.

In-frame stop codons are represented as '.' within a sequence. A trailing '.'
(terminal stop codon) is not considered an IFS and is stripped before checking.
Sequences with IFS after stripping the trailing stop are excluded entirely.
Clean sequences are written to STDOUT; summary report to STDERR.

Usage: removeIFS.py prot1.fa [prot2.fa ...]
"""

import sys


def parse_fasta(filepath):
    """Yield (header_line, seq_lines) tuples, preserving original line structure."""
    with open(filepath) as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if header is not None:
                    yield (header, seq_lines)
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield (header, seq_lines)


def process_file(filepath):
    """
    Parse FASTA, identify sequences with IFS, return (kept_records, stats).
    kept_records: list of (header, seq_lines) with trailing stop stripped.
    stats: dict with counts for reporting.
    """
    kept, removed = [], []
    trailing_total = 0

    for header, seq_lines in parse_fasta(filepath):
        name = header[1:].split()[0]
        new_lines = list(seq_lines)

        # Find last non-empty line
        last_idx = next((i for i in range(len(new_lines)-1, -1, -1) if new_lines[i]), None)

        # Strip trailing stop
        had_trailing = False
        if last_idx is not None and new_lines[last_idx].endswith('.'):
            new_lines[last_idx] = new_lines[last_idx][:-1]
            had_trailing = True

        # Count internal stops in remaining sequence
        n_internal = sum(line.count('.') for line in new_lines)

        if n_internal > 0:
            removed.append((name, n_internal))
        else:
            if had_trailing:
                trailing_total += 1
            kept.append((header, new_lines))

    stats = {
        'kept': len(kept),
        'removed': len(removed),
        'removed_seqs': removed,
        'trailing_truncated': trailing_total,
        'total_input': len(kept) + len(removed),
    }
    return kept, stats


def print_report(filepath, stats):
    n_in = stats['total_input']
    n_kept = stats['kept']
    n_removed = stats['removed']
    n_trail = stats['trailing_truncated']

    print(f"\n{filepath}", file=sys.stderr)
    print(f"  Input sequences:    {n_in}", file=sys.stderr)
    if n_trail:
        print(f"  Trailing stop ('.'): {n_trail} seqs truncated (kept)", file=sys.stderr)
    if n_removed:
        total_ifs = sum(c for _, c in stats['removed_seqs'])
        print(f"  Internal IFS:       {n_removed} seqs removed ({total_ifs} total IFS)", file=sys.stderr)
        for name, count in stats['removed_seqs']:
            print(f"CV\t{name}\t{count}", file=sys.stderr)
    else:
        print(f"  Internal IFS:       none found", file=sys.stderr)
    print(f"  Output sequences:   {n_kept}", file=sys.stderr)


def write_fasta(records):
    for header, seq_lines in records:
        print(header)
        for line in seq_lines:
            print(line)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    for filepath in sys.argv[1:]:
        try:
            kept, stats = process_file(filepath)
        except FileNotFoundError:
            print(f"ERROR: file not found: {filepath}", file=sys.stderr)
            continue
        print_report(filepath, stats)
        write_fasta(kept)
