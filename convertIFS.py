#!/usr/bin/env python3
"""
convertIFS.py — Convert in-frame stop codons in protein FASTA files.

In-frame stop codons are represented as '.' within a sequence. A trailing '.'
(terminal stop codon) is stripped. Internal '.' are replaced with 'X' (unknown
amino acid). All sequences are retained. Clean sequences are written to STDOUT;
summary report to STDERR.

Usage: convertIFS.py prot1.fa [prot2.fa ...]
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
    Parse FASTA, fix stop codons, return (converted_records, stats).
    converted_records: list of (header, seq_lines) with stops handled.
    stats: dict with counts for reporting.
    """
    records = []
    trailing_seqs, internal_seqs = [], []
    trailing_total, internal_total = 0, 0

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
            trailing_total += 1
            trailing_seqs.append(name)

        # Replace internal stops with X
        n_internal = 0
        for i, line in enumerate(new_lines):
            count = line.count('.')
            if count:
                n_internal += count
                new_lines[i] = line.replace('.', 'X')
        if n_internal:
            internal_total += n_internal
            internal_seqs.append((name, n_internal))

        records.append((header, new_lines))

    stats = {
        'total_input': len(records),
        'trailing_seqs': trailing_seqs,
        'trailing_total': trailing_total,
        'internal_seqs': internal_seqs,
        'internal_total': internal_total,
    }
    return records, stats


def print_report(filepath, stats):
    n_in = stats['total_input']
    trailing_seqs = stats['trailing_seqs']
    internal_seqs = stats['internal_seqs']

    print(f"\n{filepath}", file=sys.stderr)
    print(f"  Input sequences:    {n_in}", file=sys.stderr)
    if trailing_seqs:
        if len(trailing_seqs) == n_in:
            print(f"  Trailing stop ('.'): all {n_in} sequences terminated with '.' stop codon — truncated", file=sys.stderr)
        else:
            print(f"  Trailing stop ('.'): {len(trailing_seqs)} seqs truncated", file=sys.stderr)
            for s in trailing_seqs:
                print(f"TR\t{s}", file=sys.stderr)
    else:
        print(f"  Trailing stop ('.'): none found", file=sys.stderr)
    if internal_seqs:
        print(f"  Internal IFS:       {len(internal_seqs)} seqs converted ({stats['internal_total']} total IFS -> X)", file=sys.stderr)
        for name, count in internal_seqs:
            print(f"CV\t{name}\t{count}", file=sys.stderr)
    else:
        print(f"  Internal IFS:       none found", file=sys.stderr)


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
            records, stats = process_file(filepath)
        except FileNotFoundError:
            print(f"ERROR: file not found: {filepath}", file=sys.stderr)
            continue
        print_report(filepath, stats)
        write_fasta(records)
