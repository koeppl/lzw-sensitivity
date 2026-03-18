"""LZW compression analysis tool.

Implements the LZW (Lempel-Ziv-Welch) factorization on integer sequences and
reports the number of factors together with a detailed per-step breakdown.

Usage
-----
    python lzw_compress.py -i INPUT_PREFIX

The script reads input files named ``INPUT_PREFIX_none.txt``,
``INPUT_PREFIX_sub.txt``, ``INPUT_PREFIX_del.txt``, and
``INPUT_PREFIX_ins.txt`` (whichever exist), runs LZW factorization on each,
prints a detailed per-step table, and then emits a single summary line::

    RESULT prefix=<INPUT_PREFIX> num_none=<N> num_sub=<N> num_del=<N> num_ins=<N>

For variants whose input file is absent the count is reported as ``None``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]

# Variant suffixes processed by this script
VARIANTS: list[str] = ["none", "sub", "del", "ins"]


def load_sequence(path: Path) -> list[int]:
    """Load an integer sequence from a file in ``[a,b,c,...]`` format.

    Args:
        path: Path to the file to read.

    Returns:
        A list of integers parsed from the file.
    """
    text = path.read_text(encoding="utf-8").strip()
    # Strip surrounding brackets
    text = text[1:-1]
    if not text:
        return []
    return [int(x) for x in text.split(",")]


def format_seq(seq_tuple: tuple[int, ...] | None, label_map: dict[int, str]) -> str:
    """Format a sequence tuple as a concatenation of symbol labels.

    Args:
        seq_tuple: A tuple of integer symbols, or None.
        label_map: Mapping from integer symbol to its display label.

    Returns:
        A string of concatenated labels, or an empty string for None/empty input.
    """
    if not seq_tuple:
        return ""
    return "".join(label_map[s] for s in seq_tuple)


def lzw_compress(
    symbol_sequence: list[int],
) -> tuple[list[dict[str, str]], int]:
    """Run LZW compression on an integer symbol sequence.

    The algorithm maintains a dictionary initialised with all unique symbols.
    At each step it finds the longest dictionary match starting at the current
    position, outputs that as a factor, then registers the match extended by
    the next symbol as a new dictionary entry.

    Args:
        symbol_sequence: The input sequence of integer symbols.

    Returns:
        A tuple ``(factors_details, factor_count)`` where ``factors_details``
        is a list of per-step records suitable for display, and ``factor_count``
        is the total number of factors.
    """
    unique_symbols = sorted(set(symbol_sequence))
    label_map: dict[int, str] = {s: f"sigma_{s}" for s in unique_symbols}

    # Dictionary keyed by tuples of symbols; values are human-readable labels
    dictionary: dict[tuple[int, ...], str] = {(s,): f"sigma_{s}" for s in unique_symbols}

    factors_details: list[dict[str, str]] = []
    seq = tuple(symbol_sequence)
    n = len(seq)
    i = 0
    factor_count = 0
    dict_add_count = 0

    while i < n:
        # Find the longest dictionary match at the current position
        w_seq: tuple[int, ...] = ()
        w_lbl = ""
        for d_seq in sorted(dictionary.keys(), key=len, reverse=True):
            if seq[i : i + len(d_seq)] == d_seq:
                w_seq = d_seq
                w_lbl = dictionary[d_seq]
                break

        # Fallback: always matches at least the single symbol (guaranteed by init)
        if not w_seq:
            w_seq = (seq[i],)
            w_lbl = label_map[seq[i]]

        factor_count += 1

        # Register W+c as a new dictionary entry where c is the next symbol
        new_entry: tuple[int, ...] | None = None
        new_id = "-"
        next_pos = i + len(w_seq)
        if next_pos < n:
            c: tuple[int, ...] = (seq[next_pos],)
            combined = w_seq + c
            if combined not in dictionary:
                dict_add_count += 1
                new_id = f"D_{dict_add_count}"
                dictionary[combined] = new_id
                new_entry = combined

        factors_details.append(
            {
                "Step": str(factor_count),
                "Factor": format_seq(w_seq, label_map),
                "Dict match": w_lbl,
                "New entry (W+c)": format_seq(new_entry, label_map) if new_entry else "-",
                "New ID": new_id,
            }
        )

        i += len(w_seq)

    return factors_details, factor_count


def print_table(rows: list[dict[str, str]], headers: list[str]) -> None:
    """Print a list of row dicts as a plain-text aligned table.

    Args:
        rows: List of dicts, one per table row.
        headers: Ordered list of column header names (also used as dict keys).
    """
    widths = {h: len(h) for h in headers}
    for row in rows:
        for h in headers:
            widths[h] = max(widths[h], len(row.get(h, "")))

    header_line = " | ".join(h.ljust(widths[h]) for h in headers)
    print(header_line)
    print("-" * len(header_line))
    for row in rows:
        print(" | ".join(row.get(h, "").ljust(widths[h]) for h in headers))


def compress_file(input_path: Path) -> int:
    """Compress one sequence file with LZW and print a detailed breakdown.

    Args:
        input_path: Path to the input sequence file (``[a,b,c,...]`` format).

    Returns:
        The number of LZW factors.
    """
    seq = load_sequence(input_path)
    results, z = lzw_compress(seq)

    print(f"\n=== LZW compression: {input_path} ===")
    print(f"Input T: {''.join(f'sigma_{s}' for s in seq)}")
    print(f"Total factors z: {z}\n")

    headers = ["Step", "Factor", "Dict match", "New entry (W+c)", "New ID"]
    print("--- LZW compression process ---")
    print_table(results, headers)

    print("\n--- Factorization visualization ---")
    print("".join(f"[{item['Factor']}]" for item in results))

    return z


def main() -> None:
    """Parse arguments and run LZW compression on all available input variants."""
    parser = argparse.ArgumentParser(
        description=(
            "Run LZW compression on sequence files. "
            "Reads INPUT_none.txt, INPUT_sub.txt, INPUT_del.txt, INPUT_ins.txt "
            "(whichever exist), prints per-file details, then emits a RESULT line."
        )
    )
    parser.add_argument(
        "-i",
        required=True,
        metavar="INPUT_PREFIX",
        help="Common prefix path for input sequence files",
    )
    args = parser.parse_args()

    in_prefix = Path(args.i)

    counts: dict[str, int | None] = {v: None for v in VARIANTS}
    processed = 0
    for variant in VARIANTS:
        input_path = Path(f"{in_prefix}_{variant}.txt")
        if not input_path.exists():
            continue
        counts[variant] = compress_file(input_path)
        processed += 1

    if processed == 0:
        print(
            f"No input files found with prefix '{in_prefix}' "
            f"and variants {VARIANTS}.",
            file=sys.stderr,
        )
        sys.exit(1)

    parts = " ".join(f"num_{v}={counts[v]}" for v in VARIANTS)
    print(f"RESULT prefix={args.i} {parts}")


if __name__ == "__main__":
    main()
