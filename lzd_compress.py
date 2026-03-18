"""LZD compression analysis tool.

Implements the LZD (Lempel-Ziv-Double) factorization on integer sequences and
reports the number of factors together with a detailed per-step breakdown.

In LZD each factor is the concatenation of two dictionary matches:
  f_i = f_{i,1} · f_{i,2}
where f_{i,1} is the longest match in the current dictionary starting at the
current position, and f_{i,2} is the longest match immediately following
f_{i,1}.  The resulting phrase f_i is then added to the dictionary.

Usage
-----
    python lzd_compress.py -i INPUT_PREFIX

The script reads input files named ``INPUT_PREFIX_none.txt``,
``INPUT_PREFIX_del.txt``, etc. (whichever exist), runs LZD factorization on
each, prints a detailed per-step table, and then emits a single summary line::

    RESULT prefix=<INPUT_PREFIX> num_none=<N> num_sub=<N> num_del=<N> num_ins=<N>

For variants whose input file is absent the count is reported as ``None``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]

# Variant suffixes processed by this script.
# Note: lzd_generate.py only produces 'none' and 'del' variants; the
# 'sub' and 'ins' entries are silently skipped when their files are absent.
VARIANTS: list[str] = ["none", "sub", "del", "ins"]


def load_sequence(path: Path) -> list[int]:
    """Load an integer sequence from a file in ``[a,b,c,...]`` format.

    Args:
        path: Path to the file to read.

    Returns:
        A list of integers parsed from the file.
    """
    text = path.read_text(encoding="utf-8").strip()
    text = text[1:-1]  # remove surrounding brackets
    if not text:
        return []
    return [int(x) for x in text.split(",")]


def format_factor(
    f1_seq: tuple[int, ...],
    f2_seq: tuple[int, ...],
    label_map: dict[int, str],
) -> str:
    """Format the two sub-sequences of an LZD factor for display.

    Args:
        f1_seq:    First sub-sequence of the factor.
        f2_seq:    Second sub-sequence of the factor (may be empty).
        label_map: Mapping from integer symbol to its display label.

    Returns:
        A string ``"<f1>, <f2>"`` or just ``"<f1>"`` when f2 is empty.
    """
    left = "".join(label_map[s] for s in f1_seq)
    right = "".join(label_map[s] for s in f2_seq)
    return f"{left}, {right}" if right else left


def lzd_compress(
    symbol_sequence: list[int],
) -> tuple[list[dict[str, str]], int]:
    """Run LZD compression on an integer symbol sequence.

    The dictionary is initialised with all unique input symbols.  At each
    step the algorithm picks two consecutive longest-match phrases (f1, f2)
    from the dictionary and records their concatenation as the new factor.
    The concatenated phrase is then registered in the dictionary.

    Args:
        symbol_sequence: The input sequence of integer symbols.

    Returns:
        A tuple ``(factors_details, factor_count)`` where ``factors_details``
        is a list of per-step records for display, and ``factor_count`` is the
        total number of LZD factors.
    """
    unique_symbols = sorted(set(symbol_sequence))
    label_map: dict[int, str] = {s: f"sigma_{s}" for s in unique_symbols}

    # Dictionary keyed by symbol tuples; values are human-readable labels
    dictionary: dict[tuple[int, ...], str] = {(s,): f"sigma_{s}" for s in unique_symbols}

    factors_details: list[dict[str, str]] = []
    seq = tuple(symbol_sequence)
    n = len(seq)
    i = 0
    factor_count = 0

    while i < n:
        # Step 1: find the longest dictionary match f1 starting at position i
        f1_seq: tuple[int, ...] = ()
        f1_lbl = ""
        for d_seq in sorted(dictionary.keys(), key=len, reverse=True):
            if seq[i : i + len(d_seq)] == d_seq:
                f1_seq = d_seq
                f1_lbl = dictionary[d_seq]
                break

        i_next = i + len(f1_seq)

        # Step 2: find the longest dictionary match f2 immediately after f1
        f2_seq: tuple[int, ...] = ()
        f2_lbl = "-"
        if i_next < n:
            for d_seq in sorted(dictionary.keys(), key=len, reverse=True):
                if seq[i_next : i_next + len(d_seq)] == d_seq:
                    f2_seq = d_seq
                    f2_lbl = dictionary[d_seq]
                    break

        # Step 3: record the factor f_i = f1 · f2 and update the dictionary
        factor_count += 1
        combined_seq = f1_seq + f2_seq
        new_id = f"F_{factor_count}"
        if combined_seq not in dictionary:
            dictionary[combined_seq] = new_id

        factors_details.append(
            {
                "Step": str(factor_count),
                "Factor": format_factor(f1_seq, f2_seq, label_map),
                "G1": f1_lbl,
                "G2": f2_lbl,
                "New entry": new_id,
            }
        )

        i += len(combined_seq)

    return factors_details, factor_count


def print_table(rows: list[dict[str, str]], headers: list[str]) -> None:
    """Print a list of row dicts as a plain-text aligned table.

    Args:
        rows:    List of dicts, one per table row.
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
    """Compress one sequence file with LZD and print a detailed breakdown.

    Args:
        input_path: Path to the input sequence file (``[a,b,c,...]`` format).

    Returns:
        The number of LZD factors.
    """
    seq = load_sequence(input_path)
    results, z = lzd_compress(seq)

    print(f"\n=== LZD compression: {input_path} ===")
    print(f"Input T: {''.join(f'sigma_{s}' for s in seq)}")
    print(f"Total factors z: {z}\n")

    headers = ["Step", "Factor", "G1", "G2", "New entry"]
    print("--- LZD compression process ---")
    print_table(results, headers)

    print("\n--- Factorization visualization ---")
    print("".join(f"[{item['Factor']}]" for item in results))

    return z


def main() -> None:
    """Parse arguments and run LZD compression on all available input variants."""
    parser = argparse.ArgumentParser(
        description=(
            "Run LZD compression on sequence files. "
            "Reads INPUT_none.txt, INPUT_del.txt, etc. (whichever exist), "
            "prints per-file details, then emits a RESULT line."
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