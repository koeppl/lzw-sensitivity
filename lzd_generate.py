"""LZD text generator.

Generates the canonical LZD sensitivity test sequence T and its
1-edit (deletion) variant T', then writes them to output files.

Sequence construction follows the paper's definition:
  - addheadpattern builds the initial portion of T using parameter k
  - A middle section and a tail are appended
  - T' = T with its first element removed (1-deletion edit)
"""

from __future__ import annotations

import argparse
from pathlib import Path


def add_value(lst: list[int], j: int) -> None:
    """Append integer j twice to lst (encodes a repeated symbol pair)."""
    lst.append(j)
    lst.append(j)


def add_head_pattern(lst: list[int], k: int) -> None:
    """Append the head pattern for the LZD test sequence with parameter k.

    For i = 1..k appends the pair (k+i, i) each doubled, then appends
    the sentinel symbol 2k+1 doubled.
    """
    for i in range(1, k + 1):
        add_value(lst, k + i)
        add_value(lst, i)
    add_value(lst, 2 * k + 1)


def build_sequence(k: int) -> list[int]:
    """Construct the canonical LZD test sequence T for parameter k.

    The sequence consists of three parts:
      1. A head pattern (add_head_pattern)
      2. A structured middle section
      3. A tail that repeats symbols k..1 each four times
    """
    seq: list[int] = []

    # Build the head pattern
    add_head_pattern(seq, k)

    # Build the middle section
    for i in range(1, k + 1):
        for j in range(i + 2, k + 1):
            # Stop early when the first pair (i=1) would reach the last symbol (j=k)
            if j > k or (i == 1 and j == k):
                break
            add_value(seq, j)
            add_value(seq, i)
        if i + 1 > k:
            break
        add_value(seq, i + 1)

    # Append the repeated symbol 1
    add_value(seq, 1)

    # Build the tail: each symbol i from k down to 1, doubled twice
    for i in range(k, 0, -1):
        add_value(seq, i)
        add_value(seq, i)

    return seq


def write_pattern(seq: list[int], path: Path) -> None:
    """Write a sequence to a file in [a,b,c,...] format."""
    path.write_text("[" + ",".join(str(x) for x in seq) + "]", encoding="utf-8")


def main() -> None:
    """Parse command-line arguments, generate LZD sequences, and write output files."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate LZD test sequences T (none variant) and T' (deletion variant). "
            "Output files are written with suffixes _none.txt and _del.txt."
        )
    )
    parser.add_argument(
        "-o",
        required=True,
        metavar="PREFIX",
        help="Common prefix path for output files (e.g. 'output/lzd' -> 'output/lzd_none.txt')",
    )
    parser.add_argument(
        "-k",
        type=int,
        default=4,
        metavar="K",
        help="Sequence size parameter (default: 4)",
    )
    args = parser.parse_args()

    prefix = Path(args.o)

    # Ensure the output directory exists
    prefix.parent.mkdir(parents=True, exist_ok=True)

    T = build_sequence(args.k)

    # Write the original sequence T (no-edit variant)
    write_pattern(T, Path(str(prefix) + "_none.txt"))

    # Write T' = T with first element deleted (1-deletion variant)
    T_del = T[1:]
    write_pattern(T_del, Path(str(prefix) + "_del.txt"))


if __name__ == "__main__":
    main()
