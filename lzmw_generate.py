"""LZMW text generator.

Generates the canonical LZMW sensitivity test sequence T and its
1-edit variants (substitution, deletion, insertion), then writes
them to output files.

The sequence T is divided into three blocks:
  Block A: for i = 1..k, symbols 1..i followed by k+i
  Block B: for j = 1..k, symbols 1..k followed by k+j and alpha_j
  Block C: for j = 1..k, symbols 1..k followed by k+j

The alpha_j symbols are encoded as 1000+j to distinguish them from
ordinary symbols. The 1-edit variants modify the first element of Block C:
  - Substitution (sub): replace T[blockCStart] with the sentinel 999
  - Deletion   (del): delete T[blockCStart]
  - Insertion  (ins): replace T[blockCStart + 1] with the sentinel 999
"""

from __future__ import annotations

import argparse
from pathlib import Path


def alpha_code(j: int) -> int:
    """Return the integer code for symbol alpha_j (alpha_1 -> 1001, alpha_12 -> 1012)."""
    return 1000 + j


def make_sequence(k: int) -> tuple[list[int], int]:
    """Construct the canonical LZMW test sequence T for parameter k.

    Returns:
        A tuple (T, block_c_start) where T is the full integer sequence and
        block_c_start is the index in T where Block C begins.
    """
    seq: list[int] = []

    # Block A: for i = 1..k, emit symbols 1..i then k+i
    for i in range(1, k + 1):
        for t in range(1, i + 1):
            seq.append(t)
        seq.append(k + i)

    # Block B: for j = 1..k, emit symbols 1..k then k+j then alpha_j
    for j in range(1, k + 1):
        for t in range(1, k + 1):
            seq.append(t)
        seq.append(k + j)
        seq.append(alpha_code(j))

    # Record the start index of Block C
    block_c_start: int = len(seq)

    # Block C: for j = 1..k, emit symbols 1..k then k+j
    for j in range(1, k + 1):
        for t in range(1, k + 1):
            seq.append(t)
        seq.append(k + j)

    return seq, block_c_start


def write_sequence(seq: list[int], path: Path) -> None:
    """Write a sequence to a file in [a,b,c,...] format."""
    path.write_text("[" + ",".join(str(x) for x in seq) + "]\n", encoding="utf-8")


def main() -> None:
    """Parse command-line arguments, generate LZMW sequences, and write output files."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate LZMW test sequence T and its 1-edit variants "
            "(substitution, deletion, insertion). "
            "Output files are written with suffixes _none.txt, _sub.txt, "
            "_del.txt, and _ins.txt."
        )
    )
    parser.add_argument(
        "-o",
        required=True,
        metavar="PREFIX",
        help="Common prefix path for output files (e.g. 'output/lzmw' -> 'output/lzmw_none.txt')",
    )
    parser.add_argument(
        "-k",
        type=int,
        default=17,
        metavar="K",
        help="Sequence size parameter (default: 17)",
    )
    args = parser.parse_args()

    prefix = Path(args.o)

    # Ensure the output directory exists
    prefix.parent.mkdir(parents=True, exist_ok=True)

    T, block_c_start = make_sequence(args.k)

    # Write original sequence T (no-edit variant)
    write_sequence(T, Path(str(prefix) + "_none.txt"))

    # Substitution variant: replace the first symbol of Block C with sentinel 999
    T_sub = T.copy()
    if block_c_start < len(T_sub):
        T_sub[block_c_start] = 999
    write_sequence(T_sub, Path(str(prefix) + "_sub.txt"))

    # Deletion variant: remove the first symbol of Block C
    T_del = T.copy()
    if block_c_start < len(T_del):
        del T_del[block_c_start]
    write_sequence(T_del, Path(str(prefix) + "_del.txt"))

    # Insertion variant: replace the second symbol of Block C with sentinel 999.
    # This models the effect of inserting a new symbol immediately before the
    # second Block-C symbol: the second position now holds an alien value 999
    # instead of the expected symbol 2, which is the kind of 1-edit perturbation
    # studied in the sensitivity analysis.  Note that the sequence length is
    # unchanged (this is a substitution at offset +1, not a literal insertion).
    T_ins = T.copy()
    if block_c_start + 1 < len(T_ins):
        T_ins[block_c_start + 1] = 999
    write_sequence(T_ins, Path(str(prefix) + "_ins.txt"))


if __name__ == "__main__":
    main()
