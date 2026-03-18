"""LZW text generator.

Generates the canonical LZW max-sensitivity test sequence T and its
1-edit variants (substitution, deletion, insertion), then writes them
to output files compatible with ``lzw_compress.py``.

Sequence construction follows the definition in the max-sensitivity paper
(mirroring the ``construct_T`` function in ``lzw_ms_compress.py``) but uses
integer symbols throughout:

  - ``sigma_i`` is encoded as the integer ``i`` (1 ≤ i ≤ 3k)
  - The sentinel symbol used for edits is encoded as ``999``

Part (1): for i = 1..k emit 1..i followed by k+i.
Part (2): for j = 1..k emit core · (2k+j) · core, where
          core = 1..y_j · (k+j), and y_j is computed from L_j.
          When the edit is applied (j == 1 only) the first core is modified.

Output files follow the same convention as the other generators:

  ``PREFIX_none.txt``  — original T (no edit)
  ``PREFIX_sub.txt``   — substitution edit: replace first element of first core with 999
  ``PREFIX_del.txt``   — deletion edit: remove first element of first core
  ``PREFIX_ins.txt``   — insertion edit: insert 999 after first element of first core
"""

from __future__ import annotations

import argparse
from pathlib import Path

# Integer used as a sentinel symbol for edits (not a valid sigma_i value)
SENTINEL: int = 999


def compute_l_j(j: int) -> int:
    """Return L_j, the smallest L such that L*(L+1)/2 >= j.

    Args:
        j: Block index (1-based).

    Returns:
        The smallest positive integer L satisfying L*(L+1)//2 >= j.
    """
    L = 1
    while (L * (L + 1)) // 2 < j:
        L += 1
    return L


def compute_y_j(j: int, k: int) -> tuple[int, int]:
    """Return (y_j, L_j) for the j-th repetition block.

    Args:
        j: Block index (1-based).
        k: Sequence size parameter.

    Returns:
        A tuple ``(y_j, L_j)`` where y_j is the largest candidate in 1..k
        whose residue modulo L_j equals (2+j+L_j-1) mod L_j.
    """
    L = compute_l_j(j)
    residue = (2 + j + L - 1) % L
    candidates = [y for y in range(1, k + 1) if y % L == residue]
    return max(candidates), L


def construct_T_int(k: int, substitute: str) -> list[int]:
    """Construct the LZW test sequence T (or a 1-edit variant) as an integer list.

    ``sigma_i`` is encoded as integer ``i``.  The edit sentinel is ``SENTINEL``.

    Part (1): for i = 1..k emit [1, 2, ..., i, k+i].
    Part (2): for j = 1..k emit core + [2k+j] + core, where
              core = [1, ..., y_j, k+j].
              When ``substitute != 'none'`` the first occurrence of core
              (at j == 1) is replaced according to the edit type.

    Args:
        k:          Sequence size parameter (positive integer).
        substitute: Edit to apply to the first core of part (2).
                    One of ``'none'``, ``'insert'``, ``'delete'``,
                    ``'substitute'``.

    Returns:
        A list of integers representing T (or its edited variant).
    """
    # Part (1)
    part1: list[int] = []
    for i in range(1, k + 1):
        part1.extend(range(1, i + 1))
        part1.append(k + i)

    # Part (2)
    part2: list[int] = []
    for j in range(1, k + 1):
        yj, _Lj = compute_y_j(j, k)
        core: list[int] = list(range(1, yj + 1)) + [k + j]

        if substitute != "none" and j == 1:
            # Apply the requested 1-edit to the first core only
            if substitute == "insert":
                first_core: list[int] = [core[0], SENTINEL] + core[1:]
            elif substitute == "delete":
                first_core = core[1:]
            else:  # substitute == 'substitute'
                first_core = [SENTINEL] + core[1:]
            part2.extend(first_core)
            part2.append(2 * k + j)
            part2.extend(core)
        else:
            part2.extend(core)
            part2.append(2 * k + j)
            part2.extend(core)

    return part1 + part2


def write_sequence(seq: list[int], path: Path) -> None:
    """Write a sequence to a file in ``[a,b,c,...]`` format.

    Args:
        seq:  The integer sequence to write.
        path: Destination file path.
    """
    path.write_text("[" + ",".join(str(x) for x in seq) + "]", encoding="utf-8")


def main() -> None:
    """Parse command-line arguments, generate LZW sequences, and write output files."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate LZW max-sensitivity test sequences T and its 1-edit variants. "
            "Output files are written with suffixes _none.txt, _sub.txt, _del.txt, "
            "and _ins.txt.  Files are compatible with lzw_compress.py."
        )
    )
    parser.add_argument(
        "-o",
        required=True,
        metavar="PREFIX",
        help=(
            "Common prefix path for output files "
            "(e.g. 'output/lzw' -> 'output/lzw_none.txt', ...)"
        ),
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
    prefix.parent.mkdir(parents=True, exist_ok=True)

    # Map (substitute arg, file suffix)
    variants: list[tuple[str, str]] = [
        ("none", "none"),
        ("substitute", "sub"),
        ("delete", "del"),
        ("insert", "ins"),
    ]

    for sub, suffix in variants:
        seq = construct_T_int(args.k, sub)
        write_sequence(seq, Path(f"{prefix}_{suffix}.txt"))


if __name__ == "__main__":
    main()
