"""LZW Max-Sensitivity sequence generator and compressor.

Constructs the theoretical LZW text T (and its 1-edit variants) as defined
in the max-sensitivity paper, runs LZW factorization on each, and prints a
colour-highlighted visualization of the factorization.

Unlike the other compress scripts this file generates the text internally
(no external generate script is required).  The text consists of string
symbols (e.g. ``sigma_1``) rather than bare integers.

Usage
-----
    python lzw_ms_compress.py [-k K] [-o OUTPUT_PREFIX]

With ``-o`` the factor count for each variant is also written to
``OUTPUT_PREFIX_none.txt``, ``OUTPUT_PREFIX_sub.txt``, etc.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]


def compute_l_j(j: int) -> int:
    """Return L_j, the smallest L such that L*(L+1)/2 >= j.

    Args:
        j: The index j >= 1.

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
    candidates = [y for y in range(1, k + 1) if y % L == residue % L]
    return max(candidates), L


def construct_T(k: int, substitute: str) -> list[str]:
    """Construct the theoretical LZW text T (or a 1-edit variant).

    Part (1): for i = 1..k emit sigma_1..sigma_i followed by sigma_{k+i}.
    Part (2): for j = 1..k emit core * sigma_{2k+j} * core, where
              core = sigma_1..sigma_{y_j} * sigma_{k+j}.
              When ``substitute != 'none'`` and j == 1 the first occurrence
              of the core is replaced according to the edit type.

    Args:
        k:          Sequence size parameter.
        substitute: Edit to apply to the first core of part (2).
                    One of ``'none'``, ``'insert'``, ``'delete'``,
                    ``'substitute'``.

    Returns:
        A list of symbol strings representing T (or its variant).
    """
    sigma = [f"sigma_{i}" for i in range(1, 3 * k + 1)]

    # Part (1): prefix blocks (sigma_1..sigma_i) * sigma_{k+i}
    part1: list[str] = []
    for i in range(1, k + 1):
        part1.extend(sigma[:i])
        part1.append(sigma[k + i - 1])

    # Part (2): repeated core blocks with optional edit at j=1
    part2: list[str] = []
    y_values: list[int] = []
    for j in range(1, k + 1):
        yj, _Lj = compute_y_j(j, k)
        y_values.append(yj)
        core = sigma[:yj] + [sigma[k + j - 1]]
        if substitute != "none" and j == 1:
            # Apply the requested 1-edit to the first core only
            if substitute == "insert":
                core_sub: list[str] = [core[0]] + ["#"] + core[1:]
            elif substitute == "delete":
                core_sub = core[1:]
            else:  # substitute == 'substitute'
                core_sub = ["#"] + core[1:]
            assert core_sub, "Edited core must not be empty."
            part2.extend(core_sub)
            part2.append(sigma[2 * k + j - 1])
            part2.extend(core)
        else:
            part2.extend(core)
            part2.append(sigma[2 * k + j - 1])
            part2.extend(core)

    label = "T' (edited)" if substitute != "none" else "T (original)"
    print(f"\n=== {label} construction ===")
    print(f"y_j values = {y_values}")
    return part1 + part2


def lzw_factors(symbols: list[str]) -> list[str]:
    """Run LZW factorization on a list of string symbols.

    The dictionary is initialised with the unique symbols.  At each step
    the greedy longest match is output as a factor and the match extended by
    the next symbol is registered as a new dictionary entry.

    Args:
        symbols: The input sequence of string symbols.

    Returns:
        A list of factor strings (each factor is a concatenation of symbols).
    """
    dictionary: dict[str, int] = {
        sym: i for i, sym in enumerate(sorted(set(symbols)), start=1)
    }
    w = ""
    factors: list[str] = []
    for c in symbols:
        wc = w + c
        if wc in dictionary:
            w = wc
        else:
            factors.append(w)
            dictionary[wc] = len(dictionary) + 1
            w = c
    if w:
        factors.append(w)
    return factors


def visualize_factors_T(factors: list[str], k: int) -> str:
    """Format the LZW factorization of T for display.

    The first ``3k-1`` factors are placed on one line; subsequent factors
    are grouped in rows of four.

    Args:
        factors: The list of LZW factors.
        k:       Sequence size parameter (controls line-break positions).

    Returns:
        A multi-line string ready for printing.
    """
    lines: list[str] = []
    first_block_len = 3 * k - 1
    lines.append("|".join(factors[:first_block_len]) + "|")

    current: list[str] = []
    for f in factors[first_block_len:]:
        current.append(f)
        if len(current) == 4:
            lines.append("|".join(current) + "|")
            current = []
    if current:
        lines.append("|".join(current))

    return "\n".join(lines)


def visualize_factors_Tp(factors: list[str], k: int) -> str:
    """Format the LZW factorization of T' with ANSI colour highlighting.

    Colours indicate:
      - Cyan    : factors ending with sigma_1
      - Yellow  : factors immediately following a sigma_1-ending factor
      - Magenta : factors immediately before a large-index (>= 2k+2) factor
      - Red     : factors containing a symbol with index >= 2k+2
      - Green   : factors immediately after a large-index factor

    Args:
        factors: The list of LZW factors.
        k:       Sequence size parameter.

    Returns:
        A multi-line ANSI-coloured string ready for printing.
    """
    RESET = "\033[0m"
    CYAN = "\033[46m"
    YELLOW = "\033[43m"
    MAGENTA = "\033[45m"
    RED = "\033[41m"
    GREEN = "\033[42m"

    lines: list[str] = []
    first_block_len = 3 * k - 1
    lines.append("|".join(factors[:first_block_len]) + "|")

    # Identify positions of factors containing symbols with index >= 2k+2
    prev_big: set[int] = set()
    big: set[int] = set()
    next_big: set[int] = set()
    for idx, f in enumerate(factors):
        for part in f.split("sigma_")[1:]:
            try:
                sym_idx = int(part)
                if sym_idx >= 2 * k + 2:
                    if idx - 1 >= 0:
                        prev_big.add(idx - 1)
                    big.add(idx)
                    if idx + 1 < len(factors):
                        next_big.add(idx + 1)
                    break
            except ValueError:
                pass

    current_line = ""
    highlight_next_sigma1 = False

    for idx in range(first_block_len, len(factors)):
        f = factors[idx]
        f_out = f

        # Yellow: factor immediately after a sigma_1-ending factor
        if highlight_next_sigma1:
            f_out = f"{YELLOW}{f}{RESET}"
            highlight_next_sigma1 = False

        # Colour factors near large-index symbols (highest priority)
        if idx in prev_big:
            f_out = f"{MAGENTA}{f}{RESET}"
        elif idx in big:
            f_out = f"{RED}{f}{RESET}"
        elif idx in next_big:
            f_out = f"{GREEN}{f}{RESET}"

        # Cyan + line break: factor ending with sigma_1
        if f.endswith("sigma_1"):
            f_out = f"{CYAN}{f}{RESET}"
            current_line += f_out + "|"
            lines.append(current_line)
            current_line = ""
            highlight_next_sigma1 = True
            continue

        current_line += f_out + "|"

    if current_line:
        lines.append(current_line)

    return "\n".join(lines)


def main() -> None:
    """Parse arguments, construct T and its variants, run LZW, and print results."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate the LZW max-sensitivity sequence T and its 1-edit variants, "
            "run LZW factorization, and visualize the results."
        )
    )
    parser.add_argument(
        "-k",
        type=int,
        default=50,
        metavar="K",
        help="Sequence size parameter (default: 50)",
    )
    parser.add_argument(
        "-o",
        metavar="OUTPUT_PREFIX",
        default=None,
        help=(
            "Optional common prefix for output factor-count files. "
            "Writes OUTPUT_PREFIX_none.txt, OUTPUT_PREFIX_sub.txt, etc."
        ),
    )
    args = parser.parse_args()
    k: int = args.k
    out_prefix: str | None = args.o

    substitutions = ["none", "insert", "delete", "substitute"]
    # Map edit names to output file suffixes
    suffix_map: dict[str, str] = {
        "none": "none",
        "insert": "ins",
        "delete": "del",
        "substitute": "sub",
    }

    # Prepare output directory if needed
    if out_prefix is not None:
        Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)

    # Construct and compress T (no edit)
    T = construct_T(k, substitute="none")
    print(" ".join(T))
    factors_T = lzw_factors(T)

    print("\n=== LZW factorization: T ===")
    print(visualize_factors_T(factors_T, k))
    print(f"\nFactor count z_LZW(T) = {len(factors_T)}")

    if out_prefix is not None:
        Path(f"{out_prefix}_none.txt").write_text(str(len(factors_T)), encoding="utf-8")

    # Construct, compress, and report each 1-edit variant
    for sub in substitutions[1:]:
        T_prime = construct_T(k, substitute=sub)
        print(f"\n--- Edit type: {sub} ---")
        print(" ".join(T_prime))
        factors_Tp = lzw_factors(T_prime)

        print("\n=== LZW factorization: T' ===")
        print(visualize_factors_Tp(factors_Tp, k))
        print(f"\nFactor count z_LZW(T') = {len(factors_Tp)}")

        diff = len(factors_Tp) - len(factors_T)
        ratio = len(factors_Tp) / len(factors_T)
        print("\n=== Sensitivity comparison ===")
        print(f"Factor count difference: {diff}")
        print(f"Ratio z_LZW(T') / z_LZW(T) = {ratio:.3f}")

        if out_prefix is not None:
            suffix = suffix_map[sub]
            Path(f"{out_prefix}_{suffix}.txt").write_text(
                str(len(factors_Tp)), encoding="utf-8"
            )


if __name__ == "__main__":
    main()
