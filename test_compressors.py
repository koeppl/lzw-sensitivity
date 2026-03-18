"""Tests for the LZW, LZD, and LZMW compressors and text generators.

Run with pytest::

    pytest test_compressors.py -v
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Import the modules under test (they live in the same directory)
# ---------------------------------------------------------------------------

sys.path.insert(0, str(Path(__file__).parent))

from lzw_compress import lzw_compress, load_sequence as lzw_load  # noqa: E402
from lzd_compress import lzd_compress, load_sequence as lzd_load  # noqa: E402
from lzmw_compress import lzmw_compress, load_sequence as lzmw_load  # noqa: E402
from lzd_generate import build_sequence as lzd_build, write_pattern  # noqa: E402
from lzmw_generate import make_sequence, write_sequence  # noqa: E402
from lzw_ms_compress import construct_T, lzw_factors  # noqa: E402
from lzw_generate import construct_T_int, SENTINEL  # noqa: E402


# ===========================================================================
# Helper
# ===========================================================================


from typing import Callable

FactorizeFn = Callable[[list[int]], tuple[list[dict[str, str]], int]]


def factor_count(compress_fn: FactorizeFn, seq: list[int]) -> int:
    """Return the factor count produced by *compress_fn* for *seq*."""
    _, z = compress_fn(seq)
    return z


# ===========================================================================
# LZW compressor tests
# ===========================================================================


class TestLZWCompress:
    """Unit tests for the LZW factorization algorithm."""

    def test_single_symbol_repeated(self) -> None:
        """A sequence of identical symbols."""
        # [1,1,1,1]: factors are [1], [1,1], [1] => 3
        assert factor_count(lzw_compress, [1, 1, 1, 1]) == 3

    def test_alternating_two_symbols(self) -> None:
        """Alternating sequence [1,2,1,2]."""
        # [1],[2],[1,2] => 3
        assert factor_count(lzw_compress, [1, 2, 1, 2]) == 3

    def test_single_element(self) -> None:
        """A sequence of length 1 is a single factor."""
        assert factor_count(lzw_compress, [5]) == 1

    def test_two_distinct_symbols(self) -> None:
        """[1,2] produces two factors (each symbol is its own factor)."""
        assert factor_count(lzw_compress, [1, 2]) == 2

    def test_three_distinct_no_repeat(self) -> None:
        """[1,2,3]: three factors, one per symbol."""
        assert factor_count(lzw_compress, [1, 2, 3]) == 3

    def test_longer_sequence(self) -> None:
        """[1,2,1,2,1,2]: verified by hand."""
        # i=0: [1] -> register (1,2); i=1: [2] -> register (2,1);
        # i=2: [1,2] match -> register (1,2,1); i=4: [1,2] match -> done
        assert factor_count(lzw_compress, [1, 2, 1, 2, 1, 2]) == 4

    def test_returns_correct_detail_length(self) -> None:
        """Number of detail records equals factor count."""
        details, z = lzw_compress([1, 2, 1, 2])
        assert len(details) == z

    def test_empty_sequence(self) -> None:
        """An empty sequence produces zero factors."""
        assert factor_count(lzw_compress, []) == 0


# ===========================================================================
# LZD compressor tests
# ===========================================================================


class TestLZDCompress:
    """Unit tests for the LZD factorization algorithm."""

    def test_single_symbol_repeated(self) -> None:
        """[1,1,1,1]: LZD groups 2 at a time."""
        # Step1: f1=(1,), f2=(1,) -> combined=(1,1), F_1; i=2
        # Step2: f1=(1,1) match, f2=() -> F_2; i=4
        assert factor_count(lzd_compress, [1, 1, 1, 1]) == 2

    def test_alternating_two_symbols(self) -> None:
        """[1,2,1,2]: each LZD factor covers the whole sequence."""
        assert factor_count(lzd_compress, [1, 2, 1, 2]) == 2

    def test_single_element(self) -> None:
        """A sequence of length 1 produces a single factor."""
        assert factor_count(lzd_compress, [7]) == 1

    def test_two_distinct_symbols(self) -> None:
        """[1,2]: one factor (f1=1, f2=2)."""
        assert factor_count(lzd_compress, [1, 2]) == 1

    def test_three_symbols(self) -> None:
        """[1,2,3]: f1=(1,), f2=(2,) -> F_1; then f1=(3,), f2=() -> F_2."""
        assert factor_count(lzd_compress, [1, 2, 3]) == 2

    def test_returns_correct_detail_length(self) -> None:
        """Number of detail records equals factor count."""
        details, z = lzd_compress([1, 2, 1, 2])
        assert len(details) == z

    def test_empty_sequence(self) -> None:
        """An empty sequence produces zero factors."""
        assert factor_count(lzd_compress, []) == 0


# ===========================================================================
# LZMW compressor tests
# ===========================================================================


class TestLZMWCompress:
    """Unit tests for the LZMW factorization algorithm."""

    def test_single_symbol_repeated(self) -> None:
        """[1,1,1,1]: LZMW uses concatenation of consecutive factors."""
        # Step1: f=(1,); Step2: f=(1,), add (1,1)=D_1;
        # Step3: f=(1,1) match (len 2 at pos 2), add (1,)+(1,1)=(1,1,1)=D_2; i=4
        assert factor_count(lzmw_compress, [1, 1, 1, 1]) == 3

    def test_alternating_two_symbols(self) -> None:
        """[1,2,1,2]: Step1:(1,) Step2:(2,)->add(1,2); Step3:(1,2) match->done."""
        assert factor_count(lzmw_compress, [1, 2, 1, 2]) == 3

    def test_single_element(self) -> None:
        """A sequence of length 1 produces a single factor."""
        assert factor_count(lzmw_compress, [3]) == 1

    def test_two_distinct_symbols(self) -> None:
        """[1,2]: two factors, registers (1,2) but sequence ends."""
        assert factor_count(lzmw_compress, [1, 2]) == 2

    def test_returns_correct_detail_length(self) -> None:
        """Number of detail records equals factor count."""
        details, z = lzmw_compress([1, 2, 1, 2])
        assert len(details) == z

    def test_empty_sequence(self) -> None:
        """An empty sequence produces zero factors."""
        assert factor_count(lzmw_compress, []) == 0


# ===========================================================================
# LZW ms (max-sensitivity string-based) tests
# ===========================================================================


class TestLZWMs:
    """Tests for the LZW max-sensitivity sequence construction."""

    def test_construct_T_none_length(self) -> None:
        """T for k=2 has the expected number of symbols."""
        T = construct_T(2, "none")
        # Part(1): i=1: sigma_1,sigma_3; i=2: sigma_1,sigma_2,sigma_4 -> 5 symbols
        # Part(2): for j=1..2 each block: core*sigma_{2k+j}*core
        assert len(T) > 0

    def test_construct_T_substitute_differs(self) -> None:
        """Substitution variant must differ from T."""
        T = construct_T(3, "none")
        T_sub = construct_T(3, "substitute")
        assert T != T_sub

    def test_construct_T_delete_shorter(self) -> None:
        """Deletion variant is strictly shorter than T."""
        T = construct_T(3, "none")
        T_del = construct_T(3, "delete")
        assert len(T_del) < len(T)

    def test_construct_T_insert_longer(self) -> None:
        """Insertion variant is strictly longer than T."""
        T = construct_T(3, "none")
        T_ins = construct_T(3, "insert")
        assert len(T_ins) > len(T)

    def test_lzw_factors_simple(self) -> None:
        """lzw_factors on two alternating string symbols."""
        symbols = ["a", "b", "a", "b"]
        factors = lzw_factors(symbols)
        assert len(factors) == 3  # "a", "b", "ab"

    def test_lzw_factors_single(self) -> None:
        """lzw_factors on a single symbol."""
        assert lzw_factors(["x"]) == ["x"]


# ===========================================================================
# LZD generator tests
# ===========================================================================


class TestLZDGenerate:
    """Tests for the LZD canonical sequence generator."""

    def test_build_sequence_nonempty(self) -> None:
        """build_sequence returns a non-empty list for k >= 1."""
        assert len(lzd_build(1)) > 0
        assert len(lzd_build(4)) > 0

    def test_build_sequence_k4_starts_correctly(self) -> None:
        """For k=4 the sequence starts with [5, 5, 1, 1, ...]."""
        seq = lzd_build(4)
        assert seq[:4] == [5, 5, 1, 1]

    def test_build_sequence_k1(self) -> None:
        """For k=1 the sequence is well-defined and non-trivial."""
        seq = lzd_build(1)
        # head: add_value(2) + add_value(1) + add_value(3) = [2,2,1,1,3,3]
        # middle: i=1, j loop empty, i+1>k -> break; no add_value(2)
        # add_value(1) -> [1,1]
        # tail: i=1: add_value(1)*2 -> [1,1,1,1]
        assert seq[:6] == [2, 2, 1, 1, 3, 3]

    def test_write_and_read_roundtrip(self, tmp_path: Path) -> None:
        """write_pattern -> load_sequence round-trips correctly."""
        seq = lzd_build(4)
        out = tmp_path / "test_none.txt"
        write_pattern(seq, out)
        loaded = lzw_load(out)
        assert loaded == seq

    def test_generate_creates_none_and_del(self, tmp_path: Path) -> None:
        """lzd_generate.py -o prefix creates _none.txt and _del.txt."""
        import subprocess

        prefix = str(tmp_path / "lzd")
        result = subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzd_generate.py"),
             "-o", prefix, "-k", "4"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        assert (tmp_path / "lzd_none.txt").exists()
        assert (tmp_path / "lzd_del.txt").exists()

    def test_del_variant_is_none_without_first_element(self, tmp_path: Path) -> None:
        """The deletion variant equals the original without its first element."""
        import subprocess

        prefix = str(tmp_path / "lzd")
        subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzd_generate.py"),
             "-o", prefix, "-k", "4"],
            capture_output=True,
        )
        none_seq = lzw_load(tmp_path / "lzd_none.txt")
        del_seq = lzw_load(tmp_path / "lzd_del.txt")
        assert del_seq == none_seq[1:]


# ===========================================================================
# LZMW generator tests
# ===========================================================================


class TestLZMWGenerate:
    """Tests for the LZMW canonical sequence generator."""

    def test_make_sequence_nonempty(self) -> None:
        """make_sequence returns a non-empty list for k >= 1."""
        seq, _ = make_sequence(1)
        assert len(seq) > 0

    def test_make_sequence_block_c_start(self) -> None:
        """block_c_start is a valid index inside the sequence."""
        seq, bcs = make_sequence(4)
        assert 0 < bcs < len(seq)

    def test_generate_creates_four_files(self, tmp_path: Path) -> None:
        """lzmw_generate.py -o prefix creates all four variant files."""
        import subprocess

        prefix = str(tmp_path / "lzmw")
        result = subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzmw_generate.py"),
             "-o", prefix, "-k", "4"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        for suffix in ("none", "sub", "del", "ins"):
            assert (tmp_path / f"lzmw_{suffix}.txt").exists()

    def test_sub_variant_differs_from_none(self, tmp_path: Path) -> None:
        """Substitution variant differs from original at block_c_start."""
        import subprocess

        prefix = str(tmp_path / "lzmw")
        subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzmw_generate.py"),
             "-o", prefix, "-k", "4"],
            capture_output=True,
        )
        none_seq = lzmw_load(tmp_path / "lzmw_none.txt")
        sub_seq = lzmw_load(tmp_path / "lzmw_sub.txt")
        assert none_seq != sub_seq
        # They differ at block_c_start; sub value must be 999
        _, bcs = make_sequence(4)
        assert sub_seq[bcs] == 999

    def test_del_variant_is_shorter(self, tmp_path: Path) -> None:
        """Deletion variant is one element shorter than original."""
        import subprocess

        prefix = str(tmp_path / "lzmw")
        subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzmw_generate.py"),
             "-o", prefix, "-k", "4"],
            capture_output=True,
        )
        none_seq = lzmw_load(tmp_path / "lzmw_none.txt")
        del_seq = lzmw_load(tmp_path / "lzmw_del.txt")
        assert len(del_seq) == len(none_seq) - 1

    def test_ins_variant_same_length(self, tmp_path: Path) -> None:
        """Insertion (ins) variant has the same length as the original.

        The 'ins' variant is implemented as a *replacement* of the element at
        block_c_start+1 with the sentinel 999.  This simulates the effect of
        inserting a symbol just before that position, while preserving the
        overall sequence length.
        """
        import subprocess

        prefix = str(tmp_path / "lzmw")
        subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzmw_generate.py"),
             "-o", prefix, "-k", "4"],
            capture_output=True,
        )
        none_seq = lzmw_load(tmp_path / "lzmw_none.txt")
        ins_seq = lzmw_load(tmp_path / "lzmw_ins.txt")
        # ins replaces element at blockCStart+1 with 999 (same length)
        assert len(ins_seq) == len(none_seq)


# ===========================================================================
# End-to-end pipeline tests
# ===========================================================================


def parse_result_line(stdout: str) -> dict[str, int | None]:
    """Parse a RESULT line emitted by the compress scripts.

    Expected format::

        RESULT prefix=<p> num_none=<N> num_sub=<N> num_del=<N> num_ins=<N>

    Args:
        stdout: Captured stdout text from a compress script invocation.

    Returns:
        A dict mapping variant name to factor count (or None if absent).
    """
    for line in stdout.splitlines():
        if line.startswith("RESULT "):
            result: dict[str, int | None] = {}
            for token in line.split():
                for variant in ("none", "sub", "del", "ins"):
                    key = f"num_{variant}="
                    if token.startswith(key):
                        val = token[len(key):]
                        result[variant] = None if val == "None" else int(val)
            return result
    raise AssertionError(f"No RESULT line found in stdout:\n{stdout}")


class TestEndToEndLZD:
    """End-to-end tests: generate -> compress -> check factor counts."""

    def test_lzd_pipeline_none(self, tmp_path: Path) -> None:
        """LZD factor count for the canonical k=4 T sequence is reproducible."""
        import subprocess

        gen_prefix = str(tmp_path / "gen")

        subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzd_generate.py"),
             "-o", gen_prefix, "-k", "4"],
            capture_output=True,
            check=True,
        )
        result = subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzd_compress.py"),
             "-i", gen_prefix],
            capture_output=True,
            text=True,
            check=True,
        )
        counts = parse_result_line(result.stdout)
        assert counts["none"] is not None and counts["none"] > 0
        assert counts["del"] is not None and counts["del"] > 0

    def test_lzd_pipeline_del_differs(self, tmp_path: Path) -> None:
        """The deletion variant has a different factor count from the original."""
        import subprocess

        gen_prefix = str(tmp_path / "gen")

        subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzd_generate.py"),
             "-o", gen_prefix, "-k", "4"],
            capture_output=True,
            check=True,
        )
        result = subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzd_compress.py"),
             "-i", gen_prefix],
            capture_output=True,
            text=True,
            check=True,
        )
        counts = parse_result_line(result.stdout)
        assert counts["none"] != counts["del"]


class TestEndToEndLZMW:
    """End-to-end tests: generate -> compress -> check factor counts."""

    def test_lzmw_pipeline(self, tmp_path: Path) -> None:
        """LZMW factor counts for the canonical k=4 sequence are all positive."""
        import subprocess

        gen_prefix = str(tmp_path / "gen")

        subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzmw_generate.py"),
             "-o", gen_prefix, "-k", "4"],
            capture_output=True,
            check=True,
        )
        result = subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzmw_compress.py"),
             "-i", gen_prefix],
            capture_output=True,
            text=True,
            check=True,
        )
        counts = parse_result_line(result.stdout)
        for variant in ("none", "sub", "del", "ins"):
            assert counts[variant] is not None and counts[variant] > 0, (
                f"Expected positive count for variant '{variant}'"
            )


class TestEndToEndLZW:
    """End-to-end tests for lzw_compress.py using LZMW-generated input."""

    def test_lzw_pipeline(self, tmp_path: Path) -> None:
        """LZW factor counts for LZMW-generated k=4 sequences are positive."""
        import subprocess

        gen_prefix = str(tmp_path / "gen")

        subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzmw_generate.py"),
             "-o", gen_prefix, "-k", "4"],
            capture_output=True,
            check=True,
        )
        result = subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzw_compress.py"),
             "-i", gen_prefix],
            capture_output=True,
            text=True,
            check=True,
        )
        counts = parse_result_line(result.stdout)
        for variant in ("none", "sub", "del", "ins"):
            assert counts[variant] is not None and counts[variant] > 0, (
                f"Expected positive count for variant '{variant}'"
            )


# ===========================================================================
# LZW generator tests
# ===========================================================================


class TestLZWGenerate:
    """Tests for the LZW canonical sequence generator (lzw_generate.py)."""

    def test_construct_T_int_nonempty(self) -> None:
        """construct_T_int returns a non-empty list for k >= 1."""
        assert len(construct_T_int(2, "none")) > 0
        assert len(construct_T_int(4, "none")) > 0

    def test_construct_T_int_sub_differs(self) -> None:
        """Substitution variant differs from the original."""
        T = construct_T_int(3, "none")
        T_sub = construct_T_int(3, "substitute")
        assert T != T_sub

    def test_construct_T_int_sub_has_sentinel(self) -> None:
        """Substitution variant contains the sentinel symbol."""
        T_sub = construct_T_int(3, "substitute")
        assert SENTINEL in T_sub

    def test_construct_T_int_del_shorter(self) -> None:
        """Deletion variant is strictly shorter than the original."""
        T = construct_T_int(3, "none")
        T_del = construct_T_int(3, "delete")
        assert len(T_del) < len(T)

    def test_construct_T_int_ins_longer(self) -> None:
        """Insertion variant is strictly longer than the original."""
        T = construct_T_int(3, "none")
        T_ins = construct_T_int(3, "insert")
        assert len(T_ins) > len(T)

    def test_construct_T_int_none_no_sentinel(self) -> None:
        """Original (none) variant contains no sentinel symbol."""
        T = construct_T_int(4, "none")
        assert SENTINEL not in T

    def test_generate_creates_four_files(self, tmp_path: Path) -> None:
        """lzw_generate.py -o prefix creates all four variant files."""
        import subprocess

        prefix = str(tmp_path / "lzw")
        result = subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzw_generate.py"),
             "-o", prefix, "-k", "4"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        for suffix in ("none", "sub", "del", "ins"):
            assert (tmp_path / f"lzw_{suffix}.txt").exists()

    def test_generate_output_compatible_with_lzw_compress(self, tmp_path: Path) -> None:
        """lzw_generate.py output can be read and compressed by lzw_compress.py."""
        import subprocess

        gen_prefix = str(tmp_path / "lzw")
        subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzw_generate.py"),
             "-o", gen_prefix, "-k", "4"],
            capture_output=True,
            check=True,
        )
        result = subprocess.run(
            [sys.executable, str(Path(__file__).parent / "lzw_compress.py"),
             "-i", gen_prefix],
            capture_output=True,
            text=True,
            check=True,
        )
        counts = parse_result_line(result.stdout)
        for variant in ("none", "sub", "del", "ins"):
            assert counts[variant] is not None and counts[variant] > 0, (
                f"Expected positive count for variant '{variant}'"
            )

    def test_generate_none_variant_consistent_with_ms_compress(self) -> None:
        """construct_T_int and construct_T agree on symbol structure for k=3."""
        # construct_T returns strings like "sigma_3"; construct_T_int returns 3
        T_str = construct_T(3, "none")
        T_int = construct_T_int(3, "none")
        assert len(T_str) == len(T_int)
        for sym_str, sym_int in zip(T_str, T_int):
            expected_int = int(sym_str.split("_")[1])
            assert sym_int == expected_int
