# LZW/LZD/LZMW Sensitivity Analysis

This repository contains tools for measuring the **sensitivity** of three
classical text compressors — **LZW**, **LZD**, and **LZMW** — with respect
to a single-symbol edit of the input text.

Sensitivity here means: by how much can the number of factors (phrases) of a
factorization change when the input is modified by exactly one
*insertion*, *deletion*, or *substitution*?

---

### tl;dr

Run `benchmark.sh | grep RESULT` to generate all canonical test sequences and 
report the number of factors for all compressors.

---

## Background

Each of the three factorization algorithms builds a phrase dictionary while
parsing the input left-to-right:

| Algorithm | Dictionary update rule |
|-----------|------------------------|
| **LZW** (Lempel–Ziv–Welch) | Longest-match phrase W extended by the next symbol |
| **LZD** (Lempel–Ziv–Double) | Two consecutive longest-match phrases f₁·f₂ |
| **LZMW** (Lempel–Ziv–Miller–Wegman) | Longest previous factor concatenated with its subsequent factor |

For each algorithm the repository provides:
- A **generator** that constructs a worst-case input text T and its 1-edit
  variants (stored as integer-sequence files).
- A **compressor** that reads those files, runs the factorization, and writes
  the factor count for each variant.
- A special **LZW max-sensitivity** script that constructs and compresses a
  string-level worst-case sequence directly (no separate generator needed).

---

## Repository Structure

```
lzd_generate.py      Generate LZD canonical test sequences (replaces lzd_generate.cpp)
lzmw_generate.py     Generate LZMW canonical test sequences (replaces lzmw_generate.cpp)
lzw_generate.py      Generate LZW max-sensitivity test sequences (integer-sequence files)
lzw_compress.py      LZW compression of integer-sequence files
lzd_compress.py      LZD compression of integer-sequence files
lzmw_compress.py     LZMW compression of integer-sequence files
lzw_ms_compress.py   LZW max-sensitivity sequence generator + compressor (string-level)
test_compressors.py  pytest test suite
benchmark.sh         Benchmark script for running all generators and compressors 
```

---

## Requirements

- Python ≥ 3.10
- No third-party packages required (standard library only)
- For tests: `pytest`

```bash
pip install pytest
pytest
```

---

## Usage

### 1. Generate canonical test sequences

```bash
# LZD: produces lzd_none.txt and lzd_del.txt
python lzd_generate.py -o path/to/lzd [-k 4]

# LZMW: produces lzmw_none.txt, lzmw_sub.txt, lzmw_del.txt, lzmw_ins.txt
python lzmw_generate.py -o path/to/lzmw [-k 17]

# LZW: produces lzw_none.txt, lzw_sub.txt, lzw_del.txt, lzw_ins.txt
python lzw_generate.py -o path/to/lzw [-k 4]
```

**Arguments**

| Flag | Description |
|------|-------------|
| `-o PREFIX` | Common prefix path for output files (required) |
| `-k K` | Sequence size parameter (optional; default 4 for LZD, 17 for LZMW) |

The output files use the suffix convention:

| Suffix | Variant |
|--------|---------|
| `_none.txt` | Original sequence T |
| `_sub.txt` | T with one symbol substituted |
| `_del.txt` | T with one symbol deleted |
| `_ins.txt` | T with one symbol "inserted" (replacement at offset +1) |

### 2. Run a compressor

All three compressors share the same interface:

```bash
python lzd_compress.py  -i path/to/lzd
python lzmw_compress.py -i path/to/lzmw
python lzw_compress.py  -i path/to/lzmw
```

**Arguments**

| Flag | Description |
|------|-------------|
| `-i PREFIX` | Common prefix path for input sequence files (required) |

The script reads every `INPUT_PREFIX_<variant>.txt` file that exists, runs
the factorization, prints a detailed per-step table to stdout, and finally
emits a single summary line:

```
RESULT prefix=<PREFIX> num_none=<N> num_sub=<N> num_del=<N> num_ins=<N>
```

Variants whose input file is absent are reported as `None`.

### 3. LZW max-sensitivity (string-level)

```bash
python lzw_ms_compress.py [-k 50] [-o path/to/results]
```

This script constructs the theoretical LZW worst-case sequence T (as defined
in the sensitivity paper) at the string-symbol level, runs LZW on T and all
three 1-edit variants, and prints ANSI-colored factorization visualizations.

With `-o` it also writes factor counts to `_none.txt`, `_sub.txt`, `_del.txt`,
and `_ins.txt` files.

`lzw_generate.py` provides the same sequences in integer-file format,
making it easy to combine with the standard pipeline.

### Full example

```bash
# 1. Generate LZW sequences
python lzw_generate.py -o /tmp/demo/lzw -k 4

# 2. Compress with LZW (prints RESULT line to stdout)
python lzw_compress.py -i /tmp/demo/lzw

# Example output line:
# RESULT prefix=/tmp/demo/lzw num_none=27 num_sub=28 num_del=26 num_ins=28
```

---

## Running the Tests

```bash
pytest test_compressors.py -v
```

The test suite covers:
- Unit tests for each factorization algorithm on simple hand-verified inputs
- Generator correctness tests (file creation, sequence properties)
- End-to-end pipeline tests (generate → compress → check output)

---

## Output Format

**Sequence files** (input to compressors, output from generators):

```
[a,b,c,d,...]
```

A single line containing a comma-separated integer list enclosed in square
brackets.

**RESULT line** (stdout from compress scripts):

```
RESULT prefix=<INPUT_PREFIX> num_none=<N> num_sub=<N> num_del=<N> num_ins=<N>
```

One line printed to stdout after all variants have been processed.  Variants
whose input file does not exist are reported as ``None``.

---

## Citation and Credits

三神 摩周, クップル ドミニク:
LZW の圧縮感度.  Local Proceedings of the LA Symposium Winter 2025. (2026)
