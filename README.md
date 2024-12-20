[![Licence](https://img.shields.io/badge/License-MIT-9cf.svg?style=flat-square)](https://choosealicense.com/licenses/mit/)
[![Test](https://img.shields.io/github/actions/workflow/status/akikuno/midsv/ci.yml?branch=main&label=Test&color=brightgreen&style=flat-square)](https://github.com/akikuno/midsv/actions)
[![Python](https://img.shields.io/pypi/pyversions/midsv.svg?label=Python&color=blue&style=flat-square)](https://pypi.org/project/midsv/)
[![PyPI](https://img.shields.io/pypi/v/midsv.svg?label=PyPI&color=orange&style=flat-square)](https://pypi.org/project/midsv/)
[![Bioconda](https://img.shields.io/conda/v/bioconda/midsv?label=Bioconda&color=orange&style=flat-square)](https://anaconda.org/bioconda/midsv)


# midsv

`midsv` is a Python module to convert SAM to MIDSV format.

MIDSV (Match, Insertion, Deletion, Substitution, and inVersion) is a comma-separated format representing the difference between a reference and a query with the same length as the reference.

> ⚠️ MIDSV is for the target amplicon sequence (10-100 kbp). It may crash when whole chromosomes are used as reference due to running out of memory.

MIDSV can provides `MIDSV` and `QSCORE`.

- `MIDSV` keeps original nucleotides annotating mutations
- `QSCORE` provides Phred quality score on each nucleotide


# Installation

From [Bioconda](https://anaconda.org/bioconda/midsv) (recommended):

```bash
conda install -c bioconda midsv
```

From [PyPI](https://pypi.org/project/midsv/):

```bash
pip install midsv
```

# Usage

```python
midsv.transform(
    sam: list[list[str]] | Iterator[list[str]],
    qscore: bool = True) -> list[dict[str, str | int]],
    keep: str | list[str] = None
```

- sam: Lists or Iterator of SAM format
- qscore (bool, optional): Output QSCORE. Defaults to True.
- keep: Subset of {'FLAG', 'POS', 'SEQ', 'QUAL', 'CIGAR', 'CSTAG'} to keep the field of SAM file. Defaults to None.

- `midsv.transform()` returns a list of dictionaries incuding `QNAME`, `RNAME`, `MIDSV`, and `QSCORE`.
- `MIDSV` and `QSCORE` are comma-separated strings and have the same reference sequence length.

# Specification

## MIDSV

| Op  | Regex          | Description                  |
| --- | -------------- | ---------------------------- |
| =   | [ACGTN]        | Identical sequence           |
| +   | [ACGTN]        | Insertion to the reference   |
| -   | [ACGTN]        | Deletion from the reference  |
| *   | [ACGTN][ACGTN] | Substitution                 |
|     | [acgtn]        | Inversion                    |
| \|  |                | Separater of insertion sites |

`MIDSV` uses `|` to separate nucleotides in insertion sites.

Therefore, `+A|+C|+G|+T|=A` can be easily splited to `[+A, +C, +G, +T, =A]` by `"+A|+C|+G|+T|=A".split("|")` in Python.

## QSCORE


| Op  | Description                  |
| --- | ---------------------------- |
| -1  | Unknown                      |
| \|  | Separator at insertion sites |

`QSCORE` uses `-1` at deletion or unknown nucleotides.

As with `MIDSV`, `QSCORE` uses `|` to separate quality scores in insertion sites.


```python
import midsv

# Perfect match

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['match', '0', 'example', '1', '60', '10M', '*', '0', '0', 'ACGTACGTAC', '0123456789', 'cs:Z:=ACGTACGTAC']
    ]

midsv.transform(sam)
# [{
#   'QNAME': 'control',
#   'RNAME': 'example',
#   'MIDSV': '=A,=C,=G,=T,=A,=C,=G,=T,=A,=C',
#   'QSCORE': '15,16,17,18,19,20,21,22,23,24'
# }]

# Insertion, deletion and substitution

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['indel_sub', '0', 'example', '1', '60', '5M3I1M2D2M', '*', '0', '0', 'ACGTGTTTCGT', '01234!!!56789', 'cs:Z:=ACGT*ag+ttt=C-aa=GT']
    ]

midsv.transform(sam)
# [{
#   'QNAME': 'indel_sub',
#   'RNAME': 'example',
#   'MIDSV': '=A,=C,=G,=T,*AG,+T|+T|+T|=C,-A,-A,=G,=T',
#   'QSCORE': '15,16,17,18,19,0|0|0|20,-1,-1,21,22'
# }]

# Large deletion

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['large-deletion', '0', 'example', '1', '60', '2M', '*', '0', '0', 'AC', '01', 'cs:Z:=AC'],
    ['large-deletion', '0', 'example', '9', '60', '2M', '*', '0', '0', 'AC', '89', 'cs:Z:=AC']
    ]

midsv.transform(sam)
# [
#   {'QNAME': 'large-deletion',
#   'RNAME': 'example',
#   'MIDSV': '=A,=C,N,N,N,N,N,N,=A,=C',
#   'QSCORE': '15,16,-1,-1,-1,-1,-1,-1,23,24'}
# ]

# Inversion

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['inversion', '0', 'example', '1', '60', '5M', '*', '0', '0', 'ACGTA', '01234', 'cs:Z:=ACGTA'],
    ['inversion', '16', 'example', '6', '60', '3M', '*', '0', '0', 'CGT', '567', 'cs:Z:=CGT'],
    ['inversion', '2048', 'example', '9', '60', '2M', '*', '0', '0', 'AC', '89', 'cs:Z:=AC']
    ]

midsv.transform(sam)
# [
#   {'QNAME': 'inversion',
#   'RNAME': 'example',
#   'MIDSV': '=A,=C,=G,=T,=A,=c,=g,=t,=A,=C',
#   'QSCORE': '15,16,17,18,19,20,21,22,23,24'}
# ]

```

# Helper functions

## Read SAM file

```python
midsv.read_sam(path_of_sam: str | Path) -> list[list]
```

`midsv.read_sam` read SAM file into a list of lists.


## Read/Write JSON Line (JSONL)

```python
midsv.write_jsonl(dict: list[dict], path_of_jsonl: str | Path)
```

```python
midsv.read_jsonl(path_of_jsonl: str | Path) -> list[dict]
```

Since `midsv` returns a list of dictionaries, `midsv.write_jsonl` outputs it to a file in JSONL format.

Conversely, `midsv.read_jsonl` reads JSONL as a list of dictionaries.

