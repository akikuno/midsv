[![Licence](https://img.shields.io/badge/License-MIT-9cf.svg?style=flat-square)](https://choosealicense.com/licenses/mit/)
[![Test](https://img.shields.io/github/actions/workflow/status/akikuno/midsv/ci.yml?branch=main&label=Test&color=brightgreen&style=flat-square)](https://github.com/akikuno/midsv/actions)
[![Python](https://img.shields.io/pypi/pyversions/midsv.svg?label=Python&color=blue&style=flat-square)](https://pypi.org/project/midsv/)
[![PyPI](https://img.shields.io/pypi/v/midsv.svg?label=PyPI&color=orange&style=flat-square)](https://pypi.org/project/midsv/)
[![Bioconda](https://img.shields.io/conda/v/bioconda/midsv?label=Bioconda&color=orange&style=flat-square)](https://anaconda.org/bioconda/midsv)


# midsv

`midsv` is a Python module to convert SAM to MIDSV format.

MIDSV (Match, Insertion, Deletion, Substitution, and inVersion) is a comma-separated format representing the differences between a reference and a query with the same length as the reference.

> [!CAUTION]
> MIDSV is for target amplicon sequences (10-100 kbp). It may crash when whole chromosomes are used as references due to running out of memory.

MIDSV can provide `MIDSV` and `QSCORE`.

- `MIDSV` keeps original nucleotides while annotating mutations.
- `QSCORE` provides Phred quality scores for each nucleotide.


MIDSV (formerly named MIDS) details are described in [our paper](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001507#sec009).  

# 🛠️Installation

From [Bioconda](https://anaconda.org/bioconda/midsv) (recommended):

```bash
conda install -c bioconda midsv
```

From [PyPI](https://pypi.org/project/midsv/):

```bash
pip install midsv
```

# 📘Usage

```python
midsv.transform(
    path_sam: str | Path,
    qscore: bool = True,
    keep: str | list[str] = None
) -> list[dict[str, str | int]]
```

- path_sam: Path of SAM sfile
- qscore (bool, optional): Output QSCORE. Defaults to True.
- keep: Subset of {'FLAG', 'POS', 'SEQ', 'QUAL', 'CIGAR', 'CSTAG'} to keep fields from the SAM file. Defaults to None.

- `midsv.transform()` returns a list of dictionaries including `QNAME`, `RNAME`, `MIDSV`, and `QSCORE`.
- `MIDSV` and `QSCORE` are comma-separated strings and have the same reference sequence length.

# 📜Specification

## MIDSV

| Op  | Regex          | Description                  |
| --- | -------------- | ---------------------------- |
| =   | [ACGTN]        | Identical sequence           |
| +   | [ACGTN]        | Insertion to the reference   |
| -   | [ACGTN]        | Deletion from the reference  |
| *   | [ACGTN][ACGTN] | Substitution                 |
|     | [acgtn]        | Inversion                    |
| \|  |                | Separator of insertion sites |

`MIDSV` uses `|` to separate nucleotides in insertion sites so `+A|+C|+G|+T|=A` can be easily split into `[+A, +C, +G, +T, =A]` by `"+A|+C|+G|+T|=A".split("|")`.

## QSCORE

| Op  | Description                  |
| --- | ---------------------------- |
| -1  | Unknown                      |
| \|  | Separator at insertion sites |

`QSCORE` uses `-1` for deletions or unknown nucleotides.

As with `MIDSV`, `QSCORE` uses `|` to separate quality scores in insertion sites.


# 🧩Helper functions

## Read SAM file

```python
midsv.io.read_sam(path_sam: str | Path) -> Iterator[list[str]]
```

`midsv.io.read_sam` reads a SAM file into an iterator of lists.


## Read/Write JSON Line (JSONL)

```python
midsv.io.write_jsonl(dicts: list[dict[str, str]], path_jsonl: str | Path)
```

Since `midsv.io.transform` returns a list of dictionaries, `midsv.io.write_jsonl` outputs it to a file in JSONL format.

```python
midsv.io.read_jsonl(path_jsonl: str | Path) -> Iterator[dict[str, str]]
```

Conversely, `midsv.io.read_jsonl` reads JSONL as a list of dictionaries.


# 🖍️Examples

## Perfect match

```python
import midsv
from midsv.io import read_sam

# Perfect match

path_sam = "https://raw.githubusercontent.com/akikuno/midsv/refs/heads/main/examples/example_match.sam"
print(list(read_sam(path_sam)))
# sam = [
#     ['@SQ', 'SN:example', 'LN:10'],
#     ['match', '0', 'example', '1', '60', '10M', '*', '0', '0', 'ACGTACGTAC', '0123456789', 'cs:Z:=ACGTACGTAC']
# ]

print(midsv.transform(read_sam(path_sam), qscore=True))
# [{
#   'QNAME': 'control',
#   'RNAME': 'example',
#   'MIDSV': '=A,=C,=G,=T,=A,=C,=G,=T,=A,=C',
#   'QSCORE': '15,16,17,18,19,20,21,22,23,24'
# }]
```

## Insertion, deletion, and substitution

```python
import midsv
from midsv.io import read_sam

path_sam = "https://raw.githubusercontent.com/akikuno/midsv/refs/heads/main/examples/example_indels.sam"
print(list(read_sam(path_sam)))
# [
#     ['@SQ', 'SN:example', 'LN:10'],
#     ['indel_sub', '0', 'example', '1', '60', '5M3I1M2D2M', '*', '0', '0', 'ACGTGTTTCGT', '01234!!!56789', 'cs:Z:=ACGT*ag+ttt=C-aa=GT']
# ]

print(midsv.transform(sam, qscore=True))
# [{
#   'QNAME': 'indel_sub',
#   'RNAME': 'example',
#   'MIDSV': '=A,=C,=G,=T,*AG,+T|+T|+T|=C,-A,-A,=G,=T',
#   'QSCORE': '15,16,17,18,19,0|0|0|20,-1,-1,21,22'
# }]
```

## Large deletion

```python
import midsv
from midsv.io import read_sam

path_sam = "https://raw.githubusercontent.com/akikuno/midsv/refs/heads/main/examples/example_large_deletion.sam"
print(list(read_sam(path_sam)))
# [
#     ['@SQ', 'SN:example', 'LN:10'],
#     ['large-deletion', '0', 'example', '1', '60', '2M', '*', '0', '0', 'AC', '01', 'cs:Z:=AC'],
#     ['large-deletion', '0', 'example', '9', '60', '2M', '*', '0', '0', 'AC', '89', 'cs:Z:=AC']
# ]

print(midsv.transform(sam, qscore=True))
# [
#   {'QNAME': 'large-deletion',
#   'RNAME': 'example',
#   'MIDSV': '=A,=C,=N,=N,=N,=N,=N,=N,=A,=C',
#   'QSCORE': '15,16,-1,-1,-1,-1,-1,-1,23,24'}
# ]
```

## Inversion

```python
import midsv
from midsv.io import read_sam

path_sam = "https://raw.githubusercontent.com/akikuno/midsv/refs/heads/main/examples/example_inversion.sam"
print(list(read_sam(path_sam)))
# [
#     ['@SQ', 'SN:example', 'LN:10'],
#     ['inversion', '0', 'example', '1', '60', '5M', '*', '0', '0', 'ACGTA', '01234', 'cs:Z:=ACGTA'],
#     ['inversion', '16', 'example', '6', '60', '3M', '*', '0', '0', 'CGT', '567', 'cs:Z:=CGT'],
#     ['inversion', '2048', 'example', '9', '60', '2M', '*', '0', '0', 'AC', '89', 'cs:Z:=AC']
# ]

print(midsv.transform(sam, qscore=True))
# [
#   {'QNAME': 'inversion',
#   'RNAME': 'example',
#   'MIDSV': '=A,=C,=G,=T,=A,=c,=g,=t,=A,=C',
#   'QSCORE': '15,16,17,18,19,20,21,22,23,24'}
# ]

```

