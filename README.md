[![Licence](https://img.shields.io/badge/License-MIT-9cf.svg?style=flat-square)](https://choosealicense.com/licenses/mit/)
[![Test](https://img.shields.io/github/actions/workflow/status/akikuno/midsv/ci.yml?branch=main&label=Test&color=brightgreen&style=flat-square)](https://github.com/akikuno/midsv/actions)
[![Python](https://img.shields.io/pypi/pyversions/midsv.svg?label=Python&color=blue&style=flat-square)](https://pypi.org/project/midsv/)
[![PyPI](https://img.shields.io/pypi/v/midsv.svg?label=PyPI&color=orange&style=flat-square)](https://pypi.org/project/midsv/)
[![Bioconda](https://img.shields.io/conda/v/bioconda/midsv?label=Bioconda&color=orange&style=flat-square)](https://anaconda.org/bioconda/midsv)


# midsv

`midsv` is a Python module that converts SAM files to MIDSV format.

MIDSV (Match, Insertion, Deletion, Substitution, and inVersion) is a comma-separated format that represents differences between a reference and a query, with the same length as the reference.

> [!CAUTION]
> MIDSV is intended for targeted amplicon sequences (10-100 kbp).  
> Using whole chromosomes as references may exhaust memory and crash.  


>[!IMPORTANT]
> MIDSV requires long-format cstag tags in the SAM file.  
> Please use minimap2 with `--cs=long` option.
> or use [`cstag`](https://github.com/akikuno/cstag) tool to append long-format cstag.

The output includes `MIDSV` and, optionally, `QSCORE`.

- `MIDSV` preserves original nucleotides while annotating mutations.
- `QSCORE` provides Phred quality scores for each nucleotide.


Details of MIDSV (formerly MIDS) are described in [our paper](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001507#sec009).

# 🛠️Installation

From [Bioconda](https://anaconda.org/bioconda/midsv) (recommended):

```bash
conda install -c bioconda midsv
```

From [PyPI](https://pypi.org/project/midsv/):

```bash
pip install midsv
```


# 📜Specifications

## MIDSV

| Op  | Regex          | Description                   |
| --- | -------------- | ----------------------------- |
| =   | [ACGTN]        | Identical sequence            |
| +   | [ACGTN]        | Insertion to the reference    |
| -   | [ACGTN]        | Deletion from the reference   |
| *   | [ACGTN][ACGTN] | Substitution                  |
|     | [acgtn]        | Inversion                     |
| \|  |                | Separator for insertion sites |

`MIDSV` uses `|` to separate nucleotides in insertion sites so `+A|+C|+G|+T|=A` can be easily split into `[+A, +C, +G, +T, =A]` by `"+A|+C|+G|+T|=A".split("|")`.

## QSCORE

| Op  | Description                   |
| --- | ----------------------------- |
| -1  | Unknown                       |
| \|  | Separator for insertion sites |

`QSCORE` uses `-1` for deletions or unknown nucleotides.

As with `MIDSV`, `QSCORE` uses `|` to separate quality scores in insertion sites.



# 📘Usage

```python
midsv.transform(
    path_sam: str | Path,
    qscore: bool = False,
    keep: str | list[str] = None
) -> list[dict[str, str | int]]
```

- path_sam: Path to a SAM file on disk.
- qscore (bool, optional): Output QSCORE. Defaults to False.
- keep: Subset of {'FLAG', 'POS', 'SEQ', 'QUAL', 'CIGAR', 'CSTAG'} to include from the SAM file. Defaults to None.

- `midsv.transform()` returns a list of dictionaries containing `QNAME`, `RNAME`, `MIDSV`, and optionally `QSCORE`, plus any fields specified by `keep`.
- `MIDSV` and `QSCORE` are comma-separated strings and have the same reference sequence length.


# 🖍️Examples

## Perfect match

```python
import midsv
from midsv.io import read_sam

# Perfect match

path_sam = "examples/example_match.sam"
print(list(read_sam(path_sam)))
# sam = [
#     ['@SQ', 'SN:example', 'LN:10'],
#     ['match', '0', 'example', '1', '60', '10M', '*', '0', '0', 'ACGTACGTAC', '0123456789', 'cs:Z:=ACGTACGTAC']
# ]

print(midsv.transform(path_sam, qscore=True))
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

path_sam = "examples/example_indels.sam"
print(list(read_sam(path_sam)))
# [
#     ['@SQ', 'SN:example', 'LN:10'],
#     ['indel_sub', '0', 'example', '1', '60', '5M3I1M2D2M', '*', '0', '0', 'ACGTGTTTCGT', '01234!!!56789', 'cs:Z:=ACGT*ag+ttt=C-aa=GT']
# ]

print(midsv.transform(path_sam, qscore=True))
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

path_sam = "examples/example_large_deletion.sam"
print(list(read_sam(path_sam)))
# [
#     ['@SQ', 'SN:example', 'LN:10'],
#     ['large-deletion', '0', 'example', '1', '60', '2M', '*', '0', '0', 'AC', '01', 'cs:Z:=AC'],
#     ['large-deletion', '0', 'example', '9', '60', '2M', '*', '0', '0', 'AC', '89', 'cs:Z:=AC']
# ]

print(midsv.transform(path_sam, qscore=True))
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

path_sam = "examples/example_inversion.sam"
print(list(read_sam(path_sam)))
# [
#     ['@SQ', 'SN:example', 'LN:10'],
#     ['inversion', '0', 'example', '1', '60', '5M', '*', '0', '0', 'ACGTA', '01234', 'cs:Z:=ACGTA'],
#     ['inversion', '16', 'example', '6', '60', '3M', '*', '0', '0', 'CGT', '567', 'cs:Z:=CGT'],
#     ['inversion', '2048', 'example', '9', '60', '2M', '*', '0', '0', 'AC', '89', 'cs:Z:=AC']
# ]

print(midsv.transform(path_sam, qscore=True))
# [
#   {'QNAME': 'inversion',
#   'RNAME': 'example',
#   'MIDSV': '=A,=C,=G,=T,=A,=c,=g,=t,=A,=C',
#   'QSCORE': '15,16,17,18,19,20,21,22,23,24'}
# ]

```



# 🧩Helper functions

## Read SAM file

```python
midsv.io.read_sam(path_sam: str | Path) -> Iterator[list[str]]
```

`midsv.io.read_sam` reads a local SAM file into an iterator of string lists.


## Read/Write JSON Line (JSONL)

```python
midsv.io.write_jsonl(dicts: list[dict[str, str]], path_output: str | Path)
```

Since `midsv.transform` returns a list of dictionaries, `midsv.io.write_jsonl` outputs it to a file in JSONL format.

```python
midsv.io.read_jsonl(path_input: str | Path) -> Iterator[dict[str, str]]
```

Conversely, `midsv.io.read_jsonl` reads JSONL as an iterator of dictionaries.

## Reverse complement MIDSV

```python
from midsv import formatter

midsv_tag = "=A,=A,-G,+T|+C|=A,=A,*AG,=C"
revcomp_tag = formatter.revcomp(midsv_tag)
print(revcomp_tag)
# =G,*TC,=T,=T,+G|+A|-C,=T,=T
```

`midsv.formatter.revcomp` returns the reverse complement of a MIDSV string. Insertions are reversed and complemented with their anchor moved to the new position, following the MIDSV specification.

## Export VCF

```python
from midsv import transform
from midsv.io import write_vcf

alignments = transform("examples/example_indels.sam", qscore=False)
write_vcf(alignments, "variants.vcf", large_sv_threshold=50)
```

`midsv.io.write_vcf` writes MIDSV output to VCF and supports insertion, deletion, substitution, large insertion, large deletion, and inversion. Insertions longer than `large_sv_threshold` are emitted as symbolic `<INS>`, large deletions (or `=N` padding) use `<DEL>`, and inversions use `<INV>`. The INFO field includes `TYPE` or `SVTYPE`, `SVLEN`, `SEQ`, and `QNAME`.
