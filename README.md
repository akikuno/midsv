# MIDSV

`MIDSV` is a Python module to convert SAM to MIDSV format.

MIDSV (Match, Insertion, Deletion, Substitution, and inVersion) is a comma-separated format representing the difference between a reference and a query with the same length as the reference.

MIDSV provides `MIDSV`, `CSSPLIT`, and `QSCORE`.

- `MIDSV` is a simple representation focusing on mutations
- `CSSPLIT` keeps original nucleotides
- `QSCORE` provides Phred quality score on each nucleotide

MIDSV details are described in [our paper](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001507#sec009).  

## Installation

From [PyPI](https://pypi.org/project/MIDSV/):

```bash
pip install midsv
```

From [Bioconda](https://anaconda.org/bioconda/MIDSV):

```bash
conda install -c bioconda midsv
```

## Usage

`MIDSV.transcorm(sam: list[list])` returns a dictionary incuding `QNAME`, `RNAME`, `MIDSV`, `CSSPLIT`, and `QSCORE`.

Notice `MIDSV`, `CSSPLIT`, and `QSCORE` are comma-separated and have the same reference sequence length.

(*LN represents the length of a reference sequence).

```python
import MIDSV

# Perfect match

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['match', '0', 'example', '1', '60', '10M', '*', '0', '0', 'ACGTACGTAC', '0123456789', 'cs:Z:=ACGTACGTAC']
    ]

MIDSV.transform(sam)
# [{
#   'QNAME': 'control',
#   'RNAME': 'example',
#   'MIDSV': 'M,M,M,M,M,M,M,M,M,M',
#   'CSSPLIT': '=A,=C,=G,=T,=A,=C,=G,=T,=A,=C',
#   'QSCORE': '15,16,17,18,19,20,21,22,23,24'
# }]

# Insertion, deletion and substitution

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['indel_sub', '0', 'example', '1', '60', '5M3I1M2D2M', '*', '0', '0', 'ACGTGTTTCGT', '01234!!!56789', 'cs:Z:=ACGT*ag+ttt=C-aa=GT']
    ]

MIDSV.transform(sam)
# [{
#   'QNAME': 'indel_sub',
#   'RNAME': 'example',
#   'MIDSV': 'M,M,M,M,S,3M,D,D,M,M',
#   'CSSPLIT': '=A,=C,=G,=T,*AG,+T|+T|+T|=C,-A,-A,=G,=T',
#   'QSCORE': '15,16,17,18,19,0|0|0|20,-1,-1,21,22'
# }]

# Large deletion

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['large-deletion', '0', 'example', '1', '60', '2M', '*', '0', '0', 'AC', '01', 'cs:Z:=AC'],
    ['large-deletion', '0', 'example', '9', '60', '2M', '*', '0', '0', 'AC', '89', 'cs:Z:=AC']
    ]

MIDSV.transform(sam)
# [
#   {'QNAME': 'large-deletion',
#   'RNAME': 'example',
#   'MIDSV': 'M,M,D,D,D,D,D,D,M,M',
#   'CSSPLIT': '=A,=C,N,N,N,N,N,N,=A,=C',
#   'QSCORE': '15,16,-1,-1,-1,-1,-1,-1,23,24'}
# ]

# Inversion

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['inversion', '0', 'example', '1', '60', '5M', '*', '0', '0', 'ACGTA', '01234', 'cs:Z:=ACGTA'],
    ['inversion', '16', 'example', '6', '60', '3M', '*', '0', '0', 'CGT', '567', 'cs:Z:=CGT'],
    ['inversion', '2048', 'example', '9', '60', '2M', '*', '0', '0', 'AC', '89', 'cs:Z:=AC']
    ]

MIDSV.transform(sam)
# [
#   {'QNAME': 'inversion',
#   'RNAME': 'example',
#   'MIDSV': 'M,M,M,M,M,m,m,m,M,M',
#   'CSSPLIT': '=A,=C,=G,=T,=A,=c,=g,=t,=A,=C',
#   'QSCORE': '15,16,17,18,19,20,21,22,23,24'}
# ]

```

## Operators

### MIDSV

| Op          | Description                 |
| ----------- | --------------------------- |
| M           | Identical sequence          |
| [1-9][0-9]+ | Insertion to the reference  |
| D           | Deletion from the reference |
| S           | Substitution                |
| [mdsn]      | Inversion                   |
| N           | Unknown                     |

`MIDSV` represents insertion as an integer and append the following operators.

If five insertions follow three matches, MIDSV returns `5M,M,M` (not `5,M,M,M`) since `5M,M,M` keeps reference sequence length in a comma-separated field.

### CSSPLIT

| Op  | Regex          | Description                  |
| --- | -------------- | ---------------------------- |
| =   | [ACGTN]        | Identical sequence           |
| +   | [ACGTN]        | Insertion to the reference   |
| -   | [ACGTN]        | Deletion from the reference  |
| *   | [ACGTN][ACGTN] | Substitution                 |
|     | [acgtn]        | Inversion                    |
| N   |                | Unknown                      |
| \|  |                | Separater at insertion sites |

`CSSPLIT` uses `|` to separate nucleotides in insertion sites.

Therefore, `+A|+C|+G|+T|=A` can be easily splited to `[+A, +C, +G, +T, =A]` by `"+A|+C|+G|+T|=A".split("|")` in Python.

### QSCORE


| Op  | Description                  |
| --- | ---------------------------- |
| -1  | Unknown                      |
| \|  | Separater at insertion sites |

`QSCORE` uses `-1` at deletion sites.

As with `CSSPLIT`, `QSCORE` uses `|` to separate quality scores in insertion sites.