# mids

`mids` is a Python module to manipulate MIDS format.

MIDS is a format that represents Match/Insertion/Deletion/Substitution at each nucleotide of sequencing reads.  
The details are described in [our paper](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001507#sec002:~:text=Cas%2Dcutting%20site.-,Preprocessing,-We%20performed%20preprocessing).  

## Installation

From [PyPI](https://pypi.org/project/mids/):

```bash
pip install mids
```

From [Bioconda](https://anaconda.org/bioconda/mids):

```bash
conda install -c bioconda mids
```

## Examples

```python
import mids

# match ---------------------------------------------------

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['match', '0', 'example', '1', '60', '10M', '*', '0', '0', 'ACGTACGTAC', '0123456789', 'cs:Z:=ACGTACGTAC']
    ]

mids.transform(sam)
# [{'QNAME': 'match', 'RNAME': 'example', 'MIDS': 'M,M,M,M,M,M,M,M,M,M', 'QSCORE': '15,16,17,18,19,20,21,22,23,24'}]

# insertion -----------------------------------------------

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['insertion', '0', 'example', '1', '60', '5M3I5M', '*', '0', '0', 'ACGTATTTCGTAC', '01234!!!56789', 'cs:Z:=ACGTA+ttt=CGTAC']
    ]

mids.transform(sam)
# [{'QNAME': 'insertion', 'RNAME': 'example', 'MIDS': 'M,M,M,M,M,3M,M,M,M,M', 'CSSPLIT': 'A,C,G,T,A,tttC,G,T,A,C', 'QSCORE': '15,16,17,18,19,0|0|0|20,21,22,23,24'}]

# deleton -------------------------------------------------

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['deletion', '0', 'example', '1', '60', '5M3D2M', '*', '0', '0', 'ACGTAAC', '0123489', 'cs:Z:=ACGTA-cgt=AC']
    ]

mids.transform(sam)
# [{'QNAME': 'deletion', 'RNAME': 'example', 'MIDS': 'M,M,M,M,M,D,D,D,M,M', 'QSCORE': '15,16,17,18,19,-1,-1,-1,23,24'}]

# substitution --------------------------------------------

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['substitution', '0', 'example', '1', '60', '10M', '*', '0', '0', 'ACGTGCGTAC', '01234!6789', 'cs:Z:=ACGT*ag=CGTAC']
    ]

mids.transform(sam)
# [{'QNAME': 'substitution', 'RNAME': 'example', 'MIDS': 'M,M,M,M,S,M,M,M,M,M', 'QSCORE': '15,16,17,18,19,0,21,22,23,24'}]

# large deletion ------------------------------------------

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['large-deletion', '0', 'example', '1', '60', '2M', '*', '0', '0', 'AC', '01', 'cs:Z:=AC'],
    ['large-deletion', '0', 'example', '9', '60', '2M', '*', '0', '0', 'AC', '89', 'cs:Z:=AC']
    ]

mids.transform(sam)
# [{'QNAME': 'large-deletion', 'RNAME': 'example', 'MIDS': 'M,M,D,D,D,D,D,D,M,M', 'QSCORE': '15,16,-1,-1,-1,-1,-1,-1,23,24'}]

# inversion -----------------------------------------------

sam = [
    ['@SQ', 'SN:example', 'LN:10'],
    ['inversion', '0', 'example', '1', '60', '5M', '*', '0', '0', 'ACGTA', '01234', 'cs:Z:=ACGTA'],
    ['inversion', '16', 'example', '6', '60', '3M', '*', '0', '0', 'CGT', '567', 'cs:Z:=CGT'],
    ['inversion', '2048', 'example', '9', '60', '2M', '*', '0', '0', 'AC', '89', 'cs:Z:=AC']
    ]

mids.transform(sam)
# [{'QNAME': 'inversion', 'RNAME': 'example', 'MIDS': 'M,M,M,M,M,m,m,m,M,M', 'QSCORE': '15,16,17,18,19,20,21,22,23,24'}]

```
