# mids

`mids` is a Python module to manipulate MIDS format.

MIDS is a format that represents Match/Insertion/Deletion/Substitution at each nucleotide.  
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

## Usage

```python
import mids

# Convert a cs tag from long to short
cs = "cs:Z:=ACGT*ag=CGT"

cstag.shorten(cs)
# => cs:Z::4*ag:3


# Convert a cs tag from short to long
cs = "cs:Z::4*ag:3"
cigar = "8M"
seq = "ACGTACGT"

cstag.lengthen(cs, cigar, seq)
# => cs:Z:=ACGT*ag=CGT
```

### Call consensus

```python
import cstag

cs_list = ["cs:Z:=ACGT", "cs:Z:=AC*gt=T", "cs:Z:=C*gt=T", "cs:Z:=C*gt=T", "cs:Z:=ACT+ccc=T"]
cigar_list = ["4M", "4M", "1S3M", "3M", "3M3I1M"]
pos_list = [1, 1, 1, 2, 1]

cstag.consensus(cs_list, cigar_list, pos_list)
# => cs:Z:=AC*gt*T
```

### Mask low-quality bases in a cs tag

```python
import cstag

cs = "cs:Z:=ACGT*ac+gg-cc=T"
cigar = "5M2I2D1M"
qual = "AA!!!!AA"
phred_threshold = 10
cstag.mask(cs, cigar, qual, phred_threshold)
# => cs:Z:=ACNN*an+ng-cc=T
```

### Output HTML report

```python
import cstag

cs = "cs:Z:=AC+GGG=T-ACGT*at~gt10cg=GNNN"
output = "report"
description = "Example"

cstag.to_html(cs, output, description)
# => Output "report.html"
```
The `report.html` is :point_down:

<img width="414" alt="example_report" src="https://user-images.githubusercontent.com/15861316/158910398-67f480d2-8742-412a-b528-40e545c46513.png">
