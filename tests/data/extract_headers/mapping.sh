#!/bin/bash

wget -O - http://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/chr6.fa.gz |
    gzip -dc >tmp_chr6.fa

wget -O - http://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/chr13.fa.gz |
    gzip -dc >tmp_chr13.fa

cat tmp_*.fa |
    minimap2 -cs=long -ax map-ont - tests/data/multiple_snln/query.fq >tests/data/multiple_snln/query.sam
