#!/bin/bash

ref=tests/data/random_100bp.fa
que=tests/data/overlap/overlap.fq

minimap2 -ax map-ont --cs=long "$ref" "$que" >tests/data/overlap/overlap_cslong.sam

cat tests/data/real/tyr_cslong.sam |
    grep -e "^@" -e 09624242-03c9-43dd-9bde-8857090c6efd >tests/data/overlap/real_overlap.sam
