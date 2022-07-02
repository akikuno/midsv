#!/bin/sh

ref=tests/data/random_100bp.fa
que=tests/data/overlap/overlap.fq

minimap2 -ax map-ont --cs=long "$ref" "$que" >tests/data/overlap/overlap_cslong.sam
