#!/bin/sh

ref=tests/data/random_1500bp.fa
que=tests/data/inversion/inversion.fq

minimap2 -ax map-ont --cs "$ref" "$que" >tests/data/inversion/inv_cs.sam
minimap2 -ax map-ont --cs=long "$ref" "$que" >tests/data/inversion/inv_cslong.sam
