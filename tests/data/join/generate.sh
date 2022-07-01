#!/bin/sh

ref=tests/data/random_1500bp.fa
que=tests/data/join/inversion_deletion.fq

# minimap2 -ax map-ont --cs "$ref" "$que" >tests/data/join/inv_del_cs.sam
minimap2 -ax map-ont --cs=long "$ref" "$que" >tests/data/join/inv_del_cslong.sam
