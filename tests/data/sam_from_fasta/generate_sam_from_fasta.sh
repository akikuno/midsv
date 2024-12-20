#!/bin/sh

ref=tests/data/random_100bp.fa
que=tests/data/sam_from_fasta/control.fa

# SAM
minimap2 -ax map-ont --cs=long "$ref" "$que" >tests/data/sam_from_fasta/control.sam
