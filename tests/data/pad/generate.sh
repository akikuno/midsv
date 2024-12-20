#!/bin/sh

ref=tests/data/random_100bp.fa
que=tests/data/pad/padding.fq

minimap2 -ax map-ont --cs=long "$ref" "$que" >tests/data/pad/padding.sam
