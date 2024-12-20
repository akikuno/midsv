#!/bin/bash

###########################################################
# test
###########################################################

(
    cd midsv &&
    export PYTHONPATH="$PYTHONPATH:./src" &&
    pytest tests/ -vv ||
    exit 1
)
