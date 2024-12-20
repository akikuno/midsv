#!/bin/bash

###########################################################
# test
###########################################################

(cd midsv && pip install -qe .)
(
    cd midsv &&
    export PYTHONPATH="$PYTHONPATH:./src" &&
    pytest tests/ -vv ||
    exit 1
)
