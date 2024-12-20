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


###########################################################
# TOKEN
###########################################################

git remote set-url origin https://<PAT>@github.com/akikuno/midsv.git
