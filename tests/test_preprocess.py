import pytest
from pathlib import Path
from src.mids import preprocess


def test_check_headers():
    with pytest.raises(AttributeError):
        preprocess.check_sam_format([["no sq header"]])


def test_check_alignments():
    with pytest.raises(AttributeError):
        preprocess.check_sam_format([["@SQ\tSN:random_100bp\tLN:100"]])

