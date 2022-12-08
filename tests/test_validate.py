import pytest
from pathlib import Path
from src.midsv import validate, io
from importlib import reload

reload(validate)

###########################################################
# Validate sam format
###########################################################


def test_validate_alignments_cs_short():
    with pytest.raises(AttributeError) as excinfo:
        validate.sam_alignments([["id", "0", "test", "0", "0", "4M", "*", "0", "0", "ACGT", "0000", "cs:Z:*ga:3"]])
    assert str(excinfo.value) == "Input does not have long-formatted cs tag"


def test_validate_alignments_start_substitution():
    with pytest.raises(AssertionError):
        assert validate.sam_alignments(
            [["id", "0", "test", "0", "0", "4M", "*", "0", "0", "ACGT", "0000", "cs:Z:*ga=CGT"]]
        )


def test_validate_headers_no_header():
    with pytest.raises(AttributeError) as excinfo:
        validate.sam_headers([["no sq header"]])
        assert str(excinfo.value) == "Input does not have @SQ header"


def test_validate_alignments_no_alignment():
    with pytest.raises(AttributeError) as excinfo:
        validate.sam_alignments([["@SQ", "SN:random_100bp", "LN:100"]])
    assert str(excinfo.value) == "No alignment information"


def test_validate_alignments_no_cslong():
    path = Path("tests", "data", "splicing", "splicing_cs.sam")
    sam = io.read_sam(path)
    with pytest.raises(AttributeError) as excinfo:
        validate.sam_alignments(sam)
    assert str(excinfo.value) == "Input does not have long-formatted cs tag"

