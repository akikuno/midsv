from pathlib import Path

import pytest
from src.midsv import io, validate


@pytest.mark.parametrize(
    "input_value, expected_output",
    [
        (None, set()),
        ("FLAG", {"FLAG"}),
        (["FLAG", "POS"], {"FLAG", "POS"}),
        (["SEQ", "QUAL", "CIGAR"], {"SEQ", "QUAL", "CIGAR"}),
        ("CSTAG", {"CSTAG"}),
    ],
)
def test_keep_argument_valid(input_value, expected_output):
    assert validate.keep_argument(input_value) == expected_output


@pytest.mark.parametrize(
    "input_value",
    [
        "INVALID",
        ["FLAG", "INVALID"],
        ["SEQ", "QUAL", "INVALID"],
    ],
)
def test_keep_argument_invalid(input_value):
    with pytest.raises(
        ValueError, match=r"'keep' must be a subset of \{'FLAG', 'POS', 'SEQ', 'QUAL', 'CIGAR', 'CSTAG'\}"
    ):
        validate.keep_argument(input_value)


###########################################################
# Validate sam format
###########################################################


@pytest.mark.parametrize(
    "sam, qscore, expected_exception",
    [
        # Case 1: Valid SAM input with cs tag
        (
            [["@header"], ["read1", "0", "chr1", "100", "*", "50M", "*", "*", "0", "ACGT", "!!!!", "cs:Z:=ACGT"]],
            False,
            None,
        ),
        # Case 2: Missing cs tag
        ([["@header"], ["read1", "0", "chr1", "100", "*", "50M", "*", "*", "0", "ACGT", "!!!!"]], False, ValueError),
        # Case 3: Missing QUAL information when qscore=True
        (
            [["@header"], ["read1", "0", "chr1", "100", "*", "50M", "*", "*", "0", "ACGT", "*", "cs:Z:=ACGT"]],
            True,
            ValueError,
        ),
        # Case 4: No alignment information
        ([["@header"], ["read1", "0", "*", "*", "*", "*", "*", "*", "*", "*", "*"]], False, ValueError),
        # Case 5: Invalid SAM format (less than 10 columns)
        ([["@header"], ["read1", "0", "chr1", "100"]], False, ValueError),
    ],
)
def test_sam_alignments(sam, qscore, expected_exception):
    if expected_exception:
        with pytest.raises(expected_exception):
            validate.sam_alignments(sam, qscore)
    else:
        validate.sam_alignments(sam, qscore)


def test_validate_alignments_cs_short():
    with pytest.raises(ValueError) as excinfo:
        validate.sam_alignments([["id", "0", "test", "0", "0", "4M", "*", "0", "0", "ACGT", "0000", "cs:Z:*ga:3"]])
    assert str(excinfo.value) == "Input does not have long-formatted cs tag"


def test_validate_alignments_start_substitution():
    with pytest.raises(AssertionError):
        assert validate.sam_alignments(
            [["id", "0", "test", "0", "0", "4M", "*", "0", "0", "ACGT", "0000", "cs:Z:*ga=CGT"]]
        )


def test_validate_headers_no_header():
    with pytest.raises(ValueError) as excinfo:
        validate.sam_headers([["no sq header"]])
        assert str(excinfo.value) == "Input does not have @SQ header"


def test_validate_alignments_no_alignment():
    with pytest.raises(ValueError) as excinfo:
        validate.sam_alignments([["@SQ", "SN:random_100bp", "LN:100"]])
    assert str(excinfo.value) == "No alignment information"


def test_validate_alignments_no_cslong():
    path = Path("tests", "data", "splicing", "splicing_cs.sam")
    sam = io.read_sam(path)
    with pytest.raises(ValueError) as excinfo:
        validate.sam_alignments(sam)
    assert str(excinfo.value) == "Input does not have long-formatted cs tag"
