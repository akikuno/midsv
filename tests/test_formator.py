from pathlib import Path

import pytest

from src.midsv import formatter, io

###########################################################
# Format headers and alignments
###########################################################


def test_extract_sqheaders():
    sampath = Path("tests", "data", "extract_headers", "query.sam")
    sam = io.read_sam(str(sampath))
    test = formatter.extract_sqheaders(sam)
    answer = {"chr13": 120421639, "chr6": 149736546}
    assert test == answer


###########################################################
# Remove undesired reads
###########################################################


@pytest.mark.parametrize(
    "cigar, expected",
    [
        pytest.param("10M5I2D", ["10M", "5I", "2D"], id="case_simple_operations"),
        pytest.param("5M1S3M", ["5M", "1S", "3M"], id="case_with_softclip"),
        pytest.param("20M", ["20M"], id="case_single_operation"),
        pytest.param("", [], id="case_empty_string"),
        pytest.param("10M1I1M2S", ["10M", "1I", "1M", "2S"], id="case_complex_operations"),
    ],
)
def test_split_cigar(cigar, expected):
    assert formatter.split_cigar(cigar) == expected


def test_remove_softclips():
    path = Path("tests", "data", "softclip", "softclip_cslong.sam")
    sam = io.read_sam(path)
    samdict = formatter.alignments_to_dict(sam)
    test = formatter.remove_softclips(samdict)
    for t in test:
        assert len(t["QUAL"]) == 100


@pytest.mark.parametrize(
    "alignment, expected",
    [
        pytest.param({"POS": 1, "SEQ": "ACGT", "CIGAR": "2M1I1D1M"}, {"SEQ": "ACNT"}, id="case_simple_cigar"),
        pytest.param({"POS": 6, "SEQ": "ACGT", "CIGAR": "2H2M1I1D1M"}, {"SEQ": "NNNNNACNT"}, id="case_start_with_5nt"),
        pytest.param({"POS": 1, "SEQ": "ACGT", "CIGAR": "2H2M1I1D1M"}, {"SEQ": "ACNT"}, id="case_hardclip"),
        pytest.param({"POS": 1, "SEQ": "ACGT", "CIGAR": "2M5N2M"}, {"SEQ": "ACNNNNNGT"}, id="case_splicing"),
    ],
)
def test_padding_n_to_sequence(alignment, expected):
    result = formatter._padding_n_to_sequence(alignment)
    del result["POS"]
    del result["CIGAR"]
    assert result == expected


def test_remove_resequence():
    path = Path("tests", "data", "overlap", "overlapped.sam")
    sam = io.read_sam(path)
    samdict = formatter.organize_alignments_to_dict(sam)
    test = formatter.remove_resequence(samdict)
    count_overlap = 0
    count_nonoverlap = 0
    for t in test:
        if t["QNAME"].startswith("overlap"):
            count_overlap += 1
        else:
            count_nonoverlap += 1
    assert count_overlap == 1 and count_nonoverlap == 2


###########################################################
# alignments_to_dict
###########################################################


@pytest.mark.parametrize(
    "sam, expected",
    [
        pytest.param(
            [
                ["@HD", "VN:1.6", "SO:coordinate"],
                ["@SQ", "SN:chr1", "LN:248956422"],
                [
                    "read1",
                    "0",
                    "chr1",
                    "100",
                    "255",
                    "10M",
                    "*",
                    "0",
                    "0",
                    "ACTGACTGAA",
                    "FFFFFFFFF",
                    "cs:Z:=ACTGACTGAA",
                ],
            ],
            [
                {
                    "QNAME": "read1",
                    "FLAG": 0,
                    "RNAME": "chr1",
                    "POS": 100,
                    "CIGAR": "10M",
                    "SEQ": "ACTGACTGAA",
                    "QUAL": "FFFFFFFFF",
                    "CSTAG": "cs:Z:=ACTGACTGAA",
                }
            ],
            id="case_single_alignment",
        ),
        pytest.param(
            [
                ["@HD", "VN:1.6", "SO:coordinate"],
                ["@SQ", "SN:chr1", "LN:248956422"],
                [
                    "read1",
                    "0",
                    "chr1",
                    "100",
                    "255",
                    "10M",
                    "*",
                    "0",
                    "0",
                    "ACTGACTGAA",
                    "FFFFFFFFF",
                    "cs:Z:=ACTGACTGAA",
                ],
                ["read2", "16", "chr1", "200", "255", "8M", "*", "0", "0", "GTCAGTCA", "FFFFFFF", "cs:Z:=GTCAGTCA"],
            ],
            [
                {
                    "QNAME": "read1",
                    "FLAG": 0,
                    "RNAME": "chr1",
                    "POS": 100,
                    "CIGAR": "10M",
                    "SEQ": "ACTGACTGAA",
                    "QUAL": "FFFFFFFFF",
                    "CSTAG": "cs:Z:=ACTGACTGAA",
                },
                {
                    "QNAME": "read2",
                    "FLAG": 16,
                    "RNAME": "chr1",
                    "POS": 200,
                    "CIGAR": "8M",
                    "SEQ": "GTCAGTCA",
                    "QUAL": "FFFFFFF",
                    "CSTAG": "cs:Z:=GTCAGTCA",
                },
            ],
            id="case_multiple_alignments",
        ),
        pytest.param(
            [
                ["@HD", "VN:1.6", "SO:coordinate"],
                ["@SQ", "SN:chr1", "LN:248956422"],
                ["read1", "0", "chr1", "100", "255", "10M", "*", "0", "0", "*", "*", "cs:Z:100"],
            ],
            [],
            id="case_unmapped_read",
        ),
        pytest.param(
            [
                ["@HD", "VN:1.6", "SO:coordinate"],
                ["@SQ", "SN:chr1", "LN:248956422"],
            ],
            [],
            id="case_no_alignment",
        ),
    ],
)
def test_alignments_to_dict(sam, expected):
    assert formatter.alignments_to_dict(sam) == expected
