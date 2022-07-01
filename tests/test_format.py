import pytest
from pathlib import Path
from src.mids import format
from importlib import reload

reload(format)

###########################################################
# Check sam format
###########################################################


def test_check_headers_no_header():
    with pytest.raises(AttributeError) as excinfo:
        format.check_sam_format([["no sq header"]])
        assert str(excinfo.value) == "Input does not have @SQ header"


def test_check_alignments_no_alignment():
    with pytest.raises(AttributeError) as excinfo:
        format.check_sam_format([["@SQ", "SN:random_100bp", "LN:100"]])
    assert str(excinfo.value) == "No alignment information"


def test_check_alignments_no_cslong():
    path = Path("tests", "data", "splicing", "splicing_cs.sam")
    sam = format.read_sam(str(path))
    with pytest.raises(AttributeError) as excinfo:
        format.check_sam_format(sam)
    assert str(excinfo.value) == "Input does not have long-formatted cs tag"


def test_check_alignments_splicing():
    path = Path("tests", "data", "splicing", "splicing_cslong.sam")
    sam = format.read_sam(str(path))
    with pytest.raises(AttributeError) as excinfo:
        format.check_sam_format(sam)
    assert str(excinfo.value) == "Spliced long reads are currently not supported"


###########################################################
# Format headers and alignments
###########################################################


def test_extract_sqheaders():
    sampath = Path("tests", "data", "extract_headers", "query.sam")
    sam = format.read_sam(str(sampath))
    test = format.extract_sqheaders(sam)
    answer = {"chr13": 120421639, "chr6": 149736546}
    assert test == answer


# def test_append_reflen():
#     sampath = Path("tests", "data", "append_reflen", "query.sam")
#     sam = format.read_sam(str(sampath))
#     test = format.append_reflen(sam)
#     answer = Path("tests", "data", "append_reflen", "answer.txt").read_text()
#     answer = eval(answer)
#     assert test == answer


def test_dictionarize_sam():
    sampath = Path("tests", "data", "dictionalize_alignments", "sub_cslong.sam")
    sam = format.read_sam(str(sampath))
    test = format.dictionarize_sam(sam)
    answer = Path("tests", "data", "dictionalize_alignments", "answer.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_dictionarize_sam_inversion():
    sampath = Path("tests", "data", "inversion", "inv_cslong.sam")
    sam = format.read_sam(str(sampath))
    test = format.dictionarize_sam(sam)
    answer = Path("tests", "data", "dictionalize_alignments", "answer_inversion.txt").read_text()
    answer = eval(answer)
    assert test == answer


###########################################################
# Remove undesired reads
###########################################################


def test_remove_softclips():
    path = Path("tests", "data", "softclip", "softclip_cslong.sam")
    sam = format.read_sam(str(path))
    test = format.dictionarize_sam(sam)
    test = format.remove_softclips(test)
    for t in test:
        assert len(t["QUAL"]) == 100
