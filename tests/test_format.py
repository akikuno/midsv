import pytest
import os
from pathlib import Path
from src.midsv import format
from importlib import reload

reload(format)


###########################################################
# Read sam
###########################################################


def test_check_read_sam_str():
    path = os.path.join("tests", "data", "subindel", "subindel_cslong.sam")
    assert format.read_sam(path)


def test_check_read_sam_Path():
    path = Path("tests", "data", "subindel", "subindel_cslong.sam")
    assert format.read_sam(path)


def test_check_read_sam_TypeError():
    with pytest.raises(TypeError):
        assert format.read_sam(1)


def test_check_read_sam_FileNotFoundError():
    with pytest.raises(FileNotFoundError):
        assert format.read_sam("hoge")


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
    sam = format.read_sam(path)
    with pytest.raises(AttributeError) as excinfo:
        format.check_sam_format(sam)
    assert str(excinfo.value) == "Input does not have long-formatted cs tag"


def test_check_alignments_splicing():
    path = Path("tests", "data", "splicing", "splicing_cslong.sam")
    sam = format.read_sam(path)
    with pytest.raises(AttributeError) as excinfo:
        format.check_sam_format(sam)
    assert str(excinfo.value) == "long-read spliced alignment are currently not supported"


###########################################################
# Format headers and alignments
###########################################################


def test_extract_sqheaders():
    sampath = Path("tests", "data", "extract_headers", "query.sam")
    sam = format.read_sam(str(sampath))
    test = format.extract_sqheaders(sam)
    answer = {"chr13": 120421639, "chr6": 149736546}
    assert test == answer


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
    sam = format.read_sam(path)
    samdict = format.dictionarize_sam(sam)
    test = format.remove_softclips(samdict)
    for t in test:
        assert len(t["QUAL"]) == 100


def test_remove_overlapped():
    path = Path("tests", "data", "overlap", "overlapped.sam")
    sam = format.read_sam(path)
    samdict = format.dictionarize_sam(sam)
    test = format.remove_overlapped(samdict)
    for t in test:
        assert not t["QNAME"].startswith("overlap")

