import pytest
from pathlib import Path
from src.mids import preprocess


def test_check_headers_no_header():
    with pytest.raises(AttributeError) as excinfo:
        preprocess.check_sam_format([["no sq header"]])
        assert str(excinfo.value) == "Input does not have @SQ header"


def test_check_alignments_no_alignment():
    with pytest.raises(AttributeError) as excinfo:
        preprocess.check_sam_format([["@SQ", "SN:random_100bp", "LN:100"]])
    assert str(excinfo.value) == "No alignment information"


def test_check_alignments_no_cslong():
    path = Path("tests", "data", "splicing", "splicing_cs.sam")
    sam = preprocess.read_sam(str(path))
    with pytest.raises(AttributeError) as excinfo:
        preprocess.check_sam_format(sam)
    assert str(excinfo.value) == "Input does not have long-formatted cs tag"


def test_check_alignments_splicing():
    path = Path("tests", "data", "splicing", "splicing_cslong.sam")
    sam = preprocess.read_sam(str(path))
    with pytest.raises(AttributeError) as excinfo:
        preprocess.check_sam_format(sam)
    assert str(excinfo.value) == "Spliced long reads are currently not supported"


from importlib import reload

reload(preprocess)


def test_remove_softclips():
    path = Path("tests", "data", "softclip", "softclip_cslong.sam")
    sam = preprocess.read_sam(str(path))
    test = preprocess.remove_softclips(sam)
    answer = Path("tests", "data", "softclip", "answer_removed.txt").read_text()
    answer = eval(answer)
    assert test == answer
