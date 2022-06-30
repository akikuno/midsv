from pathlib import Path
from src.mids import preprocess, format
from importlib import reload

reload(preprocess)
reload(format)


def test_extract_sqheaders():
    sampath = Path("tests", "data", "extract_headers", "query.sam")
    sam = preprocess.read_sam(str(sampath))
    test = format.extract_sqheaders(sam)
    answer = {"chr13": 120421639, "chr6": 149736546}
    assert test == answer


# def test_append_reflen():
#     sampath = Path("tests", "data", "append_reflen", "query.sam")
#     sam = preprocess.read_sam(str(sampath))
#     test = format.append_reflen(sam)
#     answer = Path("tests", "data", "append_reflen", "answer.txt").read_text()
#     answer = eval(answer)
#     assert test == answer


def test_dictionarize_sam():
    sampath = Path("tests", "data", "dictionalize_alignments", "sub_cslong.sam")
    sam = preprocess.read_sam(str(sampath))
    test = format.dictionarize_sam(sam)
    answer = Path("tests", "data", "dictionalize_alignments", "answer.txt").read_text()
    answer = eval(answer)
    assert test == answer

