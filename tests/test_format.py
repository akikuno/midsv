from pathlib import Path
from src.midsv import format, io
from importlib import reload

reload(format)

###########################################################
# Format headers and alignments
###########################################################


def test_extract_sqheaders():
    sampath = Path("tests", "data", "extract_headers", "query.sam")
    sam = io.read_sam(str(sampath))
    test = format.extract_sqheaders(sam)
    answer = {"chr13": 120421639, "chr6": 149736546}
    assert test == answer


def test_dictionarize_sam():
    sampath = Path("tests", "data", "dictionalize_alignments", "sub_cslong.sam")
    sam = io.read_sam(str(sampath))
    test = format.dictionarize_sam(sam)
    answer = Path("tests", "data", "dictionalize_alignments", "answer.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_dictionarize_sam_inversion():
    sampath = Path("tests", "data", "inversion", "inv_cslong.sam")
    sam = io.read_sam(str(sampath))
    test = format.dictionarize_sam(sam)
    answer = Path("tests", "data", "dictionalize_alignments", "answer_inversion.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_dictionarize_sam_not_primary():
    test = format.dictionarize_sam(
        [["not-primary", "272", "test", "1", "0", "3M", "*", "0", "0", "*", "*", "cs:Z:=AGG"]]
    )
    answer = []
    assert test == answer


###########################################################
# Remove undesired reads
###########################################################


def test_remove_softclips():
    path = Path("tests", "data", "softclip", "softclip_cslong.sam")
    sam = io.read_sam(path)
    samdict = format.dictionarize_sam(sam)
    test = format.remove_softclips(samdict)
    for t in test:
        assert len(t["QUAL"]) == 100


def test_remove_overlapped():
    path = Path("tests", "data", "overlap", "overlapped.sam")
    sam = io.read_sam(path)
    samdict = format.dictionarize_sam(sam)
    test = format.remove_overlapped(samdict)
    for t in test:
        assert not t["QNAME"].startswith("overlap")
