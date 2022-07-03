from pathlib import Path
from src.mids import format
from src.mids import convert
from src.mids import proofread

from importlib import reload

reload(proofread)


def test_join_control():
    sam = Path("tests", "data", "join", "test_control.txt").read_text()
    sam = eval(sam)
    test = proofread.join(sam)
    answer = Path("tests", "data", "join", "answer_control.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_join_inversion():
    sam = Path("tests", "data", "join", "test_inv.txt").read_text()
    sam = eval(sam)
    test = proofread.join(sam)
    answer = Path("tests", "data", "join", "answer_inv.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_join_deletion():
    sam = Path("tests", "data", "join", "test_del.txt").read_text()
    sam = eval(sam)
    test = proofread.join(sam)
    answer = Path("tests", "data", "join", "answer_del.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_pad():
    sam = Path("tests", "data", "pad", "padding.sam")
    sam = format.read_sam(str(sam))
    sqheaders = format.extract_sqheaders(sam)
    samdict = format.dictionarize_sam(sam)
    for i, alignment in enumerate(samdict):
        samdict[i]["MIDS"] = convert.cstag_to_mids(alignment["CSTAG"])
        samdict[i]["CSSPLIT"] = convert.cstag_to_mids(alignment["CSTAG"])
        samdict[i]["QSCORE"] = convert.qual_to_qscore(alignment["QUAL"], alignment["MIDS"])
    test = proofread.pad(samdict, sqheaders)
    for t in test:
        mlen = len(t["MIDS"].split(","))
        clen = len(t["CSSPLIT"].split(","))
        qlen = len(t["QSCORE"].split(","))
        assert mlen == clen == qlen
    answer = Path("tests", "data", "pad", "answer.txt").read_text()
    answer = eval(answer)
    assert test == answer

