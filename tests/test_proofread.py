from pathlib import Path
from src.midsv import io
from src.midsv import format
from src.midsv import convert
from src.midsv import proofread

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


def test_join_deletion_microhomology():
    samdict = Path("tests", "data", "join", "test_del_microhomology.txt").read_text()
    samdict = eval(samdict)
    test = proofread.join(samdict)
    answer = Path("tests", "data", "join", "answer_del_microhomology.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_join_real_microhomology():
    samdict = Path("tests", "data", "join", "test_real_microhomology.txt").read_text()
    samdict = eval(samdict)
    test = proofread.join(samdict)
    test = test[0]["CSSPLIT"]
    answer = Path("tests", "data", "join", "answer_real_microhomology.txt").read_text()
    answer = eval(answer)
    assert test == answer

def test_select_keep_FLAG():
    samdict = Path("tests", "data", "join", "test_real_microhomology.txt").read_text()
    samdict = eval(samdict)
    test = proofread.select(samdict, keep={"FLAG", "SEQ"})
    test = set(list(test[0].keys()))
    answer = {'QNAME', 'RNAME', 'FLAG', 'CSSPLIT', 'QSCORE'}
    assert test == answer


def test_pad():
    sam = Path("tests", "data", "pad", "padding.sam")
    sam = io.read_sam(str(sam))
    sam = list(sam)
    sqheaders = format.extract_sqheaders(sam)
    samdict = format.dictionarize_sam(sam)
    for i, alignment in enumerate(samdict):
        samdict[i]["MIDSV"] = convert.cstag_to_midsv(alignment["CSTAG"])
        samdict[i]["CSSPLIT"] = convert.cstag_to_cssplit(alignment["CSTAG"])
        samdict[i]["QSCORE"] = convert.qual_to_qscore_midsv(alignment["QUAL"], alignment["MIDSV"])
    test = proofread.pad(samdict, sqheaders)
    for t in test:
        mlen = len(t["MIDSV"].split(","))
        clen = len(t["CSSPLIT"].split(","))
        qlen = len(t["QSCORE"].split(","))
        assert mlen == clen == qlen
    answer = Path("tests", "data", "pad", "answer_pad.txt").read_text()
    answer = eval(answer)
    assert test == answer

