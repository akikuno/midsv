from pathlib import Path
from src.midsv import transform
from src.midsv import format
from src.midsv import convert
from src.midsv import proofread


def test_integration_subindelinv():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = format.read_sam(str(sampath))
    test = transform(sam, midsv=True, cssplit=True, qscore=True)
    answer = Path("tests", "data", "integrate", "answer_integrate.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_onlymidsv():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = format.read_sam(str(sampath))
    test = transform(sam, midsv=True, cssplit=False, qscore=False)
    answer = Path("tests", "data", "integrate", "answer_integrate_onlymidsv.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_onlycssplit():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = format.read_sam(str(sampath))
    test = transform(sam, midsv=False, cssplit=True, qscore=False)
    answer = Path("tests", "data", "integrate", "answer_integrate_onlycssplit.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_real_sam():
    sampath = Path("tests", "data", "real", "tyr_cslong.sam")
    sam = format.read_sam(str(sampath))

    sqheaders = format.extract_sqheaders(sam)
    samdict_polished = transform(sam, midsv=True, cssplit=True, qscore=True)

    for alignment in samdict_polished:
        RNAME = alignment["RNAME"]
        MIDSV = alignment["MIDSV"]
        QSCORE = alignment["QSCORE"]
        RLEN = sqheaders[RNAME]
        mlen = len(MIDSV.split(","))
        qlen = len(QSCORE.split(","))
        assert mlen == qlen == RLEN


def test_integration_eachcomponent():
    sampath = Path("tests", "data", "real", "tyr_cslong.sam")
    sam = format.read_sam(str(sampath))

    format.check_sam_format(sam)

    sqheaders = format.extract_sqheaders(sam)
    samdict = format.dictionarize_sam(sam)

    samdict = format.remove_softclips(samdict)
    samdict = format.remove_overlapped(samdict)

    for i, alignment in enumerate(samdict):
        samdict[i]["MIDSV"] = convert.cstag_to_midsv(alignment["CSTAG"])

    for i, alignment in enumerate(samdict):
        samdict[i]["CSSPLIT"] = convert.cstag_to_cssplit(alignment["CSTAG"])

    for i, alignment in enumerate(samdict):
        samdict[i]["QSCORE"] = convert.qual_to_qscore(alignment["QUAL"], alignment["MIDSV"])

    samdict_polished = proofread.join(samdict)
    samdict_polished = proofread.pad(samdict_polished, sqheaders)
    samdict_polished = proofread.select(samdict_polished)

    for alignment in samdict_polished:
        RNAME = alignment["RNAME"]
        MIDSV = alignment["MIDSV"]
        CSSPLIT = alignment["CSSPLIT"]
        QSCORE = alignment["QSCORE"]
        RLEN = sqheaders[RNAME]
        mlen = len(MIDSV.split(","))
        clen = len(CSSPLIT.split(","))
        qlen = len(QSCORE.split(","))
        assert RLEN == mlen == clen == qlen
