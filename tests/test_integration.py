from pathlib import Path
from src import midsv
import pytest
from importlib import reload

reload(midsv)


def test_valueerror_midsv_cssplit():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = midsv.format.read_sam(str(sampath))
    with pytest.raises(ValueError) as e:
        _ = midsv.transform(sam, midsv=False, cssplit=False, qscore=False)
    assert str(e.value) == "Either midsv or cssplit must be True"


def test_valueerror_qscore():
    sampath = Path("tests", "data", "sam_from_fasta", "control.sam")
    sam = midsv.format.read_sam(str(sampath))
    with pytest.raises(ValueError) as e:
        _ = midsv.transform(sam, midsv=True, cssplit=True, qscore=True)
    assert str(e.value) == "qsucore must be False because the input does not contain QUAL"


def test_integration_midsv():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = midsv.format.read_sam(str(sampath))
    test = midsv.transform(sam, midsv=True, cssplit=False, qscore=False)
    answer = Path("tests", "data", "integrate", "answer_integrate_onlymidsv.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_cssplit():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = midsv.format.read_sam(str(sampath))
    test = midsv.transform(sam, midsv=False, cssplit=True, qscore=False)
    answer = Path("tests", "data", "integrate", "answer_integrate_onlycssplit.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_all_true():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = midsv.format.read_sam(str(sampath))
    test = midsv.transform(sam, midsv=True, cssplit=True, qscore=True)
    answer = Path("tests", "data", "integrate", "answer_integrate_all_true.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_cssplit_and_qual():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = midsv.format.read_sam(str(sampath))
    test = midsv.transform(sam, midsv=False, cssplit=True, qscore=True)
    answer = Path("tests", "data", "integrate", "answer_integrate_cssplit_and_qual.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_real_sam():
    sampath = Path("tests", "data", "real", "tyr_cslong.sam")
    sam = midsv.format.read_sam(str(sampath))

    sqheaders = midsv.format.extract_sqheaders(sam)
    samdict_polished = midsv.transform(sam, midsv=True, cssplit=True, qscore=True)

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
    sam = midsv.format.read_sam(str(sampath))

    midsv.format.check_sam_format(sam)

    sqheaders = midsv.format.extract_sqheaders(sam)
    samdict = midsv.format.dictionarize_sam(sam)

    samdict = midsv.format.remove_softclips(samdict)
    samdict = midsv.format.remove_overlapped(samdict)

    for i, alignment in enumerate(samdict):
        samdict[i]["MIDSV"] = midsv.convert.cstag_to_midsv(alignment["CSTAG"])

    for i, alignment in enumerate(samdict):
        samdict[i]["CSSPLIT"] = midsv.convert.cstag_to_cssplit(alignment["CSTAG"])

    for i, alignment in enumerate(samdict):
        samdict[i]["QSCORE"] = midsv.convert.qual_to_qscore_midsv(alignment["QUAL"], alignment["MIDSV"])

    samdict_polished = midsv.proofread.join(samdict)
    samdict_polished = midsv.proofread.pad(samdict_polished, sqheaders)
    samdict_polished = midsv.proofread.select(samdict_polished)

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
