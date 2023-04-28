from pathlib import Path
from src import midsv
import pytest
from importlib import reload

reload(midsv)


def test_valueerror_midsv_cssplit():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = midsv.io.read_sam(str(sampath))
    with pytest.raises(ValueError) as e:
        _ = midsv.transform(sam, midsv=False, cssplit=False, qscore=False)
    assert str(e.value) == "Either midsv or cssplit must be True"


def test_valueerror_qscore():
    sampath = Path("tests", "data", "sam_from_fasta", "control.sam")
    sam = midsv.io.read_sam(str(sampath))
    with pytest.raises(ValueError) as e:
        _ = midsv.transform(sam, midsv=True, cssplit=True, qscore=True)
    assert str(e.value) == "qscore must be False because the input does not contain QUAL"


def test_integration_midsv():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = midsv.io.read_sam(str(sampath))
    test = midsv.transform(sam, midsv=True, cssplit=False, qscore=False)
    answer = Path("tests", "data", "integrate", "answer_integrate_onlymidsv.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_cssplit():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = midsv.io.read_sam(str(sampath))
    test = midsv.transform(sam, midsv=False, cssplit=True, qscore=False)
    answer = Path("tests", "data", "integrate", "answer_integrate_onlycssplit.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_all_true():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = midsv.io.read_sam(str(sampath))
    test = midsv.transform(sam, midsv=True, cssplit=True, qscore=True)
    answer = Path("tests", "data", "integrate", "answer_integrate_all_true.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_cssplit_and_qual():
    sampath = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    sam = midsv.io.read_sam(str(sampath))
    test = midsv.transform(sam, midsv=False, cssplit=True, qscore=True)
    answer = Path("tests", "data", "integrate", "answer_integrate_cssplit_and_qual.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_real_sam():
    sampath = Path("tests", "data", "real", "tyr_cslong.sam")
    sam = midsv.io.read_sam(str(sampath))
    sam = list(sam)
    sqheaders = midsv.format.extract_sqheaders(sam)
    samdict_polished = midsv.transform(sam, midsv=True, cssplit=True, qscore=True)
    for alignment in samdict_polished:
        RNAME = alignment["RNAME"]
        MIDSV = alignment["MIDSV"]
        QSCORE = alignment["QSCORE"]
        mlen = len(MIDSV.split(","))
        qlen = len(QSCORE.split(","))
        answer = sqheaders[RNAME]
        assert mlen == qlen == answer


def test_integration_real_splicing():
    sam = midsv.io.read_sam("tests/data/splicing/real_splicing.sam")
    sam = list(sam)
    sqheaders = midsv.format.extract_sqheaders(sam)
    test = midsv.transform(sam)
    mlen = len(test[0]["MIDSV"].split(","))
    clen = len(test[0]["CSSPLIT"].split(","))
    qlen = len(test[0]["QSCORE"].split(","))
    answer = sqheaders["deletion"]
    assert mlen == clen == qlen == answer


def test_integration_eachcomponent():
    sampath = Path("tests", "data", "real", "tyr_cslong.sam")
    sam = midsv.io.read_sam(str(sampath))
    sam = list(sam)
    midsv.validate.sam_headers(sam)
    midsv.validate.sam_alignments(sam)
    sqheaders = midsv.format.extract_sqheaders(sam)
    samdict = midsv.format.dictionarize_sam(sam)
    samdict = midsv.format.remove_softclips(samdict)
    samdict = midsv.format.remove_resequence(samdict)
    for alignment in samdict:
        alignment["MIDSV"] = midsv.convert.cstag_to_midsv(alignment["CSTAG"])
    for alignment in samdict:
        alignment["CSSPLIT"] = midsv.convert.cstag_to_cssplit(alignment["CSTAG"])
    for alignment in samdict:
        alignment["QSCORE"] = midsv.convert.qual_to_qscore_cssplit(alignment["QUAL"], alignment["CSSPLIT"])
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
