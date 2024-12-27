from pathlib import Path

from src import midsv
from src.midsv import converter, formatter, io, polisher, validator


def test_integration_midsv():
    path_sam = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    test = midsv.transform(path_sam, qscore=False)
    answer = Path("tests", "data", "integrate", "answer_integrate_midsv.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_all_true():
    path_sam = Path("tests", "data", "integrate", "subindelinv_cslong_10bp.sam")
    test = midsv.transform(path_sam, qscore=True)
    answer = Path("tests", "data", "integrate", "answer_integrate_all_true.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_integration_real_sam():
    path_sam = Path("tests", "data", "real", "tyr_cslong.sam")
    sqheaders = formatter.extract_sqheaders(io.read_sam(path_sam))
    samdict_polished = midsv.transform(path_sam, qscore=True)
    for alignment in samdict_polished:
        RNAME = alignment["RNAME"]
        MIDSV = alignment["MIDSV"]
        QSCORE = alignment["QSCORE"]
        mlen = len(MIDSV.split(","))
        qlen = len(QSCORE.split(","))
        answer = sqheaders[RNAME]
        assert mlen == qlen == answer


def test_integration_real_splicing():
    path_sam = Path("tests/data/splicing/real_splicing.sam")
    sqheaders = formatter.extract_sqheaders(io.read_sam(path_sam))
    test = midsv.transform(path_sam, qscore=True)
    mlen = len(test[0]["MIDSV"].split(","))
    qlen = len(test[0]["QSCORE"].split(","))
    answer = sqheaders["deletion"]
    assert mlen == qlen == answer


def test_integration_eachcomponent():
    path_sam = Path("tests", "data", "real", "tyr_cslong.sam")
    qscore = True

    sam = midsv.io.read_sam(str(path_sam))
    sam = list(sam)
    validator.sam_headers(io.read_sam(path_sam))
    validator.sam_alignments(io.read_sam(path_sam))
    sqheaders = formatter.extract_sqheaders(io.read_sam(path_sam))
    alignments = formatter.organize_alignments_to_dict(io.read_sam(path_sam))

    alignments = converter.convert(alignments, qscore)

    alignments = polisher.polish(alignments, sqheaders)

    for alignment in alignments:
        RNAME = alignment["RNAME"]
        MIDSV = alignment["MIDSV"]
        QSCORE = alignment["QSCORE"]
        RLEN = sqheaders[RNAME]
        mlen = len(MIDSV.split(","))
        qlen = len(QSCORE.split(","))
        assert RLEN == mlen == qlen
