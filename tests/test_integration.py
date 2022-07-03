from pathlib import Path
from src.mids import integrate
from src.mids import format
from src.mids import convert
from src.mids import proofread


def test_integration_integrate():
    sampath = Path("tests", "data", "real", "tyr_cslong.sam")
    sam = format.read_sam(str(sampath))

    sqheaders = format.extract_sqheaders(sam)
    samdict_polished = integrate(sam)

    for alignment in samdict_polished:
        RNAME = alignment["RNAME"]
        MIDS = alignment["MIDS"]
        QSCORE = alignment["QSCORE"]
        RLEN = sqheaders[RNAME]
        mlen = len(MIDS.split(","))
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
        samdict[i]["MIDS"] = convert.cstag_to_mids(alignment["CSTAG"])

    for i, alignment in enumerate(samdict):
        samdict[i]["QSCORE"] = convert.qual_to_qscore(alignment["QUAL"], alignment["MIDS"])

    samdict_polished = proofread.join(samdict)
    samdict_polished = proofread.pad(samdict_polished, sqheaders)
    samdict_polished = proofread.select(samdict_polished)

    for alignment in samdict_polished:
        RNAME = alignment["RNAME"]
        MIDS = alignment["MIDS"]
        QSCORE = alignment["QSCORE"]
        RLEN = sqheaders[RNAME]
        mlen = len(MIDS.split(","))
        qlen = len(QSCORE.split(","))
        assert mlen == qlen == RLEN
