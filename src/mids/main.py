from __future__ import annotations
from src.mids import format
from src.mids import convert
from src.mids import proofread


def transform(sam: list[list]) -> list[dict]:
    """Integrated function to perform MIDS conversion

    Args:
        sam (list[list]): Lists ot SAM format

    Returns:
        list[dict]: Dictionary containing QNAME, RNAME, MIDS, and QSCORE
    """

    format.check_sam_format(sam)

    sqheaders = format.extract_sqheaders(sam)
    samdict = format.dictionarize_sam(sam)

    samdict = format.remove_softclips(samdict)
    samdict = format.remove_overlapped(samdict)

    for i, alignment in enumerate(samdict):
        samdict[i]["MIDS"] = convert.cstag_to_mids(alignment["CSTAG"])
        samdict[i]["CSSPLIT"] = convert.cstag_to_cssplit(alignment["CSTAG"])
        samdict[i]["QSCORE"] = convert.qual_to_qscore(alignment["QUAL"], alignment["MIDS"])

    samdict_polished = proofread.join(samdict)
    samdict_polished = proofread.pad(samdict_polished, sqheaders)
    samdict_polished = proofread.select(samdict_polished)
    return samdict_polished
