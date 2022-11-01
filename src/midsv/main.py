from __future__ import annotations

from . import convert
from . import format
from . import proofread


def transform(sam: list[list], midsv: bool = True, cssplit: bool = True, qscore: bool = True) -> list[dict]:
    """Integrated function to perform MIDSV conversion

    Args:
        sam (list[list]): Lists ot SAM format
        midsv (bool, optional): Output MIDSV. Defaults to True.
        cssplit (bool, optional): Output CSSPLIT. Defaults to True.
        qscore (bool, optional): Output QSCORE. Require `midsv == True` or `cssplit == True`. Defaults to True.

    Returns:
        list[dict]: Dictionary containing QNAME, RNAME, MIDSV, and QSCORE
    """

    if midsv or cssplit:
        pass
    else:
        raise ValueError("Either midsv or cssplit must be True")

    format.check_sam_format(sam)

    sqheaders = format.extract_sqheaders(sam)
    samdict = format.dictionarize_sam(sam)

    if qscore and any(s["QUAL"] == "*" for s in samdict):
        raise ValueError("qscore must be False because the input does not contain QUAL")

    samdict = format.remove_softclips(samdict)
    samdict = format.remove_overlapped(samdict)

    for alignment in samdict:
        if midsv:
            alignment["MIDSV"] = convert.cstag_to_midsv(alignment["CSTAG"])
        if cssplit:
            alignment["CSSPLIT"] = convert.cstag_to_cssplit(alignment["CSTAG"])
        if midsv and qscore:
            alignment["QSCORE"] = convert.qual_to_qscore_midsv(alignment["QUAL"], alignment["MIDSV"])
        elif cssplit and qscore:
            alignment["QSCORE"] = convert.qual_to_qscore_cssplit(alignment["QUAL"], alignment["CSSPLIT"])

    samdict_polished = proofread.join(samdict)
    samdict_polished = proofread.pad(samdict_polished, sqheaders)
    samdict_polished = proofread.select(samdict_polished)
    return samdict_polished
