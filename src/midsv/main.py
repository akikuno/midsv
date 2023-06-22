from __future__ import annotations

from typing import Generator
from . import validate
from . import format
from . import convert
from . import proofread


def transform(sam: list[list] | Generator[list], midsv: bool = True, cssplit: bool = True, qscore: bool = True, keep: set(list[str]) = set()) -> list[dict]:
    """Integrated function to perform MIDSV conversion

    Args:
        sam (list[list]): Lists ot SAM format
        midsv (bool, optional): Output MIDSV. Defaults to True.
        cssplit (bool, optional): Output CSSPLIT. Defaults to True.
        qscore (bool, optional): Output QSCORE. Require `midsv == True` or `cssplit == True`. Defaults to True.
        keep (set(str), optional): Subset of {'FLAG', 'POS', 'SEQ', 'QUAL', 'CIGAR', 'CSTAG'} to keep. Defaults to set().
    Returns:
        list[dict]: Dictionary containing QNAME, RNAME, MIDSV, and QSCORE
    """

    if midsv or cssplit:
        pass
    else:
        raise ValueError("Either midsv or cssplit must be True")

    if keep != set() and not keep.issubset({"FLAG", "POS", "SEQ", "QUAL", "CIGAR", "CSTAG"}):
        raise ValueError("'keep' must be a subset of {'FLAG', 'POS', 'SEQ', 'QUAL', 'CIGAR', 'CSTAG'}")

    sam = list(sam)
    validate.sam_headers(sam)
    validate.sam_alignments(sam)

    sqheaders = format.extract_sqheaders(sam)
    samdict = format.dictionarize_sam(sam)

    if qscore and any(s["QUAL"] == "*" for s in samdict):
        raise ValueError("qscore must be False because the input does not contain QUAL")

    samdict = format.remove_softclips(samdict)
    samdict = format.remove_resequence(samdict)

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
    samdict_polished = proofread.remove_different_length(samdict_polished, sqheaders)
    samdict_polished = proofread.select(samdict_polished, keep)
    return samdict_polished
