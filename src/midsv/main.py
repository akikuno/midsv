from __future__ import annotations

from collections.abc import Iterator

from midsv import convert, format, proofread, validate


def transform(
    sam: list[list[str]] | Iterator[list[str]],
    qscore: bool = True,
    keep: list[str] = None,
) -> list[dict]:
    """Integrated function to perform MIDSV conversion.

    Args:
        sam (list[list[str]] | Iterator[list[str]]): List or Iterator of SAM format.
        qscore (bool, optional): Output QSCORE. Require `midsv == True` or `cssplit == True`. Defaults to True.
        keep (set[str], optional): Subset of {'FLAG', 'POS', 'SEQ', 'QUAL', 'CIGAR', 'CSTAG'} to keep. Defaults to set().

    Returns:
        list[dict[str, str]]: Dictionary containing QNAME, RNAME, MIDSV, and QSCORE.
    """

    keep = set(keep) if keep else set()
    if keep != set() and not keep.issubset({"FLAG", "POS", "SEQ", "QUAL", "CIGAR", "CSTAG"}):
        raise ValueError("'keep' must be a subset of {'FLAG', 'POS', 'SEQ', 'QUAL', 'CIGAR', 'CSTAG'}")

    sam = list(sam)
    validate.sam_headers(sam)
    validate.sam_alignments(sam)

    sqheaders = format.extract_sqheaders(sam)
    samdict: list[dict[str | int]] = format.dictionarize_sam(sam)

    if qscore and any(s["QUAL"] == "*" for s in samdict):
        raise ValueError("qscore must be False because the input does not contain QUAL")

    samdict = format.remove_softclips(samdict)
    samdict = format.remove_resequence(samdict)

    for alignment in samdict:
        alignment["MIDSV"] = convert.cstag_to_midsv(alignment["CSTAG"])
        if qscore:
            alignment["QSCORE"] = convert.qual_to_qscore_midsv(alignment["QUAL"], alignment["MIDSV"])

    samdict_polished = proofread.join(samdict)
    samdict_polished = proofread.pad(samdict_polished, sqheaders)
    samdict_polished = proofread.remove_different_length(samdict_polished, sqheaders)
    samdict_polished = proofread.select(samdict_polished, keep)

    return samdict_polished
