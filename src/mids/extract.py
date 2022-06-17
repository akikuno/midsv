from __future__ import annotations

import re


def headers(sam: list[list]) -> dict[str, int]:
    """Extract SN (Reference sequence name) and LN (Reference sequence length) from SQ header

    Args:
        sam (list[list]): a list of lists of SAM format

    Returns:
        dict: a dictionary containing (multiple) SN and LN
    """
    sqheaders = [s for s in sam if "@SQ" in s]
    SNLN = {}
    for sqheader in sqheaders:
        sn_ln = [sq for sq in sqheader if re.search(("SN:|LN:"), sq)]
        sn = sn_ln[0].replace("SN:", "")
        ln = sn_ln[1].replace("LN:", "")
        SNLN.update({sn: int(ln)})
    return SNLN

def alignments(sam: list[list]) -> list[dict]:
    """Extract mapped alignments from SAM

    Args:
        sam (list[list]): a list of lists of SAM format including CS tag

    Returns:
        dict: a dictionary containing QNAME, RNAME, POS, QUAL, and CSTAG
    """
    aligns = []
    idx_cstag = -1
    for alignment in sam:
        if "@" in alignment[0]:
            continue
        if alignment[2] == "*":
            continue
        if idx_cstag == -1:
            try:
                idx_cstag = [i for i, a in enumerate(alignment) if a.startswith("cs:Z=")][0]
            except IndexError:
                print(f"IndexError: CS tag is not found")
        samdict = dict(
            QNAME=alignment[0].replace(",", "_"),
            RNAME=alignment[2],
            POS=int(alignment[3]),
            QUOL=alignment[10],
            CSTAG=alignment[idx_cstag],
        )
        aligns.append(samdict)
    aligns = sorted(aligns, key=lambda x: [x["QNAME"], x["POS"]])
    return aligns
