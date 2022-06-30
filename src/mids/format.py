from __future__ import annotations

import re


def extract_headers(sam: list[list]) -> dict[str, int]:
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


def dictionarize_alignments(sam: list[list]) -> list[dict]:
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
            for i, a in enumerate(alignment):
                if a.startswith("cs:Z:="):
                    idx_cstag = i
        samdict = dict(
            QNAME=alignment[0].replace(",", "_"),
            FLAG=int(alignment[1]),
            RNAME=alignment[2],
            POS=int(alignment[3]),
            QUAL=alignment[10],
            CSTAG=alignment[idx_cstag],
        )
        aligns.append(samdict)
    aligns = sorted(aligns, key=lambda x: [x["QNAME"], x["POS"]])
    return aligns


def append_reflen(sam: list[list]) -> list[dict]:
    """Append reference sequence length

    Args:
        sam (list[list]): a list of lists of SAM format including CS tag

    Returns:
        list[dict]: alignments appended LN named as "RLEN"
    """
    headers = extract_headers(sam)
    alignments = dictionarize_alignments(sam)
    alignments_with_reflen = []
    for alignment in alignments:
        alignment["RLEN"] = headers[alignment["RNAME"]]
        alignments_with_reflen.append(alignment)
    return alignments_with_reflen

