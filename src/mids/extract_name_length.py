from __future__ import annotations

import re

def extract_name_length(sam: list[list]) -> dict[str, int]:
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