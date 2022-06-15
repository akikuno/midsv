
import re
from itertools import groupby

def extract_name_and_length(sam: list) -> dict:
    """Extract SN (Reference sequence name) and LN (Reference sequence length) information at SQ header from SAM file

    Args:
        sam (list): alignment information in SAM format

    Returns:
        dict: a dictionary containing SN and LN
    """
    sqheaders = (s for s in sam if s.startswith("@SQ"))
    SNLN = {}
    for sqheader in sqheaders:
        sn_ln = [sq for sq in sqheader.split("\t") if re.search(("SN:|LN:"), sq)]
        sn = sn_ln[0].replace("SN:", "")
        ln = sn_ln[1].replace("LN:", "")
        SNLN.update({sn: ln})
    return SNLN

def convert_from_sam():
    pass