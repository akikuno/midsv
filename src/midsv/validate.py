from __future__ import annotations
import re


###########################################################
# Validate sam format
###########################################################


def sam_headers(sam: list[list]) -> None:
    """Check headers containing SN (Reference sequence name) and LN (Reference sequence length)

    Args:
        sam (list[list]): a list of lists of SAM format

    """
    sqheaders = [s for s in sam if "@SQ" in s]
    if not sqheaders:
        raise AttributeError("Input does not have @SQ header")


def sam_alignments(sam: list[list]) -> None:
    """Check alignments are mapped and have long-formatted cs tag

    Args:
        sam (list[list]): a list of lists of SAM format including CS tag
    """
    is_alignment = False
    for alignment in sam:
        if alignment[0].startswith("@"):
            continue
        if alignment[2] == "*":
            continue
        is_alignment = True
        if len(alignment) < 10:
            raise AttributeError("Alighment may not be SAM format")
        idx_cstag = -1
        for i, a in enumerate(alignment):
            if a.startswith("cs:Z:") and not re.search(r":[0-9]+", alignment[i]):
                idx_cstag = i
                break
        if idx_cstag == -1:
            raise AttributeError("Input does not have long-formatted cs tag")
        # if "~" in alignment[idx_cstag]:
        #     raise AttributeError("long-read spliced alignment are currently not supported")
    if not is_alignment:
        raise AttributeError("No alignment information")
