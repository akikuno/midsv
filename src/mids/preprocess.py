from __future__ import annotations
from pathlib import Path
import re

###########################################################
# Read sam
###########################################################


def read_sam(path_of_sam: str) -> list[list]:
    sam = Path(path_of_sam).read_text().strip().split("\n")
    return [s.split("\t") for s in sam]


###########################################################
# Check sam format
###########################################################


def check_headers(sam: list[list]):
    """Check headers containing SN (Reference sequence name) and LN (Reference sequence length)

    Args:
        sam (list[list]): a list of lists of SAM format

    """
    sqheaders = [s for s in sam if "@SQ" in s]
    if not sqheaders:
        raise AttributeError("Input does not have @SQ header")


def check_alignments(sam: list[list]):
    """Check alignments are mapped and have long-formatted cs tag

    Args:
        sam (list[list]): a list of lists of SAM format including CS tag
    """
    is_alignment = False
    for alignment in sam:
        if "@" in alignment[0]:
            continue
        is_alignment = True
        if len(alignment) < 10:
            raise AttributeError("Alighment may not be SAM format")
        if alignment[2] == "*":
            continue
        idx_cstag = -1
        for i, a in enumerate(alignment):
            if a.startswith("cs:Z:="):
                idx_cstag = i
                break
        if idx_cstag == -1:
            raise AttributeError("Input does not have long-formatted cs tag")
        if "~" in alignment[idx_cstag]:
            raise AttributeError("Spliced long reads are currently not supported")
    if not is_alignment:
        raise AttributeError("No alignment information")


def check_sam_format(sam: list[list]):
    check_headers(sam)
    check_alignments(sam)


###########################################################
# Remove undesired reads
###########################################################


def remove_softclips(sam: list[list]) -> list[list]:
    """_summary_

    Args:
        sam (list[list]): _description_

    Returns:
        list[list]: _description_
    """
    sam_list = []
    for alignment in sam:
        if "@" in alignment[0]:
            sam_list.append(alignment)
            continue
        cigar = alignment[5]
        if "S" not in cigar:
            sam_list.append(alignment)
            continue
        left = re.sub(r"^([0-9]+S).*", r"\1", cigar)
        if left[:-1].isdigit():
            left = int(left[:-1])
            alignment[9] = alignment[9][left:]
            alignment[10] = alignment[10][left:]
        right = re.sub(r".*([0-9]+S$)", r"\1", cigar)
        if right[:-1].isdigit():
            right = int(right[:-1])
            alignment[9] = alignment[9][:-right]
            alignment[10] = alignment[10][:-right]
        sam_list.append(alignment)
    return sam_list


def remove_overlapped(sam: list[list]) -> list[list]:
    pass
