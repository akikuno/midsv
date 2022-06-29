from __future__ import annotations
from pathlib import Path

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
    for alignment in sam:
        if "@" in alignment[0]:
            continue
        if len(alignment) < 10:
            raise AttributeError("Alighment may not be SAM format")
        if alignment[2] == "*":
            continue
        idx_cstag = [i for i, a in enumerate(alignment) if a.startswith("cs:Z:=")]
        if not idx_cstag:
            raise AttributeError("Input does not have long-formatted cs tag")


def check_sam_format(sam: list[list]):
    check_headers(sam)
    check_alignments(sam)


###########################################################
# Remove undesired reads
###########################################################

