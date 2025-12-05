from __future__ import annotations

import re
from collections.abc import Iterator
from pathlib import Path

from midsv import io

###########################################################
# Validate keep argument
###########################################################


def keep_argument(keep: str | list[str]) -> list[str]:
    if keep is None:
        keep = set()
    elif isinstance(keep, str):
        keep = [keep]

    keep_set = set(keep)

    if not keep_set.issubset({"FLAG", "POS", "SEQ", "QUAL", "CIGAR", "CSTAG"}):
        raise ValueError("'keep' must be a subset of {'FLAG', 'POS', 'SEQ', 'QUAL', 'CIGAR', 'CSTAG'}")

    return keep


###########################################################
# Validate sam format
###########################################################


def sam_headers(sam: list[list[str]] | Iterator[list[str]]) -> None:
    """Check headers containing SN (Reference sequence name) and LN (Reference sequence length)

    Args:
        sam (list[list[str]] | Iterator[list[str]]): a list of lists of SAM format

    """
    sqheaders = [s for s in sam if "@SQ" in s]
    if not sqheaders:
        raise ValueError("Input does not have @SQ header")


def sam_alignments(sam: list[list[str]] | Iterator[list[str]], qscore: bool = False) -> None:
    """Check alignments are mapped and have long-formatted cs tag

    Args:
        sam (list[list[str]]): a list of lists of SAM format including CS tag
    """
    has_alignment = False
    for alignment in sam:
        if alignment[0].startswith("@"):
            continue

        if len(alignment) < 10:
            raise ValueError("Alignment may not be SAM format because it has less than 10 columns")

        if alignment[2] == "*" or alignment[9] == "*":  # No alignment of reads
            continue
        has_alignment = True

        if qscore and alignment[10] == "*":
            raise ValueError("Input does not have QUAL information")

        idx_cstag = None
        for i, a in enumerate(alignment):
            if a.startswith("cs:Z:") and not re.search(r":[0-9]+", alignment[i]):
                idx_cstag = i
                break
        if idx_cstag is None:
            raise ValueError("Input does not have long-formatted cs tag")

    if not has_alignment:
        raise ValueError("No alignment information")


###########################################################
# Validate sam format
###########################################################


def validate_sam(path_sam: str | Path, qscore: bool = False) -> None:
    """Check headers containing SN (Reference sequence name) and LN (Reference sequence length)

    Args:
        sam (list[list[str]] | Iterator[list[str]]): a list of lists of SAM format

    """
    path_sam = Path(path_sam)
    if not path_sam.exists():
        raise FileNotFoundError(f"{path_sam} does not exist")
    sam_headers(io.read_sam(path_sam))
    sam_alignments(io.read_sam(path_sam), qscore)
