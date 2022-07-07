from __future__ import annotations
import re
from pathlib import Path
from typing import Union
from itertools import groupby

###########################################################
# Read sam
###########################################################


def read_sam(path_of_sam: Union[str, Path]) -> list[list]:
    if "str" in str(type(path_of_sam)):
        path_of_sam = Path(path_of_sam)
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
        if "~" in alignment[idx_cstag]:
            raise AttributeError("long-read spliced alignment are currently not supported")
    if not is_alignment:
        raise AttributeError("No alignment information")


def check_sam_format(sam: list[list]):
    check_headers(sam)
    check_alignments(sam)


###########################################################
# Format headers and alignments
###########################################################


def extract_sqheaders(sam: list[list]) -> dict[str, int]:
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


def dictionarize_sam(sam: list[list]) -> list[dict]:
    """Extract mapped alignments from SAM

    Args:
        sam (list[list]): a list of lists of SAM format including CS tag

    Returns:
        dict: a dictionary containing QNAME, RNAME, POS, QUAL, CSTAG and RLEN
    """
    aligns = []
    for alignment in sam:
        if alignment[0].startswith("@"):
            continue
        if alignment[2] == "*":
            continue
        if alignment[10] == "*":
            continue
        for i, a in enumerate(alignment):
            if a.startswith("cs:Z:") and not re.search(r":[0-9]+", alignment[i]):
                idx_cstag = i
        samdict = dict(
            QNAME=alignment[0].replace(",", "_"),
            FLAG=int(alignment[1]),
            RNAME=alignment[2],
            POS=int(alignment[3]),
            CIGAR=alignment[5],
            QUAL=alignment[10],
            CSTAG=alignment[idx_cstag],
        )
        aligns.append(samdict)
    aligns = sorted(aligns, key=lambda x: [x["QNAME"], x["POS"]])
    return aligns


###########################################################
# Remove undesired reads
###########################################################


def split_cigar(cigar: str) -> list[str]:
    cigar_iter = iter(re.split(r"([MIDNSHPX=])", cigar))
    cigar_splitted = [i + op for i, op in zip(cigar_iter, cigar_iter)]
    return cigar_splitted


def remove_softclips(sam: list[dict]) -> list[dict]:
    """Remove softclip quality score from QUAL.

    Args:
        sam (list[list]): disctionalized SAM

    Returns:
        list[list]: disctionalized SAM with trimmed softclips in QUAL
    """
    sam_list = []
    for alignment in sam:
        cigar = alignment["CIGAR"]
        if "S" not in cigar:
            sam_list.append(alignment)
            continue
        cigar_split = split_cigar(cigar)
        left, right = cigar_split[0], cigar_split[-1]
        if "S" in left:
            left = int(left[:-1])
            alignment["QUAL"] = alignment["QUAL"][left:]
        if "S" in right:
            right = int(right[:-1])
            alignment["QUAL"] = alignment["QUAL"][:-right]
        sam_list.append(alignment)
    return sam_list


def remove_overlapped(samdict: list[list]) -> list[list]:
    """Remove overlapped reads within the same QNAME.

    Args:
        sam (list[list]): disctionalized SAM

    Returns:
        list[list]: disctionalized SAM with removed overlaped reads
    """
    samdict.sort(key=lambda x: x["QNAME"])
    sam_groupby = groupby(samdict, lambda x: x["QNAME"])
    sam_nonoverlapped = []
    for _, alignments in sam_groupby:
        alignments = sorted(alignments, key=lambda x: x["POS"])
        is_overraped = False
        end_of_previous_read = -1
        for alignment in alignments:
            start_of_current_read = alignment["POS"]
            if end_of_previous_read > start_of_current_read:
                is_overraped = True
                break
            alignment_length = 0
            cigar = alignment["CIGAR"]
            cigar_split = split_cigar(cigar)
            for cig in cigar_split:
                if "M" in cig or "D" in cig:
                    alignment_length += int(cig[:-1])
            end_of_previous_read = start_of_current_read + alignment_length - 1
        if not is_overraped:
            sam_nonoverlapped += alignments
    return sam_nonoverlapped
