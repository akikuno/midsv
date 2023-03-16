from __future__ import annotations
import re
from itertools import groupby
from copy import deepcopy

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
        if alignment[9] == "*":
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
            SEQ=alignment[9],
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
    """Remove softclip information from SEQ and QUAL.

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
            alignment["SEQ"] = alignment["SEQ"][left:]
            alignment["QUAL"] = alignment["QUAL"][left:]
        if "S" in right:
            right = int(right[:-1])
            alignment["SEQ"] = alignment["SEQ"][:-right]
            alignment["QUAL"] = alignment["QUAL"][:-right]
        sam_list.append(alignment)
    return sam_list


def return_end_of_current_read(alignment:dict) -> int:
    start_of_current_read = alignment["POS"]
    cigar = alignment["CIGAR"]
    cigar_split = split_cigar(cigar)
    alignment_length = 0
    for cig in cigar_split:
        if "M" in cig or "D" in cig or "N" in cig:
            alignment_length += int(cig[:-1])
    return start_of_current_read + alignment_length - 1


def realign_sequence(alignment: dict) -> dict:
    """Discard insertion, and add deletion and spliced nucreotide to unify sequence length
    """
    cigar = alignment["CIGAR"]
    cigar_split = split_cigar(cigar)
    sequence = alignment["SEQ"]
    sequence_ignored = ["N"] * alignment["POS"]
    start = 0
    for cig in cigar_split:
        if "M" in cig:
            end = start + int(cig[:-1])
            sequence_ignored.append(sequence[start: end])
            start = end
        elif "I" in cig:
            start += int(cig[:-1])
        elif any(x in cig for x in ["D", "N"]):
            sequence_ignored.append("N" * int(cig[:-1]))
    realignment = deepcopy(alignment)
    realignment["SEQ"] = "".join(sequence_ignored)
    return realignment



def remove_resequence(samdict: list[list]) -> list[list]:
    """Remove non-microhomologic overlapped reads within the same QNAME.
    The overlapped sequences can be (1) realignments by microhomology or (2) resequence by sequencing error.
    The 'realignments' is not sequencing errors, and it preserves the same sequence.
    In contrast, the 'resequence' is a sequencing error with the following characteristics:
    (1) The shorter reads are completely included in the longer reads
    (2) Overlapped but not the same DNA sequence
    The resequenced fragments will be discarded and the longest alignment will be retain.
    Example reads are in `tests/data/overlap/real_overlap.sam` and `tests/data/overlap/real_overlap2.sam`

    Args:
        sam (list[list]): disctionalized SAM

    Returns:
        list[list]: disctionalized SAM with removed overlaped reads
    """
    samdict.sort(key=lambda x: x["QNAME"])
    sam_groupby = groupby(samdict, lambda x: x["QNAME"])
    sam_nonoverlapped = []
    for _, alignments in sam_groupby:
        alignments = list(alignments)
        if len(alignments) == 1:
            sam_nonoverlapped += alignments
            continue
        alignments = [realign_sequence(alignment) for alignment in alignments]
        alignments = sorted(alignments, key=lambda x: [x["POS"], -len(x["SEQ"])])
        is_overraped = False
        end_of_previous_read = -1
        previous_read = alignments[0]["SEQ"]
        for i, alignment in enumerate(alignments):
            if i == 0:
                start_of_previous_read = alignment["POS"] - 1
                end_of_previous_read = return_end_of_current_read(alignment)
                continue
            start_of_current_read = alignment["POS"] - 1
            end_of_current_read = return_end_of_current_read(alignment)
            # (1) The shorter reads are completely included in the longer reads
            if start_of_previous_read <= start_of_current_read and end_of_previous_read >= end_of_current_read:
                is_overraped = True
                break
            else:
                start_overlap = max(start_of_previous_read, start_of_current_read)
                end_overlap = min(end_of_previous_read, end_of_current_read)
                for prev, curr in zip(previous_read[start_overlap: end_overlap], alignment["SEQ"][start_overlap: end_overlap]):
                    if prev == "N" or curr == "N":
                        continue
                    # (2) Overlapped but not the same DNA sequence
                    if prev != curr:
                        is_overraped = True
                        break
            start_of_previous_read = start_of_current_read
            end_of_previous_read = return_end_of_current_read(alignment)
        if is_overraped:
            # The longest alignment will be retain
            sam_nonoverlapped.append(alignments[0])
        else:
            sam_nonoverlapped += alignments
    return sam_nonoverlapped
