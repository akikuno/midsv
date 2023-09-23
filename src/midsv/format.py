from __future__ import annotations
import re
from itertools import groupby
from midsv.utils.cigar_handler import split_cigar

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
# Remove softclips
###########################################################


def process_softclip(cigar: str, seq: str, qual: str) -> tuple[str, str]:
    cigar_split = split_cigar(cigar)
    left, right = cigar_split[0], cigar_split[-1]
    left_clip = int(left[:-1]) if "S" in left else 0
    right_clip = int(right[:-1]) if "S" in right else 0
    return seq[left_clip : -right_clip or None], qual[left_clip : -right_clip or None]


def remove_softclips(samdict: list[dict]) -> list[dict]:
    """Remove softclip information from SEQ and QUAL.

    Args:
        samdict (list[dict]): disctionalized SAM

    Returns:
        list[dict]: disctionalized SAM with trimmed softclips in SEQ and QUAL
    """
    samdict_trimmed_softclips = []
    for alignment in samdict:
        if "S" in alignment["CIGAR"]:
            new_seq, new_qual = process_softclip(alignment["CIGAR"], alignment["SEQ"], alignment["QUAL"])
            updated_alignment = {**alignment, "SEQ": new_seq, "QUAL": new_qual}
            samdict_trimmed_softclips.append(updated_alignment)
        else:
            samdict_trimmed_softclips.append(alignment)
    return samdict_trimmed_softclips


###########################################################
# Remove undesired reads
###########################################################


def return_end_of_current_read(alignment: dict) -> int:
    start_of_current_read = alignment["POS"]
    cigar_split = split_cigar(alignment["CIGAR"])
    alignment_length = sum(int(cig[:-1]) for cig in cigar_split if cig[-1] in "MDN")
    return start_of_current_read + alignment_length - 1


def get_min_position(samdict: list[dict]) -> int:
    return min(alignments["POS"] for alignments in samdict)


def realign_sequence(alignment: dict, min_position: int = 0) -> dict:
    """Discard insertion, and add deletion and spliced nucreotide to unify sequence length"""

    cigar = alignment["CIGAR"]
    sequence = alignment["SEQ"]
    cigar_split = split_cigar(cigar)
    sequence_unified_length = ["N"] * (alignment["POS"] - min_position)
    start = 0
    for cig in cigar_split:
        if "M" in cig:
            end = start + int(cig[:-1])
            sequence_unified_length.append(sequence[start:end])
            start = end
        elif "I" in cig:
            start += int(cig[:-1])
        elif any(x in cig for x in ["D", "N"]):
            sequence_unified_length.append("N" * int(cig[:-1]))
    alignment["SEQ"] = "".join(sequence_unified_length)
    return alignment


def check_overlap(start_prev: int, end_prev: int, start_curr: int, end_curr: int) -> bool:
    # Check if shorter read is completely included in longer read
    return start_prev <= start_curr and end_prev >= end_curr


def has_mismatched_overlap(seq_prev: str, seq_curr: str, start: int, end: int) -> bool:
    # Check for DNA sequence overlap but not the same DNA sequence
    for prev, curr in zip(seq_prev[start:end], seq_curr[start:end]):
        if prev != "N" and curr != "N" and prev != curr:
            return True
    return False


def remove_resequence(samdict: list[dict]) -> list[dict]:
    """Remove non-microhomologic overlapped reads within the same QNAME.
    The overlapped sequences can be (1) realignments by microhomology or (2) resequence by sequencing error.
    The 'realignments' is not sequencing errors, and it preserves the same sequence.
    In contrast, the 'resequence' is a sequencing error with the following characteristics:
    (1) The shorter reads that are completely included in the longer reads
    (2) Overlapped but not the same DNA sequence
    The resequenced fragments will be discarded and the longest alignment will be retain.
    Example reads are in `tests/data/overlap/real_overlap.sam` and `tests/data/overlap/real_overlap2.sam`

    Args:
        samdict (list[dict]): disctionalized SAM

    Returns:
        list[dict]: disctionalized SAM with removed overlaped reads
    """
    samdict.sort(key=lambda x: x["QNAME"])
    sam_groupby = groupby(samdict, lambda x: x["QNAME"])
    sam_nonoverlapped = []
    for _, alignments in sam_groupby:
        alignments = list(alignments)
        if len(alignments) == 1:
            sam_nonoverlapped += alignments
            continue

        min_position = get_min_position(alignments)
        alignments = [realign_sequence(alignment, min_position) for alignment in alignments]
        alignments.sort(key=lambda x: x["POS"])

        is_overlapped = False
        prev_alignment = alignments[0]
        start_prev, end_prev = prev_alignment["POS"] - 1, return_end_of_current_read(prev_alignment)

        for curr_alignment in alignments[1:]:
            start_curr, end_curr = curr_alignment["POS"] - 1, return_end_of_current_read(curr_alignment)

            if check_overlap(start_prev, end_prev, start_curr, end_curr):
                is_overlapped = True
                break

            start_overlap, end_overlap = max(start_prev, start_curr), min(end_prev, end_curr)
            if has_mismatched_overlap(prev_alignment["SEQ"], curr_alignment["SEQ"], start_overlap, end_overlap):
                is_overlapped = True
                break

            start_prev, end_prev = start_curr, end_curr
            prev_alignment = curr_alignment

        if is_overlapped:
            # Retain the longest alignment
            longest_alignment = max(alignments, key=lambda x: len(x["SEQ"]))
            sam_nonoverlapped.append(longest_alignment)
        else:
            sam_nonoverlapped += alignments

    return sam_nonoverlapped
