from __future__ import annotations

import re
from collections.abc import Iterator
from itertools import groupby

###########################################################
# Format headers and alignments
###########################################################


def extract_sqheaders(sam: list[list[str]] | Iterator[list[str]]) -> dict[str, int]:
    """Extract SN (Reference sequence name) and LN (Reference sequence length) from SQ header

    Args:
        sam (list[list[str]] | Iterator[list[str]]): a list of lists of SAM format

    Returns:
        dict[str, int]: a dictionary containing (multiple) SN and LN
    """
    sqheaders = [s for s in sam if "@SQ" in s]
    header_snln = {}
    for sqheader in sqheaders:
        snln = [sq for sq in sqheader if re.search(("SN:|LN:"), sq)]
        sn = snln[0].replace("SN:", "")
        ln = snln[1].replace("LN:", "")
        header_snln.update({sn: int(ln)})
    return header_snln


###########################################################
# Remove undesired reads
###########################################################


def split_cigar(cigar: str) -> list[str]:
    cigar_iter = iter(re.split(r"([MIDNSHPX=])", cigar))
    cigar_splitted = [i + op for i, op in zip(cigar_iter, cigar_iter)]
    return cigar_splitted


def remove_softclips(alignments: list[dict[str, str | int]]) -> list[dict[str, str | int]]:
    """Remove softclip information from SEQ and QUAL.

    Args:
        alignments (list[dict[str, str | int]]): disctionalized alignments

    Returns:
        list[dict[str, str | int]]: disctionalized SAM with trimmed softclips in QUAL
    """
    alignments_softclips_removed = []
    for alignment in alignments:
        cigar = alignment["CIGAR"]
        if "S" not in cigar:
            alignments_softclips_removed.append(alignment)
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
        alignments_softclips_removed.append(alignment)
    return alignments_softclips_removed


def _return_end_of_current_read(alignment: dict[str, str | int]) -> int:
    start_of_current_read = alignment["POS"]
    cigar = alignment["CIGAR"]
    cigar_split = split_cigar(cigar)
    alignment_length = 0
    for cig in cigar_split:
        if "M" in cig or "D" in cig or "N" in cig:
            alignment_length += int(cig[:-1])
    return start_of_current_read + alignment_length - 1


def _padding_n_to_sequence(alignment: dict[str, str | int]) -> dict[str, str | int]:
    """Discard insertion, and add deletion and spliced nucleotides to unify sequence length"""
    cigar = alignment["CIGAR"]
    cigar_operations = split_cigar(cigar)
    sequence = alignment["SEQ"]

    unified_sequence = []
    unified_sequence.append("N" * (alignment["POS"] - 1))
    start_position = 0

    def process_match(operation: str):
        nonlocal start_position
        length = int(operation[:-1])
        end_position = start_position + length
        unified_sequence.append(sequence[start_position:end_position])
        start_position = end_position

    def process_insertion(operation: str):
        nonlocal start_position
        length = int(operation[:-1])
        start_position += length

    def process_deletion_or_splice(operation: str):
        length = int(operation[:-1])
        unified_sequence.append("N" * length)

    for operation in cigar_operations:
        if "M" in operation:
            process_match(operation)
        elif "I" in operation:
            process_insertion(operation)
        elif "D" in operation or "N" in operation:
            process_deletion_or_splice(operation)

    alignment["SEQ"] = "".join(unified_sequence)
    return alignment


def remove_resequence(alignments: list[dict[str, str | int]]) -> list[dict[str, str | int]]:
    """Remove non-microhomologic overlapped reads within the same QNAME.
    The overlapped sequences can be (1) realignments by microhomology or (2) resequence by sequencing error.
    The 'realignments' is not sequencing errors, and it preserves the same sequence.
    In contrast, the 'resequence' is a sequencing error with the following characteristics:
    (1) The shorter reads that are completely included in the longer reads
    (2) Overlapped but not the same DNA sequence
    The resequenced fragments will be discarded and the longest alignment will be retain.
    Example reads are in `tests/data/overlap/real_overlap.sam` and `tests/data/overlap/real_overlap2.sam`

    Args:
        alignments (list[dict[str, str | int]]): disctionalized alignments

    Returns:
        list[dict[str, str | int]]: disctionalized SAM with removed overlaped reads
    """

    def is_resequence(prev_read: dict[str, str | int], curr_read: dict[str, str | int]) -> bool:
        """Check if the current read is a resequence of the previous read."""
        start_overlap = max(prev_read["POS"], curr_read["POS"])
        end_overlap = min(_return_end_of_current_read(prev_read), _return_end_of_current_read(curr_read))

        if prev_read["POS"] <= curr_read["POS"] and _return_end_of_current_read(
            prev_read
        ) >= _return_end_of_current_read(curr_read):
            return True  # Completely contained

        overlap_prev = prev_read["SEQ"][start_overlap - prev_read["POS"] : end_overlap - prev_read["POS"]]
        overlap_curr = curr_read["SEQ"][start_overlap - curr_read["POS"] : end_overlap - curr_read["POS"]]

        for prev_base, curr_base in zip(overlap_prev, overlap_curr):
            if prev_base != "=N" and curr_base != "=N" and prev_base != curr_base:
                return True  # Overlapped but different sequences

        return False

    alignments.sort(key=lambda x: x["QNAME"])
    grouped_alignments = groupby(alignments, key=lambda x: x["QNAME"])
    filtered_alignments = []

    for _, group in grouped_alignments:
        group = sorted((_padding_n_to_sequence(alignment) for alignment in group), key=lambda x: x["POS"])

        retained_alignments = []
        longest_alignment = max(group, key=lambda x: len(x["SEQ"]))

        for i, alignment in enumerate(group):
            if i == 0:
                retained_alignments.append(alignment)
                continue

            if is_resequence(retained_alignments[-1], alignment):
                retained_alignments = [longest_alignment]
                break

            retained_alignments.append(alignment)

        filtered_alignments.extend(retained_alignments)

    return filtered_alignments


###########################################################
# alignments_to_dict
###########################################################


def alignments_to_dict(sam: list[list[str]] | Iterator[list[str]]) -> list[dict[str, str | int]]:
    """Extract mapped alignments from SAM

    Args:
        sam (list[list[str]] | Iterator[list[str]]): a list of lists of SAM format including CS tag

    Returns:
        list[dict[str, str | int]]: a dictionary containing QNAME, RNAME, POS, QUAL, CSTAG and RLEN
    """
    aligns = []
    for alignment in sam:
        if alignment[0].startswith("@"):
            continue
        if alignment[2] == "*" or alignment[9] == "*":
            continue

        idx_cstag = None
        for i, a in enumerate(alignment):
            if a.startswith("cs:Z:") and not re.search(r":[0-9]+", alignment[i]):
                idx_cstag = i
        if idx_cstag is None:
            continue

        alignments = dict(
            QNAME=alignment[0].replace(",", "_"),
            FLAG=int(alignment[1]),
            RNAME=alignment[2],
            POS=int(alignment[3]),
            CIGAR=alignment[5],
            SEQ=alignment[9],
            QUAL=alignment[10],
            CSTAG=alignment[idx_cstag],
        )
        aligns.append(alignments)
    return aligns


def organize_alignments_to_dict(sam: list[list[str]] | Iterator[list[str]]) -> list[dict[str, str | int]]:
    """Extract mapped alignments from SAM

    Args:
        sam (list[list[str]] | Iterator[list[str]]): a list of lists of SAM format including CS tag

    Returns:
        list[dict[str, str | int]]: a dictionary containing QNAME, RNAME, POS, QUAL, CSTAG and RLEN
    """
    aligns = alignments_to_dict(sam)
    aligns = remove_softclips(aligns)
    aligns = remove_resequence(aligns)

    return sorted(aligns, key=lambda x: [x["QNAME"], x["POS"]])
