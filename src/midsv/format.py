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


def dictionarize_sam(sam: list[list[str]] | Iterator[list[str]]) -> list[dict[str, str | int]]:
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

    return sorted(aligns, key=lambda x: [x["QNAME"], x["POS"]])


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
    unified_sequence = ["N"] * alignment["POS"]
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
    alignments.sort(key=lambda x: x["QNAME"])
    sam_groupby = groupby(alignments, lambda x: x["QNAME"])
    sam_nonoverlapped = []
    for _, alignments in sam_groupby:
        alignments = list(alignments)
        if len(alignments) == 1:
            sam_nonoverlapped += alignments
            continue
        alignments = [_padding_n_to_sequence(alignment) for alignment in alignments]
        alignments = sorted(alignments, key=lambda x: [x["POS"]])
        is_overraped = False
        end_of_previous_read = -1
        previous_read = alignments[0]["SEQ"]
        for i, alignment in enumerate(alignments):
            if i == 0:
                start_of_previous_read = alignment["POS"] - 1
                end_of_previous_read = _return_end_of_current_read(alignment)
                continue
            start_of_current_read = alignment["POS"] - 1
            end_of_current_read = _return_end_of_current_read(alignment)
            # (1) The shorter reads that are completely included in the longer reads
            if start_of_previous_read <= start_of_current_read and end_of_previous_read >= end_of_current_read:
                is_overraped = True
                break
            else:
                start_overlap = max(start_of_previous_read, start_of_current_read)
                end_overlap = min(end_of_previous_read, end_of_current_read)
                for prev, curr in zip(
                    previous_read[start_overlap:end_overlap], alignment["SEQ"][start_overlap:end_overlap]
                ):
                    if prev == "N" or curr == "N":
                        continue
                    # (2) Overlapped but not the same DNA sequence
                    if prev != curr:
                        is_overraped = True
                        break
            start_of_previous_read = start_of_current_read
            end_of_previous_read = _return_end_of_current_read(alignment)
        if is_overraped:
            # The longest alignment will be retain
            alignment = sorted(alignments, key=lambda x: -len(x["SEQ"]))[0]
            sam_nonoverlapped.append(alignment)
        else:
            sam_nonoverlapped += alignments
    return sam_nonoverlapped
