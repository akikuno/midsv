from __future__ import annotations

from copy import deepcopy
from itertools import groupby


def is_forward_strand(flag: int) -> bool:
    """
    Determines if the read is mapped to the forward strand.

    Args:
    flag (int): The SAM flag value.

    Returns:
    bool: True if the read is mapped to the forward strand, False otherwise.
    """
    # Check if the bit 4 (0-based) is NOT set
    return not (flag & 16)


def process_inversion(current_alignment: dict[str, int | str], first_strand: bool) -> None:
    """Detect and mark inversion in the current alignment."""
    current_read_is_forward = is_forward_strand(current_alignment["FLAG"])
    if first_strand is not current_read_is_forward:
        current_alignment["MIDSV"] = current_alignment["MIDSV"].lower()


def calculate_microhomology(previous_alignment: dict[str, int | str], current_alignment: dict[str, int | str]) -> int:
    """Calculate the length of microhomology between two alignments."""
    previous_midsv = previous_alignment["MIDSV"].split(",")
    current_midsv = current_alignment["MIDSV"].split(",")
    if "QSCORE" in current_alignment:
        previous_qscore = previous_alignment["QSCORE"].split(",")
        current_qscore = current_alignment["QSCORE"].split(",")

    num_microhomology = 0
    min_length = min(len(previous_midsv), len(current_midsv))
    for i in range(1, len(current_midsv) + 1):
        if i == min_length + 1:
            break
        prev_index = len(previous_midsv) - i
        if previous_midsv[prev_index:] == current_midsv[:i]:
            if "QSCORE" in current_alignment and previous_qscore[prev_index:] != current_qscore[:i]:
                break
            num_microhomology = i
    return num_microhomology


def remove_microhomology(current_alignment: dict[str, str], num_microhomology: int) -> None:
    """Remove microhomology and update the current alignment."""
    current_midsv = current_alignment["MIDSV"].split(",")
    current_alignment["MIDSV"] = ",".join(current_midsv[num_microhomology:])

    if "QSCORE" in current_alignment:
        current_qscore = current_alignment["QSCORE"].split(",")
        current_alignment["QSCORE"] = ",".join(current_qscore[num_microhomology:])

    current_alignment["POS"] += num_microhomology


def fill_gap(sam_template: dict[str, int | str], gap: int) -> None:
    """Fill the gap between alignments with unknown nucleotides."""
    sam_template["MIDSV"] += ",=N" * gap
    if "QSCORE" in sam_template:
        sam_template["QSCORE"] += ",-1" * gap


def merge(alignments: list[dict[str, int | str]]) -> list[dict[str, int | str]]:
    """Merge splitted reads including large deletion or inversion.

    Args:
        alignments (list[dict[str, int | str]]): dictionarized SAM

    Returns:
        list[dict[str, int | str]]: SAM with joined splitted reads to single read
    """
    sam_sorted = sorted(alignments, key=lambda x: [x["QNAME"], x["POS"]])
    sam_groupby = groupby(sam_sorted, key=lambda x: x["QNAME"])
    sam_merged = []

    for *_, records in sam_groupby:
        records = list(records)
        if len(records) == 1:
            sam_merged.append(records[0])
            continue

        sam_template = deepcopy(records[0])
        first_strand = is_forward_strand(sam_template["FLAG"])

        for i, current_alignment in enumerate(records[1:], start=1):
            process_inversion(current_alignment, first_strand)

            previous_alignment = records[i - 1]
            num_microhomology = calculate_microhomology(previous_alignment, current_alignment)
            remove_microhomology(current_alignment, num_microhomology)

            previous_end = previous_alignment["POS"] + len(previous_alignment["MIDSV"].split(",")) - 1
            current_start = current_alignment["POS"] - 1

            fill_gap(sam_template, current_start - previous_end)

            sam_template["MIDSV"] += "," + current_alignment["MIDSV"]
            if "QSCORE" in sam_template:
                sam_template["QSCORE"] += "," + current_alignment["QSCORE"]

        sam_merged.append(sam_template)

    return sam_merged


def pad(alignments: list[dict[str, int | str]], sqheaders: dict[str, int]) -> list[dict[str, int | str]]:
    """Padding left and right flanks as "=" in MIDSV, "-1" in QUAL

    Args:
        sam (list[dict[str, int | str]]): dictionarized SAM
        sqheaders (dict[str, int]): dictionary as {SQ:LN}

    Returns:
        list[dict[str, int | str]]: dictionarized SAM with padding as "=N" in MIDSV and CSSPLIT, and "-1" in QUAL
    """
    alignments_padding = []
    for alignment in alignments:
        ref_length = sqheaders[alignment["RNAME"]]
        left_pad = max(0, alignment["POS"] - 1)
        right_pad = max(0, ref_length - (len(alignment["MIDSV"].split(",")) + left_pad))
        left_pad_midsv, right_pad_midsv = "=N," * left_pad, ",=N" * right_pad
        left_pad_qscore, right_pad_qscore = "-1," * left_pad, ",-1" * right_pad

        alignment["MIDSV"] = left_pad_midsv + alignment["MIDSV"] + right_pad_midsv
        if "QSCORE" in alignment:
            alignment["QSCORE"] = left_pad_qscore + alignment["QSCORE"] + right_pad_qscore

        alignments_padding.append(alignment)

    return alignments_padding


def remove_different_length(
    alignments: list[dict[str, int | str]], sqheaders: dict[str, int]
) -> list[dict[str, int | str]]:
    """remove different sequence length of the reference

    Args:
        sam (list[dict[str, int | str]]): dictionarized SAM
        sqheaders (dict[str, int]): dictionary as {SQ:LN}

    Returns:
        list[dict[str, int | str]]: filtered SAM by different sequence length of the reference
    """
    alignments_filtered = []
    for alignment in alignments:
        ref_length = sqheaders[alignment["RNAME"]]
        if len(alignment["MIDSV"].split(",")) != ref_length:
            continue
        alignments_filtered.append(alignment)
    return alignments_filtered


def select(alignments: list[dict[str, int | str]], keep: set[str] = None) -> list[dict[str, int | str]]:
    """Select QNAME, RNAME, MIDSV, CSSPLIT and QSCORE

    Args:
        alignments (list[dict[str, int | str]]): dictionarized SAM
        keep (set(str), optional): Subset of {'FLAG', 'POS', 'SEQ', 'QUAL', 'CIGAR', 'CSTAG'} to keep. Defaults to set().
    Returns:
        list[dict[str, int | str]]: dictionarized SAM of QNAME, RNAME, MIDSV, CSSPLIT and QSCORE
    """
    keep = set(keep) if keep else set()
    keys_to_delete = {"FLAG", "POS", "SEQ", "QUAL", "CIGAR", "CSTAG"} - keep
    selected = []
    for record in alignments:
        for key in keys_to_delete:
            record.pop(key)
        selected.append(record)
    return selected


###############################################################################
# main
###############################################################################


def polish(
    alignments: list[dict[str, int | str]], sqheaders: dict[str, int], keep: set[str] = None
) -> list[dict[str, int | str]]:
    """Polish SAM by merging splitted reads, padding, removing different length, and selecting fields
    Args:
        alignments (list[dict[str, int | str]]): dictionarized SAM
        sqheaders (dict[str, int]): dictionary as {SQ:LN}
        keep (set(str), optional): Subset of {'FLAG', 'POS', 'SEQ', 'QUAL', 'CIGAR', 'CSTAG'} to keep. Defaults to set().

    Returns:
        list[dict[str, int | str]]: polished SAM
    """
    alignments_polished = merge(alignments)
    alignments_polished = pad(alignments_polished, sqheaders)
    alignments_polished = remove_different_length(alignments_polished, sqheaders)
    return select(alignments_polished, keep)
