from __future__ import annotations
from itertools import groupby
from copy import deepcopy


def join(sam: list[dict]) -> list[dict]:
    """Join splitted reads including large deletion or inversion.

    Args:
        sam (list[dict]): dictionarized SAM

    Returns:
        list[dict]: SAM with joined splitted reads to single read
    """
    sam_sorted = sorted(sam, key=lambda x: [x["QNAME"], x["POS"]])
    sam_groupby = groupby(sam_sorted, key=lambda x: x["QNAME"])
    sam_joined = []
    for _, alignments in sam_groupby:
        alignments = list(alignments)
        if len(alignments) == 1:
            sam_joined.append(alignments)
            continue
        for i, alignment in enumerate(alignments):
            # 1. Determine the strand (strand_first) of the first read
            if i == 0:
                sam_dict = deepcopy(alignment)
                if alignment["FLAG"] == 0 or alignment["FLAG"] == 2048:
                    strand_first = 0
                else:
                    strand_first = 1
                continue
            # 2. If the strand of the next read is different from strand_first, lowercase it as an Inversion.
            if alignment["FLAG"] == 0 or alignment["FLAG"] == 2048:
                strand = 0
            else:
                strand = 1
            if strand_first != strand:
                alignment["MIDS"] = alignment["MIDS"].lower()
            # 3. Fill in the gap between the first read and the next read with a D (deletion)
            previous_end = alignments[i - 1]["POS"] - 1
            previous_end += len(alignments[i - 1]["MIDS"].split(","))
            current_start = alignments[i]["POS"] - 1
            gap = current_start - previous_end
            sam_dict["MIDS"] += ",D" * gap
            sam_dict["QSCORE"] += ",-1" * gap
            # 4. Update sam_dict
            sam_dict["MIDS"] += "," + alignment["MIDS"]
            sam_dict["QSCORE"] += "," + alignment["QSCORE"]
        sam_joined.append(sam_dict)
    return sam_joined


def pad(sam: list[dict], sqheaders: dict) -> list[dict]:
    pass