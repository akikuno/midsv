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
            sam_joined.append(alignments[0])
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
                alignment["MIDSV"] = alignment["MIDSV"].lower()
                alignment["CSSPLIT"] = alignment["CSSPLIT"].lower()
            # 3. Fill in the gap between the first read and the next read with a D (deletion)
            previous_end = alignments[i - 1]["POS"] - 1
            previous_end += len(alignments[i - 1]["MIDSV"].split(","))
            current_start = alignments[i]["POS"] - 1
            gap = current_start - previous_end
            sam_dict["MIDSV"] += ",D" * gap
            sam_dict["CSSPLIT"] += ",N" * gap
            sam_dict["QSCORE"] += ",-1" * gap
            # 4. Update sam_dict
            sam_dict["MIDSV"] += "," + alignment["MIDSV"]
            sam_dict["CSSPLIT"] += "," + alignment["CSSPLIT"]
            sam_dict["QSCORE"] += "," + alignment["QSCORE"]
        sam_joined.append(sam_dict)
    return sam_joined


def pad(samdict: list[dict], sqheaders: dict) -> list[dict]:
    """Padding left and right flanks as "=" in MIDSV, "-1" in QUAL

    Args:
        sam (list[dict]): dictionarized SAM
        sqheaders (dict): dictionary as {SQ:LN}

    Returns:
        list[dict]: dictionarized SAM with padding as "N" in MIDSV and CSSPLIT, and "-1" in QUAL
    """
    samdict_padding = []
    for alignment in samdict:
        reflength = sqheaders[alignment["RNAME"]]
        leftpad = max(0, alignment["POS"] - 1)
        rightpad = reflength - (len(alignment["MIDSV"].split(",")) + leftpad)
        rightpad = max(0, rightpad)
        leftpad_midsv, rightpad_midsv = "N," * leftpad, ",N" * rightpad
        leftpad_cssplit, rightpad_cssplit = "N," * leftpad, ",N" * rightpad
        leftpad_qscore, rightpad_qscore = "-1," * leftpad, ",-1" * rightpad
        alignment["MIDSV"] = leftpad_midsv + alignment["MIDSV"] + rightpad_midsv
        alignment["CSSPLIT"] = leftpad_cssplit + alignment["CSSPLIT"] + rightpad_cssplit
        alignment["QSCORE"] = leftpad_qscore + alignment["QSCORE"] + rightpad_qscore
        samdict_padding.append(alignment)
    return samdict_padding


def select(samdict: list[dict]) -> list[dict]:
    """Select QNAME, RNAME, MIDSV, CSSPLIT and QSCORE

    Args:
        sam (list[dict]): dictionarized SAM

    Returns:
        list[dict]: dictionarized SAM of QNAME, RNAME, MIDSV, CSSPLIT and QSCORE
    """
    return [
        {"QNAME": m["QNAME"], "RNAME": m["RNAME"], "MIDSV": m["MIDSV"], "CSSPLIT": m["CSSPLIT"], "QSCORE": m["QSCORE"]}
        for m in samdict
    ]

