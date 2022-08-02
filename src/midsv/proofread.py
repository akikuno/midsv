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
                if "MIDSV" in alignment:
                    alignment["MIDSV"] = alignment["MIDSV"].lower()
                if "CSSPLIT" in alignment:
                    alignment["CSSPLIT"] = alignment["CSSPLIT"].lower()
            # 3. Fill in the gap between the first read and the next read with a D (deletion)
            previous_end = alignments[i - 1]["POS"] - 1
            if "MIDSV" in alignment:
                previous_end += len(alignments[i - 1]["MIDSV"].split(","))
            else:
                previous_end += len(alignments[i - 1]["CSSPLIT"].split(","))
            current_start = alignments[i]["POS"] - 1
            gap = current_start - previous_end
            if "MIDSV" in sam_dict:
                sam_dict["MIDSV"] += ",D" * gap
            if "CSSPLIT" in sam_dict:
                sam_dict["CSSPLIT"] += ",N" * gap
            if "QSCORE" in sam_dict:
                sam_dict["QSCORE"] += ",-1" * gap
            # 4. Update sam_dict
            if "MIDSV" in sam_dict:
                sam_dict["MIDSV"] += "," + alignment["MIDSV"]
            if "CSSPLIT" in sam_dict:
                sam_dict["CSSPLIT"] += "," + alignment["CSSPLIT"]
            if "QSCORE" in sam_dict:
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
        if "MIDSV" in alignment:
            rightpad = reflength - (len(alignment["MIDSV"].split(",")) + leftpad)
        else:
            rightpad = reflength - (len(alignment["CSSPLIT"].split(",")) + leftpad)
        rightpad = max(0, rightpad)
        leftpad_midsv, rightpad_midsv = "N," * leftpad, ",N" * rightpad
        leftpad_cssplit, rightpad_cssplit = "N," * leftpad, ",N" * rightpad
        leftpad_qscore, rightpad_qscore = "-1," * leftpad, ",-1" * rightpad
        if "MIDSV" in alignment:
            alignment["MIDSV"] = leftpad_midsv + alignment["MIDSV"] + rightpad_midsv
        if "CSSPLIT" in alignment:
            alignment["CSSPLIT"] = leftpad_cssplit + alignment["CSSPLIT"] + rightpad_cssplit
        if "QSCORE" in alignment:
            alignment["QSCORE"] = leftpad_qscore + alignment["QSCORE"] + rightpad_qscore
        samdict_padding.append(alignment)
    return samdict_padding


def filter_length(samdict: list[dict]) -> list[dict]:
    """Extract full length reads

    Args:
        samdict (list[dict]): dictionarized SAM

    Returns:
        List[dict]: dictionarized SAM filtered by length
    """
    if "MIDSV" in samdict[0]:
        reflength = len(samdict[0]["MIDSV"].split(","))
    else:
        reflength = len(samdict[0]["CSSPLIT"].split(","))
    if reflength < 100:
        return samdict
    samdict_filtered = []
    for alignment in samdict:
        if reflength > 1000:
            threshold = 50
        else:
            threshold = 10
        is_filter = False
        filterN = ["N"] * threshold
        if "MIDSV" in alignment:
            leftN = alignment["MIDSV"].split(",")[:threshold]
            rightN = alignment["MIDSV"].split(",")[-threshold:]
        else:
            leftN = alignment["CSSPLIT"].split(",")[:threshold]
            rightN = alignment["CSSPLIT"].split(",")[-threshold:]
        if filterN == leftN or filterN == rightN:
            is_filter = True
        if not is_filter:
            samdict_filtered.append(alignment)
    return samdict_filtered


def select(samdict: list[dict]) -> list[dict]:
    """Select QNAME, RNAME, MIDSV, CSSPLIT and QSCORE

    Args:
        sam (list[dict]): dictionarized SAM

    Returns:
        list[dict]: dictionarized SAM of QNAME, RNAME, MIDSV, CSSPLIT and QSCORE
    """
    selected = []
    for m in samdict:
        for delkeys in ["FLAG", "POS", "QUAL", "CIGAR", "CSTAG"]:
            m.pop(delkeys)
        selected.append(m)
    return selected
