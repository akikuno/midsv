from __future__ import annotations
from itertools import groupby
from copy import deepcopy


def join(samdict: list[dict]) -> list[dict]:
    """Join splitted reads including large deletion or inversion.

    Args:
        samdict (list[dict]): dictionarized SAM

    Returns:
        list[dict]: SAM with joined splitted reads to single read
    """
    sam_sorted = sorted(samdict, key=lambda x: [x["QNAME"], x["POS"]])
    sam_groupby = groupby(sam_sorted, key=lambda x: x["QNAME"])
    sam_joined = []
    for _, alignments in sam_groupby:
        alignments = list(alignments)
        if len(alignments) == 1:
            sam_joined.append(alignments[0])
            continue
        for i, alignment in enumerate(alignments):
            # Determine the strand (strand_first) of the first read
            if i == 0:
                sam_template = deepcopy(alignment)
                if alignment["FLAG"] == 0 or alignment["FLAG"] == 2048:
                    strand_first = 0
                else:
                    strand_first = 1
                continue
            # If the strand of the next read is different from strand_first, lowercase it as an Inversion.
            if alignment["FLAG"] == 0 or alignment["FLAG"] == 2048:
                strand = 0
            else:
                strand = 1
            if strand_first != strand:
                if "MIDSV" in alignment:
                    alignment["MIDSV"] = alignment["MIDSV"].lower()
                if "CSSPLIT" in alignment:
                    alignment["CSSPLIT"] = alignment["CSSPLIT"].lower()
            previous_alignment = alignments[i - 1]
            previous_end = previous_alignment["POS"] - 1
            if "MIDSV" in alignment:
                previous_end += len(previous_alignment["MIDSV"].split(","))
            else:
                previous_end += len(previous_alignment["CSSPLIT"].split(","))
            current_start = alignment["POS"] - 1
            # Remove microhomology
            if "CSSPLIT" in alignment:
                previous_cssplit = previous_alignment["CSSPLIT"].split(",")
                current_cssplit = alignment["CSSPLIT"].split(",")
                if "QSCORE" in alignment:
                    previous_qscore = previous_alignment["QSCORE"].split(",")
                    current_qscore = alignment["QSCORE"].split(",")
                num_microhomology = 0
                for i in range(min(len(previous_cssplit), len(current_cssplit))):
                    if previous_cssplit[-i:] == current_cssplit[:i]:
                        if "QSCORE" in alignment and previous_qscore[-i:] == current_qscore[:i]:
                            num_microhomology = i
                # Update
                alignment["CSSPLIT"] = ",".join(current_cssplit[num_microhomology:])
                if "QSCORE" in alignment:
                    alignment["QSCORE"] = ",".join(current_qscore[num_microhomology:])
                if "MIDSV" in alignment:
                    current_midsv = alignment["MIDSV"].split(",")
                    alignment["MIDSV"] = ",".join(current_midsv[num_microhomology:])
                current_start += num_microhomology
                alignment["POS"] = current_start + 1
            # Fill in the gap between the first read and the next read with a D (deletion)
            gap = current_start - previous_end
            if "MIDSV" in sam_template:
                sam_template["MIDSV"] += ",D" * gap
            if "CSSPLIT" in sam_template:
                sam_template["CSSPLIT"] += ",N" * gap
            if "QSCORE" in sam_template:
                sam_template["QSCORE"] += ",-1" * gap
            # Update sam_template
            if "MIDSV" in sam_template:
                sam_template["MIDSV"] += "," + alignment["MIDSV"]
            if "CSSPLIT" in sam_template:
                sam_template["CSSPLIT"] += "," + alignment["CSSPLIT"]
            if "QSCORE" in sam_template:
                sam_template["QSCORE"] += "," + alignment["QSCORE"]
        sam_joined.append(sam_template)
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
