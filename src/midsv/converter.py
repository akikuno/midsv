from __future__ import annotations

import re

###########################################################
# MIDSV conversion (from CS tag to MIDSV)
###########################################################


def split_cstag(cstag: str) -> list[str]:
    """Split cstag
    Args:
        cstag (str): a string of CS tags
    Returns:
        list[str]: splitted CS tags
    Example:
        >>> from midsv import convert
        >>> cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
        >>> convert.split_cstag(cstag)
        "['=ACGT', '*ag', '=C', '-g', '=T', '+t', '=ACGT']"
    """
    cstag_splitted = []
    for cs in re.split(r"(=|\*|\-|\+|\~)", cstag):
        if not cs or cs == "cs:Z:":
            continue
        if re.match(r"(=|\*|\-|\+|\~)", cs):
            cstag_splitted.append(cs)
        else:
            cstag_splitted[-1] += cs
    return cstag_splitted


###########################################################
# Convert to MIDSV
###########################################################


def _process_insertion(cstag_splitted: list[str], i: int, midsv_tags: list[str]) -> None:
    """Process insertion operation in cstag.

    Args:
        cstag_splitted (list[str]): Split cstag components.
        i (int): Current index.
        midsv_tags (list[str]): List to append processed results.
    """
    insertion: str = "+" + "|+".join(cstag_splitted[i][1:])
    if i + 1 == len(cstag_splitted):
        midsv_tags.append(insertion)
        return

    next_cstag: str = cstag_splitted[i + 1]
    next_op: str = next_cstag[0]

    if next_op == "*":
        midsv_tags.append(insertion + "|" + next_cstag)
        cstag_splitted[i + 1] = "*"
    elif next_op == "~":
        insertion += "|=N"
        midsv_tags.append(insertion)
        splice_match = re.match(r"([a-z]+)([0-9]+)([a-z]+)", next_cstag[1:])
        if splice_match:
            left, splice, right = splice_match.groups()
            cstag_splitted[i + 1] = f"{next_op}{left}{int(splice) - 1}{right}"
    else:
        insertion += "|" + next_cstag[:2]
        midsv_tags.append(insertion)
        cstag_splitted[i + 1] = next_op + next_cstag[2:]


def _process_splice(cs: str, midsv_tags: list[str]) -> None:
    """Process splice operation in cstag.

    Args:
        cs (str): Current cstag component.
        midsv_tags (list[str]): List to append processed results.
    """
    splice_match = re.match(r"([a-z]+)([0-9]+)([a-z]+)", cs[1:])
    if splice_match:
        _, splice, _ = splice_match.groups()
        midsv_tags.extend(["=N"] * int(splice))


def _process_match(cs: str, midsv_tags: list[str]) -> None:
    """Process default operation in cstag.

    Args:
        cs (str): Current cstag component.
        midsv_tags (list[str]): List to append processed results.
    """
    cs_list: list[str] = list(cs[1:])
    midsv_tags.append(cs[0] + f",{cs[0]}".join(cs_list))


def cstag_to_midsv(cstag: str) -> str:
    """Generate MIDSV, a comma-separated nucreotide sequence

    Args:
        cstag (str): a long format cstag

    Returns:
        str: MIDSV

    Examples:
        >>> cstag = "cs:Z:=A+ttt=CC-aa=T*ag=TT"
        >>> convert.cstag_to_midsv(cstag)
        "=A,+T|+T|+T|=C,=C,-A,-A,=T,*AG,=T,=T"

        >>> cstag = "cs:Z:=A~ta10cg=T"
        >>> convert.cstag_to_midsv(cstag)
        "=A,=N,=N,=N,=N,=N,=N,=N,=N,=N,=N,=T"
    """
    cstag_splitted: list[str] = split_cstag(cstag)

    midsv_converted: list[str] = []

    for i, cs in enumerate(cstag_splitted):
        if len(cs) == 1:
            continue

        op: str = cs[0]

        if op == "+":
            _process_insertion(cstag_splitted, i, midsv_converted)
        elif op == "*":
            midsv_converted.append(cs)
        elif op == "~":
            _process_splice(cs, midsv_converted)
        else:
            _process_match(cs, midsv_converted)

    return ",".join(cs.upper() for cs in midsv_converted)


###########################################################
# Phred score
###########################################################


def ascii_to_phred(ascii: str) -> str:
    return str(ord(ascii) - 33)


def qual_to_qscore(qual: str, midsv_tag: str) -> str:
    """Convert ascii quality to phred score.
    To adjust the same length as midsv tags, insertion is discarded and deletion is interpolated as -1.

    Args:
        qual (str): QUAL in SAM format
        midsv_tag (str): midsv_tag

    Returns:
        str: Phred quality score with indel compensation ('-1' is assigned to the deletion loci.)
    """
    qscore = []
    idx = 0
    for tag in midsv_tag.split(","):
        if tag.startswith("-") or tag == "=N":
            qscore.append("-1")
            idx -= 1
        elif tag.startswith("+"):
            num_insertion = len(tag.split("|")) - 1
            insertion = []
            for j in range(idx, idx + num_insertion):
                insertion.append(ascii_to_phred(qual[j]) + "|")
            if tag.split("|")[-1].startswith("-") or tag.split("|")[-1] == "=N":
                insertion.append("-1")
                idx -= 1
            else:
                insertion.append(ascii_to_phred(qual[j + 1]))
            qscore.append("".join(insertion))
            idx += num_insertion
        else:
            qscore.append(ascii_to_phred(qual[idx]))
        idx += 1
    return ",".join(qscore)


###########################################################
# main
###########################################################


def convert(samdict: list[dict[str, str | int]], qscore: bool = False) -> list[dict[str, str | int]]:
    for alignment in samdict:
        alignment["MIDSV"] = cstag_to_midsv(alignment["CSTAG"])
        if qscore:
            alignment["QSCORE"] = qual_to_qscore(alignment["QUAL"], alignment["MIDSV"])
    return samdict
