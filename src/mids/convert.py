from __future__ import annotations
import re

###########################################################
# MIDS conversion (from CS tag to MIDS)
###########################################################


def split(cstag: str) -> list[str]:
    """Split cstag
    Args:
        cstag (str): a string of CS tags
    Returns:
        list[str]: splitted CS tags
    Example:
        >>> from mids import convert
        >>> cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
        >>> convert.split(cstag)
        "['=ACGT', '*ag', '=C', '-g', '=T', '+t', '=ACGT']"
    """
    cstag_splitted = []
    for cs in re.split(r"(=|\*|\-|\+)", cstag):
        if not cs or cs == "cs:Z:":
            continue
        if re.match(r"(=|\*|\-|\+)", cs):
            cstag_splitted.append(cs)
        else:
            cstag_splitted[-1] += cs
    return cstag_splitted


def to_mids(cstag_splitted: list[str]) -> list[str]:
    """Convert to MIDS
    Args:
        cstag_splitted (list[str]): splitted CS tags
    Returns:
        list[str]: splitted MIDS
    Example:
        >>> from mids import convert
        >>> cstag_splitted = ['=ACGT', '*ag', '=C', '-gg', '=T', '+t', '=ACGT']
        >>> convert.to_mids(cstag_splitted)
        "['MMMM', 'S', 'M', 'DD', 'M', 'I', 'MMMM']"
    """
    mids = []
    for cs in cstag_splitted:
        num = len(cs) - 1
        if cs.startswith("="):
            mids.append("M" * num)
        elif cs.startswith("+"):
            mids.append("I" * num)
        elif cs.startswith("-"):
            mids.append("D" * num)
        else:
            mids.append("S")
    return mids


def numerize_insertion(mids: list[str]) -> list[str]:
    """Convert insertion to numeric numbers
    Example:
        >>> from mids import convert
        >>> mids = ['MMM', 'III', 'D', 'S']
        >>> convert.numerize_insertion(mids)
        "['MMM', 3, 'D', 'S']"
    """
    for i, m in enumerate(mids):
        if m.startswith("I"):
            mids[i] = len(m)
    return mids


def slide_insertion(mids: list[str]) -> list[str]:
    """Append one base from the next index at an inserted base.

    Args:
        mids (list[str]): numerized MIDS

    Returns:
        list[str]: slided MIDS

    Example:
        >>> from mids import convert
        >>> mids = [3, 'M', 4, "S", "MM"]
        >>> convert.slide_insertion(mids)
        "['3M', '4S', "MM"]"
    """
    for i, m in enumerate(mids):
        if not m:
            continue
        if isinstance(m, int):
            mids[i] = str(m) + mids[i + 1][0]
            mids[i + 1] = mids[i + 1][1:]
    return [m for m in mids if m]


def to_string(mids: list[str]) -> str:
    """Convert to string in CSV format

    Args:
        mids (list[str]): a list of MIDS

    Returns:
        str: MIDS

    Example:
        >>> from mids import convert
        >>> mids = ['3M', '4S', "MM", "DDD"]
        >>> convert.to_string(mids)
        "3M,4S,M,M,D,D,D"
    """
    mids_csv = []
    for m in mids:
        if m[0].isdigit():
            mids_csv.append(m)
            continue
        for mm in m:
            mids_csv.append(mm)
    return ",".join(mids_csv)


def cstag_to_mids(cstag: str) -> str:
    cstag_splited = split(cstag)
    mids_splited = to_mids(cstag_splited)
    mids_numerized = numerize_insertion(mids_splited)
    mids_slided = slide_insertion(mids_numerized)
    return to_string(mids_slided)


###########################################################
# Phred score
###########################################################


# def ascii_to_qscore(ascii_quality: str) -> str:
#     qscore = []
#     for ascii in ascii_quality:
#         qscore.append(str(ord(ascii) - 33))
#     return ",".join(qscore)


def ascii_to_phred(ascii: str) -> str:
    return str(ord(ascii) - 33)


def to_qscore_with_indel_compensation(qual: str, mids: str) -> str:
    """Convert ascii quality to phred score.
    To adjust the same length as MIDS, insertion is discarded and deletion is interpolated as -1.

    Args:
        qual (str): QUAL in SAM format
        mids (str): MIDS

    Returns:
        str: Phred quality score with indel compensation
    """
    qscore = []
    idx = 0
    for m in mids.split(","):
        if m == "D":
            qscore.append(str(-1))
            idx -= 1
        elif m[0].isdigit():
            idx += int(m[:-1])
            qscore.append(ascii_to_phred(qual[idx]))
        else:
            qscore.append(ascii_to_phred(qual[idx]))
        idx += 1
    return ",".join(qscore)