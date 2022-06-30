from __future__ import annotations

import re


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


def slide_insertion(CSTAGS: list[str]) -> list[str]:
    """Append one base from the next index at an inserted base.

    Args:
        CSTAGS (list[str]): a list of CS tags

    Returns:
        list[str]: a list of slided CS tags

    Example:
        >>> mids = ['Iacgt', 'MMMM']
        >>> slide_insertion(mids)
        "['IacgtM', 'MMM']"
    """
    for i, cs in enumerate(CSTAGS):
        if "+" in cs:
            CSTAGS[i] = cs + CSTAGS[i + 1][0]
            CSTAGS[i + 1] = CSTAGS[i + 1][1:]
    return CSTAGS


# def to_mids(CSTAGS: list[str]) -> str:
#     """Convert to MIDS format

#     Args:
#         CSTAGS (list[str]): a list of CS tags

#     Returns:
#         str: a MIDS sequence

#     Example:
#         >>> to_mids(["IacgtacgtM"])
#         "8M"
#         >>> to_mids(["Dacgtacgt"])
#         "D,D,D,D,D,D,D,D"
#         >>> to_mids(["MMMM"])
#         "M,M,M,M"
#     """
#     cstags = []
#     for cs in CSTAGS:
#         if "I" == cs[0]:
#             # "IacgtacgtA" -> "8M"
#             cstags.append(f"{len(cs)-2}{cs[-1]}")
#         elif "D" == cs[0]:
#             # "Dacgtacgt" -> "D,D,D,D,D,D,D,D"
#             cstags.append(re.sub("[acgt]", "D,", cs[1:]).rstrip(","))
#         elif "M" == cs[0]:
#             # "MMM" -> "M,M,M"
#             cstags.append(cs.replace("M", "M,").rstrip(","))
#         else:
#             cstags.append(cs)
#     return ",".join(cstags)


def cstag_to_mids(cstag: str) -> str:
    """
    Input:  "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
    Output: "M,M,M,M,S,M,D,M,1M,M,M,M"
    """
    # cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
    _cstag_split = split(cstag)
    _cstags_slide = slide_insertion(_cstag_split)
    return to_mids(_cstags_slide)


def ascii_to_qscore(ascii_quality: str) -> list[int]:
    qscore = []
    for ascii in ascii_quality:
        qscore.append(ord(ascii) - 33)
    return qscore

