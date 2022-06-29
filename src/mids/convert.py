from __future__ import annotations

import re


def split(cstag: str) -> list[str]:
    """Split cstag

    Args:
        sctag (str): a string of CS tags

    Returns:
        list[str]: a list of splitted CS tags

    Example:
        >>> cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
        >>> split(cstag)
        "['MMMM', 'S', 'M', 'Dg', 'M', 'It', 'MMMM']"
        =ACGT *ag =C -g =T +t =ACGT
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
        if "I" in cs:
            CSTAGS[i] = cs + CSTAGS[i + 1][0]
            CSTAGS[i + 1] = CSTAGS[i + 1][1:]
    return CSTAGS


def to_mids(CSTAGS: list[str]) -> str:
    """Convert to MIDS format

    Args:
        CSTAGS (list[str]): a list of CS tags

    Returns:
        str: a MIDS sequence

    Example:
        >>> to_mids(["IacgtacgtM"])
        "8M"
        >>> to_mids(["Dacgtacgt"])
        "D,D,D,D,D,D,D,D"
        >>> to_mids(["MMMM"])
        "M,M,M,M"
    """
    cstags = []
    for cs in CSTAGS:
        if "I" == cs[0]:
            # "IacgtacgtA" -> "8M"
            cstags.append(f"{len(cs)-2}{cs[-1]}")
        elif "D" == cs[0]:
            # "Dacgtacgt" -> "D,D,D,D,D,D,D,D"
            cstags.append(re.sub("[acgt]", "D,", cs[1:]).rstrip(","))
        elif "M" == cs[0]:
            # "MMM" -> "M,M,M"
            cstags.append(cs.replace("M", "M,").rstrip(","))
        else:
            cstags.append(cs)
    return ",".join(cstags)


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

