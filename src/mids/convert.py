from __future__ import annotations

import re

def split(cstag: str) -> list[str]:
    """Split cstag along with annotating MIDS

    Args:
        CSTAGS (list[str]): a list of CS tags

    Returns:
        list[str]: a list of slided CS tags

    Example:
        >>> cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
        >>> split(cstag)
        "['MMMM', 'S', 'M', 'Dg', 'M', 'It', 'MMMM']"
    """
    _cstag = cstag.replace("cs:Z:", "")
    _cstag = _cstag.lstrip("=")
    _cstag = _cstag.replace("-", "=D")
    _cstag = _cstag.replace("+", "=I")
    _cstag = re.sub("[ACGT]", "M", _cstag)
    _cstag = re.sub(r"\*[acgt][acgt]", "=S", _cstag)
    return _cstag.split("=")

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
