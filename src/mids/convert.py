import re

def slide_insertion(CSTAGS: [str]) -> list[str]:
    """To append one base from the next index at an inserted base.

    Args:
        alignments (list[str]): _description_

    Returns:
        list[str]: _description_
    Example:
        >>> _input = ['Iacgt', 'MMMM']
        >>> slide_insertion(_input)
        "['IacgtM', 'MMM']"
    """
    for i, cs in enumerate(CSTAGS):
        if "I" in cs:
            CSTAGS[i] = cs + CSTAGS[i + 1][0]
            CSTAGS[i + 1] = CSTAGS[i + 1][1:]
    return CSTAGS

def to_csv_with_fixed_length(CSTAG: str) -> str:
    if "I" == CSTAG[0]:
        # "IacgtacgtA" -> "8M"
        return f"{len(CSTAG)-2}{CSTAG[-1]}"
    elif "D" == CSTAG[0]:
        # "Dacgtacgt" -> "D,D,D,D,D,D,D,D"
        return re.sub("[acgt]", "D,", CSTAG[1:]).rstrip(",")
    elif "M" == CSTAG[0]:
        # "MMM" -> "M,M,M"
        return CSTAG.replace("M", "M,").rstrip(",")
    else:
        return CSTAG


def cstag_to_mids(CSTAG: str) -> str:
    """
    Input:  "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
    Output: "M,M,M,M,S,M,D,M,1M,M,M,M"
    """
    # CSTAG = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
    _cstag = CSTAG.replace("cs:Z:", "")
    _cstag = _cstag.lstrip("=")
    _cstag = _cstag.replace("-", "=D")
    _cstag = _cstag.replace("+", "=I")
    _cstag = re.sub("[ACGT]", "M", _cstag)
    _cstag = re.sub(r"\*[acgt][acgt]", "=S", _cstag)
    _cstag_split = _cstag.split("=")
    _cstags_slide = slide_insertion(_cstag_split)
    _mids = [to_csv_with_fixed_length(cs) for cs in _cstags_slide]
    return ",".join(_mids)