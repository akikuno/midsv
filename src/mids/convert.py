from __future__ import annotations
import re

###########################################################
# MIDS conversion (from CS tag to MIDS)
###########################################################


def split_cstag(cstag: str) -> list[str]:
    """Split cstag
    Args:
        cstag (str): a string of CS tags
    Returns:
        list[str]: splitted CS tags
    Example:
        >>> from mids import convert
        >>> cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
        >>> convert.split_cstag(cstag)
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


def slide_insertion(mids_numerized: list[str]) -> list[str]:
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
    for i, m in enumerate(mids_numerized):
        if not m:
            continue
        if i + 1 == len(mids_numerized):
            mids_numerized[i] = str(mids_numerized[i])
            continue
        if isinstance(m, int):
            mids_numerized[i] = str(m) + mids_numerized[i + 1][0]
            mids_numerized[i + 1] = mids_numerized[i + 1][1:]
    return [m for m in mids_numerized if m]


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
    """Integrated function to convert cstag to MIDS

    Args:
        cstag (str): a long format cstag

    Returns:
        str: MIDS format
    """
    cstag_splited = split_cstag(cstag)
    mids_splited = to_mids(cstag_splited)
    mids_numerized = numerize_insertion(mids_splited)
    mids_slided = slide_insertion(mids_numerized)
    return to_string(mids_slided)


###########################################################
# Split CS tags
###########################################################


def cstag_to_cssplit(cstag: str) -> str:
    """Generate CS SPLIT, a comma-separated nucreotide sequence corresponding to MIDS

    Args:
        cstag (str): a long format cstag

    Returns:
        str: cs split
    """
    cstag_list = split_cstag(cstag)
    cssplit = []
    for i, cs in enumerate(cstag_list):
        if len(cs) == 1:
            continue
        if cs[0] == "+":
            insertion = list(cs[1:])
            if i + 1 == len(cstag_list):
                cssplit.append("".join(insertion))
                break
            next_cstag = cstag_list[i + 1]
            next_op = next_cstag[0]
            if next_op == "*":
                next_cs = next_cstag[-1]
                cstag_list[i + 1] = next_op
            else:
                next_cs = next_cstag[1]
                cstag_list[i + 1] = next_op + next_cstag[2:]
            insertion.append(next_cs)
            cssplit.append("".join(insertion))
        elif cs[0] == "*":
            cssplit.append(cs[-1])
        else:
            cssplit.append(",".join(cs[1:]))
    return ",".join(cssplit)


###########################################################
# Phred score
###########################################################


def ascii_to_phred(ascii: str) -> str:
    return str(ord(ascii) - 33)


def qual_to_qscore(qual: str, mids: str) -> str:
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
            num_insertion = int(m[:-1])
            insertion = []
            for j in range(idx, idx + num_insertion):
                insertion.append(ascii_to_phred(qual[j]) + "|")
            insertion.append(ascii_to_phred(qual[j + 1]))
            qscore.append("".join(insertion))
            idx += num_insertion
        else:
            qscore.append(ascii_to_phred(qual[idx]))
        idx += 1
    return ",".join(qscore)
