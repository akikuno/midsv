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
    for cs in re.split(r"(=|\*|\-|\+)", cstag):
        if not cs or cs == "cs:Z:":
            continue
        if re.match(r"(=|\*|\-|\+)", cs):
            cstag_splitted.append(cs)
        else:
            cstag_splitted[-1] += cs
    return cstag_splitted


def to_midsv(cstag_splitted: list[str]) -> list[str]:
    """Convert to MIDSV
    Args:
        cstag_splitted (list[str]): splitted CS tags
    Returns:
        list[str]: splitted MIDSV
    Example:
        >>> from midsv import convert
        >>> cstag_splitted = ['=ACGT', '*ag', '=C', '-gg', '=T', '+t', '=ACGT']
        >>> convert.to_midsv(cstag_splitted)
        "['MMMM', 'S', 'M', 'DD', 'M', 'I', 'MMMM']"
    """
    midsv = []
    for cs in cstag_splitted:
        num = len(cs) - 1
        if cs.startswith("="):
            midsv.append("M" * num)
        elif cs.startswith("+"):
            midsv.append("I" * num)
        elif cs.startswith("-"):
            midsv.append("D" * num)
        else:
            midsv.append("S")
    return midsv


def numerize_insertion(midsv: list[str]) -> list[str]:
    """Convert insertion to numeric numbers
    Example:
        >>> from midsv import convert
        >>> midsv = ['MMM', 'III', 'D', 'S']
        >>> convert.numerize_insertion(midsv)
        "['MMM', 3, 'D', 'S']"
    """
    for i, m in enumerate(midsv):
        if m.startswith("I"):
            midsv[i] = len(m)
    return midsv


def slide_insertion(midsv_numerized: list[str]) -> list[str]:
    """Append one base from the next index at an inserted base.

    Args:
        midsv (list[str]): numerized MIDSV

    Returns:
        list[str]: slided MIDSV

    Example:
        >>> from midsv import convert
        >>> midsv = [3, 'M', 4, "S", "MM"]
        >>> convert.slide_insertion(midsv)
        "['3M', '4S', "MM"]"
    """
    for i, m in enumerate(midsv_numerized):
        if not m:
            continue
        if i + 1 == len(midsv_numerized):
            midsv_numerized[i] = str(midsv_numerized[i])
            continue
        if isinstance(m, int):
            midsv_numerized[i] = str(m) + midsv_numerized[i + 1][0]
            midsv_numerized[i + 1] = midsv_numerized[i + 1][1:]
    return [m for m in midsv_numerized if m]


def to_string(midsv: list[str]) -> str:
    """Convert to string in CSV format

    Args:
        midsv (list[str]): a list of MIDSV

    Returns:
        str: MIDSV

    Example:
        >>> from midsv import convert
        >>> midsv = ['3M', '4S', "MM", "DDD"]
        >>> convert.to_string(midsv)
        "3M,4S,M,M,D,D,D"
    """
    midsv_csv = []
    for m in midsv:
        if m[0].isdigit():
            midsv_csv.append(m)
            continue
        for mm in m:
            midsv_csv.append(mm)
    return ",".join(midsv_csv)


def cstag_to_midsv(cstag: str) -> str:
    """Integrated function to convert cstag to MIDSV

    Args:
        cstag (str): a long format cstag

    Returns:
        str: MIDSV format
    """
    cstag_splited = split_cstag(cstag)
    midsv_splited = to_midsv(cstag_splited)
    midsv_numerized = numerize_insertion(midsv_splited)
    midsv_slided = slide_insertion(midsv_numerized)
    return to_string(midsv_slided)


###########################################################
# Split CS tags
###########################################################


def cstag_to_cssplit(cstag: str) -> str:
    """Generate CS SPLIT, a comma-separated nucreotide sequence corresponding to MIDSV

    Args:
        cstag (str): a long format cstag

    Returns:
        str: cs split

    Examples:
        >>> cstag = "cs:Z:=A+ttt=CC-aa=T*ag=TT"
        >>> test = convert.cstag_to_cssplit(cstag)
        "=A,+T|+T|+T|=C,=C,-A,-A,=T,*AG,=T,=T"
    """
    cstag_list = split_cstag(cstag)
    cssplit = []
    for i, cs in enumerate(cstag_list):
        op = cs[0]
        if len(cs) == 1:
            continue
        if op == "+":
            insertion = list(cs[1:])
            insertion = "+" + "|+".join(insertion)
            if i + 1 == len(cstag_list):
                cssplit.append(insertion)
                break
            next_cstag = cstag_list[i + 1]
            next_op = next_cstag[0]
            if next_op == "*":
                cssplit.append(insertion + "|" + next_cstag)
                cstag_list[i + 1] = "*"
            else:
                cstag_list[i + 1] = next_op + next_cstag[2:]
                insertion = insertion + "|" + next_cstag[:2]
                cssplit.append(insertion)
        elif op == "*":
            cssplit.append(cs)
        else:
            cs = list(cs[1:])
            cs = op + f",{op}".join(cs)
            cssplit.append(cs)
    return ",".join([cs.upper() for cs in cssplit])


###########################################################
# Phred score
###########################################################


def ascii_to_phred(ascii: str) -> str:
    return str(ord(ascii) - 33)


def qual_to_qscore(qual: str, midsv: str) -> str:
    """Convert ascii quality to phred score.
    To adjust the same length as MIDSV, insertion is discarded and deletion is interpolated as -1.

    Args:
        qual (str): QUAL in SAM format
        midsv (str): MIDSV

    Returns:
        str: Phred quality score with indel compensation
    """
    qscore = []
    idx = 0
    for m in midsv.split(","):
        if m == "D":
            qscore.append(str(-1))
            idx -= 1
        elif m[0].isdigit():
            num_insertion = int(m[:-1])
            insertion = []
            for j in range(idx, idx + num_insertion):
                insertion.append(ascii_to_phred(qual[j]) + "|")
            if m[-1] == "D":
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

