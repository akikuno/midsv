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


def to_midsv(cstag_splitted: list[str]) -> list[str]:
    """Convert to MIDSV
    Args:
        cstag_splitted (list[str]): splitted CS tags
    Returns:
        list[str]: splitted MIDSV
    Examples:
        >>> from midsv import convert
        >>> cstag_splitted = ['=ACGT', '*ag', '=C', '-gg', '=T', '+t', '=ACGT']
        >>> convert.to_midsv(cstag_splitted)
        "['MMMM', 'S', 'M', 'DD', 'M', 'I', 'MMMM']"

        >>> cstag_splitted = ['=ACGT', '~tg5ca', '=ACGT']
        >>> convert.to_midsv(cstag_splitted)
        "['MMMM', 'DDDDD', 'MMMM']"
    """
    midsv_splitted = []
    for cs in cstag_splitted:
        num = len(cs) - 1
        if cs.startswith("="):
            midsv_splitted.append("M" * num)
        elif cs.startswith("+"):
            midsv_splitted.append("I" * num)
        elif cs.startswith("-"):
            midsv_splitted.append("D" * num)
        elif cs.startswith("~"):
            match = re.match(r"([a-z]+)([0-9]+)([a-z]+)", cs[1:])
            _, splice, _ = match.groups()
            midsv_splitted.append("D" * (int(splice)))
        else:
            midsv_splitted.append("S")
    return midsv_splitted

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
        >>> convert.cstag_to_cssplit(cstag)
        "=A,+T|+T|+T|=C,=C,-A,-A,=T,*AG,=T,=T"

        >>> cstag = "cs:Z:=A~ta10cg=T"
        >>> convert.cstag_to_cssplit(cstag)
        "=A,N,N,N,N,N,N,N,N,N,N,=T"
    """
    cstag_splitted = split_cstag(cstag)
    cssplits = []
    for i, cs in enumerate(cstag_splitted):
        if len(cs) == 1:
            continue
        op = cs[0]
        if op == "+":
            insertion = list(cs[1:])
            insertion = "+" + "|+".join(insertion)
            if i + 1 == len(cstag_splitted):
                cssplits.append(insertion)
                break
            next_cstag = cstag_splitted[i + 1]
            next_op = next_cstag[0]
            if next_op == "*":
                cssplits.append(insertion + "|" + next_cstag)
                cstag_splitted[i + 1] = "*"
            elif next_op == "~":
                insertion = insertion + "|" + "N"
                cssplits.append(insertion)
                match = re.match(r"([a-z]+)([0-9]+)([a-z]+)", next_cstag[1:])
                left, splice, right = match.groups()
                cstag_splitted[i + 1] = f"{next_op}{left}{int(splice)-1}{right}"
            else:
                insertion = insertion + "|" + next_cstag[:2]
                cssplits.append(insertion)
                cstag_splitted[i + 1] = next_op + next_cstag[2:]
        elif op == "*":
            cssplits.append(cs)
        elif op == "~":
            match = re.match(r"([a-z]+)([0-9]+)([a-z]+)", cs[1:])
            _, splice, _ = match.groups()
            cssplits.extend(["N"] * int(splice))
        else:
            cs = list(cs[1:])
            cs = op + f",{op}".join(cs)
            cssplits.append(cs)
    return ",".join([cs.upper() for cs in cssplits])


###########################################################
# Phred score
###########################################################


def ascii_to_phred(ascii: str) -> str:
    return str(ord(ascii) - 33)


def qual_to_qscore_midsv(qual: str, midsv: str) -> str:
    """Convert ascii quality to phred score.
    To adjust the same length as MIDSV, insertion is discarded and deletion is interpolated as -1.

    Args:
        qual (str): QUAL in SAM format
        midsv (str): MIDSV

    Returns:
        str: Phred quality score with indel compensation ('-1' is assigned to the deletion loci.)
    """
    qscore = []
    idx = 0
    for m in midsv.split(","):
        if m == "D":
            qscore.append("-1")
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


def qual_to_qscore_cssplit(qual: str, cssplit: str) -> str:
    """Convert ascii quality to phred score.
    To adjust the same length as cssplit, insertion is discarded and deletion is interpolated as -1.

    Args:
        qual (str): QUAL in SAM format
        cssplit (str): CSSPLIT

    Returns:
        str: Phred quality score with indel compensation ('-1' is assigned to the deletion loci.)
    """
    qscore = []
    idx = 0
    for cs in cssplit.split(","):
        if cs.startswith("-") or cs == "N": # N = spliced loci
            qscore.append("-1")
            idx -= 1
        elif cs.startswith("+"):
            num_insertion = len(cs.split("|")) - 1
            insertion = []
            for j in range(idx, idx + num_insertion):
                insertion.append(ascii_to_phred(qual[j]) + "|")
            if cs.split("|")[-1].startswith("-") or cs.split("|")[-1] == "N":
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
