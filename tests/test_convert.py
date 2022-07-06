from pathlib import Path
from src.mids import convert, format
from importlib import reload

reload(convert)

###########################################################
# MIDS conversion
###########################################################


def test_split():
    cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
    test = convert.split_cstag(cstag)
    answer = ["=ACGT", "*ag", "=C", "-g", "=T", "+t", "=ACGT"]
    assert test == answer


def test_to_mids():
    cstag_splitted = ["=ACGT", "*ag", "=C", "-gg", "=T", "+t", "=ACGT"]
    test = convert.to_mids(cstag_splitted)
    answer = ["MMMM", "S", "M", "DD", "M", "I", "MMMM"]
    assert test == answer


def test_numerize_insertion():
    mids = ["MMM", "III", "D", "S"]
    test = convert.numerize_insertion(mids)
    answer = ["MMM", 3, "D", "S"]
    assert test == answer


def test_slide_insertion():
    mids = [3, "M", 4, "S", "MM"]
    test = convert.slide_insertion(mids)
    answer = ["3M", "4S", "MM"]
    assert test == answer


def test_cstag_to_mids():
    cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
    test = convert.cstag_to_mids(cstag)
    answer = "M,M,M,M,S,M,D,M,1M,M,M,M"
    assert test == answer


def test_cstag_to_mids_insertion_at_last():
    cstag = "cs:Z:+tt*ag-aa+tt"
    test = convert.cstag_to_mids(cstag)
    answer = "2S,D,D,2"
    assert test == answer


def test_to_string():
    mids = ["3M", "4S", "MM", "DDD"]
    test = convert.to_string(mids)
    answer = "3M,4S,M,M,D,D,D"
    assert test == answer


###########################################################
# Split CS tags
###########################################################


def test_cstag_to_cssplit_repeat_match():
    cstag = "cs:Z:=ACGTACGT"
    test = convert.cstag_to_cssplit(cstag)
    answer = "=A,=C,=G,=T,=A,=C,=G,=T"
    assert test == answer


def test_cstag_to_cssplit():
    cstag = "cs:Z:=A+ttt=CC-aa=T*ag=TT"
    test = convert.cstag_to_cssplit(cstag)
    answer = "=A,+T|+T|+T|=C,=C,-A,-A,=T,*AG,=T,=T"
    assert test == answer


def test_cstag_to_cssplit_insertion_substitution():
    cstag = "cs:Z:=A+tt*ag-aa"
    test = convert.cstag_to_cssplit(cstag)
    answer = "=A,+T|+T|*AG,-A,-A"
    assert test == answer


def test_cstag_to_cssplit_insertion_at_last():
    cstag = "cs:Z:+tt*ag-aa+tt"
    test = convert.cstag_to_cssplit(cstag)
    answer = "+T|+T|*AG,-A,-A,+T|+T"
    assert test == answer


###########################################################
# Phred score
###########################################################


def test_cstag_to_qscore_insertion():
    qual = "@!!!@@"
    mids = "M,3M,D,M"
    test = convert.qual_to_qscore(qual, mids)
    answer = "31,0|0|0|31,-1,31"
    assert test == answer


def test_cstag_to_qscore_deletion():
    qual = "@0"
    mids = "M,D,D,D,D,D,M"
    test = convert.qual_to_qscore(qual, mids)
    answer = "31,-1,-1,-1,-1,-1,15"
    assert test == answer


def test_cstag_to_qscore_substitution():
    qual = "@012@"
    mids = "M,S,S,S,M"
    test = convert.qual_to_qscore(qual, mids)
    answer = "31,15,16,17,31"
    assert test == answer


def test_qual_to_qscore_real():
    sampath = Path("tests", "data", "phredscore", "subindel_cslong.sam")
    sam = format.read_sam(str(sampath))
    samdict = format.dictionarize_sam(sam)
    for i, alignment in enumerate(samdict):
        samdict[i]["MIDS"] = convert.cstag_to_mids(alignment["CSTAG"])
        samdict[i]["QSCORE"] = convert.qual_to_qscore(alignment["QUAL"], alignment["MIDS"])
    test = samdict
    answer = Path("tests", "data", "phredscore", "answer.txt").read_text()
    answer = eval(answer)
    assert test == answer

