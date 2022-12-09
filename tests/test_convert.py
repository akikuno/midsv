from pathlib import Path
from src.midsv import io, convert, format
from importlib import reload

reload(convert)

###########################################################
# MIDSVV conversion
###########################################################


def test_split():
    cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
    test = convert.split_cstag(cstag)
    answer = ["=ACGT", "*ag", "=C", "-g", "=T", "+t", "=ACGT"]
    assert test == answer


def test_to_midsv():
    cstag_splitted = ["=ACGT", "*ag", "=C", "-gg", "=T", "+t", "=ACGT"]
    test = convert.to_midsv(cstag_splitted)
    answer = ["MMMM", "S", "M", "DD", "M", "I", "MMMM"]
    assert test == answer


def test_to_midsv_splice():
    cstag_splitted = ['=ACGT', '~tg5ca', '=ACGT']
    test = convert.to_midsv(cstag_splitted)
    answer = ['MMMM', 'DDDDD', 'MMMM']
    assert test == answer


def test_numerize_insertion():
    midsv = ["MMM", "III", "D", "S"]
    test = convert.numerize_insertion(midsv)
    answer = ["MMM", 3, "D", "S"]
    assert test == answer


def test_slide_insertion():
    midsv = [3, "M", 4, "S", "MM"]
    test = convert.slide_insertion(midsv)
    answer = ["3M", "4S", "MM"]
    assert test == answer


def test_cstag_to_midsv():
    cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
    test = convert.cstag_to_midsv(cstag)
    answer = "M,M,M,M,S,M,D,M,1M,M,M,M"
    assert test == answer


def test_cstag_to_midsv_insertion_at_last():
    cstag = "cs:Z:+tt*ag-aa+tt"
    test = convert.cstag_to_midsv(cstag)
    answer = "2S,D,D,2"
    assert test == answer


def test_to_string():
    midsv = ["3M", "4S", "MM", "DDD"]
    test = convert.to_string(midsv)
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


def test_cstag_to_cssplit_splice():
    cstag = "cs:Z:=A~ta10cg=T"
    test = convert.cstag_to_cssplit(cstag)
    answer = "=A,N,N,N,N,N,N,N,N,N,N,=T"
    assert test == answer


###########################################################
# Phred score
###########################################################

# ----------------------------------------------------------
# qual_to_qscore_midsv
# ----------------------------------------------------------
def test_qual_to_qscore_midsv_insertion():
    qual = "@!!!@@"
    midsv = "M,3M,D,M"
    test = convert.qual_to_qscore_midsv(qual, midsv)
    answer = "31,0|0|0|31,-1,31"
    assert test == answer


def test_qual_to_qscore_midsv_deletion():
    qual = "@0"
    midsv = "M,D,D,D,D,D,M"
    test = convert.qual_to_qscore_midsv(qual, midsv)
    answer = "31,-1,-1,-1,-1,-1,15"
    assert test == answer


def test_qual_to_qscore_midsv_substitution():
    qual = "@012@"
    midsv = "M,S,S,S,M"
    test = convert.qual_to_qscore_midsv(qual, midsv)
    answer = "31,15,16,17,31"
    assert test == answer


def test_qual_to_qscore_midsv_indel():
    qual = "@0@"
    midsv = "M,1D,D,M"
    test = convert.qual_to_qscore_midsv(qual, midsv)
    answer = "31,15|-1,-1,31"
    assert test == answer


def test_qual_to_qscore_midsv_ins_sub():
    qual = "@012!"
    midsv = "M,3S"
    test = convert.qual_to_qscore_midsv(qual, midsv)
    answer = "31,15|16|17|0"
    assert test == answer


def test_qual_to_qscore_midsv_real():
    sampath = Path("tests", "data", "phredscore", "subindel_cslong.sam")
    sam = io.read_sam(str(sampath))
    samdict = format.dictionarize_sam(sam)
    for i, alignment in enumerate(samdict):
        samdict[i]["MIDSV"] = convert.cstag_to_midsv(alignment["CSTAG"])
        samdict[i]["QSCORE"] = convert.qual_to_qscore_midsv(alignment["QUAL"], alignment["MIDSV"])
    test = samdict
    answer = Path("tests", "data", "phredscore", "answer.txt").read_text()
    answer = eval(answer)
    assert test == answer

# ----------------------------------------------------------
# qual_to_qscore_cssplit
# ----------------------------------------------------------

def test_qual_to_qscore_cssplit_insertion():
    qual = "@!!!@@"
    cssplit = "=A,+T|+T|+T|=A,-A,=A"
    test = convert.qual_to_qscore_cssplit(qual, cssplit)
    answer = "31,0|0|0|31,-1,31"
    assert test == answer


def test_qual_to_qscore_cssplit_deletion():
    qual = "@0"
    cssplit = "=A,-A,-A,-A,-A,-A,=A"
    test = convert.qual_to_qscore_cssplit(qual, cssplit)
    answer = "31,-1,-1,-1,-1,-1,15"
    assert test == answer


def test_qual_to_qscore_cssplit_substitution():
    qual = "@012@"
    cssplit = "=A,*AG,*AG,*AG,=A"
    test = convert.qual_to_qscore_cssplit(qual, cssplit)
    answer = "31,15,16,17,31"
    assert test == answer


def test_qual_to_qscore_cssplit_indel():
    qual = "@0@"
    cssplit = "=A,+T|-A,-A,=A"
    test = convert.qual_to_qscore_cssplit(qual, cssplit)
    answer = "31,15|-1,-1,31"
    assert test == answer


def test_qual_to_qscore_cssplit_ins_sub():
    qual = "@012!"
    cssplit = "=A,+T|+T|+T|*AG"
    test = convert.qual_to_qscore_cssplit(qual, cssplit)
    answer = "31,15|16|17|0"
    assert test == answer

def test_qual_to_qscore_cssplit_splicing():
    qual = "!@"
    cssplit = '=A,N,N,N,N,N,N,N,N,N,N,=T'
    test = convert.qual_to_qscore_cssplit(qual, cssplit)
    answer = '0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,31'
    assert test == answer

def test_qual_to_qscore_cssplit_splicing_inversion():
    qual = "!!!@@"
    cssplit = '+C|+A|+G|N,N,N,N,N,N,N,N,N,N,=C,=C'
    test = convert.qual_to_qscore_cssplit(qual, cssplit)
    answer = '0|0|0|-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,31,31'
    assert test == answer
