import pytest
from src.midsv import converter

###########################################################
# MIDSV conversion
###########################################################


def test_split():
    cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
    test = converter.split_cstag(cstag)
    answer = ["=ACGT", "*ag", "=C", "-g", "=T", "+t", "=ACGT"]
    assert test == answer


@pytest.mark.parametrize(
    "cstag_splitted, i, expected",
    [
        (["+ttt", "*ag"], 0, ["+t|+t|+t|*ag"]),
        (["+a", "~ta10cg"], 0, ["+a|=N"]),
        (["+a", "=T"], 0, ["+a|=T"]),
    ],
)
def test_process_insertion(cstag_splitted, i, expected):
    cssplits = []
    converter._process_insertion(cstag_splitted, i, cssplits)
    assert cssplits == expected


@pytest.mark.parametrize(
    "cs, expected",
    [
        ("~ta10cg", ["=N"] * 10),
        ("~cg5ta", ["=N"] * 5),
    ],
)
def test_process_splice(cs, expected):
    cssplits = []
    converter._process_splice(cs, cssplits)
    assert cssplits == expected


@pytest.mark.parametrize(
    "cs, expected",
    [
        ("=ACGT", ["=A,=C,=G,=T"]),
        ("-GTC", ["-G,-T,-C"]),
    ],
)
def test_process_match(cs, expected):
    cssplits = []
    converter._process_match(cs, cssplits)
    assert cssplits == expected


@pytest.mark.parametrize(
    "cstag, expected",
    [
        ("cs:Z:=A+ttt=CC-aa=T*ag=TT", "=A,+T|+T|+T|=C,=C,-A,-A,=T,*AG,=T,=T"),
        ("cs:Z:=A~ta10cg=T", "=A,=N,=N,=N,=N,=N,=N,=N,=N,=N,=N,=T"),
        ("cs:Z:=ACGT*ag=TT", "=A,=C,=G,=T,*AG,=T,=T"),
        ("cs:Z:=A+g=TT", "=A,+G|=T,=T"),
        ("cs:Z:=A~cg3ta=TT", "=A,=N,=N,=N,=T,=T"),
    ],
)
def test_cstag_to_midsv(cstag, expected):
    assert converter.cstag_to_midsv(cstag) == expected


###########################################################
# qual_to_qscore
###########################################################


def test_qual_to_qscore_insertion():
    qual = "@!!!@@"
    cssplit = "=A,+T|+T|+T|=A,-A,=A"
    test = converter.qual_to_qscore(qual, cssplit)
    answer = "31,0|0|0|31,-1,31"
    assert test == answer


def test_qual_to_qscore_deletion():
    qual = "@0"
    cssplit = "=A,-A,-A,-A,-A,-A,=A"
    test = converter.qual_to_qscore(qual, cssplit)
    answer = "31,-1,-1,-1,-1,-1,15"
    assert test == answer


def test_qual_to_qscore_substitution():
    qual = "@012@"
    cssplit = "=A,*AG,*AG,*AG,=A"
    test = converter.qual_to_qscore(qual, cssplit)
    answer = "31,15,16,17,31"
    assert test == answer


def test_qual_to_qscore_indel():
    qual = "@0@"
    cssplit = "=A,+T|-A,-A,=A"
    test = converter.qual_to_qscore(qual, cssplit)
    answer = "31,15|-1,-1,31"
    assert test == answer


def test_qual_to_qscore_ins_sub():
    qual = "@012!"
    cssplit = "=A,+T|+T|+T|*AG"
    test = converter.qual_to_qscore(qual, cssplit)
    answer = "31,15|16|17|0"
    assert test == answer


def test_qual_to_qscore_splicing():
    qual = "!@"
    cssplit = "=A,=N,=N,=N,=N,=N,=N,=N,=N,=N,=N,=T"
    test = converter.qual_to_qscore(qual, cssplit)
    answer = "0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,31"
    assert test == answer


def test_qual_to_qscore_splicing_inversion():
    qual = "!!!@@"
    cssplit = "+C|+A|+G|=N,=N,=N,=N,=N,=N,=N,=N,=N,=N,=C,=C"
    test = converter.qual_to_qscore(qual, cssplit)
    answer = "0|0|0|-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,31,31"
    assert test == answer
