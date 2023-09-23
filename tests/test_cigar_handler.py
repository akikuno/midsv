from src.midsv.utils.cigar_handler import split_cigar


def test_split_cigar():
    assert split_cigar("10M5D5M") == ["10M", "5D", "5M"]
    assert split_cigar("5M10N5M") == ["5M", "10N", "5M"]
    assert split_cigar("") == []
