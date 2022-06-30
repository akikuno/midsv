from pathlib import Path
from src.mids import convert


# substitution = Path("tests", "data", "substitution", "sub_cslong.sam").read_text().strip().split("\n")
# substitution = [s.split("\t") for s in substitution]

# headers = extract.headers(substitution)
# alignments = extract.alignments(substitution)
# alignments = format.append_reference_sequence_length(headers, alignments)

# mids_list = []
# for alignment in alignments:
#     mids_list.append(convert.cstag_to_mids(alignment["CSTAG"]))

# alignments_with_mids = []
# for alignment, mids in zip(alignments, mids_list):
#     alignment["MIDS"] = mids
#     alignments_with_mids.append(alignment)


def test_split():
    cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
    test = convert.split(cstag)
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


def test_to_string():
    mids = ["3M", "4S", "MM", "DDD"]
    test = convert.to_string(mids)
    answer = "3M,4S,M,M,D,D,D"
    assert test == answer


def test_ascii_to_qscore():
    ascii = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJK"
    qscore = convert.ascii_to_qscore(ascii)
    for i in range(43):
        assert qscore[i] == i
