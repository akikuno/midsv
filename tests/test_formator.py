from pathlib import Path

from src.midsv import formatter, io

###########################################################
# Format headers and alignments
###########################################################


def test_extract_sqheaders():
    sampath = Path("tests", "data", "extract_headers", "query.sam")
    sam = io.read_sam(str(sampath))
    test = formatter.extract_sqheaders(sam)
    answer = {"chr13": 120421639, "chr6": 149736546}
    assert test == answer


# TODO

# def test_organize_alignments_to_dict():
#     sampath = Path("tests", "data", "dictionalize_alignments", "sub_cslong.sam")
#     sam = io.read_sam(str(sampath))
#     test = formatter.organize_alignments_to_dict(sam)
#     answer = Path("tests", "data", "dictionalize_alignments", "answer.txt").read_text()
#     answer = eval(answer)
#     assert test == answer


# def test_organize_alignments_to_dict_inversion():
#     sampath = Path("tests", "data", "inversion", "inv_cslong.sam")
#     sam = io.read_sam(str(sampath))
#     test = formatter.organize_alignments_to_dict(sam)
#     answer = Path("tests", "data", "dictionalize_alignments", "answer_inversion.txt").read_text()
#     answer = eval(answer)
#     assert test == answer


# def test_organize_alignments_to_dict_not_primary():
#     test = formatter.organize_alignments_to_dict(
#         [["not-primary", "272", "test", "1", "0", "3M", "*", "0", "0", "*", "*", "cs:Z:=AGG"]]
#     )
#     answer = []
#     assert test == answer


###########################################################
# Remove undesired reads
###########################################################


def test_remove_softclips():
    path = Path("tests", "data", "softclip", "softclip_cslong.sam")
    sam = io.read_sam(path)
    samdict = formatter.alignments_to_dict(sam)
    test = formatter.remove_softclips(samdict)
    for t in test:
        assert len(t["QUAL"]) == 100


def test_padding_n_to_sequence():
    alignment = {"POS": 0, "SEQ": "ACGT", "CIGAR": "2M1I1D1M"}
    test = formatter._padding_n_to_sequence(alignment)
    del test["POS"]
    del test["CIGAR"]
    answer = {"SEQ": "ACNT"}
    assert test == answer


def test_padding_n_to_sequence_start_5nt():
    alignment = {"POS": 5, "SEQ": "ACGT", "CIGAR": "2H2M1I1D1M"}
    test = formatter._padding_n_to_sequence(alignment)
    del test["POS"]
    del test["CIGAR"]
    answer = {"SEQ": "NNNNNACNT"}
    assert test == answer


def test_padding_n_to_sequence_hardclip():
    alignment = {"POS": 0, "SEQ": "ACGT", "CIGAR": "2H2M1I1D1M"}
    test = formatter._padding_n_to_sequence(alignment)
    del test["POS"]
    del test["CIGAR"]
    answer = {"SEQ": "ACNT"}
    assert test == answer


def test_padding_n_to_sequence_splicing():
    alignment = {"POS": 0, "SEQ": "ACGT", "CIGAR": "2M5N2M"}
    test = formatter._padding_n_to_sequence(alignment)
    del test["POS"]
    del test["CIGAR"]
    answer = {"SEQ": "ACNNNNNGT"}
    assert test == answer


def test_remove_resequence():
    path = Path("tests", "data", "overlap", "overlapped.sam")
    sam = io.read_sam(path)
    samdict = formatter.organize_alignments_to_dict(sam)
    test = formatter.remove_resequence(samdict)
    count_overlap = 0
    count_nonoverlap = 0
    for t in test:
        if t["QNAME"].startswith("overlap"):
            count_overlap += 1
        else:
            count_nonoverlap += 1
    assert count_overlap == 1 and count_nonoverlap == 2
