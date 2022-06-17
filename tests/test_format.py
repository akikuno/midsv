from pathlib import Path
from src.mids import format


def test_extract_headers():
    query = Path("tests", "data", "extract_headers", "query.sam").read_text().strip().split("\n")
    query = [q.split("\t") for q in query]
    value = format.extract_headers(query)
    answer = {"chr13": 120421639, "chr6": 149736546}
    assert value == answer


def test_extract_alignments():
    query = Path("tests", "data", "extract_alignments", "sub_cslong.sam").read_text().strip().split("\n")
    query = [q.split("\t") for q in query]
    list_eval = [
        {
            "QNAME": "control",
            "RNAME": "random_100bp",
            "POS": 1,
            "QUOL": "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
            "CSTAG": "cs:Z:=ACTGTGCGGCATACTTAATTATACATTTGAAACGCGCCCAAGTGACGCTAGGCAAGTCAGAGCAGGTTCCCGTGTTAGCTTAAGGGTAAACATACAAGTC",
        },
        {
            "QNAME": "sub-1nt",
            "RNAME": "random_100bp",
            "POS": 1,
            "QUOL": "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
            "CSTAG": "cs:Z:=ACTGTGCGGCATACTTAATTATACATTTGAAACGCGCCCAAGTGACGCT*ag=GGCAAGTCAGAGCAGGTTCCCGTGTTAGCTTAAGGGTAAACATACAAGTC",
        },
    ]
    assert format.extract_alignments(query) == list_eval


def test_append_reference_sequence_length():
    query = Path("tests", "data", "extract_alignments", "sub_cslong.sam").read_text().strip().split("\n")
    query = [q.split("\t") for q in query]
    answer = format.append_reference_sequence_length(query)
    correct_answer = [
        {
            "QNAME": "control",
            "RNAME": "random_100bp",
            "POS": 1,
            "QUOL": "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
            "CSTAG": "cs:Z:=ACTGTGCGGCATACTTAATTATACATTTGAAACGCGCCCAAGTGACGCTAGGCAAGTCAGAGCAGGTTCCCGTGTTAGCTTAAGGGTAAACATACAAGTC",
            "RLEN": 100,
        },
        {
            "QNAME": "sub-1nt",
            "RNAME": "random_100bp",
            "POS": 1,
            "QUOL": "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
            "CSTAG": "cs:Z:=ACTGTGCGGCATACTTAATTATACATTTGAAACGCGCCCAAGTGACGCT*ag=GGCAAGTCAGAGCAGGTTCCCGTGTTAGCTTAAGGGTAAACATACAAGTC",
            "RLEN": 100,
        },
    ]
    assert answer == correct_answer
