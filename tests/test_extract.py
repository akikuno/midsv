from pathlib import Path
from src.mids import extract


def test_extract_headers():
    query = Path("tests", "data","extract_headers", "query.sam").read_text().strip().split("\n")
    query = [q.split("\t") for q in query]
    dict_eval = {
        "chr13": 120421639,
        "chr6": 149736546
    }
    assert extract.headers(query) == dict_eval

def test_extract_alignments():
    query = Path("tests", "data","extract_alignments", "sub_cslong.sam").read_text().strip().split("\n")
    query = [q.split("\t") for q in query]
    list_eval = [
        {'QNAME': 'control', 'RNAME': 'random_100bp', 'POS': 1, 'QUOL': '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', 'CSTAG': 'cs:Z:=ACTGTGCGGCATACTTAATTATACATTTGAAACGCGCCCAAGTGACGCTAGGCAAGTCAGAGCAGGTTCCCGTGTTAGCTTAAGGGTAAACATACAAGTC'}, {'QNAME': 'sub-1nt', 'RNAME': 'random_100bp', 'POS': 1, 'QUOL': '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', 'CSTAG': 'cs:Z:=ACTGTGCGGCATACTTAATTATACATTTGAAACGCGCCCAAGTGACGCT*ag=GGCAAGTCAGAGCAGGTTCCCGTGTTAGCTTAAGGGTAAACATACAAGTC'}
    ]
    assert extract.alignments(query) == list_eval