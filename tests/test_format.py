from pathlib import Path
from src.mids import extract
from src.mids import format

def test_append_reference_sequence_length():
    query = Path("tests", "data","extract_alignments", "sub_cslong.sam").read_text().strip().split("\n")
    query = [q.split("\t") for q in query]
    headers = extract.headers(query)
    alignments = extract.alignments(query)
    answer = format.append_reference_sequence_length(headers, alignments)
    correct_answer = [
        {'QNAME': 'control', 'RNAME': 'random_100bp', 'POS': 1, 'QUOL': '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', 'CSTAG': 'cs:Z:=ACTGTGCGGCATACTTAATTATACATTTGAAACGCGCCCAAGTGACGCTAGGCAAGTCAGAGCAGGTTCCCGTGTTAGCTTAAGGGTAAACATACAAGTC', 'RLEN': 100}, {'QNAME': 'sub-1nt', 'RNAME': 'random_100bp', 'POS': 1, 'QUOL': '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', 'CSTAG': 'cs:Z:=ACTGTGCGGCATACTTAATTATACATTTGAAACGCGCCCAAGTGACGCT*ag=GGCAAGTCAGAGCAGGTTCCCGTGTTAGCTTAAGGGTAAACATACAAGTC', 'RLEN': 100}
        ]
    assert answer == correct_answer