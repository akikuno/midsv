from pathlib import Path
from src.mids import extract_name_length


def test_extract_name_length():
    query = Path("tests", "data","multiple_snln", "query.sam").read_text().strip().split("\n")
    query = [q.split("\t") for q in query]
    dict_eval = {
        "chr13": 120421639,
        "chr6": 149736546
    }
    assert extract_name_length(query) == dict_eval