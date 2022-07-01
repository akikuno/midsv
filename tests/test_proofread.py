from pathlib import Path
from src.mids import format
from src.mids import convert
from src.mids import proofread

from importlib import reload

reload(proofread)


def test_join_inversion():
    sam = Path("tests", "data", "join", "test_inv.txt").read_text()
    sam = eval(sam)
    test = proofread.join(sam)
    answer = Path("tests", "data", "join", "answer_inv.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_join_deletion():
    sam = Path("tests", "data", "join", "test_del.txt").read_text()
    sam = eval(sam)
    test = proofread.join(sam)
    answer = Path("tests", "data", "join", "answer_del.txt").read_text()
    answer = eval(answer)
    assert test == answer
