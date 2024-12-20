import pytest
from pathlib import Path
from src.midsv import io
from importlib import reload

reload(io)

###########################################################
# Read sam
###########################################################


def test_check_read_sam_str():
    path = str(Path("tests", "data", "subindel", "subindel_cslong.sam"))
    test = io.read_sam(path)
    test = list(test)
    answer = eval(Path("tests", "data", "read_sam", "answer.txt").read_text())
    assert test == answer


def test_check_read_sam_Path():
    path = Path("tests", "data", "subindel", "subindel_cslong.sam")
    test = io.read_sam(path)
    test = list(test)
    answer = eval(Path("tests", "data", "read_sam", "answer.txt").read_text())
    assert test == answer


def test_check_read_sam_TypeError():
    with pytest.raises(TypeError):
        assert io.read_sam(1)


def test_check_read_sam_FileNotFoundError():
    with pytest.raises(FileNotFoundError):
        assert io.read_sam("hoge")

###########################################################
# Read / Write jsonl
###########################################################


def test_read_jsonl():
    path_jsonl = Path("tests", "data", "read_jsonl", "test.jsonl")
    test = io.read_jsonl(path_jsonl)
    test = list(test)
    answer = [{"hoge": 1, "fuga": 2}, {"foo": "3", "bar": "4"}]
    assert test == answer


def test_write_jsonl(tmp_path):
    dicts = [{"hoge": 1, "fuga": 2}, {"foo": "3", "bar": "4"}]
    output_path = Path(tmp_path, "tmp.jsonl")
    io.write_jsonl(dicts, output_path)
    assert output_path.read_text() == '{"hoge": 1, "fuga": 2}\n{"foo": "3", "bar": "4"}\n'
