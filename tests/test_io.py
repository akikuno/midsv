from pathlib import Path

import pytest

from src.midsv import io

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


def test_write_vcf(tmp_path):
    alignments = [
        {"QNAME": "large-deletion", "RNAME": "example", "MIDSV": "=A,=C,=N,=N,=N,=N,=N,=N,=G,=T"},
        {"QNAME": "indel_sub", "RNAME": "example", "MIDSV": "=A,=C,=G,=T,*AG,+T|+T|+T|=C,-A,-A,=G,=T"},
        {"QNAME": "inversion", "RNAME": "example", "MIDSV": "=A,=C,=G,=T,=A,=c,=g,=t,=A,=C"},
        {"QNAME": "long-ins", "RNAME": "longins", "MIDSV": "+AAAAAA|=G,=T"},
    ]
    output_path = Path(tmp_path, "variants.vcf")
    io.write_vcf(alignments, output_path, large_sv_threshold=5)
    content = output_path.read_text().strip().split("\n")
    expected = [
        "##fileformat=VCFv4.3",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "example\t3\t.\tN\t<DEL>\t.\tPASS\tTYPE=DEL;SVLEN=-6;SEQ=NNNNNN;QNAME=large-deletion",
        "example\t5\t.\tA\tG\t.\tPASS\tTYPE=SUB;QNAME=indel_sub",
        "example\t6\t.\tC\tCTTT\t.\tPASS\tTYPE=INS;SVLEN=3;SEQ=TTT;QNAME=indel_sub",
        "example\t6\t.\tC\t<INV>\t.\tPASS\tSVTYPE=INV;SVLEN=3;SEQ=CGT;QNAME=inversion",
        "example\t7\t.\tA\t<DEL>\t.\tPASS\tTYPE=DEL;SVLEN=-2;SEQ=AA;QNAME=indel_sub",
        "longins\t1\t.\tG\t<INS>\t.\tPASS\tTYPE=INS;SVLEN=6;SEQ=AAAAAA;QNAME=long-ins",
    ]
    assert content == expected
