import pytest

from src.midsv import polisher

###########################################################
# merge
###########################################################


@pytest.mark.parametrize(
    "current_alignment, first_read_is_forward, expected_mid_sv",
    [
        ({"FLAG": 0, "MIDSV": "=A,=C,=G,=T"}, True, "=A,=C,=G,=T"),
        ({"FLAG": 16, "MIDSV": "=A,=C,=G,=T"}, True, "=a,=c,=g,=t"),
        ({"FLAG": 16, "MIDSV": "=A,=C,=G,=T"}, False, "=A,=C,=G,=T"),
        ({"FLAG": 0, "MIDSV": "=A,=C,=G,=T"}, False, "=a,=c,=g,=t"),
    ],
)
def test_process_inversion(current_alignment, first_read_is_forward, expected_mid_sv):
    polisher.process_inversion(current_alignment, first_read_is_forward)
    assert current_alignment["MIDSV"] == expected_mid_sv


@pytest.mark.parametrize(
    "previous_alignment, current_alignment, expected",
    [
        (
            {"MIDSV": "=A,=C,=G,=T", "QSCORE": "10,20,30,40"},
            {"MIDSV": "=G,=T,=A,=C", "QSCORE": "30,40,50,60"},
            2,  # Microhomology length for "=G,=T"
        ),
        (
            {"MIDSV": "=A,=C,=G,=T", "QSCORE": "10,20,30,40"},
            {"MIDSV": "=C,=A,=C,=G", "QSCORE": "40,50,60,70"},
            0,  # No microhomology
        ),
        (
            {"MIDSV": "=A,=C,=G,=T", "QSCORE": "10,20,30,40"},
            {"MIDSV": "=A,=C,=G,=T", "QSCORE": "10,20,30,40"},
            4,  # Full match
        ),
        (
            {"MIDSV": "=A,=C,=G,=T", "QSCORE": "10,20,30,40"},
            {"MIDSV": "=C,=G,=T,=A", "QSCORE": "20,30,40,50"},
            3,  # Microhomology length for "=C,=G,=T"
        ),
        (
            {"MIDSV": "=A,=C,=G,=T", "QSCORE": "10,20,30,40"},
            {"MIDSV": "=C,=G,=T,=A", "QSCORE": "25,30,40,50"},  # Slight QScore mismatch
            0,  # Microhomology stops at QScore mismatch
        ),
        (
            {"MIDSV": "=A,=C,=G", "QSCORE": "10,20,30"},
            {"MIDSV": "=G,=A,=C", "QSCORE": "30,40,50"},
            1,  # Microhomology length for "=G"
        ),
        (
            {"MIDSV": "=A,=C,=G,=T,=A", "QSCORE": "10,20,30,40,50"},
            {"MIDSV": "=C,=A,=C,=G", "QSCORE": "40,50,60,70"},
            0,  # No microhomology for odd length
        ),
        (
            {"MIDSV": "=A,=C,=G", "QSCORE": "10,20,30"},
            {"MIDSV": "=C,=G,=C,=G", "QSCORE": "20,30,60,70"},
            2,  # prevのほうが短い
        ),
    ],
    ids=[
        "partial_match",  # 部分一致がある場合
        "no_match",  # 一致がない場合
        "full_match",  # 完全一致
        "multi_base_match",  # 複数の塩基にわたる一致
        "qscore_mismatch",  # QSCOREの不一致で一致が途切れる
        "odd_length_match",  # 奇数長で部分一致
        "odd_length_no_match",  # 奇数長で一致なし
        "odd_length_2_match",  # prevのほうが短い
    ],
)
def test_calculate_microhomology(previous_alignment, current_alignment, expected):
    assert polisher.calculate_microhomology(previous_alignment, current_alignment) == expected


@pytest.mark.parametrize(
    "current_alignment, num_microhomology, expected_alignment",
    [
        ({"MIDSV": "A,B,C,D", "POS": 1}, 2, {"MIDSV": "C,D", "POS": 3}),
        (
            {"MIDSV": "A,B,C,D", "QSCORE": "10,20,30,40", "POS": 1},
            1,
            {"MIDSV": "B,C,D", "QSCORE": "20,30,40", "POS": 2},
        ),
        ({"MIDSV": "A,B,C,D", "QSCORE": "10,20,30,40", "POS": 1}, 3, {"MIDSV": "D", "QSCORE": "40", "POS": 4}),
        ({"MIDSV": "A,B,C,D", "POS": 1}, 0, {"MIDSV": "A,B,C,D", "POS": 1}),
    ],
)
def test_remove_microhomology(
    current_alignment: dict[str, str], num_microhomology: int, expected_alignment: dict[str, str]
) -> None:
    polisher.remove_microhomology(current_alignment, num_microhomology)
    assert current_alignment == expected_alignment


@pytest.mark.parametrize(
    "sam_template, gap, expected_template",
    [
        ({"MIDSV": "A,B,C", "POS": 1}, 2, {"MIDSV": "A,B,C,=N,=N", "POS": 1}),
        (
            {"MIDSV": "A,B,C", "QSCORE": "10,20,30", "POS": 1},
            1,
            {"MIDSV": "A,B,C,=N", "QSCORE": "10,20,30,-1", "POS": 1},
        ),
        (
            {"MIDSV": "A,B,C", "QSCORE": "10,20,30", "POS": 1},
            3,
            {"MIDSV": "A,B,C,=N,=N,=N", "QSCORE": "10,20,30,-1,-1,-1", "POS": 1},
        ),
        ({"MIDSV": "A,B,C", "POS": 1}, 0, {"MIDSV": "A,B,C", "POS": 1}),
    ],
)
def test_fill_gap(sam_template: dict[str, int | str], gap: int, expected_template: dict[str, int | str]) -> None:
    polisher.fill_gap(sam_template, gap)
    assert sam_template == expected_template


# def test_join_control():
#     sam = Path("tests", "data", "join", "test_control.txt").read_text()
#     sam = eval(sam)
#     test = polisher.merge(sam)
#     answer = Path("tests", "data", "join", "answer_control.txt").read_text()
#     answer = eval(answer)
#     assert test == answer


# def test_join_inversion():
#     sam = Path("tests", "data", "join", "test_inv.txt").read_text()
#     sam = eval(sam)
#     test = polisher.merge(sam)
#     answer = Path("tests", "data", "join", "answer_inv.txt").read_text()
#     answer = eval(answer)
#     assert test == answer


# def test_join_deletion():
#     sam = Path("tests", "data", "join", "test_del.txt").read_text()
#     sam = eval(sam)
#     test = polisher.merge(sam)
#     answer = Path("tests", "data", "join", "answer_del.txt").read_text()
#     answer = eval(answer)
#     assert test == answer


# def test_join_deletion_microhomology():
#     samdict = Path("tests", "data", "join", "test_del_microhomology.txt").read_text()
#     samdict = eval(samdict)
#     test = polisher.merge(samdict)
#     answer = Path("tests", "data", "join", "answer_del_microhomology.txt").read_text()
#     answer = eval(answer)
#     assert test == answer


# def test_join_real_microhomology():
#     samdict = Path("tests", "data", "join", "test_real_microhomology.txt").read_text()
#     samdict = eval(samdict)
#     test = polisher.merge(samdict)
#     test = test[0]["CSSPLIT"]
#     answer = Path("tests", "data", "join", "answer_real_microhomology.txt").read_text()
#     answer = eval(answer)
#     assert test == answer


# def test_select_keep_FLAG():
#     samdict = Path("tests", "data", "join", "test_real_microhomology.txt").read_text()
#     samdict = eval(samdict)
#     test = polisher.select(samdict, keep={"FLAG", "SEQ"})
#     test = set(list(test[0].keys()))
#     answer = {"QNAME", "RNAME", "FLAG", "CSSPLIT", "QSCORE"}
#     assert test == answer


# def test_pad():
#     sam = Path("tests", "data", "pad", "padding.sam")
#     sam = io.read_sam(str(sam))
#     sam = list(sam)
#     sqheaders = format.extract_sqheaders(sam)
#     samdict = format.dictionarize_sam(sam)
#     for i, alignment in enumerate(samdict):
#         samdict[i]["MIDSV"] = convert.cstag_to_midsv(alignment["CSTAG"])
#         samdict[i]["CSSPLIT"] = convert.cstag_to_cssplit(alignment["CSTAG"])
#         samdict[i]["QSCORE"] = convert.qual_to_qscore_midsv(alignment["QUAL"], alignment["MIDSV"])
#     test = polisher.pad(samdict, sqheaders)
#     for t in test:
#         mlen = len(t["MIDSV"].split(","))
#         clen = len(t["CSSPLIT"].split(","))
#         qlen = len(t["QSCORE"].split(","))
#         assert mlen == clen == qlen
#     answer = Path("tests", "data", "pad", "answer_pad.txt").read_text()
#     answer = eval(answer)
#     assert test == answer
