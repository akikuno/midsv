from __future__ import annotations

from pathlib import Path

from midsv import converter, formatter, io, polisher, validator


def transform(
    path_sam: Path | str,
    qscore: bool = False,
    keep: str | list[str] = None,
) -> list[dict[str, str | int]]:
    """Integrated function to perform MIDSV conversion.

    Args:
        path_sam (str | Path): Path of a SAM file.
        qscore (bool, optional): Output QSCORE. Defaults to False.
        keep (str | list[str], optional): Subset of 'FLAG', 'POS', 'CIGAR', 'SEQ', 'QUAL', 'CSTAG' to keep. Defaults to None.

    Returns:
        list[dict[str, str]]: Dictionary containing QNAME, RNAME, MIDSV, QSCORE, and fields specified by the keep argument.
    """

    keep = validator.keep_argument(keep)

    path_sam = Path(path_sam)
    validator.sam_headers(io.read_sam(path_sam))
    validator.sam_alignments(io.read_sam(path_sam), qscore)

    sqheaders: dict[str, str | int] = formatter.extract_sqheaders(io.read_sam(path_sam))
    samdict: list[dict[str, str | int]] = formatter.organize_alignments_to_dict(io.read_sam(path_sam))

    samdict = converter.convert(samdict, qscore)

    samdict_polished = polisher.polish(samdict, sqheaders, keep)

    return samdict_polished
