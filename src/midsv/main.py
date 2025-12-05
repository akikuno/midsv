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
    # Validation
    keep = validator.keep_argument(keep)
    validator.validate_sam(path_sam, qscore)

    # Formatting
    sqheaders: dict[str, int] = formatter.extract_sqheaders(io.read_sam(path_sam))
    alignments: list[dict[str, str | int]] = formatter.organize_alignments_to_dict(io.read_sam(path_sam))

    # Conversion to MIDSV
    alignments = converter.convert(alignments, qscore)

    # Polishing
    alignments = polisher.polish(alignments, sqheaders)

    return alignments
