from __future__ import annotations

import re


def split_cigar(cigar: str) -> list[str]:
    return re.findall(r"\d+[MIDNSHP=X]", cigar)
