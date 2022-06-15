from __future__ import annotations

import re

def append_reference_sequence_length(headers: list[dict], alignments: list[dict]) -> list[dict]:
    """Append reference sequence length

    Args:
        headers (list[dict]): headers including SQ (reference name) and LN (reference sequence length)from SAM
        alignments (list[dict]): alignments from SAM

    Returns:
        list[dict]: alignments appended LN named as "RLEN"
    """
    alingment_with_reference_sequence_length = []
    for alignment in alignments:
        alignment["RLEN"] = headers[alignment["RNAME"]]
        alingment_with_reference_sequence_length.append(alignment)
    return alingment_with_reference_sequence_length

