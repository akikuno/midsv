from __future__ import annotations

import re

def append_reference_sequence_length(headers: list[dict], alignments: list[dict]):
    alingment_with_reference_sequence_length = []
    for alignment in alignments:
        alignment["RLEN"] = headers[alignment["RNAME"]]
        alingment_with_reference_sequence_length.append(alignment)
    return alingment_with_reference_sequence_length

