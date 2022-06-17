# from pathlib import Path
# from src.mids import extract, format, convert


# substitution = Path("tests", "data", "substitution", "sub_cslong.sam").read_text().strip().split("\n")
# substitution = [s.split("\t") for s in substitution]

# headers = extract.headers(substitution)
# alignments = extract.alignments(substitution)
# alignments = format.append_reference_sequence_length(headers, alignments)

# mids_list = []
# for alignment in alignments:
#     mids_list.append(convert.cstag_to_mids(alignment["CSTAG"]))

# alignments_with_mids = []
# for alignment, mids in zip(alignments, mids_list):
#     alignment["MIDS"] = mids
#     alignments_with_mids.append(alignment)

