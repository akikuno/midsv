from pathlib import Path
from src.mids.convert_from_sam import *


substitution = Path("tests", "data", "substitution", "sub_cslong.sam").read_text().strip().split("\n")

substitution = [s.split("\t") for s in substitution]

