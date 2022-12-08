from __future__ import annotations
import json
from pathlib import Path


###########################################################
# Read sam
###########################################################


def read_sam(path_of_sam: str | Path) -> list[list]:
    sam = Path(path_of_sam).read_text().strip().split("\n")
    return [s.split("\t") for s in sam]


###########################################################
# Read / Write jsonl
###########################################################


def read_jsonl(filepath: str | Path) -> list[dict]:
    dicts = []
    with open(filepath, "r") as f:
        for line in f:
            dicts.append(json.JSONDecoder().decode(line))
    return dicts


def write_jsonl(dicts: list[dict], filepath: str | Path):
    with open(filepath, "w") as output:
        for d in dicts:
            json.dump(d, output)
            output.write("\n")
