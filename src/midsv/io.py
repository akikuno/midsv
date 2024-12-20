from __future__ import annotations

import json
from pathlib import Path
from collections.abc import Iterator

###########################################################
# Read sam
###########################################################


def read_sam(path_of_sam: str | Path) -> Iterator[list]:
    sam = Path(path_of_sam).read_text().strip().split("\n")
    return (s.split("\t") for s in sam)


###########################################################
# Read / Write jsonl
###########################################################


def read_jsonl(filepath: str | Path) -> Iterator[dict]:
    with open(filepath, "r") as f:
        for line in f:
            yield json.JSONDecoder(strict=False).decode(line)


def write_jsonl(dicts: list[dict], filepath: str | Path):
    with open(filepath, "w") as output:
        for d in dicts:
            json.dump(d, output)
            output.write("\n")
