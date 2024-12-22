from __future__ import annotations

import json
from collections.abc import Iterator
from pathlib import Path

###########################################################
# Read sam
###########################################################


def read_sam(path_of_sam: str | Path) -> Iterator[list[str]]:
    sam = Path(path_of_sam).read_text().strip().split("\n")
    return (s.split("\t") for s in sam)


###########################################################
# Read / Write jsonl
###########################################################


def read_jsonl(filepath: str | Path) -> Iterator[dict[str, str]]:
    with open(filepath, "r") as f:
        for line in f:
            yield json.JSONDecoder(strict=False).decode(line)


def write_jsonl(dicts: list[dict[str, str]], filepath: str | Path):
    with open(filepath, "w") as output:
        for line in dicts:
            json.dump(line, output)
            output.write("\n")
