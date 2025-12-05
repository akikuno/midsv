from __future__ import annotations

import json
from collections.abc import Iterator
from pathlib import Path

###########################################################
# Read sam
###########################################################


def read_sam(path_sam: str | Path) -> Iterator[list[str]]:
    sam = Path(path_sam).read_text().strip().split("\n")
    return (s.split("\t") for s in sam)


###########################################################
# Read / Write jsonl
###########################################################


def read_jsonl(path_input: str | Path) -> Iterator[dict[str, str]]:
    with open(path_input, "r") as f:
        for line in f:
            yield json.JSONDecoder(strict=False).decode(line)


def write_jsonl(dicts: list[dict[str, str]], path_output: str | Path):
    with open(path_output, "w") as f:
        for line in dicts:
            json.dump(line, f)
            f.write("\n")
