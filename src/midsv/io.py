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


###########################################################
# MIDSV to VCF
###########################################################


def _format_info(info: dict[str, str | int]) -> str:
    if not info:
        return "."
    priority = ["SVTYPE", "TYPE", "SVLEN", "SEQ", "QNAME"]
    seen = set()
    info_items: list[tuple[str, str | int]] = []
    for key in priority:
        if key in info:
            info_items.append((key, info[key]))
            seen.add(key)
    for key in sorted(k for k in info.keys() if k not in seen):
        info_items.append((key, info[key]))
    return ";".join(f"{key}={value}" for key, value in info_items)


def _reference_base(token: str) -> str:
    if token.startswith("+"):
        anchor = next((p for p in reversed(token.split("|")) if not p.startswith("+")), "")
        token = anchor or "=N"
    if len(token) < 2:
        return "N"
    ref_base = token[1]
    return ref_base.upper()


def _is_inversion_token(token: str) -> bool:
    return any(char.islower() for char in token if char.isalpha())


def _parse_insertion_token(token: str) -> tuple[str, str, str, str]:
    parts = token.split("|")
    inserted = "".join(part[1:].upper() for part in parts if part.startswith("+"))
    anchor = next((part for part in reversed(parts) if not part.startswith("+")), "")
    anchor_op = anchor[0] if anchor else "="
    anchor_body = anchor[1:].upper() if len(anchor) > 1 else "N"
    anchor_ref = anchor_body[0] if anchor_body else "N"
    anchor_alt = anchor_ref
    if anchor_op == "*" and len(anchor_body) >= 2:
        anchor_alt = anchor_body[1:]
    elif anchor_op == "-":
        anchor_alt = ""
    return inserted, anchor_ref, anchor_alt, anchor_op


def _alignment_to_vcf_records(
    alignment: dict[str, str | int], large_sv_threshold: int
) -> list[dict[str, object]]:
    chrom = str(alignment["RNAME"])
    qname = str(alignment.get("QNAME", ""))
    tokens = str(alignment["MIDSV"]).split(",")
    records: list[dict[str, str | int]] = []
    pos = 1
    inv_start = None
    inv_bases: list[str] = []
    n_start = None
    n_len = 0

    def flush_inversion() -> None:
        nonlocal inv_start, inv_bases
        if inv_start is None:
            return
        ref_seq = "".join(inv_bases).upper() or "N"
        info = {"SVTYPE": "INV", "SVLEN": len(inv_bases), "SEQ": ref_seq}
        if qname:
            info["QNAME"] = qname
        records.append({"CHROM": chrom, "POS": inv_start, "REF": ref_seq[0], "ALT": "<INV>", "INFO": info})
        inv_start = None
        inv_bases = []

    def flush_unknown_run() -> None:
        nonlocal n_start, n_len
        if n_start is None:
            return
        ref_seq = "N" * n_len
        info = {"TYPE": "DEL", "SVLEN": -n_len, "SEQ": ref_seq}
        if qname:
            info["QNAME"] = qname
        records.append({"CHROM": chrom, "POS": n_start, "REF": ref_seq[0], "ALT": "<DEL>", "INFO": info})
        n_start = None
        n_len = 0

    idx = 0
    while idx < len(tokens):
        token = tokens[idx]

        if _is_inversion_token(token):
            flush_unknown_run()
            if inv_start is None:
                inv_start = pos
            inv_bases.append(_reference_base(token))
            pos += 1
            idx += 1
            continue
        flush_inversion()

        if token.upper() == "=N":
            if n_start is None:
                n_start = pos
                n_len = 0
            n_len += 1
            pos += 1
            idx += 1
            continue
        flush_unknown_run()

        if token.startswith("="):
            pos += 1
            idx += 1
            continue

        if token.startswith("*"):
            ref_base = _reference_base(token)
            alt_base = token[2:].upper() if len(token) >= 3 else ref_base
            info = {"TYPE": "SUB"}
            if qname:
                info["QNAME"] = qname
            records.append({"CHROM": chrom, "POS": pos, "REF": ref_base, "ALT": alt_base, "INFO": info})
            pos += 1
            idx += 1
            continue

        if token.startswith("-"):
            start_pos = pos
            deleted_seq = token[1:].upper()
            consumed = 1
            for look_ahead in tokens[idx + 1 :]:
                if not look_ahead.startswith("-") or _is_inversion_token(look_ahead):
                    break
                deleted_seq += look_ahead[1:].upper()
                consumed += 1
            info = {"TYPE": "DEL", "SVLEN": -len(deleted_seq), "SEQ": deleted_seq}
            if qname:
                info["QNAME"] = qname
            records.append({"CHROM": chrom, "POS": start_pos, "REF": deleted_seq[0], "ALT": "<DEL>", "INFO": info})
            pos += consumed
            idx += consumed
            continue

        if token.startswith("+"):
            inserted_seq, anchor_ref, anchor_alt, anchor_op = _parse_insertion_token(token)
            info = {"TYPE": "INS", "SVLEN": len(inserted_seq)}
            if inserted_seq:
                info["SEQ"] = inserted_seq
            if qname:
                info["QNAME"] = qname
            alt = "<INS>" if len(inserted_seq) > large_sv_threshold else anchor_ref + inserted_seq
            records.append({"CHROM": chrom, "POS": pos, "REF": anchor_ref, "ALT": alt, "INFO": info})
            if anchor_op == "*" and anchor_alt and anchor_alt != anchor_ref:
                sub_info = {"TYPE": "SUB"}
                if qname:
                    sub_info["QNAME"] = qname
                records.append({"CHROM": chrom, "POS": pos, "REF": anchor_ref, "ALT": anchor_alt, "INFO": sub_info})
            pos += 1
            idx += 1
            continue

        pos += 1
        idx += 1

    flush_inversion()
    flush_unknown_run()

    return records


def write_vcf(alignments: list[dict[str, str | int]], path_output: str | Path, large_sv_threshold: int = 50) -> None:
    """Export MIDSV alignments to VCF format.

    Args:
        alignments (list[dict[str, str | int]]): Output of midsv.transform including MIDSV.
        path_output (str | Path): Destination VCF path.
        large_sv_threshold (int, optional): Insertions longer than this use symbolic ALT. Defaults to 50.
    """
    records: list[dict[str, object]] = []
    for alignment in alignments:
        records.extend(_alignment_to_vcf_records(alignment, large_sv_threshold))

    for order, record in enumerate(records):
        record["_order"] = order

    records.sort(key=lambda r: (r["CHROM"], r["POS"], r["_order"]))

    with open(path_output, "w") as f:
        f.write("##fileformat=VCFv4.3\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for record in records:
            info = record.get("INFO", {})
            if isinstance(info, dict) and "END" not in info:
                svtype = info.get("SVTYPE") or info.get("TYPE")
                if svtype == "DEL":
                    svlen = info.get("SVLEN")
                    if svlen is not None:
                        try:
                            svlen_int = int(svlen)
                        except (TypeError, ValueError):
                            svlen_int = None
                        if svlen_int is not None:
                            info["END"] = int(record["POS"]) + abs(svlen_int) - 1
            info_str = _format_info(record["INFO"])
            f.write(
                f"{record['CHROM']}\t{record['POS']}\t.\t{record['REF']}\t{record['ALT']}\t.\tPASS\t{info_str}\n"
            )
