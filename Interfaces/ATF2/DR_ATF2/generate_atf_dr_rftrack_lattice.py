#!/usr/bin/env python3
"""Generate the checked-in RF-Track DR lattice data from a SAD daihon.

The generator intentionally understands only the small SAD syntax subset used
by ``atfdr-design-20111111b.sad``.  It does not execute SAD.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re


REFERENCE_SAD_DAIHON = "atfdr-design-20111111b.sad"
REFERENCE_SAD_RELATIVE_PATH = "operation/daihon/atfdr-design-20111111b.sad"
OUTPUT_FILENAME = "ATF_DR_RFTrack_lattice.json"
SUPPORTED_TYPES = ("DRIFT", "BEND", "QUAD", "SEXT", "CAVI", "MONI", "MARK")
ATTRIBUTE_NAMES = ("L", "ANGLE", "K1", "K2", "VOLT", "FREQ", "ROTATE")


def _without_comments(source: str) -> str:
    return "\n".join(line.split("!", 1)[0] for line in source.splitlines())


def _parse_attributes(body: str) -> dict[str, float]:
    names = "|".join(ATTRIBUTE_NAMES)
    pattern = rf"\b({names})\s*=\s*([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+\-]?\d+)?)"
    return {name: float(value) for name, value in re.findall(pattern, body, re.I)}


def parse_sad_daihon(path: Path) -> dict[str, object]:
    raw = path.read_bytes()
    source = raw.decode("utf-8")
    uncommented = _without_comments(source)

    line_match = re.search(
        r"\bLINE\s+RING0\s*=\s*\((.*?)\)\s*;",
        uncommented,
        flags=re.I | re.S,
    )
    if line_match is None:
        raise ValueError("LINE RING0 was not found")
    sequence = re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", line_match.group(1))

    definitions: dict[str, dict[str, object]] = {}
    type_pattern = "|".join(SUPPORTED_TYPES)
    for group in re.finditer(
        rf"^\s*({type_pattern})\s+(.*?);",
        uncommented,
        flags=re.I | re.M | re.S,
    ):
        element_type = group.group(1).upper()
        for definition in re.finditer(
            r"([A-Za-z_][A-Za-z0-9_.]*)\s*=\s*\((.*?)\)",
            group.group(2),
            flags=re.S,
        ):
            name = definition.group(1)
            definitions[name] = {
                "type": element_type,
                "attributes": _parse_attributes(definition.group(2)),
            }

    missing = sorted(set(sequence) - set(definitions))
    if missing:
        raise ValueError(f"Undefined elements in LINE RING0: {missing}")

    circumference = sum(
        float(definitions[name]["attributes"].get("L", 0.0))  # type: ignore[union-attr]
        for name in sequence
    )
    return {
        "metadata": {
            "reference_sad_daihon": REFERENCE_SAD_DAIHON,
            "reference_sad_relative_path": REFERENCE_SAD_RELATIVE_PATH,
            "source_sha256": hashlib.sha256(raw).hexdigest(),
            "ring_line": "RING0",
            "circumference_m": circumference,
            "note": "Generated without executing SAD; numeric values and order are parsed from the daihon.",
        },
        "sequence": sequence,
        "definitions": definitions,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("sad_daihon", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name(OUTPUT_FILENAME),
    )
    args = parser.parse_args()

    if args.sad_daihon.name != REFERENCE_SAD_DAIHON:
        raise ValueError(
            f"Expected {REFERENCE_SAD_DAIHON}, got {args.sad_daihon.name}"
        )
    data = parse_sad_daihon(args.sad_daihon)
    args.output.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
