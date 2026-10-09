#!/usr/bin/env python3
"""Convert a SAD daihon to an auditable RF-Track lattice.

The default interpretation follows the current SAD manual: K0, K1, K2 ... are
integrated field components and BEND face angles are ``E1 * ANGLE + AE1`` and
``E2 * ANGLE + AE2``.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping


HERE = Path(__file__).resolve().parent


_KNOWN_DAIHON_DEFAULTS: dict[str, tuple[str, float, Mapping[str, str]]] = {
    "LINAC_200311_daihon.sad": ("linac", 80.0, {}),
    "atfbt199912_daihon.sad": ("bt", 1542.282, {"L2422T": "L2423T"}),
    "atfdr-design-20111111b.sad": ("RING0", 1299.9999, {}),
    "atfexff20081015.sad": ("ATF2", 1300.0, {}),
}


def source_defaults(source: Path) -> tuple[str, float, Mapping[str, str]]:
    """Return source-specific facts only for the archived ATF daihons."""
    return _KNOWN_DAIHON_DEFAULTS.get(source.name, ("", 0.0, {}))


def select_line(document: SadDocument, requested: str | None, source_default: str) -> str:
    """Use a known source default or select the unique LINE without guessing."""
    if requested:
        selected = requested.upper()
        if selected not in document.lines:
            raise ValueError(f"LINE {requested!r} is not defined; candidates: {', '.join(sorted(document.lines))}")
        return selected
    if source_default and source_default.upper() in document.lines:
        return source_default.upper()
    candidates = sorted(document.lines)
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise ValueError("No SAD LINE definition was found in the daihon.")
    raise ValueError(
        "This daihon has multiple beam lines, so a physical route cannot be inferred safely. "
        f"Choose one with --line. Candidates: {', '.join(candidates)}"
    )

SUPPORTED_TYPES = {
    "DRIFT", "BEND", "QUAD", "SEXT", "OCT", "DECA", "DODECA", "MULT",
    "SOL", "CAVI", "CAV", "TCAVI", "COORD", "MONI", "MARK", "APERT",
}

MAX_MULTIPOLE_ORDER = 21
MAX_TFS_MULTIPOLE_ORDER = 20
RFTRACK_MRAD_PER_RAD = 1.0e3
_DEFINITION_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_.$]*)\s*=\s*\((.*?)\)", re.S)
_PARAM_RE = re.compile(
    r"([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.*?)(?=\s+[A-Za-z][A-Za-z0-9_]*\s*=|$)",
    re.S,
)
_NUMBER_WITH_UNIT_RE = re.compile(
    r"(?<![A-Za-z_])([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+\-]?\d+)?)\s*"
    r"(MM|CM|M|MRAD|DEG|RAD|EV|KEV|MEV|GEV|V|KV|MV|HZ|KHZ|MHZ|GHZ)\b",
    re.I,
)
_ASSIGNMENT_RE = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$", re.S)
_LINE_ITEM_RE = re.compile(
    r"[\s,]*(?:(?P<count>\d+)\s*\*\s*)?(?P<reverse>-)?(?P<name>[A-Za-z_][A-Za-z0-9_.$]*)"
)


_UNIT_FACTORS = {
    "MM": 1e-3,
    "CM": 1e-2,
    "M": 1.0,
    "MRAD": 1e-3,
    "DEG": math.pi / 180.0,
    "RAD": 1.0,
    "EV": 1e-6,
    "KEV": 1e-3,
    "MEV": 1.0,
    "GEV": 1e3,
    "V": 1.0,
    "KV": 1e3,
    "MV": 1e6,
    "HZ": 1.0,
    "KHZ": 1e3,
    "MHZ": 1e6,
    "GHZ": 1e9,
}


@dataclass
class Element:
    name: str
    sad_type: str
    attributes: dict[str, float]


@dataclass
class SadDocument:
    sha256: str
    momentum_mev_c: float | None
    definitions: dict[str, Element]
    lines: dict[str, list[str]]
    warnings: list[str]
    unsupported_element_types: list[str]


def _strip_comments(source: str) -> str:
    return "\n".join(line.split("!", 1)[0] for line in source.splitlines())


def _evaluate_expression(expression: str, variables: Mapping[str, float]) -> float:
    """Evaluate literal SAD arithmetic without evaluating SAD code or calls."""
    def replace_unit(match: re.Match[str]) -> str:
        value = match.group(1).replace("D", "E").replace("d", "e")
        return f"({value}*{_UNIT_FACTORS[match.group(2).upper()]:.17g})"

    converted = _NUMBER_WITH_UNIT_RE.sub(replace_unit, expression.strip())

    converted = re.sub(r"(?<=\d)[dD](?=[+\-]?\d)", "E", converted)
    tree = ast.parse(converted, mode="eval")

    def visit(node: ast.AST) -> float:
        if isinstance(node, ast.Expression):
            return visit(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.Name) and node.id.upper() in variables:
            return float(variables[node.id.upper()])
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = visit(node.operand)
            return value if isinstance(node.op, ast.UAdd) else -value
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
            left, right = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
        raise ValueError(f"Unsupported SAD numeric expression: {expression!r}")

    return float(visit(tree))


def _parse_attributes(body: str, variables: Mapping[str, float]) -> dict[str, float]:
    return {
        key.upper(): _evaluate_expression(value, variables)
        for key, value in _PARAM_RE.findall(body.strip())
    }


def _parse_line_items(raw_items: str) -> list[str]:
    items: list[str] = []
    cursor = 0
    while cursor < len(raw_items):
        match = _LINE_ITEM_RE.match(raw_items, cursor)
        if not match:
            tail = raw_items[cursor:].strip()
            if tail:
                raise ValueError(f"Unsupported SAD LINE item syntax: {tail!r}")
            break
        count = int(match.group("count") or 1)
        name = match.group("name").upper()
        if match.group("reverse"):
            name = "-" + name
        items.extend([name] * count)
        cursor = match.end()
    return items


def _expand_line(lines: Mapping[str, list[str]], name: str, stack: tuple[str, ...] = ()) -> list[str]:
    key = name.upper()
    if key.startswith("-"):
        raise ValueError(
            f"Reverse SAD LINE traversal ({key}) is not implemented; "
            "a validated inverse RF-Track map is required."
        )
    if key in stack:
        raise ValueError(f"Recursive SAD LINE: {' -> '.join((*stack, key))}")
    if key not in lines:
        return [key]
    result: list[str] = []
    for child in lines[key]:
        result.extend(_expand_line(lines, child, (*stack, key)))
    return result


def parse_sad(path: Path, aliases_on_redefinition: Mapping[str, str] | None = None) -> SadDocument:
    """Parse the shared, intentionally bounded SAD subset used by ATF daihons."""
    raw = path.read_bytes()
    source = raw.decode("utf-8", errors="replace")
    cleaned = _strip_comments(source)
    variables: dict[str, float] = {}
    definitions: dict[str, Element] = {}
    lines: dict[str, list[str]] = {}
    warnings: list[str] = []
    unsupported_element_types: list[str] = []
    for raw_line in cleaned.splitlines():
        line = raw_line.strip()
        if not line or re.match(r"^(?:LINE|" + "|".join(SUPPORTED_TYPES) + r")\b", line, re.I):
            continue
        for assignment in re.finditer(r"\b([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([^;\n]+)", line):
            name, expression = assignment.groups()
            try:
                variables[name.upper()] = _evaluate_expression(expression, variables)
            except (SyntaxError, ValueError, ZeroDivisionError):
                # This also encounters continuation-line element definitions;
                # they are not scalar variables and are handled below.
                continue

    for block in (item.strip() for item in cleaned.split(";") if item.strip()):
        line_match = re.match(r"^\s*LINE\b(.*)$", block, re.I | re.S)
        if line_match:
            for raw_name, raw_items in _DEFINITION_RE.findall(line_match.group(1)):
                lines[raw_name.upper()] = _parse_line_items(raw_items)
            continue

        typed = re.match(r"^\s*([A-Za-z]+)\s+(.*)$", block, re.I | re.S)
        if typed and typed.group(1).upper() in SUPPORTED_TYPES:
            sad_type, body = typed.group(1).upper(), typed.group(2)
            for raw_name, raw_attributes in _DEFINITION_RE.findall(body):
                name = raw_name.upper()
                attributes = _parse_attributes(raw_attributes, variables)
                if name in definitions and definitions[name] != Element(name, sad_type, attributes):
                    alias = (aliases_on_redefinition or {}).get(name)
                    if alias and alias not in definitions:
                        definitions[alias] = Element(alias, sad_type, attributes)
                        warnings.append(f"{name}: second definition retained as source alias {alias}.")
                        continue
                    warnings.append(f"{name}: later {sad_type} definition overrides earlier definition.")
                definitions[name] = Element(name, sad_type, attributes)
            continue
        if typed and _DEFINITION_RE.search(typed.group(2)):
            unsupported_element_types.append(typed.group(1).upper())
            continue

        assignment = _ASSIGNMENT_RE.match(block)
        if assignment:
            name, expression = assignment.groups()
            try:
                variables[name.upper()] = _evaluate_expression(expression, variables)
            except (SyntaxError, ValueError, ZeroDivisionError):
                continue

    momentum_match = re.search(r"\bMOMENTUM\s*=\s*(.*?);", cleaned, re.I | re.S)
    momentum: float | None = None
    if momentum_match:
        momentum = _evaluate_expression(momentum_match.group(1), variables)
    return SadDocument(
        sha256=hashlib.sha256(raw).hexdigest(),
        momentum_mev_c=momentum,
        definitions=definitions,
        lines=lines,
        warnings=warnings,
        unsupported_element_types=sorted(set(unsupported_element_types)),
    )


def _keyword(element: Element) -> str:
    if element.sad_type == "DRIFT":
        return "DRIFT"
    if element.sad_type == "QUAD":
        return "QUADRUPOLE"
    if element.sad_type == "SEXT":
        return "SEXTUPOLE"
    if element.sad_type in {"OCT", "DECA", "DODECA", "MULT"}:
        return "MULTIPOLE"
    if element.sad_type == "SOL":
        return "SOLENOID"
    if element.sad_type in {"CAVI", "CAV"}:
        return "RFCAVITY"
    if element.sad_type == "TCAVI":
        return "TRANSVERSE_RFCAVITY"
    if element.sad_type == "MONI":
        return "MONITOR"
    if element.sad_type in {"MARK", "COORD", "APERT"}:
        return "MARKER"
    if element.sad_type == "BEND":
        if abs(element.attributes.get("ANGLE", 0.0)) > 0.0:
            return "SBEND"
        if element.name.startswith(("ZV", "ZY")):
            return "VKICKER"
        if element.name.startswith(("ZH", "ZX")):
            return "HKICKER"
        if abs(element.attributes.get("K0", 0.0)) > 0.0 or abs(element.attributes.get("K1", 0.0)) > 0.0:
            raise ValueError(
                f"{element.name}: zero-angle BEND with field strength has ambiguous corrector orientation; "
                "use an explicit ZH/ZX or ZV/ZY name."
            )
        return "HKICKER"
    raise ValueError(f"Unsupported SAD element type {element.sad_type}")


def _record(element: Element, s_exit: float) -> dict[str, Any]:
    attrs = element.attributes
    length = float(attrs.get("L", 0.0))
    keyword = _keyword(element)
    angle = float(attrs.get("ANGLE", 0.0)) if keyword == "SBEND" else 0.0
    k0l = float(attrs.get("K0", 0.0))
    k1l = float(attrs.get("K1", 0.0))
    k2l = float(attrs.get("K2", 0.0))
    e1 = float(attrs.get("E1", 0.0)) * angle + float(attrs.get("AE1", 0.0))
    e2 = float(attrs.get("E2", 0.0)) * angle + float(attrs.get("AE2", 0.0))
    phase = float(attrs.get("PHI", -math.pi / 2.0))
    knl = [0.0] * (MAX_MULTIPOLE_ORDER + 1)
    ksl = [0.0] * (MAX_MULTIPOLE_ORDER + 1)
    for key, value in attrs.items():
        normal = re.fullmatch(r"K(\d+)", key)
        skew = re.fullmatch(r"SK(\d+)", key)
        if normal and int(normal.group(1)) < len(knl):
            knl[int(normal.group(1))] = float(value)
        if skew and int(skew.group(1)) < len(ksl):
            ksl[int(skew.group(1))] = float(value)
    order_for_type = {"OCT": 3, "DECA": 4, "DODECA": 5}
    if element.sad_type in order_for_type:
        order = order_for_type[element.sad_type]
        knl[order] = float(attrs.get(f"K{order}", 0.0))
        ksl[order] = float(attrs.get(f"SK{order}", 0.0))
    if element.sad_type == "QUAD":
        knl[1] = k1l
    if element.sad_type == "SEXT":
        knl[2] = k2l
    mapping_status = "exact-static"
    mapping_note = ""
    if element.sad_type == "COORD":
        mapping_status = "preserved-unmodelled-coordinate-transform"
        mapping_note = "SAD COORD has no one-to-one RF-Track lattice element; exported as a zero-length marker."
    elif element.sad_type == "APERT":
        mapping_status = "requires-explicit-aperture-geometry"
        mapping_note = "SAD APERT boundary parameters are retained, but are not reduced to an RF-Track aperture shape automatically."
    elif element.sad_type == "TCAVI":
        mapping_status = "requires-external-rf-field-model"
        mapping_note = "SAD transverse cavity requires a calibrated RF field map or an explicit RFTrack model."
    elif element.sad_type == "SOL" and any(key in attrs for key in ("GEO", "BOUND", "DX", "DY", "DPX", "DPY")):
        mapping_status = "solenoid-field-supported-geometry-review-required"
        mapping_note = "SAD SOL boundary/overlap geometry is retained but cannot be reproduced by one isolated RFTrack Solenoid."
    return {
        "name": element.name,
        "sad_type": element.sad_type,
        "keyword": keyword,
        "s": s_exit,
        "length": length,
        "angle": angle,
        # SAD BEND K0 is an additional integrated dipole kick.  It is K0L
        # in TFS; RF-Track SBend.set_K0() instead receives K0L / L.
        "k0l": k0l,
        "k1l": k1l if keyword == "QUADRUPOLE" or keyword == "SBEND" else 0.0,
        "k2l": k2l if keyword == "SEXTUPOLE" else 0.0,
        # A zero-angle SAD BEND used as a corrector is powered through K0,
        # not through its geometrical ANGLE.
        "hkick": k0l if keyword == "HKICKER" else 0.0,
        "vkick": k0l if keyword == "VKICKER" else 0.0,
        "tilt": float(attrs.get("ROTATE", 0.0)),
        "e1": e1,
        "e2": e2,
        "knl": knl,
        "ksl": ksl,
        "solenoid_bz_t": float(attrs.get("BZ", 0.0)),
        "radius_m": float(attrs.get("RADIUS", 0.0)),
        "offset_x_m": float(attrs.get("DX", 0.0)),
        "offset_y_m": float(attrs.get("DY", 0.0)),
        "volt_v": float(attrs.get("VOLT", 0.0)) if element.sad_type in {"CAVI", "CAV", "TCAVI"} else 0.0,
        "lag": (phase + math.pi / 2.0) / (2.0 * math.pi) if element.sad_type in {"CAVI", "CAV", "TCAVI"} else 0.0,
        "freq_hz": float(attrs.get("FREQ", 0.0)) if element.sad_type in {"CAVI", "CAV", "TCAVI"} else 0.0,
        "mapping_status": mapping_status,
        "mapping_note": mapping_note,
        "source_attributes": attrs,
    }


def convert(source: Path, line: str | None = None) -> dict[str, Any]:
    """Return a common structural lattice manifest using standard SAD semantics."""
    source = source.resolve()
    default_line, default_momentum, aliases = source_defaults(source)
    document = parse_sad(source, aliases)
    if document.unsupported_element_types:
        raise ValueError(
            f"{source.name}: unsupported SAD element types: {document.unsupported_element_types}. "
            "Add an explicit mapping; the converter will not silently drop them."
        )
    selected_line = select_line(document, line, default_line)
    reference_momentum = default_momentum or document.momentum_mev_c
    if not reference_momentum:
        raise ValueError(
            f"{source}: no SAD MOMENTUM and no known reference momentum; "
            "supply a daihon with MOMENTUM before creating an RF-Track lattice."
        )
    sequence = _expand_line(document.lines, selected_line)
    missing = sorted(set(sequence) - set(document.definitions))
    if missing:
        raise ValueError(f"{source.name}: undefined elements in LINE {selected_line}: {missing}")
    s = 0.0
    records: list[dict[str, Any]] = []
    for name in sequence:
        element = document.definitions[name]
        s += float(element.attributes.get("L", 0.0))
        records.append(_record(element, s))
    return {
        "format": "atf-sad-rftrack-lattice/v1",
        "source_id": source.stem,
        "source": {
            # The input directory is local context.  File name plus digest
            # identifies the source without making a shared export local-path dependent.
            "filename": source.name, "sha256": document.sha256,
            "line": selected_line, "sad_momentum_mev_c": document.momentum_mev_c,
        },
        "conventions": {
            "strengths": "SAD integrated K0, K1, K2, ...",
            "bend_faces": "E1 * ANGLE + AE1; E2 * ANGLE + AE2",
            "reference_momentum_mev_c": float(reference_momentum),
        },
        "summary": {
            "occurrences": len(records), "length_m": s,
            "keywords": dict(sorted(Counter(record["keyword"] for record in records).items())),
            "mapping_status": dict(sorted(Counter(record["mapping_status"] for record in records).items())),
        },
        "warnings": document.warnings,
        "records": records,
    }


def build_rftrack_lattice(
    manifest: Mapping[str, Any],
    *,
    charge: float = -1.0,
    rf_mode: str = "disabled",
) -> Any:
    if rf_mode not in {"disabled", "travelling_wave"}:
        raise ValueError("rf_mode must be 'disabled' or 'travelling_wave'")
    if charge == 0.0:
        raise ValueError("charge must be non-zero")
    import numpy as np
    import RF_Track as rft

    momentum = float(manifest["conventions"]["reference_momentum_mev_c"])
    p_over_q = momentum / float(charge)
    lattice = rft.Lattice()
    for record in manifest["records"]:
        kind = str(record["keyword"])
        length = float(record["length"])
        if kind in {"DRIFT", "MARKER", "TRANSVERSE_RFCAVITY"}:
            element = rft.Drift(length)
        elif kind == "MONITOR":
            element = rft.Bpm(length)
        elif kind in {"HKICKER", "VKICKER"}:
            element = rft.Corrector(length)
            element.set_kick(
                p_over_q,
                RFTRACK_MRAD_PER_RAD * float(record["hkick"]),
                RFTRACK_MRAD_PER_RAD * float(record["vkick"]),
            )
        elif kind == "QUADRUPOLE":
            k1 = float(record["k1l"]) / length if length else 0.0
            element = rft.Quadrupole(length, p_over_q, k1)
        elif kind == "SEXTUPOLE":
            k2 = float(record["k2l"]) / length if length else 0.0
            element = rft.Sextupole(length, p_over_q, k2) if length else rft.Multipole(0.0)
        elif kind == "MULTIPOLE":
            strengths = np.asarray(record["knl"], dtype=float) + 1j * np.asarray(record["ksl"], dtype=float)
            element = rft.Multipole(length)
            element.set_KnL(p_over_q, strengths)
        elif kind == "SOLENOID":
            element = rft.Solenoid(length, float(record["solenoid_bz_t"]), float(record["radius_m"]))
        elif kind == "SBEND":
            element = rft.SBend(
                length, float(record["angle"]), p_over_q,
                float(record["e1"]), float(record["e2"]),
            )
            if length:
                element.set_K0(float(record["k0l"]) / length)
            elif record["k0l"]:
                raise ValueError(f"{record['name']}: nonzero BEND K0 requires nonzero length")
            element.set_K1L(float(record["k1l"]))
        elif kind == "RFCAVITY":
            if rf_mode == "disabled":
                element = rft.Drift(length)
            else:
                frequency = float(record["freq_hz"])
                voltage = float(record["volt_v"])
                if length <= 0.0 or frequency <= 0.0:
                    raise ValueError(
                        f"{record['name']}: a thin SAD CAVI cannot be converted to a travelling-wave structure"
                    )
                cell_length = rft.clight / (3.0 * frequency)
                cells = int(round(length / cell_length))
                if cells < 1 or not math.isclose(length, cells * cell_length, abs_tol=2e-9):
                    raise ValueError(
                        f"{record['name']}: length/frequency do not define an integral 2pi/3 travelling-wave structure"
                    )
                element = rft.TW_Structure(voltage / length, 0.0, frequency, 2.0 * math.pi / 3.0, cells)
                element.set_phid(-90.0 + 360.0 * float(record["lag"]))
        else:  # guarded by _keyword; keep future additions fail-closed
            raise ValueError(f"{record['name']}: unsupported RF-Track mapping {kind}")
        element.set_name(str(record["name"]))
        # RF-Track permits placement offsets only after the element belongs to
        # a Lattice/Volume.
        lattice.append(element)
        dx, dy, roll = float(record["offset_x_m"]), float(record["offset_y_m"]), float(record["tilt"])
        if dx or dy or roll:
            lattice[int(lattice.size()) - 1].set_offsets(dx, dy, roll)
    return lattice


_TFS_COLUMNS = ["NAME", "KEYWORD", "S", "L", "ANGLE"]
for _order in range(MAX_TFS_MULTIPOLE_ORDER + 1):
    _TFS_COLUMNS.extend((f"K{_order}L", f"K{_order}SL"))
_TFS_COLUMNS.extend(("HKICK", "VKICK", "TILT", "E1", "E2", "VOLT", "LAG", "FREQ"))


def write_tfs(manifest: Mapping[str, Any], output: Path) -> None:
    """Write the static magnetic/RF input representation understood by RF-Track."""
    source = manifest["source"]
    conventions = manifest["conventions"]
    summary = manifest["summary"]
    high_order = [
        record["name"] for record in manifest["records"]
        if any(abs(float(value)) > 0.0 for value in (
            *record["knl"][MAX_TFS_MULTIPOLE_ORDER + 1:],
            *record["ksl"][MAX_TFS_MULTIPOLE_ORDER + 1:],
        ))
    ]
    if high_order:
        raise ValueError(
            "RF-Track TFS import supports multipoles through K20/SK20; "
            f"use build_rftrack_lattice() for SAD K21/SK21 ({', '.join(high_order)})."
        )
    powered_correctors = [
        record["name"] for record in manifest["records"]
        if abs(float(record["hkick"])) > 0.0 or abs(float(record["vkick"])) > 0.0
    ]
    if powered_correctors:
        raise ValueError(
            "RF-Track's TFS importer does not retain HKICK/VKICK values; "
            f"use build_rftrack_lattice() for powered correctors ({', '.join(powered_correctors)})."
        )
    pc_gev = float(conventions["reference_momentum_mev_c"]) / 1000.0
    rows: list[str] = []
    for record in manifest["records"]:
        values: list[str] = []
        for column in _TFS_COLUMNS:
            if column == "NAME":
                values.append(f'"{record["name"]}"')
            elif column == "KEYWORD":
                keyword = "MARKER" if record["keyword"] == "TRANSVERSE_RFCAVITY" else record["keyword"]
                values.append(f'"{keyword}"')
            elif column == "S":
                values.append(f'{record["s"]:.15g}')
            elif column == "L":
                values.append(f'{record["length"]:.15g}')
            elif column == "ANGLE":
                values.append(f'{record["angle"]:.15g}')
            elif match := re.fullmatch(r"K(\d+)L", column):
                order = int(match.group(1))
                if record["keyword"] == "MULTIPOLE":
                    value = record["knl"][order]
                elif order == 0:
                    value = record["k0l"]
                elif order == 1:
                    value = record["k1l"]
                elif order == 2:
                    value = record["k2l"]
                else:
                    value = 0.0
                values.append(f"{value:.15g}")
            elif match := re.fullmatch(r"K(\d+)SL", column):
                values.append(f'{record["ksl"][int(match.group(1))]:.15g}')
            elif column == "HKICK":
                values.append(f'{record["hkick"]:.15g}')
            elif column == "VKICK":
                values.append(f'{record["vkick"]:.15g}')
            elif column == "TILT":
                values.append(f'{record["tilt"]:.15g}')
            elif column == "E1":
                values.append(f'{record["e1"]:.15g}')
            elif column == "E2":
                values.append(f'{record["e2"]:.15g}')
            elif column == "VOLT":
                values.append(f'{record["volt_v"] / 1e6:.15g}')
            elif column == "LAG":
                values.append(f'{record["lag"]:.15g}')
            elif column == "FREQ":
                values.append(f'{record["freq_hz"] / 1e6:.15g}')
            else:
                values.append("0")
        rows.append(" ".join(values))
    output.write_text(
        "\n".join((
            '@ NAME %s "ATF_SAD_RFTRACK"', '@ TYPE %s "TWISS"',
            f'@ SEQUENCE %s "{source["line"]}"', '@ PARTICLE %s "ELECTRON"',
            '@ MASS %le 0.00051099895', '@ CHARGE %le -1',
            f'@ ENERGY %le {math.hypot(pc_gev, 0.00051099895):.15g}',
            f'@ PC %le {pc_gev:.15g}', f'@ LENGTH %le {summary["length_m"]:.15g}',
            f'@ N_ELEMENTS %d {summary["occurrences"]}',
            f'@ SAD_SOURCE_SHA256 %s "{source["sha256"]}"',
            '@ SOURCE_ROLE %s "Direct SAD daihon translation; not an optics calculation"',
            '* ' + ' '.join(_TFS_COLUMNS),
            '$ ' + ' '.join('%s' if name in {'NAME', 'KEYWORD'} else '%le' for name in _TFS_COLUMNS),
            *rows, '',
        )), encoding="utf-8"
    )


def write_outputs(manifest: Mapping[str, Any], output_dir: Path) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = str(manifest["source_id"])
    json_path = output_dir / f"atf_{stem}_sad_lattice.json"
    tfs_path = output_dir / f"atf_{stem}_sad_lattice.tfs"
    write_tfs(manifest, tfs_path)
    json_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return json_path, tfs_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "generated_sad_lattices")
    parser.add_argument("--sad", type=Path, required=True, help="SAD daihon file.")
    parser.add_argument("--line", help="Beam-line name only when the daihon has multiple possible routes.")
    args = parser.parse_args()
    manifest = convert(args.sad, args.line)
    json_path, tfs_path = write_outputs(manifest, args.output_dir)
    print(json.dumps({
        "json": str(json_path), "tfs": str(tfs_path), **manifest["summary"],
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
