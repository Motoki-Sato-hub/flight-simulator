"""Read-only helpers for archived ATF DR SAD operation files.

The legacy SAD server writes magnet tables to ``operation/output`` and reads
its correction settings from ``operation/input/RING_ORBIT.DAT``.  This module
does not interpret or run SAD.  It makes the archived files explicit Python
inputs for RFTrack validation and keeps their provenance visible.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re
from typing import Iterable


_HEADER_TIME = re.compile(r"!\s*(\d{4})\.([A-Za-z]+)\.(\d{2})\s+(\d{2}:\d{2}:\d{2})")
_MAGNET_LINE = re.compile(
    r"^\s*(\S+)\s+([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)"
    r"\s+([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s*$"
)
_FUDGE_LINE = re.compile(r"^\s*(\S+)\s+([-+0-9.eEdD]+)\s*$")


def _number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


@dataclass(frozen=True)
class SadMagnetSetting:
    """One six-column record in a SAD ``.DAT``/``.SAV`` table."""

    name: str
    current_a: float
    model_k0: float
    momentum_gev_c: float
    parameter4: float
    parameter5: float


@dataclass(frozen=True)
class SadMagnetTable:
    """An archived SAD output table with its timestamp and records."""

    path: Path
    timestamp: datetime | None
    records: tuple[SadMagnetSetting, ...]

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(record.name for record in self.records)

    def by_name(self) -> dict[str, SadMagnetSetting]:
        return {record.name: record for record in self.records}


@dataclass(frozen=True)
class SadRingInput:
    """The small, operation-relevant subset of ``RING_ORBIT.DAT``."""

    path: Path
    tune_x: float | None
    tune_y: float | None
    cod_steer_x: int | None
    cod_steer_y: int | None
    dispersion_frequency_hz: float | None
    momentum_compaction: float | None
    dispersion_steer_x: int | None
    dispersion_steer_y: int | None
    skew_family_sd: str | None
    skew_family_sf: str | None
    skew_calibration_sd: float | None
    skew_calibration_sf: float | None
    skew_probe_steerers: tuple[str | None, str | None]
    skew_probe_changes: tuple[float | None, float | None]
    fixed_steerers: tuple[str, ...]
    has_bpm_measurements: bool


def read_magnet_table(path: str | Path) -> SadMagnetTable:
    """Read a six-column SAD output table without applying any calibration."""
    table_path = Path(path)
    timestamp = None
    records: list[SadMagnetSetting] = []
    months = {month: index for index, month in enumerate(
        ("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"),
        start=1,
    )}
    for line in table_path.read_text(encoding="utf-8").splitlines():
        header = _HEADER_TIME.match(line)
        if header:
            year, month, day, clock = header.groups()
            if month in months:
                timestamp = datetime.strptime(
                    f"{year}-{months[month]:02d}-{day} {clock}", "%Y-%m-%d %H:%M:%S"
                )
            continue
        matched = _MAGNET_LINE.match(line)
        if matched:
            name, *values = matched.groups()
            records.append(SadMagnetSetting(name, *(_number(value) for value in values)))
    if not records:
        raise ValueError(f"No SAD magnet records found in {table_path}")
    return SadMagnetTable(table_path, timestamp, tuple(records))


def read_fudge_factors(path: str | Path) -> dict[str, float]:
    """Read the static multiplicative quadrupole calibration factors."""
    factors: dict[str, float] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("!"):
            continue
        matched = _FUDGE_LINE.match(line)
        if matched:
            name, value = matched.groups()
            factors[name] = _number(value)
    if not factors:
        raise ValueError(f"No quadrupole fudge factors found in {path}")
    return factors


def _first_number(text: str, key: str) -> float | None:
    matched = re.search(rf"\b{re.escape(key)}\s*->\s*([-+0-9.eEdD]+)", text)
    return _number(matched.group(1)) if matched else None


def _first_string(text: str, key: str) -> str | None:
    matched = re.search(rf'\b{re.escape(key)}\s*->\s*"([^"]+)"', text)
    return matched.group(1) if matched else None


def _first_integer(text: str, key: str) -> int | None:
    value = _first_number(text, key)
    return int(value) if value is not None else None


def read_ring_input(path: str | Path) -> SadRingInput:
    """Extract correction configuration, without evaluating SAD syntax."""
    input_path = Path(path)
    text = input_path.read_text(encoding="utf-8")
    tune = re.search(
        r"ATFTune\s*->\s*\{\s*([-+0-9.eEdD]+)\s*,\s*([-+0-9.eEdD]+)\s*\}", text
    )
    fixed = re.search(r"FixSteer\s*->\s*\{(?P<body>.*?)\}", text, re.DOTALL)
    fixed_steerers = tuple(re.findall(r'"([^"]+)"', fixed.group("body"))) if fixed else ()
    return SadRingInput(
        path=input_path,
        tune_x=_number(tune.group(1)) if tune else None,
        tune_y=_number(tune.group(2)) if tune else None,
        cod_steer_x=_first_integer(text, "NSteerX"),
        cod_steer_y=_first_integer(text, "NSteerY"),
        dispersion_frequency_hz=_first_number(text, "FrequencyChange"),
        momentum_compaction=_first_number(text, "MomentumCompaction"),
        dispersion_steer_x=next(
            iter(reversed([int(item) for item in re.findall(r"NSteerX\s*->\s*(\d+)", text)])), None
        ),
        dispersion_steer_y=next(
            iter(reversed([int(item) for item in re.findall(r"NSteerY\s*->\s*(\d+)", text)])), None
        ),
        skew_family_sd=_first_string(text, "SD1R"),
        skew_family_sf=_first_string(text, "SF1R"),
        skew_calibration_sd=_first_number(text, "sdcalib"),
        skew_calibration_sf=_first_number(text, "sfcalib"),
        skew_probe_steerers=(_first_string(text, "steer1"), _first_string(text, "steer2")),
        skew_probe_changes=(_first_number(text, "steerchange1"), _first_number(text, "steerchange2")),
        fixed_steerers=fixed_steerers,
        has_bpm_measurements=bool(re.search(r"ATFBPM\s*->\s*\{\s*\{", text)),
    )


def output_summary(table: SadMagnetTable) -> dict[str, float | int | str | None]:
    """Return compact, display-ready statistics for one archived SAD output."""
    return {
        "file": table.path.name,
        "timestamp": table.timestamp.isoformat(sep=" ") if table.timestamp else None,
        "records": len(table.records),
        "current_abs_max_a": max(abs(record.current_a) for record in table.records),
        "model_k0_abs_max": max(abs(record.model_k0) for record in table.records),
    }


def common_records(*tables: SadMagnetTable) -> Iterable[tuple[str, tuple[SadMagnetSetting, ...]]]:
    """Yield identically named records, preserving the order of the first table."""
    if not tables:
        return
    lookup = [table.by_name() for table in tables]
    for name in tables[0].names:
        if all(name in mapping for mapping in lookup):
            yield name, tuple(mapping[name] for mapping in lookup)


def rftrack_skew_name_to_sad(name: str) -> str:
    """Map RFTrack's internal skew-element suffix to the SAD magnet name."""
    return name.removesuffix("$SKEW")


def sad_skew_name_to_rftrack(name: str, rftrack_names: Iterable[str]) -> str:
    """Return the unique RFTrack skew actuator corresponding to one SAD name."""
    candidates = [item for item in rftrack_names if rftrack_skew_name_to_sad(item) == name]
    if len(candidates) != 1:
        raise KeyError(f"Expected one RFTrack skew actuator for {name}, found {candidates}")
    return candidates[0]
