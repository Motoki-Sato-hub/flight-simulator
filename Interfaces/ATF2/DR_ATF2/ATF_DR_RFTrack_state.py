"""Serializable, known machine state for the ATF DR RF-Track model.

This module represents the part of an accelerator state which can normally be
read from controls: normal magnet integrated strengths, corrector kicks, and
the skew-corrector settings.  It intentionally excludes survey/alignment
errors.  Those are not known from ordinary magnet PVs and must not silently be
treated as a state-synchronised digital twin input.

The state format is JSON-compatible so a future controls adapter can write the
same schema without depending on RF-Track.  In an offline simulation it is
also useful for copying the known state of a simulated machine into a separate
model lattice before calculating its response matrix.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .ATF_DR_RFTrack_lattice import NOMINAL_MOMENTUM_MEV_C


STATE_FORMAT = "atf-dr-rftrack-known-machine-state"
STATE_VERSION = 1
RFTRACK_MRAD_PER_RAD = 1.0e3


def _single_element(lattice, name: str):
    element = lattice[name]
    if isinstance(element, list):
        if len(element) != 1:
            raise ValueError(f"Expected one element named {name}, got {len(element)}")
        return element[0]
    return element


def _p_over_q(momentum_mev_c: float, charge: float) -> float:
    if charge == 0.0:
        raise ValueError("charge must be non-zero")
    return float(momentum_mev_c) / float(charge)


def capture_known_machine_state(
    lattice,
    *,
    momentum_mev_c: float = NOMINAL_MOMENTUM_MEV_C,
    charge: float = -1.0,
) -> dict[str, Any]:
    """Capture controls-visible settings from an RF-Track ATF DR lattice.

    Integrated strengths are retained in SAD-compatible normalised units:
    ``K1L`` [m^-1] for quadrupoles and bend gradients, and ``K2L`` [m^-2]
    for sextupoles.  Corrector values use physical radians at this public API
    boundary, rather than RF-Track's internal mrad coordinate convention.
    """
    p_over_q = _p_over_q(momentum_mev_c, charge)
    normal_magnets: dict[str, dict[str, Any]] = {}

    for element in lattice.get_quadrupoles():
        normal_magnets[element.get_name()] = {
            "kind": "quadrupole",
            "k1l_m_inv": float(element.get_K1L(p_over_q)),
        }
    for element in lattice.get_sextupoles():
        normal_magnets[element.get_name()] = {
            "kind": "sextupole",
            "k2l_m_inv2": float(element.get_K2L(p_over_q)),
        }
    for element in lattice.get_sbends():
        normal_magnets[element.get_name()] = {
            "kind": "sbend",
            "k1l_m_inv": float(element.get_K1L()),
        }

    correctors: dict[str, dict[str, float]] = {}
    for element in lattice.get_correctors():
        kick_mrad = np.asarray(element.get_kick(p_over_q), dtype=float)
        correctors[element.get_name()] = {
            "x_kick_rad": float(kick_mrad[0] / RFTRACK_MRAD_PER_RAD),
            "y_kick_rad": float(kick_mrad[1] / RFTRACK_MRAD_PER_RAD),
        }

    # The RF-Track lattice contains thin ``$SKEW`` companions at SD1R/SF1R.
    # They are independently powered correction channels and therefore known
    # settings, unlike the separate ``$ROLL`` companions used for hidden
    # alignment errors.
    skew_correctors: dict[str, dict[str, float]] = {}
    for name in normal_magnets:
        skew_name = f"{name}$SKEW"
        try:
            element = _single_element(lattice, skew_name)
        except Exception:
            continue
        strengths = np.asarray(element.get_KnL(p_over_q), dtype=complex).reshape(-1)
        if strengths.size >= 2:
            skew_correctors[skew_name] = {
                "skew_k1l_m_inv": float(strengths[1].imag),
            }

    return {
        "format": STATE_FORMAT,
        "version": STATE_VERSION,
        "scope": (
            "Controls-visible strengths and corrector settings only; "
            "alignment/survey errors are intentionally excluded."
        ),
        "reference_momentum_mev_c": float(momentum_mev_c),
        "charge_e": float(charge),
        "normal_magnets": normal_magnets,
        "correctors": correctors,
        "skew_correctors": skew_correctors,
    }


def apply_known_machine_state(
    lattice,
    state: dict[str, Any],
    *,
    strict: bool = True,
) -> None:
    """Apply a captured known state to an independently built lattice.

    ``strict=True`` protects against applying a snapshot to a different
    daihon/device naming scheme.  It is appropriate for correction studies.
    ``strict=False`` is reserved for a deliberate partial-state import.
    """
    if state.get("format") != STATE_FORMAT or state.get("version") != STATE_VERSION:
        raise ValueError("Unsupported ATF DR RF-Track machine-state format")
    p_over_q = _p_over_q(state["reference_momentum_mev_c"], state["charge_e"])

    def element_or_skip(name: str):
        try:
            return _single_element(lattice, name)
        except Exception:
            if strict:
                raise ValueError(f"State contains unavailable lattice element {name}")
            return None

    for name, values in state.get("normal_magnets", {}).items():
        element = element_or_skip(name)
        if element is None:
            continue
        kind = values.get("kind")
        if kind == "quadrupole":
            element.set_K1L(p_over_q, float(values["k1l_m_inv"]))
        elif kind == "sextupole":
            element.set_K2L(p_over_q, float(values["k2l_m_inv2"]))
        elif kind == "sbend":
            element.set_K1L(float(values["k1l_m_inv"]))
        else:
            raise ValueError(f"Unsupported state magnet kind {kind!r} for {name}")

    for name, values in state.get("correctors", {}).items():
        element = element_or_skip(name)
        if element is None:
            continue
        element.set_kick(
            p_over_q,
            RFTRACK_MRAD_PER_RAD * float(values["x_kick_rad"]),
            RFTRACK_MRAD_PER_RAD * float(values["y_kick_rad"]),
        )

    for name, values in state.get("skew_correctors", {}).items():
        element = element_or_skip(name)
        if element is None:
            continue
        strengths = np.asarray(element.get_KnL(p_over_q), dtype=complex).copy()
        flat = strengths.reshape(-1)
        if flat.size < 2:
            raise ValueError(f"Skew-state element {name} has no K1L component")
        flat[1] = 1j * float(values["skew_k1l_m_inv"])
        element.set_KnL(p_over_q, strengths)


def write_known_machine_state(path: str | Path, state: dict[str, Any]) -> Path:
    """Write a human-readable JSON snapshot after validating its schema."""
    if state.get("format") != STATE_FORMAT or state.get("version") != STATE_VERSION:
        raise ValueError("Unsupported ATF DR RF-Track machine-state format")
    destination = Path(path)
    destination.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")
    return destination


def read_known_machine_state(path: str | Path) -> dict[str, Any]:
    """Read and minimally validate a JSON machine-state snapshot."""
    state = json.loads(Path(path).read_text(encoding="utf-8"))
    if state.get("format") != STATE_FORMAT or state.get("version") != STATE_VERSION:
        raise ValueError("Unsupported ATF DR RF-Track machine-state format")
    return state


__all__ = [
    "STATE_FORMAT",
    "STATE_VERSION",
    "capture_known_machine_state",
    "apply_known_machine_state",
    "write_known_machine_state",
    "read_known_machine_state",
]
