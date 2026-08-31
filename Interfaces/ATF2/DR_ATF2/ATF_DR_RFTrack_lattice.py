"""RF-Track lattice builder for the ATF damping ring.

Reference SAD daihon: ``atfdr-design-20111111b.sad``

The checked-in JSON is generated from the SAD element definitions and
``LINE RING0`` by ``generate_atf_dr_rftrack_lattice.py``.  SAD itself is not
required to build or track this lattice.

This first model is deliberately transverse-only: synchrotron radiation and
the thin SAD RF cavity are both disabled.  Consequently the reference energy
is constant and the lattice is suitable for one-turn optics, closed-orbit,
dispersion and response/correction studies.  Longitudinal capture and damping
studies require a separately validated radiation/RF model.
"""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
from typing import Any

import numpy as np


REFERENCE_SAD_DAIHON = "atfdr-design-20111111b.sad"
REFERENCE_SAD_RELATIVE_PATH = "operation/daihon/atfdr-design-20111111b.sad"
LATTICE_DATA_FILENAME = "ATF_DR_RFTrack_lattice.json"
NOMINAL_MOMENTUM_MEV_C = 1299.9999


def load_lattice_data() -> dict[str, Any]:
    path = Path(__file__).with_name(LATTICE_DATA_FILENAME)
    with path.open(encoding="utf-8") as stream:
        data = json.load(stream)
    reference = data.get("metadata", {}).get("reference_sad_daihon")
    if reference != REFERENCE_SAD_DAIHON:
        raise ValueError(
            f"Unexpected SAD reference in {path}: {reference!r}"
        )
    return data


def get_lattice_metadata() -> dict[str, Any]:
    """Return a copy of the source and geometry metadata."""
    return dict(load_lattice_data()["metadata"])


def _instance_name(
    source_name: str,
    element_type: str,
    occurrences: Counter[str],
    bpm_index: int,
    drift_index: int,
) -> str:
    if source_name == "M":
        return f"MB{bpm_index}R"
    if element_type == "DRIFT":
        return f"DRIFT_{drift_index}"
    if source_name.upper().startswith(("ZH", "ZV")):
        return source_name.upper()
    if source_name == "IDUM":
        return f"IDUMB{occurrences[source_name]}R"
    return f"{source_name.upper()}.{occurrences[source_name]}"


def build_atf_dr_lattice(
    momentum_mev_c: float = NOMINAL_MOMENTUM_MEV_C,
    charge: float = -1.0,
    *,
    rf_mode: str = "disabled",
):
    """Build the ATF DR ``RING0`` lattice using RF-Track elements.

    Parameters
    ----------
    momentum_mev_c:
        Reference momentum used to convert normalized multipole strengths.
    charge:
        Particle charge in elementary-charge units.  The default is electron.
    rf_mode:
        Only ``"disabled"`` is currently accepted.  The zero-length SAD
        cavity is represented by a named, zero-length drift.  This avoids the
        invalid zero-length Pillbox cavity produced by the old TFS importer.
    """
    if rf_mode != "disabled":
        raise NotImplementedError(
            "The ATF DR radiation/RF longitudinal model has not been validated; "
            "use rf_mode='disabled' for transverse correction studies."
        )
    if charge == 0:
        raise ValueError("charge must be non-zero")

    import RF_Track as rft

    data = load_lattice_data()
    definitions = data["definitions"]
    p_over_q = float(momentum_mev_c) / float(charge)
    lattice = rft.Lattice()

    start = rft.Drift(0.0)
    start.set_name("RING0$START")
    lattice.append(start)

    occurrences: Counter[str] = Counter()
    bpm_index = 0
    drift_index = 0

    for source_name in data["sequence"]:
        definition = definitions[source_name]
        element_type = definition["type"]
        attributes = definition["attributes"]
        length = float(attributes.get("L", 0.0))
        occurrences[source_name] += 1

        if element_type == "DRIFT":
            element = rft.Drift(length)
            drift_index += 1
        elif element_type == "BEND":
            angle = float(attributes.get("ANGLE", 0.0))
            if source_name.upper().startswith(("ZH", "ZV")):
                element = rft.Corrector(length)
            elif abs(angle) > 0.0:
                # SAD BEND is represented as a rectangular bend.  RF-Track's
                # face rotations use +ANGLE/2 for the equivalent map.
                element = rft.SBend(
                    length,
                    angle,
                    p_over_q,
                    angle / 2.0,
                    angle / 2.0,
                )
                element.set_K1(float(attributes.get("K1", 0.0)))
            else:
                # Zero-angle injection/extraction bends and BHE are passive in
                # this model; operational DR steerers are the ZH/ZV elements.
                element = rft.Drift(length)
                drift_index += 1
        elif element_type == "QUAD":
            # SAD's K1 in this daihon is the integrated strength K1*L.
            # RF-Track's Quadrupole constructor expects K1 [m^-2].
            k1 = float(attributes.get("K1", 0.0)) / length if length else 0.0
            element = rft.Quadrupole(
                length,
                p_over_q,
                k1,
            )
        elif element_type == "SEXT":
            # Likewise SAD K2 is integrated K2*L, while RF-Track expects
            # normalized K2 [m^-3].
            k2 = float(attributes.get("K2", 0.0)) / length if length else 0.0
            if length:
                element = rft.Sextupole(length, p_over_q, k2)
            else:
                # Preserve SAD's zero-length sextupole-compensation markers
                # without producing NaN from Sextupole(K2L/L) at L=0.
                element = rft.Multipole()
                element.set_strengths(np.zeros(3, dtype=complex))
        elif element_type == "MONI":
            element = rft.Bpm(length)
            if source_name == "M":
                bpm_index += 1
        elif element_type == "MARK":
            element = rft.Drift(length)
            drift_index += 1
        elif element_type == "CAVI":
            element = rft.Drift(0.0)
            drift_index += 1
        else:  # pragma: no cover - protected by the generator
            raise ValueError(f"Unsupported SAD element type: {element_type}")

        instance_name = _instance_name(
            source_name,
            element_type,
            occurrences,
            bpm_index,
            drift_index,
        )
        element.set_name(instance_name)
        lattice.append(element)

        if source_name in {"SD1R", "SF1R"}:
            # SAD's skew-coupling correction drives the skew-quadrupole
            # component of these sextupole locations (skewcor.n).  Keep the
            # physical normal sextupole as above and add a thin, initially
            # zero skew K1L component.  This preserves the daihon geometry
            # while making the same 34 SD1R and 34 SF1R actuators available
            # to the RF-Track correction model.
            skew = rft.Multipole(0.0)
            skew.set_KnL(p_over_q, np.zeros(2, dtype=complex))
            skew.set_name(f"{instance_name}$SKEW")
            lattice.append(skew)

    end = rft.Drift(0.0)
    end.set_name("RING0$END")
    lattice.append(end)
    return lattice


__all__ = [
    "REFERENCE_SAD_DAIHON",
    "REFERENCE_SAD_RELATIVE_PATH",
    "NOMINAL_MOMENTUM_MEV_C",
    "build_atf_dr_lattice",
    "get_lattice_metadata",
    "load_lattice_data",
]
