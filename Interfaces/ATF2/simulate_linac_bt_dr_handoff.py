"""Simulation-only Linac--BT--DR reference handoff audit.

The script deliberately does *not* claim a physical DR injection model.  The
historical SAD Linac+BT line and the DR ring start at different coordinate
origins, and the septum, injection kicker, longitudinal synchronization, and
apertures between them are not represented.  It records the affine transverse
handoff required to map the Linac+BT design particle to the DR periodic closed
orbit, then uses that explicitly labelled provisional handoff for a one-turn
DR survival check.

The SAD conversion project is found beside ``flight-simulator`` by default;
set ``ATF2_SAD_RFT_ROOT`` to override it.  No control-system or real-machine
connection is made.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

import numpy as np
import RF_Track as rft

from Interfaces.ATF2.InterfaceATF2_DR_RFTrack import InterfaceATF2_DR_RFTrack


DEFAULT_INPUT_MOMENTUM_MEV_C = 82.36546122
# Calibrated in the SAD-derived *combined* RFTrack model, not copied from the
# BT-only 1.542 GeV design.  It gives 1300.0023 MeV/c for the reference probe.
DEFAULT_CAVITY_VOLTAGE_MV = 76.1400


@dataclass(frozen=True)
class TransverseHandoff:
    """A calibratable 4D map from BT exit to the DR ring reference point.

    Coordinates are ``[x_mm, xp_mrad, y_mm, yp_mrad]``.  The default created
    by :meth:`reference_anchored` is intentionally only an identity map plus
    a design-particle offset.  A surveyed linear transport map and measured
    injection-kicker/septum settings must replace it before this becomes an
    injection model.
    """

    matrix: np.ndarray
    offset_mm_mrad: np.ndarray
    provenance: str

    def __post_init__(self) -> None:
        if np.asarray(self.matrix).shape != (4, 4):
            raise ValueError("handoff matrix must have shape (4, 4)")
        if np.asarray(self.offset_mm_mrad).shape != (4,):
            raise ValueError("handoff offset must have shape (4,)")

    @classmethod
    def reference_anchored(
        cls,
        source_coordinates: np.ndarray,
        target_coordinates: np.ndarray,
    ) -> "TransverseHandoff":
        return cls(
            matrix=np.eye(4),
            offset_mm_mrad=np.asarray(target_coordinates) - np.asarray(source_coordinates),
            provenance="reference-anchored identity map; not a surveyed injection map",
        )

    def apply(self, coordinates: np.ndarray) -> np.ndarray:
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.shape != (4,):
            raise ValueError("handoff coordinates must have shape (4,)")
        return np.asarray(self.matrix, dtype=float) @ coordinates + np.asarray(
            self.offset_mm_mrad, dtype=float
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "matrix": np.asarray(self.matrix, dtype=float).tolist(),
            "offset_mm_mrad": np.asarray(self.offset_mm_mrad, dtype=float).tolist(),
            "provenance": self.provenance,
        }


def _sad_project_root() -> Path:
    override = os.environ.get("ATF2_SAD_RFT_ROOT")
    if override:
        return Path(override).expanduser()
    return Path(__file__).resolve().parents[3] / "ATF2_LinacBT_RFTrack"


def _load_rftrack_model(root: Path) -> ModuleType:
    source = root / "rftrack_model.py"
    if not source.is_file():
        raise FileNotFoundError(
            f"Cannot find SAD RFTrack conversion helper: {source}. "
            "Set ATF2_SAD_RFT_ROOT to the project directory."
        )
    spec = importlib.util.spec_from_file_location("atf2_sad_rftrack_model", source)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {source}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _track_dr_one_turn(lattice, coordinates: np.ndarray, momentum_mev_c: float) -> dict[str, object]:
    phase_space = np.r_[coordinates, 0.0, momentum_mev_c]
    bunch = rft.Bunch6d(rft.electronmass, 1.0, -1.0, phase_space)
    result = lattice.track(bunch)
    survived = result.size() == 1
    return {
        "survived": survived,
        "exit_coordinates": (
            np.asarray(result.get_phase_space()[0, :4], dtype=float).tolist()
            if survived
            else None
        ),
    }


def run_handoff_audit(
    *,
    input_momentum_mev_c: float = DEFAULT_INPUT_MOMENTUM_MEV_C,
    cavity_voltage_mv: float = DEFAULT_CAVITY_VOLTAGE_MV,
    handoff: TransverseHandoff | None = None,
) -> dict[str, object]:
    """Return an auditable, explicitly provisional Linac--BT--DR handoff."""
    helper = _load_rftrack_model(_sad_project_root())
    combined_tfs = helper.GENERATED / "atf2_linac_bt_RFTrack.twiss"
    linac_bt, metadata = helper.prepare_accelerating_lattice(
        combined_tfs,
        input_momentum_mev_c=input_momentum_mev_c,
        cavity_voltage_mv=cavity_voltage_mv,
    )
    linac_bt_exit = helper.track_reference(linac_bt, input_momentum_mev_c)
    source_coordinates = np.array(
        [
            linac_bt_exit["x_mm"],
            linac_bt_exit["xp_mrad"],
            linac_bt_exit["y_mm"],
            linac_bt_exit["yp_mrad"],
        ],
        dtype=float,
    )

    dr = InterfaceATF2_DR_RFTrack(population=1.0)
    closed_orbit = dr.ring_correction.find_closed_orbit()
    target_coordinates = np.asarray(closed_orbit.initial_coordinates, dtype=float)
    if handoff is None:
        handoff = TransverseHandoff.reference_anchored(
            source_coordinates, target_coordinates
        )
    # The absolute flight time from the transport is intentionally discarded:
    # longitudinal capture/RF synchronization is not in the current DR model.
    mapped_coordinates = handoff.apply(source_coordinates)
    dr_turn = _track_dr_one_turn(dr.lattice, mapped_coordinates, dr.Pref)

    return {
        "simulation_only": True,
        "status": "provisional_transverse_handoff_only",
        "linac_bt": {
            "input_momentum_mev_c": input_momentum_mev_c,
            "powered_cavity_voltage_mv": cavity_voltage_mv,
            "reference_exit": linac_bt_exit,
            "reference_transmission": 1.0,
        },
        "dr": {
            "reference_momentum_mev_c": dr.Pref,
            "periodic_closed_orbit": target_coordinates.tolist(),
            "one_turn_after_provisional_handoff": dr_turn,
        },
        "handoff": {
            "source_coordinates_mm_mrad": source_coordinates.tolist(),
            "target_coordinates_mm_mrad": target_coordinates.tolist(),
            "mapped_coordinates_mm_mrad": mapped_coordinates.tolist(),
            "transform": handoff.as_dict(),
            "longitudinal_coordinate": "not connected; transport flight time is discarded",
            "missing_physics": [
                "BT-to-DR coordinate survey",
                "septum and injection-kicker maps",
                "apertures and loss monitors",
                "longitudinal RF synchronization and capture",
                "measured injection Twiss and dispersion",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-momentum-mev", type=float, default=DEFAULT_INPUT_MOMENTUM_MEV_C)
    parser.add_argument("--cavity-voltage-mv", type=float, default=DEFAULT_CAVITY_VOLTAGE_MV)
    args = parser.parse_args()
    print(
        json.dumps(
            run_handoff_audit(
                input_momentum_mev_c=args.input_momentum_mev,
                cavity_voltage_mv=args.cavity_voltage_mv,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
