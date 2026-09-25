"""Simulation-only Linac--BT--DR reference handoff audit.

The script deliberately does *not* claim a physical DR injection model.  The
historical SAD Linac+BT line and the DR ring start at different coordinate
origins.  The BT line itself includes the historical septa and nominal BK1R
kicker through IPZT, but their survey/pulse calibration into the DR periodic
coordinates, longitudinal synchronization, and apertures are not represented.
It records the affine transverse handoff required to map the Linac+BT design
particle to the DR periodic closed orbit, then uses that explicitly labelled
provisional handoff for a one-turn DR survival check.

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
    a design-particle offset.  The SAD BT line already contains its historical
    septa and nominal BK1R kicker, but a surveyed IPZT-to-RING0 map and the
    measured pulsed settings must replace this handoff before it becomes an
    injection model.
    """

    matrix: np.ndarray
    offset_mm_mrad: np.ndarray
    provenance: str
    source_dispersion_mm_mrad: np.ndarray | None = None
    target_dispersion_mm_mrad: np.ndarray | None = None
    reference_momentum_mev_c: float | None = None

    def __post_init__(self) -> None:
        if np.asarray(self.matrix).shape != (4, 4):
            raise ValueError("handoff matrix must have shape (4, 4)")
        if np.asarray(self.offset_mm_mrad).shape != (4,):
            raise ValueError("handoff offset must have shape (4,)")
        dispersions = (self.source_dispersion_mm_mrad, self.target_dispersion_mm_mrad)
        if any(item is not None for item in dispersions):
            if any(item is None for item in dispersions):
                raise ValueError("source and target dispersion must be supplied together")
            if any(np.asarray(item).shape != (4,) for item in dispersions):
                raise ValueError("handoff dispersions must each have shape (4,)")
            if self.reference_momentum_mev_c is None or self.reference_momentum_mev_c <= 0.0:
                raise ValueError("a positive reference momentum is required with dispersion")

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

    @staticmethod
    def _normalising_matrix(beta_m: float, alpha: float) -> np.ndarray:
        """Return the unit-determinant Courant--Snyder normalising matrix."""
        if beta_m <= 0.0:
            raise ValueError("Twiss beta must be positive")
        root_beta = float(np.sqrt(beta_m))
        return np.array(
            ((root_beta, 0.0), (-alpha / root_beta, 1.0 / root_beta)),
            dtype=float,
        )

    @classmethod
    def twiss_dispersion_matched(
        cls,
        source_coordinates: np.ndarray,
        target_coordinates: np.ndarray,
        *,
        source_twiss: tuple[float, float, float, float],
        target_twiss: tuple[float, float, float, float],
        source_dispersion_mm_mrad: np.ndarray,
        target_dispersion_mm_mrad: np.ndarray,
        reference_momentum_mev_c: float,
        provenance: str | None = None,
    ) -> "TransverseHandoff":
        """Make a zero-phase 4D symplectic Twiss/dispersion matching map.

        This maps the *design covariance* at the historical BT endpoint to
        the chosen DR RING0 reference point.  It is useful as a reproducible
        optics baseline, but is explicitly not a surveyed IPZT-to-RING0 map:
        longitudinal path length, septum/kicker pulse calibration, coupling,
        and an arbitrary betatron phase remain outside this construction.
        """
        bx_s, ax_s, by_s, ay_s = source_twiss
        bx_t, ax_t, by_t, ay_t = target_twiss
        x_map = cls._normalising_matrix(bx_t, ax_t) @ np.linalg.inv(
            cls._normalising_matrix(bx_s, ax_s)
        )
        y_map = cls._normalising_matrix(by_t, ay_t) @ np.linalg.inv(
            cls._normalising_matrix(by_s, ay_s)
        )
        matrix = np.zeros((4, 4), dtype=float)
        matrix[:2, :2] = x_map
        matrix[2:, 2:] = y_map
        source = np.asarray(source_coordinates, dtype=float)
        target = np.asarray(target_coordinates, dtype=float)
        return cls(
            matrix=matrix,
            offset_mm_mrad=target - matrix @ source,
            provenance=(
                provenance
                or "zero-phase symplectic Twiss/dispersion match from SAD BT IPZT "
                "to RFTrack DR RING0; design-optics baseline, not surveyed injection map"
            ),
            source_dispersion_mm_mrad=np.asarray(source_dispersion_mm_mrad, dtype=float),
            target_dispersion_mm_mrad=np.asarray(target_dispersion_mm_mrad, dtype=float),
            reference_momentum_mev_c=float(reference_momentum_mev_c),
        )

    def apply(self, coordinates: np.ndarray) -> np.ndarray:
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.shape != (4,):
            raise ValueError("handoff coordinates must have shape (4,)")
        return np.asarray(self.matrix, dtype=float) @ coordinates + np.asarray(
            self.offset_mm_mrad, dtype=float
        )

    def apply_phase_space(self, phase_space: np.ndarray) -> np.ndarray:
        """Apply the 4D affine map and its optional first-order dispersion.

        ``phase_space`` follows RF-Track's ``[x, xp, y, yp, ct, p]`` native
        units.  With a matched-dispersion map, betatron coordinates are first
        separated from the source dispersion and then reconstructed with the
        target dispersion.  Time and momentum remain unchanged.
        """
        phase_space = np.asarray(phase_space, dtype=float)
        if phase_space.ndim != 2 or phase_space.shape[1] != 6:
            raise ValueError("phase_space must have shape (particles, 6)")
        mapped = phase_space.copy()
        matrix = np.asarray(self.matrix, dtype=float)
        mapped[:, :4] = (
            matrix @ phase_space[:, :4].T
        ).T + np.asarray(self.offset_mm_mrad, dtype=float)
        if self.source_dispersion_mm_mrad is not None:
            delta = (
                phase_space[:, 5] - float(self.reference_momentum_mev_c)
            ) / float(self.reference_momentum_mev_c)
            source = np.asarray(self.source_dispersion_mm_mrad, dtype=float)
            target = np.asarray(self.target_dispersion_mm_mrad, dtype=float)
            mapped[:, :4] += np.outer(delta, target - matrix @ source)
        return mapped

    def as_dict(self) -> dict[str, object]:
        return {
            "matrix": np.asarray(self.matrix, dtype=float).tolist(),
            "offset_mm_mrad": np.asarray(self.offset_mm_mrad, dtype=float).tolist(),
            "provenance": self.provenance,
            "source_dispersion_mm_mrad": (
                None if self.source_dispersion_mm_mrad is None
                else np.asarray(self.source_dispersion_mm_mrad, dtype=float).tolist()
            ),
            "target_dispersion_mm_mrad": (
                None if self.target_dispersion_mm_mrad is None
                else np.asarray(self.target_dispersion_mm_mrad, dtype=float).tolist()
            ),
            "reference_momentum_mev_c": self.reference_momentum_mev_c,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> "TransverseHandoff":
        """Restore an auditable handoff exported by :meth:`as_dict`.

        This is intentionally a small data-only adapter: measurement fitting,
        file ownership, and any controls connection remain outside the lattice
        model.  Unknown keys are rejected so an incomplete calibration is not
        silently treated as a valid handoff.
        """
        expected = {
            "matrix", "offset_mm_mrad", "provenance",
            "source_dispersion_mm_mrad", "target_dispersion_mm_mrad",
            "reference_momentum_mev_c",
        }
        unknown = set(payload) - expected
        missing = {"matrix", "offset_mm_mrad", "provenance"} - set(payload)
        if unknown or missing:
            raise ValueError(
                "handoff JSON has unexpected/missing keys: "
                f"unknown={sorted(unknown)}, missing={sorted(missing)}"
            )
        return cls(
            matrix=np.asarray(payload["matrix"], dtype=float),
            offset_mm_mrad=np.asarray(payload["offset_mm_mrad"], dtype=float),
            provenance=str(payload["provenance"]),
            source_dispersion_mm_mrad=(
                None if payload.get("source_dispersion_mm_mrad") is None else
                np.asarray(payload["source_dispersion_mm_mrad"], dtype=float)
            ),
            target_dispersion_mm_mrad=(
                None if payload.get("target_dispersion_mm_mrad") is None else
                np.asarray(payload["target_dispersion_mm_mrad"], dtype=float)
            ),
            reference_momentum_mev_c=(
                None if payload.get("reference_momentum_mev_c") is None else
                float(payload["reference_momentum_mev_c"])
            ),
        )


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
                "surveyed IPZT-to-RING0 coordinate map and pulsed kicker settings",
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
