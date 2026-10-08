"""Simulation-only 6D pipeline for SAD Linac -> BT -> DR RF-Track studies.

The pipeline deliberately keeps the three model boundaries visible:

* the SAD-derived accelerating Linac+BT lattice is prepared at 1.3 GeV/c;
* a calibratable transverse map transfers its exit bunch to the DR reference;
* the DR ``RING0`` model is tracked turn-by-turn for a requested short horizon.

It is an error-free *software* connection, not yet a physical injection model.
The SAD BT ``CELLST`` includes the historical septa and nominal ``BK1R``
kicker through ``IPZT``.  What remains unavailable is the surveyed/calibrated
map from that endpoint into the periodic DR ``RING0`` coordinates, including
the pulsed kicker setting and aperture/loss model.  Reports therefore call the
multi-turn quantity ``ring_survival`` rather than measured/final transmission.
"""

from __future__ import annotations

import importlib.util
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from types import ModuleType
from typing import Any, Mapping

import numpy as np
import RF_Track as rft

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import ATFDRRingCorrection
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import (
    NOMINAL_MOMENTUM_MEV_C,
    build_atf_dr_lattice,
)
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_emittance import (
    find_synchronous_orbit,
)
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_twiss import load_atf_dr_twiss
from Interfaces.ATF2.simulate_linac_bt_dr_handoff import (
    DEFAULT_CAVITY_VOLTAGE_MV,
    DEFAULT_INPUT_MOMENTUM_MEV_C,
    TransverseHandoff,
)


# The SAD helper is dynamically imported so that this repository remains
# self-contained.  Preserve the static *entrance* Twiss after its first
# successful read.  The IPZT optics must instead be calculated from the
# actual RF-Track lattice: changing the accelerating voltage / final energy
# changes its optics, so the historical 1.542 GeV TFS endpoint is not valid
# for the present 1.3 GeV model.
_SAD_ENTRANCE_TWISS_CACHE: dict[str, tuple[float, float, float, float]] = {}
_DR_START_OPTICS_CACHE: tuple[tuple[float, float, float, float], np.ndarray] | None = None

# The online optimisation interface labels eight Linac RF phase setpoints
# CM1L..CM8L, while the SAD Linac has sixteen powered CA structures.  Until a
# hardware RF distribution table is supplied, use the transparent contiguous
# two-structure assumption below.  Every report exports this assumption so it
# cannot be confused with a validated klystron-to-structure wiring map.
KLYSTRON_CAVITY_GROUPS: dict[str, tuple[str, str]] = {
    f"L{group}": (f"CA{2 * group - 1}L", f"CA{2 * group}L")
    for group in range(1, 9)
}


def _dr_start_optics() -> tuple[tuple[float, float, float, float], np.ndarray]:
    """Load the RF-Track RING0-start Twiss and native-coordinate dispersion."""
    global _DR_START_OPTICS_CACHE
    if _DR_START_OPTICS_CACHE is None:
        target_table = load_atf_dr_twiss()
        target = next(
            (row for row in target_table if str(row["NAME"]) == "RING0$START"),
            None,
        )
        if target is None:
            raise ValueError("RFTrack DR Twiss table has no RING0$START row")
        target_twiss = tuple(
            float(target[name]) for name in ("BETX", "ALFX", "BETY", "ALFY")
        )
        target_dispersion = np.array(
            [float(target["DX"]), float(target["DPX"]),
             float(target["DY"]), float(target["DPY"])]
        )
        _DR_START_OPTICS_CACHE = (target_twiss, target_dispersion)
    return _DR_START_OPTICS_CACHE


def _sad_entrance_twiss(
    helper: ModuleType, combined_tfs: Path
) -> tuple[float, float, float, float]:
    """Read the SAD Linac entrance Twiss which seeds the design baseline."""
    key = str(combined_tfs.resolve())
    if key not in _SAD_ENTRANCE_TWISS_CACHE:
        _, rows = helper.read_tfs(combined_tfs)
        entrance = next((row for row in rows if row["NAME"] == "IPP1L"), None)
        if entrance is None:
            raise ValueError(f"{combined_tfs} has no IPP1L Linac entrance marker")
        _SAD_ENTRANCE_TWISS_CACHE[key] = tuple(
            float(entrance[name]) for name in ("BETX", "ALFX", "BETY", "ALFY")
        )
    return _SAD_ENTRANCE_TWISS_CACHE[key]


def _propagate_uncoupled_twiss(
    matrix: np.ndarray,
    entrance_twiss: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    """Propagate two uncoupled projected Twiss pairs through a 4D map.

    RF acceleration changes the geometric emittance; normalising each output
    covariance by its own determinant therefore matters.  The formula stays
    valid in RF-Track's native ``mm, mrad`` coordinates, for which mm/mrad is
    numerically metres.
    """
    transfer = np.asarray(matrix, dtype=float)
    if transfer.shape != (4, 4):
        raise ValueError("Linac+BT transverse transport matrix must be 4 by 4")
    bx, ax, by, ay = entrance_twiss

    def propagate(block: np.ndarray, beta: float, alpha: float) -> tuple[float, float]:
        gamma = (1.0 + alpha**2) / beta
        covariance = np.array(((beta, -alpha), (-alpha, gamma)), dtype=float)
        output = block @ covariance @ block.T
        emittance = float(np.sqrt(np.linalg.det(output)))
        if not np.isfinite(emittance) or emittance <= 0.0:
            raise ValueError("Linac+BT map produced a non-positive projected emittance")
        return float(output[0, 0] / emittance), float(-output[0, 1] / emittance)

    bx_out, ax_out = propagate(transfer[:2, :2], bx, ax)
    by_out, ay_out = propagate(transfer[2:, 2:], by, ay)
    return bx_out, ax_out, by_out, ay_out


def _track_single_particle(lattice: Any, phase_space: np.ndarray) -> np.ndarray:
    """Track one particle and fail rather than hiding a linearisation loss."""
    bunch = rft.Bunch6d(rft.electronmass, 1.0, -1.0, np.asarray(phase_space, dtype=float))
    result = lattice.track(bunch)
    if result.size() != 1:
        raise RuntimeError("Linac+BT reference particle was lost while calculating its map")
    return np.asarray(result.get_phase_space(), dtype=float)[0]


def _rftrack_linac_bt_exit_optics(
    lattice: Any,
    *,
    input_momentum_mev_c: float,
    entrance_twiss: tuple[float, float, float, float],
) -> tuple[tuple[float, float, float, float], np.ndarray]:
    """Compute configured Linac+BT IPZT optics from a centred RF-Track map.

    The output dispersion is differentiated against *exit* relative momentum.
    This matches the convention used by :class:`TransverseHandoff` and avoids
    treating an input-energy offset as though it had the same relative size
    after acceleration.
    """
    reference = np.zeros(6, dtype=float)
    reference[5] = input_momentum_mev_c
    reference_exit = _track_single_particle(lattice, reference)
    steps = np.array((1.0e-3, 1.0e-4, 1.0e-3, 1.0e-4), dtype=float)
    transfer = np.empty((4, 4), dtype=float)
    for column, step in enumerate(steps):
        positive = reference.copy()
        negative = reference.copy()
        positive[column] += step
        negative[column] -= step
        transfer[:, column] = (
            _track_single_particle(lattice, positive)[:4]
            - _track_single_particle(lattice, negative)[:4]
        ) / (2.0 * step)

    relative_input_step = 1.0e-4
    positive = reference.copy()
    negative = reference.copy()
    positive[5] *= 1.0 + relative_input_step
    negative[5] *= 1.0 - relative_input_step
    positive_exit = _track_single_particle(lattice, positive)
    negative_exit = _track_single_particle(lattice, negative)
    relative_exit_step = (
        positive_exit[5] - negative_exit[5]
    ) / (2.0 * reference_exit[5])
    if not np.isfinite(relative_exit_step) or abs(relative_exit_step) < 1.0e-14:
        raise RuntimeError("Could not resolve Linac+BT output relative-momentum step")
    dispersion = (positive_exit[:4] - negative_exit[:4]) / (2.0 * relative_exit_step)
    return _propagate_uncoupled_twiss(transfer, entrance_twiss), dispersion


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
    spec = importlib.util.spec_from_file_location("atf2_sad_rftrack_model_pipeline", source)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {source}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@dataclass(frozen=True)
class BunchSummary:
    particles: int
    charge_e: float
    survival_fraction_from_input: float
    mean_x_mm: float
    mean_xp_mrad: float
    mean_y_mm: float
    mean_yp_mrad: float
    mean_p_mev_c: float
    rms_x_mm: float
    rms_y_mm: float
    rms_p_mev_c: float


@dataclass(frozen=True)
class ProjectedOptics:
    """Projected, dispersion-subtracted two-dimensional bunch optics."""

    particles: int
    geometric_emittance_mm_mrad: float
    beta_m: float
    alpha: float
    mismatch_to_design: float


@dataclass(frozen=True)
class InjectionOptics:
    """DR-start projected optics used as a fast transmission proxy."""

    x: ProjectedOptics
    y: ProjectedOptics
    reference: str


@dataclass(frozen=True)
class EntranceBunchTwiss:
    """Explicit RF-Track 6D parameters at the post-accelerator Linac entrance.

    RF-Track uses mm, mrad and mm/c in ``Bunch6d`` phase space.  Its
    ``Bunch6d_twiss.emitt_*`` inputs are *normalised* mm mrad emittances.
    ``sigma_p_mev_c`` is passed to RF-Track as ``sigma_pt`` and is therefore
    an absolute momentum spread in MeV/c.  The values must come from a
    measured or independently validated ATF operating point: the SAD MARK
    values are preserved in the conversion output but their statistical and
    emittance conventions are not specified by the daihon.
    """

    emitt_x_norm_mm_mrad: float
    emitt_y_norm_mm_mrad: float
    beta_x_m: float
    beta_y_m: float
    alpha_x: float = 0.0
    alpha_y: float = 0.0
    sigma_t_mm_c: float = 0.0
    sigma_p_mev_c: float = 0.0
    mean_x_mm: float = 0.0
    mean_xp_mrad: float = 0.0
    mean_y_mm: float = 0.0
    mean_yp_mrad: float = 0.0

    def __post_init__(self) -> None:
        for name in (
            "emitt_x_norm_mm_mrad", "emitt_y_norm_mm_mrad", "beta_x_m",
            "beta_y_m", "sigma_t_mm_c", "sigma_p_mev_c",
        ):
            if getattr(self, name) < 0.0:
                raise ValueError(f"{name} must be non-negative")
        if self.beta_x_m == 0.0 or self.beta_y_m == 0.0:
            raise ValueError("beta_x_m and beta_y_m must be positive")


@dataclass(frozen=True)
class EntranceBunchPhaseSpace:
    """A directly supplied 6D bunch at the pipeline's ``IPP1L`` boundary.

    ``coordinates_mm_mrad_mm_c_mev_c`` has rows in the exact RF-Track order
    ``[x, xp, y, yp, t, p]`` with units ``[mm, mrad, mm, mrad, mm/c, MeV/c]``.
    Unlike :class:`EntranceBunchTwiss`, this preserves correlations, tails,
    and non-Gaussian structure in a reconstructed or measured distribution.
    It is valid only at the SAD combined-lattice entrance marker ``IPP1L``;
    a distribution from another diagnostic must first be transported to that
    boundary using a separately documented reconstruction.
    """

    coordinates_mm_mrad_mm_c_mev_c: np.ndarray
    charge_e: float
    provenance: str
    location: str = "IPP1L"

    def __post_init__(self) -> None:
        coordinates = np.asarray(self.coordinates_mm_mrad_mm_c_mev_c, dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[1:] != (6,):
            raise ValueError("entrance phase space must have shape (particles, 6)")
        if coordinates.shape[0] < 1:
            raise ValueError("entrance phase space must contain at least one particle")
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("entrance phase-space coordinates must all be finite")
        if np.any(coordinates[:, 5] <= 0.0):
            raise ValueError("entrance phase-space momenta must be positive MeV/c")
        if not np.isfinite(self.charge_e) or self.charge_e < 0.0:
            raise ValueError("entrance bunch charge_e must be finite and non-negative")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("entrance bunch provenance must be a non-empty string")
        if self.location != "IPP1L":
            raise ValueError(
                "this pipeline starts at SAD IPP1L; input location must be exactly 'IPP1L'"
            )
        # Make the stored data independent of a caller-owned mutable array.
        object.__setattr__(self, "coordinates_mm_mrad_mm_c_mev_c", coordinates.copy())

    @property
    def particles(self) -> int:
        return int(self.coordinates_mm_mrad_mm_c_mev_c.shape[0])


_ENTRANCE_BUNCH_JSON_KEYS = {
    "coordinate_order", "coordinates_mm_mrad_mm_c_mev_c", "charge_e",
    "provenance", "location",
}
_ENTRANCE_BUNCH_COORDINATE_ORDER = "x_mm,xp_mrad,y_mm,yp_mrad,t_mm_c,p_mev_c"


def load_entrance_bunch_json(path: str | Path) -> EntranceBunchPhaseSpace:
    """Load a validated, direct 6D IPP1L bunch without controls access.

    The strict JSON schema prevents a silent mm/metre, slope/radian, or
    longitudinal-coordinate mix-up.  Required keys are ``coordinate_order``,
    ``coordinates_mm_mrad_mm_c_mev_c``, ``charge_e``, ``provenance``, and
    ``location``.  The coordinate order must be exactly
    ``x_mm,xp_mrad,y_mm,yp_mrad,t_mm_c,p_mev_c`` and location exactly
    ``IPP1L``.  A file can therefore be versioned alongside a reconstructed
    bunch while remaining entirely offline.
    """
    input_path = Path(path)
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("entrance bunch JSON must be an object")
    unknown = set(payload) - _ENTRANCE_BUNCH_JSON_KEYS
    missing = _ENTRANCE_BUNCH_JSON_KEYS - set(payload)
    if unknown or missing:
        raise ValueError(
            f"entrance bunch JSON keys invalid: unknown={sorted(unknown)}, "
            f"missing={sorted(missing)}"
        )
    if payload["coordinate_order"] != _ENTRANCE_BUNCH_COORDINATE_ORDER:
        raise ValueError(
            "entrance bunch coordinate_order must be "
            f"'{_ENTRANCE_BUNCH_COORDINATE_ORDER}'"
        )
    return EntranceBunchPhaseSpace(
        coordinates_mm_mrad_mm_c_mev_c=np.asarray(
            payload["coordinates_mm_mrad_mm_c_mev_c"], dtype=float
        ),
        charge_e=float(payload["charge_e"]),
        provenance=payload["provenance"],
        location=payload["location"],
    )


@dataclass(frozen=True)
class PipelineResult:
    input: BunchSummary
    linac_bt_exit: BunchSummary
    dr_injection: BunchSummary
    dr_injection_optics: InjectionOptics
    dr_after_turns: BunchSummary
    requested_dr_turns: int
    completed_dr_turns: int
    turn_history_sample_every: int | None
    dr_turn_history_turns: tuple[int, ...]
    dr_turn_history: tuple[BunchSummary, ...]
    handoff: dict[str, object]
    apertures: dict[str, dict[str, tuple[float, float, str | None]]]
    scope: str
    limitations: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "input": asdict(self.input),
            "linac_bt_exit": asdict(self.linac_bt_exit),
            "dr_injection": asdict(self.dr_injection),
            "dr_injection_optics": asdict(self.dr_injection_optics),
            "dr_after_turns": asdict(self.dr_after_turns),
            "requested_dr_turns": self.requested_dr_turns,
            "completed_dr_turns": self.completed_dr_turns,
            "turn_history_sample_every": self.turn_history_sample_every,
            "dr_turn_history_turns": list(self.dr_turn_history_turns),
            "dr_turn_history": [asdict(summary) for summary in self.dr_turn_history],
            "handoff": self.handoff,
            "apertures": self.apertures,
            "scope": self.scope,
            "limitations": list(self.limitations),
        }


def _summary(bunch: Any, input_charge_e: float) -> BunchSummary:
    phase_space = np.asarray(bunch.get_phase_space(), dtype=float)
    charge_e = abs(float(bunch.get_total_charge()))
    if phase_space.size == 0:
        return BunchSummary(
            particles=0, charge_e=charge_e,
            survival_fraction_from_input=0.0 if input_charge_e else 1.0,
            mean_x_mm=np.nan, mean_xp_mrad=np.nan, mean_y_mm=np.nan,
            mean_yp_mrad=np.nan, mean_p_mev_c=np.nan, rms_x_mm=np.nan,
            rms_y_mm=np.nan, rms_p_mev_c=np.nan,
        )
    means = np.mean(phase_space, axis=0)
    rms = np.std(phase_space, axis=0)
    return BunchSummary(
        particles=int(phase_space.shape[0]),
        charge_e=charge_e,
        survival_fraction_from_input=charge_e / input_charge_e if input_charge_e else 1.0,
        mean_x_mm=float(means[0]), mean_xp_mrad=float(means[1]),
        mean_y_mm=float(means[2]), mean_yp_mrad=float(means[3]),
        mean_p_mev_c=float(means[5]), rms_x_mm=float(rms[0]),
        rms_y_mm=float(rms[2]), rms_p_mev_c=float(rms[5]),
    )


def _projected_optics(
    coordinates: np.ndarray,
    *,
    target_beta_m: float,
    target_alpha: float,
) -> ProjectedOptics:
    """Calculate projected Courant--Snyder optics in RF-Track native units."""
    particles = int(coordinates.shape[0])
    if particles < 2:
        return ProjectedOptics(
            particles=particles,
            geometric_emittance_mm_mrad=np.nan,
            beta_m=np.nan,
            alpha=np.nan,
            mismatch_to_design=np.nan,
        )
    covariance = np.cov(coordinates, rowvar=False, ddof=1)
    emittance = float(np.sqrt(max(np.linalg.det(covariance), 0.0)))
    if emittance == 0.0:
        return ProjectedOptics(
            particles=particles,
            geometric_emittance_mm_mrad=0.0,
            beta_m=np.nan,
            alpha=np.nan,
            mismatch_to_design=np.nan,
        )
    beta = float(covariance[0, 0] / emittance)
    alpha = float(-covariance[0, 1] / emittance)
    gamma = float(covariance[1, 1] / emittance)
    target_gamma = (1.0 + target_alpha**2) / target_beta_m
    mismatch = float(0.5 * (
        beta * target_gamma + target_beta_m * gamma - 2.0 * alpha * target_alpha
    ))
    return ProjectedOptics(
        particles=particles,
        geometric_emittance_mm_mrad=emittance,
        beta_m=beta,
        alpha=alpha,
        mismatch_to_design=mismatch,
    )


def _apply_apertures(
    lattice: Any,
    apertures_mm: Mapping[str, tuple[float, float] | tuple[float, float, str]] | None,
    *,
    section: str,
) -> dict[str, tuple[float, float, str | None]]:
    """Attach explicitly supplied RF-Track half-apertures by element name.

    The SAD conversions contain optics but no chamber/survey table.  This
    routine never invents an aperture: an empty mapping means that no physical
    loss boundary is present in that section.
    """
    configured: dict[str, tuple[float, float, str | None]] = {}
    for name, specification in (apertures_mm or {}).items():
        if len(specification) not in {2, 3}:
            raise ValueError(f"{section} {name}: aperture must be (x_mm, y_mm[, shape])")
        x_mm, y_mm = float(specification[0]), float(specification[1])
        if not np.isfinite(x_mm) or not np.isfinite(y_mm) or x_mm <= 0.0 or y_mm <= 0.0:
            raise ValueError(f"{section} {name}: aperture dimensions must be positive finite mm values")
        shape = str(specification[2]) if len(specification) == 3 else None
        elements = lattice.get_elements_by_name(name)
        elements = elements if isinstance(elements, list) else [elements]
        if not elements or elements[0] is None:
            raise ValueError(f"{section} {name}: no lattice element exists for aperture assignment")
        for element in elements:
            if shape is None:
                element.set_aperture(x_mm, y_mm)
            else:
                element.set_aperture(x_mm, y_mm, shape)
        configured[name] = (x_mm, y_mm, shape)
    return configured


class ATF2LinacBTDRRFTrack:
    """SAD Linac+BT and DR tracking with explicit, replaceable handoff."""

    def __init__(
        self,
        *,
        input_momentum_mev_c: float = DEFAULT_INPUT_MOMENTUM_MEV_C,
        cavity_voltage_mv: float = DEFAULT_CAVITY_VOLTAGE_MV,
        klystron_phase_offsets_deg: Mapping[str, float] | None = None,
        magnet_reference_klystron_phase_offsets_deg: Mapping[str, float] | None = None,
        handoff: TransverseHandoff | None = None,
        handoff_mode: str = "reference_anchored",
        linac_bt_apertures_mm: Mapping[str, tuple[float, float] | tuple[float, float, str]] | None = None,
        dr_apertures_mm: Mapping[str, tuple[float, float] | tuple[float, float, str]] | None = None,
        dr_rf_mode: str = "disabled",
        dr_quantum_radiation: bool = False,
        dr_rf_phase_deg: float | None = None,
    ):
        if dr_rf_mode not in {"disabled", "equilibrium"}:
            raise ValueError("dr_rf_mode must be 'disabled' or 'equilibrium'")
        if handoff_mode not in {"reference_anchored", "sad_optics_matched"}:
            raise ValueError(
                "handoff_mode must be 'reference_anchored' or 'sad_optics_matched'"
            )
        self.input_momentum_mev_c = float(input_momentum_mev_c)
        self.cavity_voltage_mv = float(cavity_voltage_mv)
        supplied_klystron_phases = dict(klystron_phase_offsets_deg or {})
        unknown_klystrons = set(supplied_klystron_phases) - set(KLYSTRON_CAVITY_GROUPS)
        if unknown_klystrons:
            raise ValueError(
                f"Unknown Linac klystron phase setpoint(s): {sorted(unknown_klystrons)}"
            )
        self.klystron_phase_offsets_deg = {
            name: float(supplied_klystron_phases.get(name, 0.0))
            for name in KLYSTRON_CAVITY_GROUPS
        }
        if not all(np.isfinite(value) for value in self.klystron_phase_offsets_deg.values()):
            raise ValueError("Linac klystron phase offsets must be finite degrees")
        cavity_phase_adjustments_deg = {
            cavity: self.klystron_phase_offsets_deg[klystron]
            for klystron, cavities in KLYSTRON_CAVITY_GROUPS.items()
            for cavity in cavities
        }
        supplied_magnet_reference_phases = (
            supplied_klystron_phases
            if magnet_reference_klystron_phase_offsets_deg is None
            else dict(magnet_reference_klystron_phase_offsets_deg)
        )
        unknown_magnet_reference_klystrons = (
            set(supplied_magnet_reference_phases) - set(KLYSTRON_CAVITY_GROUPS)
        )
        if unknown_magnet_reference_klystrons:
            raise ValueError(
                "Unknown Linac magnet-reference klystron phase setpoint(s): "
                f"{sorted(unknown_magnet_reference_klystrons)}"
            )
        self.magnet_reference_klystron_phase_offsets_deg = {
            name: float(supplied_magnet_reference_phases.get(name, 0.0))
            for name in KLYSTRON_CAVITY_GROUPS
        }
        if not all(
            np.isfinite(value)
            for value in self.magnet_reference_klystron_phase_offsets_deg.values()
        ):
            raise ValueError("Linac magnet-reference klystron phase offsets must be finite degrees")
        magnet_cavity_phase_adjustments_deg = {
            cavity: self.magnet_reference_klystron_phase_offsets_deg[klystron]
            for klystron, cavities in KLYSTRON_CAVITY_GROUPS.items()
            for cavity in cavities
        }
        helper = _load_rftrack_model(_sad_project_root())
        combined_tfs = helper.GENERATED / "atf2_linac_bt_RFTrack.twiss"
        self.linac_bt_lattice, self.linac_bt_metadata = helper.prepare_accelerating_lattice(
            combined_tfs,
            input_momentum_mev_c=self.input_momentum_mev_c,
            cavity_voltage_mv=self.cavity_voltage_mv,
            cavity_phase_adjustments_deg=cavity_phase_adjustments_deg,
            magnet_cavity_phase_adjustments_deg=magnet_cavity_phase_adjustments_deg,
        )
        self.linac_bt_metadata["klystron_cavity_groups_assumption"] = {
            name: list(cavities) for name, cavities in KLYSTRON_CAVITY_GROUPS.items()
        }
        self.linac_bt_metadata["klystron_phase_offsets_deg"] = dict(
            self.klystron_phase_offsets_deg
        )
        self.linac_bt_metadata["magnet_reference_klystron_phase_offsets_deg"] = dict(
            self.magnet_reference_klystron_phase_offsets_deg
        )
        self.linac_bt_apertures = _apply_apertures(
            self.linac_bt_lattice, linac_bt_apertures_mm, section="Linac+BT"
        )
        reference_exit = helper.track_reference(
            self.linac_bt_lattice, self.input_momentum_mev_c
        )
        self.reference_exit = np.array(
            [reference_exit["x_mm"], reference_exit["xp_mrad"],
             reference_exit["y_mm"], reference_exit["yp_mrad"]],
            dtype=float,
        )
        self.reference_exit_momentum_mev_c = float(reference_exit["p_MeV_c"])
        self.dr_rf_mode = dr_rf_mode
        self.dr_lattice = build_atf_dr_lattice(
            rf_mode=dr_rf_mode,
            radiation_quantum=dr_quantum_radiation,
            rf_phase_deg=dr_rf_phase_deg,
        )
        self.dr_apertures = _apply_apertures(
            self.dr_lattice, dr_apertures_mm, section="DR"
        )
        self.dr_momentum_mev_c = NOMINAL_MOMENTUM_MEV_C
        self.dr_start_twiss, self.dr_start_dispersion_mm_mrad = _dr_start_optics()
        self.dr_synchronous_orbit: np.ndarray | None = None
        if dr_rf_mode == "disabled":
            closed_orbit = ATFDRRingCorrection(
                self.dr_lattice, self.dr_momentum_mev_c
            ).find_closed_orbit()
            self.dr_closed_orbit = np.asarray(
                closed_orbit.initial_coordinates, dtype=float
            )
        else:
            # This aligns only the RFTrack ring's coordinate origin and RF
            # phase.  It does not infer an ATF injection timing or energy trim.
            synchronous = find_synchronous_orbit(self.dr_lattice, max_iterations=50)
            self.dr_synchronous_orbit = np.asarray(synchronous.coordinates, dtype=float)
            self.dr_closed_orbit = self.dr_synchronous_orbit[:4].copy()
        self.handoff_mode = handoff_mode
        if handoff is not None:
            self.handoff = handoff
        elif handoff_mode == "reference_anchored":
            self.handoff = TransverseHandoff.reference_anchored(
                self.reference_exit, self.dr_closed_orbit
            )
        else:
            self.handoff = self._sad_optics_matched_handoff(helper, combined_tfs)
        self.linac_bt_corrector_names = tuple(
            element.get_name() for element in self.linac_bt_lattice.get_correctors()
        )

    def _sad_optics_matched_handoff(
        self, helper: ModuleType, combined_tfs: Path
    ) -> TransverseHandoff:
        """Build a configured IPZT-to-RING0 covariance-matching baseline.

        SAD provides the entrance design Twiss, while the exit Twiss and
        dispersion are calculated from a centred RF-Track map at this
        instance's input energy and cavity voltage.  This is essential after
        scaling the historical 1.542 GeV lattice to 1.3 GeV.  The resulting
        map is still an optics baseline, not a surveyed IPZT-to-RING0 map.
        """
        entrance_twiss = _sad_entrance_twiss(helper, combined_tfs)
        source_twiss, source_dispersion = _rftrack_linac_bt_exit_optics(
            self.linac_bt_lattice,
            input_momentum_mev_c=self.input_momentum_mev_c,
            entrance_twiss=entrance_twiss,
        )
        target_twiss, target_dispersion = _dr_start_optics()
        return TransverseHandoff.twiss_dispersion_matched(
            self.reference_exit,
            self.dr_closed_orbit,
            source_twiss=source_twiss,
            target_twiss=target_twiss,
            source_dispersion_mm_mrad=source_dispersion,
            target_dispersion_mm_mrad=target_dispersion,
            reference_momentum_mev_c=self.reference_exit_momentum_mev_c,
            provenance=(
                "zero-phase symplectic Twiss/dispersion match from the SAD IPP1L "
                "design Twiss propagated by a centred RF-Track Linac+BT map at the "
                "configured energy/RF voltage, to RFTrack DR RING0; design-optics "
                "baseline, not surveyed injection map"
            ),
        )

    def make_reference_bunch(self) -> Any:
        """Return the zero-amplitude SAD-Linac entrance reference particle."""
        return rft.Bunch6d(
            rft.electronmass,
            1.0,
            -1.0,
            np.array([[0.0, 0.0, 0.0, 0.0, 0.0, self.input_momentum_mev_c]]),
        )

    def set_klystron_phase_offsets_deg(self, values: Mapping[str, float]) -> None:
        """Change RF phase offsets while retaining the currently loaded magnets.

        ``values`` is a partial mapping over ``L1`` through ``L8``; omitted
        settings retain their current value.  A group changes both of its
        assumed CA structures by the same relative phase.  No quadrupole or
        bend is re-scaled here, which is the appropriate mode for evaluating
        an RF phase scan against fixed magnet settings.  If a new reference
        optics handoff is wanted instead, construct a fresh machine instance.
        """
        unknown = set(values) - set(KLYSTRON_CAVITY_GROUPS)
        if unknown:
            raise ValueError(f"Unknown Linac klystron phase setpoint(s): {sorted(unknown)}")
        proposed = dict(self.klystron_phase_offsets_deg)
        for name, value in values.items():
            if not np.isfinite(value):
                raise ValueError(f"{name}: klystron phase offset must be finite degrees")
            proposed[name] = float(value)
        for name, cavities in KLYSTRON_CAVITY_GROUPS.items():
            delta = proposed[name] - self.klystron_phase_offsets_deg[name]
            if delta == 0.0:
                continue
            for cavity_name in cavities:
                elements = self.linac_bt_lattice.get_elements_by_name(cavity_name)
                elements = elements if isinstance(elements, list) else [elements]
                if len(elements) != 1:
                    raise ValueError(
                        f"{cavity_name}: expected one RF structure, got {len(elements)}"
                    )
                element = elements[0]
                element.set_phid(element.get_phid() + delta)
        self.klystron_phase_offsets_deg = proposed
        self.linac_bt_metadata["klystron_phase_offsets_deg"] = dict(proposed)

    @staticmethod
    def _corrector_plane(name: str) -> int:
        """Return the active local-strength component of SAD correctors.

        The SAD import preserves vertical devices through the element rotation;
        both horizontal and vertical named devices use component zero locally.
        The global direction is intentionally determined by finite-difference
        response, not inferred from the name at this API boundary.
        """
        return 0

    def _corrector_element(self, name: str) -> Any:
        if name not in self.linac_bt_corrector_names:
            raise ValueError(f"Unknown SAD Linac+BT corrector: {name}")
        elements = self.linac_bt_lattice.get_elements_by_name(name)
        elements = elements if isinstance(elements, list) else [elements]
        if len(elements) != 1:
            raise ValueError(f"Expected exactly one corrector named {name}, got {len(elements)}")
        return elements[0]

    def get_linac_bt_correctors(
        self, names: list[str] | tuple[str, ...] | None = None
    ) -> dict[str, float]:
        """Return RF-Track native kick strengths for selected correctors.

        These values are deliberately *not* labelled as supply current or
        calibrated radians.  A controls/supply calibration belongs in a
        digital-twin adapter; keeping the native quantity here makes offline
        response-matrix and optimiser studies reproducible.
        """
        selected = self.linac_bt_corrector_names if names is None else tuple(names)
        values: dict[str, float] = {}
        for name in selected:
            strength = np.asarray(self._corrector_element(name).get_strength(), dtype=float)
            values[name] = float(strength[self._corrector_plane(name)])
        return values

    def set_linac_bt_correctors(
        self,
        values: Mapping[str, float],
    ) -> None:
        """Set RF-Track native kick strengths without touching real controls."""
        for name, value in values.items():
            if not np.isfinite(value):
                raise ValueError(f"{name}: corrector strength must be finite")
            element = self._corrector_element(name)
            strength = np.asarray(element.get_strength(), dtype=float)
            strength[self._corrector_plane(name)] = float(value)
            element.set_strength(float(strength[0]), float(strength[1]))

    def linac_bt_dr_orbit_response_matrix(
        self,
        names: list[str] | tuple[str, ...] | None = None,
        *,
        step: float = 1.0e-4,
    ) -> tuple[np.ndarray, tuple[str, ...]]:
        """Finite-difference response from correctors to DR injection orbit.

        Rows are ``[x_mm, xp_mrad, y_mm, yp_mrad]`` at the DR handoff; columns
        follow the returned corrector names.  This supports a response-matrix
        baseline against which RL can be compared.  It has no real-machine
        side effects and restores every setting before returning.
        """
        if not np.isfinite(step) or step == 0.0:
            raise ValueError("response-matrix step must be finite and non-zero")
        selected = self.linac_bt_corrector_names if names is None else tuple(names)
        baseline = self.get_linac_bt_correctors(selected)
        reference = self.track(self.make_reference_bunch(), dr_turns=0)
        baseline_orbit = np.array((
            reference.dr_injection.mean_x_mm,
            reference.dr_injection.mean_xp_mrad,
            reference.dr_injection.mean_y_mm,
            reference.dr_injection.mean_yp_mrad,
        ))
        response = np.empty((4, len(selected)), dtype=float)
        try:
            for column, name in enumerate(selected):
                self.set_linac_bt_correctors({name: baseline[name] + step})
                shifted = self.track(self.make_reference_bunch(), dr_turns=0)
                shifted_orbit = np.array((
                    shifted.dr_injection.mean_x_mm,
                    shifted.dr_injection.mean_xp_mrad,
                    shifted.dr_injection.mean_y_mm,
                    shifted.dr_injection.mean_yp_mrad,
                ))
                response[:, column] = (shifted_orbit - baseline_orbit) / step
                self.set_linac_bt_correctors({name: baseline[name]})
        finally:
            self.set_linac_bt_correctors(baseline)
        return response, selected

    def make_entrance_bunch(
        self,
        twiss: EntranceBunchTwiss,
        *,
        particles: int,
        charge_e: float = 1.0e9,
    ) -> Any:
        """Create an explicit finite 6D entrance bunch for the SAD Linac.

        This deliberately has no ATF-default optics.  Supplying ``twiss`` is
        an acknowledgement that the distribution is a stated experimental or
        design assumption, rather than an unverified interpretation of a SAD
        MARK.  It is a convenient bridge for both measured fitted Twiss data
        and controlled sensitivity scans.
        """
        if particles < 1:
            raise ValueError("particles must be positive")
        if charge_e < 0.0:
            raise ValueError("charge_e must be non-negative")
        distribution = rft.Bunch6d_twiss()
        distribution.emitt_x = twiss.emitt_x_norm_mm_mrad
        distribution.emitt_y = twiss.emitt_y_norm_mm_mrad
        distribution.beta_x = twiss.beta_x_m
        distribution.beta_y = twiss.beta_y_m
        distribution.alpha_x = twiss.alpha_x
        distribution.alpha_y = twiss.alpha_y
        distribution.sigma_t = twiss.sigma_t_mm_c
        distribution.sigma_pt = twiss.sigma_p_mev_c
        distribution.mean_x = twiss.mean_x_mm
        distribution.mean_xp = twiss.mean_xp_mrad
        distribution.mean_y = twiss.mean_y_mm
        distribution.mean_yp = twiss.mean_yp_mrad
        return rft.Bunch6d_QR(
            rft.electronmass,
            charge_e,
            -1.0,
            self.input_momentum_mev_c,
            distribution,
            particles,
        )

    def make_entrance_bunch_from_phase_space(
        self, entrance_bunch: EntranceBunchPhaseSpace
    ) -> Any:
        """Build a fresh RF-Track bunch from validated direct 6D input.

        A fresh object is returned on every call because ``Lattice.track``
        consumes/changes its bunch.  This makes endpoint, short-turn, and
        long-turn candidates start from exactly the same supplied distribution
        in the benchmark, rather than accidentally tracking one copy twice.
        """
        if not isinstance(entrance_bunch, EntranceBunchPhaseSpace):
            raise TypeError("entrance_bunch must be an EntranceBunchPhaseSpace")
        return rft.Bunch6d(
            rft.electronmass,
            entrance_bunch.charge_e,
            -1.0,
            entrance_bunch.coordinates_mm_mrad_mm_c_mev_c.copy(),
        )

    def handoff_bunch(self, bunch: Any) -> Any:
        """Apply the transverse handoff while preserving relative time/energy."""
        phase_space = np.asarray(bunch.get_phase_space(), dtype=float)
        if phase_space.size == 0:
            return rft.Bunch6d(
                rft.electronmass, 0.0, -1.0, np.empty((0, 6), dtype=float)
            )
        mapped = self.handoff.apply_phase_space(phase_space)
        # Keep bunch length but remove Linac+BT absolute time of flight.  In
        # RF mode the relative bunch coordinate is placed at the model's
        # synchronous phase.  Momentum is deliberately not shifted: forcing
        # it to the ring value would hide a real injection-energy mismatch.
        mapped[:, 4] -= np.mean(mapped[:, 4])
        if self.dr_synchronous_orbit is not None:
            mapped[:, 4] += self.dr_synchronous_orbit[4]
        return rft.Bunch6d(
            rft.electronmass,
            abs(float(bunch.get_total_charge())),
            -1.0,
            mapped,
        )

    def injection_optics(self, bunch: Any) -> InjectionOptics:
        """Return a fast DR-start optics mismatch proxy for an injected bunch."""
        phase_space = np.asarray(bunch.get_phase_space(), dtype=float)
        if phase_space.size == 0:
            phase_space = np.empty((0, 6), dtype=float)
        delta = (
            phase_space[:, 5] - self.reference_exit_momentum_mev_c
        ) / self.reference_exit_momentum_mev_c
        betatron = phase_space[:, :4] - np.outer(
            delta, self.dr_start_dispersion_mm_mrad
        )
        bx, ax, by, ay = self.dr_start_twiss
        return InjectionOptics(
            x=_projected_optics(
                betatron[:, :2], target_beta_m=bx, target_alpha=ax
            ),
            y=_projected_optics(
                betatron[:, 2:], target_beta_m=by, target_alpha=ay
            ),
            reference=(
                "RFTrack RING0$START Twiss; first-order target dispersion "
                "subtracted using Linac+BT design reference momentum"
            ),
        )

    def track(
        self,
        bunch: Any,
        *,
        dr_turns: int = 1,
        record_turn_history: bool = False,
        turn_history_sample_every: int | None = None,
    ) -> PipelineResult:
        """Track a supplied 6D bunch through Linac+BT and then DR turns.

        ``turn_history_sample_every`` retains every Nth turn (and the final
        completed turn) when ``record_turn_history`` is true.  It keeps a
        10^4--10^5-turn validation auditable without storing an unnecessary
        summary object for every turn.  It reduces memory/report size, not
        RF-Track's per-turn computation cost.
        """
        if dr_turns < 0:
            raise ValueError("dr_turns must be non-negative")
        if turn_history_sample_every is not None and turn_history_sample_every < 1:
            raise ValueError("turn_history_sample_every must be positive")
        history_every = (
            int(turn_history_sample_every)
            if turn_history_sample_every is not None else 1
        )
        input_charge_e = abs(float(bunch.get_total_charge()))
        input_summary = _summary(bunch, input_charge_e)
        linac_bt_exit = self.linac_bt_lattice.track(bunch)
        linac_bt_summary = _summary(linac_bt_exit, input_charge_e)
        dr_bunch = self.handoff_bunch(linac_bt_exit)
        injection_summary = _summary(dr_bunch, input_charge_e)
        injection_optics = self.injection_optics(dr_bunch)
        completed_turns = 0
        turn_history: list[BunchSummary] = []
        turn_history_turns: list[int] = []
        for turn in range(1, dr_turns + 1):
            if dr_bunch.size() == 0:
                break
            dr_bunch = self.dr_lattice.track(dr_bunch)
            completed_turns += 1
            if record_turn_history and (
                turn % history_every == 0 or turn == dr_turns or dr_bunch.size() == 0
            ):
                turn_history_turns.append(turn)
                turn_history.append(_summary(dr_bunch, input_charge_e))
        after_summary = _summary(dr_bunch, input_charge_e)
        return PipelineResult(
            input=input_summary,
            linac_bt_exit=linac_bt_summary,
            dr_injection=injection_summary,
            dr_injection_optics=injection_optics,
            dr_after_turns=after_summary,
            requested_dr_turns=dr_turns,
            completed_dr_turns=completed_turns,
            turn_history_sample_every=(history_every if record_turn_history else None),
            dr_turn_history_turns=tuple(turn_history_turns),
            dr_turn_history=tuple(turn_history),
            handoff=self.handoff.as_dict(),
            apertures={
                "linac_bt": self.linac_bt_apertures,
                "dr": self.dr_apertures,
            },
            scope=(
                "SAD Linac+BT 6D tracking plus "
                + ("RF/radiation DR ring survival" if self.dr_rf_mode == "equilibrium"
                   else "transverse-only DR ring survival")
            ),
            limitations=(
                "The SAD BT line includes historical septa and a nominal BK1R kicker through IPZT.  The selected IPZT-to-RING0 handoff is a design baseline, not a surveyed/calibrated injection map.",
                "The periodic RING0 lattice does not model injection pulse timing or measured kicker settings.",
                "A directly supplied entrance 6D bunch is valid only at IPP1L and does not by itself calibrate its diagnostic reconstruction or transport to that point.",
                *(('No aperture/loss-monitor model is loaded, so ring_survival is not physical final transmission.',)
                  if not (self.linac_bt_apertures or self.dr_apertures) else (
                    'Configured apertures cover only explicitly named elements; their completeness and values must be checked against survey/engineering data before interpreting losses physically.',
                )),
                *(() if self.dr_rf_mode == "equilibrium" else (
                    "DR RF/radiation longitudinal capture is disabled in this pipeline.",
                )),
                *(('DR RF phase is aligned to the model synchronous orbit, but BT time-of-flight and injection-energy calibration are not yet established.',)
                  if self.dr_rf_mode == "equilibrium" else ()),
            ),
        )
