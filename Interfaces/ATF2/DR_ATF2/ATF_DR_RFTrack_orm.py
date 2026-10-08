"""ORM measurement and small, regularised model-update tools for the ATF DR.

The functions in this module deliberately separate a *virtual machine* from
its RFTrack model.  The virtual machine can have hidden quadrupole-family
scales, BPM rolls and corrector calibration gains.  An ORM is measured from
the virtual-machine BPM readout; a small set of physical model parameters can
then be fitted to it.

This is a diagnostic/model-update layer.  It neither accesses control-system
PVs nor applies a correction to a real machine.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

import numpy as np

from .ATF_DR_RFTrack_correction import ATFDRRingCorrection
from .ATF_DR_RFTrack_lattice import NOMINAL_MOMENTUM_MEV_C


RFTRACK_MRAD_PER_RAD = 1.0e3


@dataclass(frozen=True)
class OrbitResponseMatrix:
    """A full two-plane ORM measured at BPM readout coordinates.

    Rows are ``[x_BPM_1, ..., x_BPM_N, y_BPM_1, ..., y_BPM_N]``.  Each
    column is the response to one physical horizontal or vertical corrector
    kick in radians.
    """

    bpm_names: tuple[str, ...]
    corrector_names: tuple[str, ...]
    corrector_planes: tuple[str, ...]
    matrix: np.ndarray
    perturbation_rad: float

    @property
    def shape(self) -> tuple[int, int]:
        return tuple(self.matrix.shape)


@dataclass(frozen=True)
class RegularisedFitResult:
    """Linearised ORM-fit result and its in-sample residual reduction."""

    parameter_names: tuple[str, ...]
    parameter_delta: np.ndarray
    singular_values: np.ndarray
    rank: int
    residual_norm_before: float
    residual_norm_after: float

    @property
    def residual_reduction(self) -> float:
        if self.residual_norm_before == 0.0:
            return 0.0
        return 1.0 - self.residual_norm_after / self.residual_norm_before


def _single_element(lattice, name: str):
    element = lattice[name]
    if isinstance(element, list):
        if len(element) != 1:
            raise ValueError(f"Expected one element named {name}, got {len(element)}")
        return element[0]
    return element


def quadrupole_family(name: str) -> str:
    """Return the SAD family token from an RFTrack instance name."""
    return name.split(".", maxsplit=1)[0]


def capture_quadrupole_k1l(
    lattice,
    *,
    momentum_mev_c: float = NOMINAL_MOMENTUM_MEV_C,
    charge: float = -1.0,
) -> dict[str, float]:
    """Capture normalised integrated strengths before applying family errors."""
    p_over_q = float(momentum_mev_c) / float(charge)
    return {
        element.get_name(): float(element.get_K1L(p_over_q))
        for element in lattice.get_quadrupoles()
    }


def apply_quadrupole_family_scale(
    lattice,
    reference_k1l: Mapping[str, float],
    relative_scale: Mapping[str, float],
    *,
    momentum_mev_c: float = NOMINAL_MOMENTUM_MEV_C,
    charge: float = -1.0,
) -> None:
    """Apply ``K1L = K1L_reference * (1 + family_scale)``.

    ``reference_k1l`` is explicit so repeated calls do not compound an error.
    Families absent from ``relative_scale`` remain at their reference values.
    """
    p_over_q = float(momentum_mev_c) / float(charge)
    present = set()
    for element in lattice.get_quadrupoles():
        name = element.get_name()
        if name not in reference_k1l:
            raise ValueError(f"No reference K1L for quadrupole {name}")
        family = quadrupole_family(name)
        present.add(family)
        element.set_K1L(
            p_over_q,
            float(reference_k1l[name]) * (1.0 + float(relative_scale.get(family, 0.0))),
        )
    unknown = set(relative_scale) - present
    if unknown:
        raise ValueError(f"Unknown quadrupole families: {sorted(unknown)}")


def apply_quadrupole_individual_scale(
    lattice,
    reference_k1l: Mapping[str, float],
    relative_scale: Mapping[str, float],
    *,
    momentum_mev_c: float = NOMINAL_MOMENTUM_MEV_C,
    charge: float = -1.0,
) -> None:
    """Apply per-quadrupole relative K1L errors from an explicit reference.

    This is intentionally separate from :func:`apply_quadrupole_family_scale`.
    A virtual machine can therefore first receive coherent family errors and
    then independent component errors without accidentally compounding either
    error source.
    """
    p_over_q = float(momentum_mev_c) / float(charge)
    present = set()
    for element in lattice.get_quadrupoles():
        name = element.get_name()
        if name not in reference_k1l:
            raise ValueError(f"No reference K1L for quadrupole {name}")
        present.add(name)
        element.set_K1L(
            p_over_q,
            float(reference_k1l[name]) * (1.0 + float(relative_scale.get(name, 0.0))),
        )
    unknown = set(relative_scale) - present
    if unknown:
        raise ValueError(f"Unknown quadrupoles: {sorted(unknown)}")


def bpm_roll_readout(
    positions_mm: np.ndarray,
    rolls_rad: Sequence[float],
    offsets_mm: np.ndarray | None = None,
) -> np.ndarray:
    """Apply the ATF BPM-roll convention to true ``[x, y]`` positions.

    This is a diagnostics error only: it rotates/translates the readout, not
    the beam coordinates in the lattice.  Additive offsets cancel exactly in
    an ideal central-difference ORM, but retaining them here makes that fact
    explicit when the same virtual BPM is used for COD measurements.
    """
    positions = np.asarray(positions_mm, dtype=float)
    rolls = np.asarray(rolls_rad, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 2 or rolls.shape != (positions.shape[0],):
        raise ValueError("positions must have shape (n_bpm, 2) and rolls shape (n_bpm,)")
    offsets = (
        np.zeros_like(positions)
        if offsets_mm is None
        else np.asarray(offsets_mm, dtype=float)
    )
    if offsets.shape != positions.shape:
        raise ValueError("offsets_mm must have shape (n_bpm, 2)")
    cosine = np.cos(rolls)
    sine = np.sin(rolls)
    return np.column_stack((
        cosine * positions[:, 0] + sine * positions[:, 1] + offsets[:, 0],
        -sine * positions[:, 0] + cosine * positions[:, 1] + offsets[:, 1],
    ))


def flatten_bpm_planes(positions_mm: np.ndarray) -> np.ndarray:
    """Convert ``(BPM, x/y)`` positions into the canonical full-ORM row order."""
    positions = np.asarray(positions_mm, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 2:
        raise ValueError("positions must have shape (n_bpm, 2)")
    return np.concatenate((positions[:, 0], positions[:, 1]))


def measure_full_orbit_response(
    machine: ATFDRRingCorrection,
    corrector_names: Sequence[str],
    corrector_planes: Sequence[str],
    *,
    bpm_names: Sequence[str] | None = None,
    bpm_rolls_rad: Sequence[float] | None = None,
    bpm_offsets_mm: np.ndarray | None = None,
    corrector_gains: Mapping[str, float] | None = None,
    perturbation_rad: float = 1.0e-5,
) -> OrbitResponseMatrix:
    """Measure a two-plane ORM of a virtual machine by central differences.

    ``corrector_gains`` are physical multiplicative gains.  A value of 1.01
    means that a requested 10-urad kick changes the virtual machine by
    10.1 urad.  The model normally uses a gain of one until an ORM fit updates
    it.
    """
    names = tuple(corrector_names)
    planes = tuple(plane.lower() for plane in corrector_planes)
    if not names or len(names) != len(planes):
        raise ValueError("corrector_names and corrector_planes must be non-empty and equal length")
    if any(plane not in {"x", "y"} for plane in planes):
        raise ValueError("corrector planes must be 'x' or 'y'")
    if perturbation_rad <= 0.0:
        raise ValueError("perturbation_rad must be positive")

    baseline = machine.find_closed_orbit(bpm_names=bpm_names)
    readout_rolls = (
        np.zeros(len(baseline.bpm_names), dtype=float)
        if bpm_rolls_rad is None
        else np.asarray(bpm_rolls_rad, dtype=float)
    )
    if readout_rolls.shape != (len(baseline.bpm_names),):
        raise ValueError("bpm_rolls_rad must have one value per selected BPM")
    gains = dict(corrector_gains or {})
    matrix = np.empty((2 * len(baseline.bpm_names), len(names)), dtype=float)

    for column, (name, plane) in enumerate(zip(names, planes)):
        element = machine._single_element(name)
        original = machine._get_corrector_kick(element)
        physical_probe = perturbation_rad * float(gains.get(name, 1.0))
        index = 0 if plane == "x" else 1
        plus = original.copy()
        minus = original.copy()
        plus[index] += physical_probe
        minus[index] -= physical_probe
        try:
            machine._set_corrector_kick(element, plus)
            positive = machine.find_closed_orbit(
                initial_coordinates=baseline.initial_coordinates,
                bpm_names=baseline.bpm_names,
            )
            machine._set_corrector_kick(element, minus)
            negative = machine.find_closed_orbit(
                initial_coordinates=baseline.initial_coordinates,
                bpm_names=baseline.bpm_names,
            )
        finally:
            machine._set_corrector_kick(element, original)
        matrix[:, column] = (
            flatten_bpm_planes(bpm_roll_readout(
                positive.bpm_positions, readout_rolls, bpm_offsets_mm
            ))
            - flatten_bpm_planes(bpm_roll_readout(
                negative.bpm_positions, readout_rolls, bpm_offsets_mm
            ))
        ) / (2.0 * perturbation_rad)

    return OrbitResponseMatrix(
        bpm_names=baseline.bpm_names,
        corrector_names=names,
        corrector_planes=planes,
        matrix=matrix,
        perturbation_rad=float(perturbation_rad),
    )


def fit_linearised_orm(
    measured: np.ndarray,
    model: np.ndarray,
    sensitivity: np.ndarray,
    parameter_sigma: Sequence[float],
    parameter_names: Sequence[str],
    *,
    row_weights: np.ndarray | None = None,
    rcond: float = 1.0e-8,
) -> RegularisedFitResult:
    r"""Fit a small model update to an ORM with Gaussian parameter priors.

    The linearisation is ``vec(R_meas - R_model) = J delta_p``.  The prior
    term is ``sum((delta_p / sigma_p)^2)``.  This makes unobservable parameter
    combinations revert to their nominal value rather than acquire arbitrary
    large values.
    """
    measured = np.asarray(measured, dtype=float)
    model = np.asarray(model, dtype=float)
    sensitivity = np.asarray(sensitivity, dtype=float)
    sigma = np.asarray(parameter_sigma, dtype=float)
    names = tuple(parameter_names)
    if measured.shape != model.shape:
        raise ValueError("measured and model ORM shapes must match")
    observations = measured.size
    if sensitivity.shape != (observations, len(names)):
        raise ValueError("sensitivity must have shape (ORM.size, parameter count)")
    if sigma.shape != (len(names),) or np.any(sigma <= 0.0):
        raise ValueError("parameter_sigma must be positive and match parameter_names")
    if rcond < 0.0:
        raise ValueError("rcond must be non-negative")

    residual = (measured - model).reshape(-1)
    if row_weights is None:
        weights = np.ones(observations, dtype=float)
    else:
        weights = np.asarray(row_weights, dtype=float).reshape(-1)
        if weights.shape != (observations,) or np.any(weights < 0.0):
            raise ValueError("row_weights must be non-negative and match ORM.size")

    weighted_sensitivity = sensitivity * weights[:, None]
    weighted_residual = residual * weights
    # Augment least squares with the normalised parameter prior.
    augmented_matrix = np.vstack((weighted_sensitivity, np.diag(1.0 / sigma)))
    augmented_target = np.concatenate((weighted_residual, np.zeros(len(names))))
    u, singular_values, vt = np.linalg.svd(augmented_matrix, full_matrices=False)
    cutoff = rcond * singular_values[0] if singular_values.size else np.inf
    keep = singular_values > cutoff
    delta = np.zeros(len(names), dtype=float)
    if np.any(keep):
        delta = vt[keep].T @ ((u[:, keep].T @ augmented_target) / singular_values[keep])
    fitted_residual = residual - sensitivity @ delta
    return RegularisedFitResult(
        parameter_names=names,
        parameter_delta=delta,
        singular_values=singular_values,
        rank=int(np.count_nonzero(keep)),
        residual_norm_before=float(np.linalg.norm(residual)),
        residual_norm_after=float(np.linalg.norm(fitted_residual)),
    )


__all__ = [
    "OrbitResponseMatrix",
    "RegularisedFitResult",
    "apply_quadrupole_family_scale",
    "apply_quadrupole_individual_scale",
    "bpm_roll_readout",
    "capture_quadrupole_k1l",
    "fit_linearised_orm",
    "flatten_bpm_planes",
    "measure_full_orbit_response",
    "quadrupole_family",
]
