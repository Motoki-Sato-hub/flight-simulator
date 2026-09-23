"""Radiation/RF pilot for ATF DR equilibrium-emittance studies.

This module is intentionally separate from the transverse closed-orbit
correction module.  It uses RF-Track's incoherent synchrotron radiation and
the SAD CAV voltage/frequency represented by
``build_atf_dr_lattice(rf_mode='equilibrium')``.

The reported values are simulated projected and normal-mode emittances.  They
become a Kubo-2003 comparison only after the RF phase, radiation model and
actuator/error conventions have been benchmarked against SAD or measurement.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import RF_Track as rft

from .ATF_DR_RFTrack_lattice import (
    NOMINAL_MOMENTUM_MEV_C,
    build_atf_dr_lattice,
)


SAD_RF_FREQUENCY_HZ = 714e6
RF_WAVELENGTH_MM = rft.clight / SAD_RF_FREQUENCY_HZ * 1e3


@dataclass(frozen=True)
class SynchronousOrbitResult:
    """Six-dimensional fixed point, with longitudinal position modulo RF."""

    coordinates: np.ndarray
    residual: np.ndarray
    iterations: int


@dataclass(frozen=True)
class EmittanceResult:
    """Geometric emittances calculated from a tracked bunch covariance."""

    turns: int
    survived_particles: int
    projected_x_m_rad: float
    projected_y_m_rad: float
    eigen_emittances_m_rad: tuple[float, float]
    covariance_m_rad: np.ndarray


@dataclass(frozen=True)
class EquilibriumEmittanceResult:
    """Radiation equilibrium from a numerical one-turn envelope equation."""

    synchronous_orbit: SynchronousOrbitResult
    spectral_radius: float
    damping_turns_one_over_e: float
    quantum_particles: int
    projected_x_m_rad: float
    projected_y_m_rad: float
    eigen_emittances_m_rad: tuple[float, float]
    one_turn_map: np.ndarray
    diffusion_native: np.ndarray
    covariance_native: np.ndarray


@dataclass(frozen=True)
class TransverseEnvelopeResult:
    """Four-dimensional radiation envelope at a six-dimensional fixed point.

    This diagnostic intentionally retains the RF/deterministic synchronous
    orbit but solves only ``[x, px, y, py]``.  It is useful when a historical
    lattice has a longitudinal RF-map mismatch in the thin-pillbox model.
    It includes transverse diffusion accumulated during a turn, but it is not
    a replacement for the full 6D normal-mode emittance required for a final
    Table-II claim.
    """

    synchronous_orbit: SynchronousOrbitResult
    spectral_radius: float
    damping_turns_one_over_e: float
    quantum_particles: int
    projected_x_m_rad: float
    projected_y_m_rad: float
    eigen_emittances_m_rad: tuple[float, float]
    one_turn_map: np.ndarray
    diffusion_native: np.ndarray
    covariance_native: np.ndarray


def build_equilibrium_lattice(
    *,
    quantum: bool,
    radiation_steps: int = 10,
    rf_phase_deg: float | None = None,
    rf_voltage_scale: float = 1.0,
    lattice_data_path: str | Path | None = None,
) -> object:
    """Build the radiation/RF ATF DR lattice for deterministic or quantum runs."""
    return build_atf_dr_lattice(
        rf_mode="equilibrium",
        radiation_quantum=quantum,
        radiation_steps=radiation_steps,
        rf_phase_deg=rf_phase_deg,
        rf_voltage_scale=rf_voltage_scale,
        lattice_data_path=lattice_data_path,
    )


def _track_one(lattice, coordinates: Sequence[float]) -> np.ndarray:
    bunch = rft.Bunch6d(
        rft.electronmass,
        0.0,
        -1.0,
        np.asarray(coordinates, dtype=float),
    )
    tracked = lattice.track(bunch)
    if tracked.size() != 1:
        raise RuntimeError("Synchronous-orbit probe was lost in one turn")
    return np.asarray(tracked.get_phase_space()[0], dtype=float)


def _rf_periodic_residual(tracked: np.ndarray, initial: np.ndarray) -> np.ndarray:
    """Return one-turn residual with ct reduced modulo the RF wavelength."""
    residual = tracked - initial
    residual[4] = (residual[4] + 0.5 * RF_WAVELENGTH_MM) % RF_WAVELENGTH_MM - 0.5 * RF_WAVELENGTH_MM
    return residual


def find_synchronous_orbit(
    lattice,
    *,
    initial_coordinates: Sequence[float] | None = None,
    max_iterations: int = 20,
) -> SynchronousOrbitResult:
    """Find a deterministic radiation/RF 6D fixed point by damped Newton steps."""
    z = (
        np.array((0.0, 0.0, 0.0, 0.0, 0.0, NOMINAL_MOMENTUM_MEV_C))
        if initial_coordinates is None
        else np.asarray(initial_coordinates, dtype=float).copy()
    )
    if z.shape != (6,):
        raise ValueError("initial_coordinates must contain [x, px, y, py, ct, p]")
    steps = np.array((1e-4, 1e-7, 1e-4, 1e-7, 1e-3, 1e-4))
    scales = np.array((1e-3, 1e-5, 1e-3, 1e-5, 1e-2, 1e-3))

    for iteration in range(max_iterations + 1):
        residual = _rf_periodic_residual(_track_one(lattice, z), z)
        merit = float(np.max(np.abs(residual / scales)))
        if merit < 1.0:
            return SynchronousOrbitResult(z.copy(), residual, iteration)
        if iteration == max_iterations:
            break
        jacobian = np.empty((6, 6))
        for column, step in enumerate(steps):
            offset = np.zeros(6)
            offset[column] = step
            plus = _rf_periodic_residual(_track_one(lattice, z + offset), z + offset)
            minus = _rf_periodic_residual(_track_one(lattice, z - offset), z - offset)
            jacobian[:, column] = (plus - minus) / (2.0 * step)
        update = np.linalg.lstsq(jacobian, residual, rcond=None)[0]
        accepted = False
        # The initial RF-off transverse orbit can be far from the 6D
        # radiation/RF fixed point.  Keep trying smaller steps rather than
        # declaring a physically valid fixed point absent after one lost probe.
        for damping in (
            1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625,
            0.0078125, 0.00390625, 0.001953125, 0.0009765625,
        ):
            candidate = z - damping * update
            candidate[4] = (
                (candidate[4] + 0.5 * RF_WAVELENGTH_MM) % RF_WAVELENGTH_MM
                - 0.5 * RF_WAVELENGTH_MM
            )
            try:
                candidate_residual = _rf_periodic_residual(
                    _track_one(lattice, candidate), candidate
                )
            except RuntimeError:
                continue
            if np.max(np.abs(candidate_residual / scales)) < merit:
                z = candidate
                accepted = True
                break
        if not accepted:
            raise RuntimeError("Synchronous-orbit Newton iteration did not reduce residual")
    raise RuntimeError("Synchronous orbit did not converge")


def gaussian_bunch(
    synchronous_orbit: SynchronousOrbitResult,
    *,
    particles: int,
    emit_x_m_rad: float,
    emit_y_m_rad: float,
    beta_x_m: float,
    beta_y_m: float,
    sigma_delta: float,
    sigma_ct_mm: float,
    seed: int = 1,
):
    """Create an uncoupled Gaussian bunch around the synchronous orbit.

    Input emittances are geometric.  The initial bunch is a controllable seed;
    quantum radiation, not this initial choice, must set the long-turn limit.
    """
    if particles < 2:
        raise ValueError("particles must be at least 2")
    if min(emit_x_m_rad, emit_y_m_rad, beta_x_m, beta_y_m, sigma_delta, sigma_ct_mm) <= 0:
        raise ValueError("beam sizes and emittances must be positive")
    rng = np.random.default_rng(seed)
    samples = np.zeros((particles, 6))
    # RF-Track ATF coordinates are [mm, mrad, mm, mrad, ct(mm), p(MeV/c)].
    samples[:, 0] = rng.normal(0.0, np.sqrt(emit_x_m_rad * beta_x_m) * 1e3, particles)
    samples[:, 1] = rng.normal(0.0, np.sqrt(emit_x_m_rad / beta_x_m) * 1e3, particles)
    samples[:, 2] = rng.normal(0.0, np.sqrt(emit_y_m_rad * beta_y_m) * 1e3, particles)
    samples[:, 3] = rng.normal(0.0, np.sqrt(emit_y_m_rad / beta_y_m) * 1e3, particles)
    samples[:, 4] = rng.normal(0.0, sigma_ct_mm, particles)
    samples[:, 5] = rng.normal(
        0.0, NOMINAL_MOMENTUM_MEV_C * sigma_delta, particles
    )
    samples += synchronous_orbit.coordinates
    return rft.Bunch6d(rft.electronmass, 0.0, -1.0, samples)


def emittance_from_bunch(bunch, *, turns: int) -> EmittanceResult:
    """Calculate projected and 4D normal-mode geometric emittances."""
    phase = np.asarray(bunch.get_phase_space(), dtype=float)
    if phase.shape[0] < 10:
        raise RuntimeError("Too few survived particles for an emittance estimate")
    # Transform [x(mm), px(mrad), y(mm), py(mrad)] to [m, rad, m, rad].
    transverse = phase[:, :4].copy()
    transverse *= 1e-3
    covariance = np.cov(transverse, rowvar=False, ddof=1)
    projected_x = float(np.sqrt(max(np.linalg.det(covariance[:2, :2]), 0.0)))
    projected_y = float(np.sqrt(max(np.linalg.det(covariance[2:, 2:]), 0.0)))
    symplectic = np.array(
        ((0.0, 1.0, 0.0, 0.0), (-1.0, 0.0, 0.0, 0.0),
         (0.0, 0.0, 0.0, 1.0), (0.0, 0.0, -1.0, 0.0))
    )
    values = np.linalg.eigvals(symplectic @ covariance)
    eigen = np.sort(np.abs(np.imag(values)))[::2]
    return EmittanceResult(
        turns=turns,
        survived_particles=int(phase.shape[0]),
        projected_x_m_rad=projected_x,
        projected_y_m_rad=projected_y,
        eigen_emittances_m_rad=(float(eigen[0]), float(eigen[1])),
        covariance_m_rad=covariance,
    )


def _native_to_transverse_covariance(covariance: np.ndarray) -> np.ndarray:
    """Convert [mm, mrad, mm, mrad] covariance to [m, rad, m, rad]."""
    return np.asarray(covariance[:4, :4], dtype=float) * 1e-6


def _emittances_from_transverse_covariance(covariance: np.ndarray):
    projected_x = float(np.sqrt(max(np.linalg.det(covariance[:2, :2]), 0.0)))
    projected_y = float(np.sqrt(max(np.linalg.det(covariance[2:, 2:]), 0.0)))
    symplectic = np.array(
        ((0.0, 1.0, 0.0, 0.0), (-1.0, 0.0, 0.0, 0.0),
         (0.0, 0.0, 0.0, 1.0), (0.0, 0.0, -1.0, 0.0))
    )
    eigenvalues = np.linalg.eigvals(symplectic @ covariance)
    eigen = np.sort(np.abs(np.imag(eigenvalues)))[::2]
    return projected_x, projected_y, (float(eigen[0]), float(eigen[1]))


def one_turn_map(lattice, synchronous_orbit: SynchronousOrbitResult) -> np.ndarray:
    """Numerically linearize the deterministic radiation/RF map at its fixed point."""
    z = synchronous_orbit.coordinates
    steps = np.array((1e-4, 1e-5, 1e-4, 1e-5, 1e-3, 1e-4))
    matrix = np.empty((6, 6))
    for column, step in enumerate(steps):
        offset = np.zeros(6)
        offset[column] = step
        plus = _rf_periodic_residual(_track_one(lattice, z + offset), z + offset)
        minus = _rf_periodic_residual(_track_one(lattice, z - offset), z - offset)
        matrix[:, column] = (plus - minus) / (2.0 * step)
    return matrix + np.eye(6)


def quantum_diffusion(
    lattice,
    synchronous_orbit: SynchronousOrbitResult,
    *,
    particles: int,
) -> np.ndarray:
    """Estimate one-turn quantum diffusion in native RF-Track coordinates."""
    if particles < 10:
        raise ValueError("particles must be at least 10 for quantum diffusion")
    coordinates = np.repeat(
        synchronous_orbit.coordinates[None, :], particles, axis=0
    )
    bunch = rft.Bunch6d(rft.electronmass, 0.0, -1.0, coordinates)
    tracked = lattice.track(bunch)
    if tracked.size() < 10:
        raise RuntimeError("Too few particles survived the quantum-diffusion turn")
    residual = np.asarray(tracked.get_phase_space(), dtype=float)
    residual -= synchronous_orbit.coordinates
    residual[:, 4] = (
        (residual[:, 4] + 0.5 * RF_WAVELENGTH_MM) % RF_WAVELENGTH_MM
        - 0.5 * RF_WAVELENGTH_MM
    )
    return np.cov(residual, rowvar=False, ddof=1)


def equilibrium_emittance_for_lattices(
    deterministic,
    quantum,
    *,
    quantum_particles: int = 64,
    initial_coordinates: Sequence[float] | None = None,
) -> EquilibriumEmittanceResult:
    """Solve equilibrium for matched deterministic and quantum lattices.

    A deterministic radiation/RF lattice supplies the synchronous orbit and
    damping map ``M``.  A separate quantum-radiation one-turn ensemble supplies
    ``Q``.  The discrete Lyapunov equation ``Sigma=M Sigma M^T+Q`` then gives
    the long-turn covariance directly.
    """
    try:
        from scipy.linalg import solve_discrete_lyapunov
    except ImportError as error:  # pragma: no cover - dependency is in rftrack-env
        raise RuntimeError("equilibrium_emittance requires scipy") from error

    synchronous = find_synchronous_orbit(
        deterministic,
        initial_coordinates=initial_coordinates,
        max_iterations=50,
    )
    matrix = one_turn_map(deterministic, synchronous)
    spectral_radius = float(np.max(np.abs(np.linalg.eigvals(matrix))))
    if spectral_radius >= 1.0:
        raise RuntimeError(
            f"Radiation/RF one-turn map is not stable (spectral radius={spectral_radius:.8f})"
        )
    diffusion = quantum_diffusion(
        quantum, synchronous, particles=quantum_particles
    )
    covariance = solve_discrete_lyapunov(matrix, diffusion)
    transverse = _native_to_transverse_covariance(covariance)
    projected_x, projected_y, eigen = _emittances_from_transverse_covariance(
        transverse
    )
    return EquilibriumEmittanceResult(
        synchronous_orbit=synchronous,
        spectral_radius=spectral_radius,
        damping_turns_one_over_e=float(-1.0 / np.log(spectral_radius)),
        quantum_particles=quantum_particles,
        projected_x_m_rad=projected_x,
        projected_y_m_rad=projected_y,
        eigen_emittances_m_rad=eigen,
        one_turn_map=matrix,
        diffusion_native=diffusion,
        covariance_native=covariance,
    )


def transverse_envelope_for_lattices(
    deterministic,
    quantum,
    *,
    quantum_particles: int = 64,
    initial_coordinates: Sequence[float] | None = None,
    reference_coordinates: Sequence[float] | None = None,
) -> TransverseEnvelopeResult:
    """Evaluate the stable 4D transverse radiation envelope.

    The input lattices must use identical deterministic and quantum radiation
    settings.  A longitudinally unstable full 6D RF map is deliberately not
    rejected here; callers must label this as a transverse diagnostic.
    """
    try:
        from scipy.linalg import solve_discrete_lyapunov
    except ImportError as error:  # pragma: no cover
        raise RuntimeError("transverse_envelope requires scipy") from error

    if reference_coordinates is None:
        synchronous = find_synchronous_orbit(
            deterministic,
            initial_coordinates=initial_coordinates,
            max_iterations=50,
        )
    else:
        coordinates = np.asarray(reference_coordinates, dtype=float).copy()
        if coordinates.shape != (6,):
            raise ValueError("reference_coordinates must contain [x, px, y, py, ct, p]")
        synchronous = SynchronousOrbitResult(
            coordinates=coordinates,
            residual=np.full(6, np.nan),
            iterations=0,
        )
    matrix = one_turn_map(deterministic, synchronous)[:4, :4]
    spectral_radius = float(np.max(np.abs(np.linalg.eigvals(matrix))))
    if spectral_radius >= 1.0:
        raise RuntimeError(
            "Transverse radiation one-turn map is not stable "
            f"(spectral radius={spectral_radius:.8f})"
        )
    diffusion = quantum_diffusion(
        quantum, synchronous, particles=quantum_particles
    )[:4, :4]
    covariance = solve_discrete_lyapunov(matrix, diffusion)
    transverse = _native_to_transverse_covariance(covariance)
    projected_x, projected_y, eigen = _emittances_from_transverse_covariance(
        transverse
    )
    return TransverseEnvelopeResult(
        synchronous_orbit=synchronous,
        spectral_radius=spectral_radius,
        damping_turns_one_over_e=float(-1.0 / np.log(spectral_radius)),
        quantum_particles=quantum_particles,
        projected_x_m_rad=projected_x,
        projected_y_m_rad=projected_y,
        eigen_emittances_m_rad=eigen,
        one_turn_map=matrix,
        diffusion_native=diffusion,
        covariance_native=covariance,
    )


def equilibrium_emittance(
    *,
    quantum_particles: int = 64,
    rf_phase_deg: float | None = None,
    rf_voltage_scale: float = 1.0,
    lattice_data_path: str | Path | None = None,
) -> EquilibriumEmittanceResult:
    """Solve the nominal ATF DR radiation/quantum equilibrium envelope."""
    return equilibrium_emittance_for_lattices(
        build_equilibrium_lattice(
            quantum=False, rf_phase_deg=rf_phase_deg,
            rf_voltage_scale=rf_voltage_scale,
            lattice_data_path=lattice_data_path,
        ),
        build_equilibrium_lattice(
            quantum=True, rf_phase_deg=rf_phase_deg,
            rf_voltage_scale=rf_voltage_scale,
            lattice_data_path=lattice_data_path,
        ),
        quantum_particles=quantum_particles,
    )


def track_emittance(
    lattice,
    bunch,
    *,
    turns: int,
    sample_every: int | None = None,
) -> tuple[object, list[EmittanceResult]]:
    """Track a quantum-radiation bunch and return covariance samples."""
    if turns < 1:
        raise ValueError("turns must be positive")
    every = turns if sample_every is None else int(sample_every)
    if every < 1:
        raise ValueError("sample_every must be positive")
    samples = []
    tracked = bunch
    for turn in range(1, turns + 1):
        tracked = lattice.track(tracked)
        if tracked.size() < 10:
            raise RuntimeError(f"Only {tracked.size()} particles survived at turn {turn}")
        if turn % every == 0 or turn == turns:
            samples.append(emittance_from_bunch(tracked, turns=turn))
    return tracked, samples
