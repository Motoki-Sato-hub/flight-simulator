"""Periodic closed-orbit and response correction for the ATF DR RF-Track model.

The SAD reference implementation (``operation/lib/cod.n``) evaluates the
closed orbit around the ring, builds steerer responses and solves a constrained
least-squares problem.  This module provides the corresponding model-side
operations without requiring SAD.

Coordinates returned by RF-Track are in mm and mrad.  Corrector values are in
caller units; ``actuator_scale`` converts one caller unit to the value passed to
``RF_Track.Corrector.set_strength``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .ATF_DR_RFTrack_lattice import NOMINAL_MOMENTUM_MEV_C


class ClosedOrbitError(RuntimeError):
    """Raised when the periodic fixed point cannot be found."""


class ParticleLostError(ClosedOrbitError):
    """Raised when RF-Track loses the probe particle."""


@dataclass(frozen=True)
class ClosedOrbitResult:
    initial_coordinates: np.ndarray
    bpm_names: tuple[str, ...]
    bpm_positions: np.ndarray
    momentum_mev_c: float
    iterations: int
    residual_norm: float

    @property
    def x(self) -> np.ndarray:
        return self.bpm_positions[:, 0]

    @property
    def y(self) -> np.ndarray:
        return self.bpm_positions[:, 1]


@dataclass(frozen=True)
class DispersionResult:
    bpm_names: tuple[str, ...]
    values: np.ndarray
    relative_momentum_step: float

    @property
    def x(self) -> np.ndarray:
        return self.values[:, 0]

    @property
    def y(self) -> np.ndarray:
        return self.values[:, 1]


@dataclass(frozen=True)
class OrbitResponseResult:
    plane: str
    bpm_names: tuple[str, ...]
    corrector_names: tuple[str, ...]
    matrix: np.ndarray
    baseline_orbit: np.ndarray
    perturbation: float
    actuator_scale: float


@dataclass(frozen=True)
class DispersionResponseResult:
    plane: str
    bpm_names: tuple[str, ...]
    corrector_names: tuple[str, ...]
    dispersion_matrix: np.ndarray
    orbit_matrix: np.ndarray
    baseline_dispersion: np.ndarray
    baseline_orbit: np.ndarray
    relative_momentum_step: float
    perturbation: float
    actuator_scale: float


@dataclass(frozen=True)
class OrbitCorrectionResult:
    plane: str
    bpm_names: tuple[str, ...]
    corrector_names: tuple[str, ...]
    delta_correctors: np.ndarray
    measured_orbit: np.ndarray
    target_orbit: np.ndarray
    predicted_orbit: np.ndarray
    selected_correctors: tuple[str, ...]
    singular_values: np.ndarray
    rank: int

    @property
    def rms_before(self) -> float:
        residual = self.measured_orbit - self.target_orbit
        return float(np.sqrt(np.mean(residual**2)))

    @property
    def rms_predicted(self) -> float:
        residual = self.predicted_orbit - self.target_orbit
        return float(np.sqrt(np.mean(residual**2)))


@dataclass(frozen=True)
class DispersionCorrectionResult:
    plane: str
    bpm_names: tuple[str, ...]
    corrector_names: tuple[str, ...]
    delta_correctors: np.ndarray
    measured_dispersion: np.ndarray
    target_dispersion: np.ndarray
    predicted_dispersion: np.ndarray
    measured_orbit: np.ndarray
    target_orbit: np.ndarray
    predicted_orbit: np.ndarray
    selected_correctors: tuple[str, ...]
    singular_values: np.ndarray
    rank: int
    relative_momentum_step: float

    @staticmethod
    def _rms(values: np.ndarray) -> float:
        return float(np.sqrt(np.mean(values**2)))

    @property
    def dispersion_rms_before(self) -> float:
        return self._rms(self.measured_dispersion - self.target_dispersion)

    @property
    def dispersion_rms_predicted(self) -> float:
        return self._rms(self.predicted_dispersion - self.target_dispersion)

    @property
    def orbit_rms_before(self) -> float:
        return self._rms(self.measured_orbit - self.target_orbit)

    @property
    def orbit_rms_predicted(self) -> float:
        return self._rms(self.predicted_orbit - self.target_orbit)


@dataclass(frozen=True)
class CouplingResponseResult:
    """Vertical response to horizontal probe kicks versus skew K1L.

    Rows are ordered by probe corrector, then BPM.  The signal is therefore
    the same observable as SAD ``skewcor.n``: the vertical closed orbit seen
    after each of two (or more) horizontal-steerer changes.
    """

    bpm_names: tuple[str, ...]
    probe_corrector_names: tuple[str, ...]
    skew_corrector_names: tuple[str, ...]
    matrix: np.ndarray
    baseline_signal: np.ndarray
    probe_perturbation: float
    skew_perturbation: float
    actuator_scale: float
    skew_actuator_scale: float


@dataclass(frozen=True)
class CouplingCorrectionResult:
    """A proposed skew-quadrupole correction of the coupling observable."""

    bpm_names: tuple[str, ...]
    probe_corrector_names: tuple[str, ...]
    skew_corrector_names: tuple[str, ...]
    probe_perturbation: float
    delta_skew_correctors: np.ndarray
    measured_signal: np.ndarray
    target_signal: np.ndarray
    predicted_signal: np.ndarray
    selected_skew_correctors: tuple[str, ...]
    singular_values: np.ndarray
    rank: int

    @staticmethod
    def _rms(values: np.ndarray) -> float:
        return float(np.sqrt(np.mean(values**2)))

    @property
    def rms_before(self) -> float:
        return self._rms(self.measured_signal - self.target_signal)

    @property
    def rms_predicted(self) -> float:
        return self._rms(self.predicted_signal - self.target_signal)


class ATFDRRingCorrection:
    """Closed-orbit solver and corrector-response model for an RF-Track ring."""

    def __init__(
        self,
        lattice,
        momentum_mev_c: float = NOMINAL_MOMENTUM_MEV_C,
        charge: float = -1.0,
        *,
        population: float = 0.0,
        actuator_scale: float = 1.0,
        skew_actuator_scale: float = 1.0,
    ):
        if charge == 0:
            raise ValueError("charge must be non-zero")
        if actuator_scale == 0:
            raise ValueError("actuator_scale must be non-zero")
        if skew_actuator_scale == 0:
            raise ValueError("skew_actuator_scale must be non-zero")
        self.lattice = lattice
        self.momentum_mev_c = float(momentum_mev_c)
        self.charge = float(charge)
        self.population = float(population)
        self.actuator_scale = float(actuator_scale)
        self.skew_actuator_scale = float(skew_actuator_scale)
        self.bpm_names = tuple(bpm.get_name() for bpm in lattice.get_bpms())

    @staticmethod
    def _plane_index(plane: str) -> int:
        normalized = plane.lower()
        if normalized == "x":
            return 0
        if normalized == "y":
            return 1
        raise ValueError("plane must be 'x' or 'y'")

    def _single_element(self, name: str):
        element = self.lattice[name]
        if isinstance(element, list):
            if len(element) != 1:
                raise ValueError(f"Expected one element named {name}, got {len(element)}")
            return element[0]
        return element

    def _track_coordinates(
        self,
        coordinates: Sequence[float],
        momentum_mev_c: float,
    ) -> np.ndarray:
        import RF_Track as rft

        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.shape != (4,):
            raise ValueError("coordinates must contain [x, px, y, py]")
        phase_space = np.concatenate((coordinates, [0.0, float(momentum_mev_c)]))
        bunch = rft.Bunch6d(
            rft.electronmass,
            self.population,
            self.charge,
            phase_space,
        )
        tracked = self.lattice.track(bunch)
        if tracked.size() != 1:
            raise ParticleLostError(
                "RF-Track lost the closed-orbit probe particle in one turn"
            )
        return np.asarray(tracked.get_phase_space()[0, :4], dtype=float)

    def _bpm_positions(
        self, bpm_names: Sequence[str] | None
    ) -> tuple[tuple[str, ...], np.ndarray]:
        names = self.bpm_names if bpm_names is None else tuple(bpm_names)
        readings = []
        for name in names:
            if name not in self.bpm_names:
                raise ValueError(f"Unknown BPM: {name}")
            reading = np.asarray(self._single_element(name).get_reading(), dtype=float)
            readings.append(reading[:2])
        return names, np.asarray(readings, dtype=float)

    def find_closed_orbit(
        self,
        *,
        momentum_mev_c: float | None = None,
        initial_coordinates: Sequence[float] | None = None,
        bpm_names: Sequence[str] | None = None,
        finite_difference_step: float = 1e-5,
        tolerance: float = 1e-8,
        max_iterations: int = 8,
    ) -> ClosedOrbitResult:
        """Solve ``one_turn(z) - z = 0`` with a numerical Newton method."""
        momentum = (
            self.momentum_mev_c
            if momentum_mev_c is None
            else float(momentum_mev_c)
        )
        z = (
            np.zeros(4, dtype=float)
            if initial_coordinates is None
            else np.asarray(initial_coordinates, dtype=float).copy()
        )
        if z.shape != (4,):
            raise ValueError("initial_coordinates must contain [x, px, y, py]")
        if finite_difference_step <= 0:
            raise ValueError("finite_difference_step must be positive")

        residual_norm = float("inf")
        iterations = 0
        for iterations in range(max_iterations + 1):
            residual = self._track_coordinates(z, momentum) - z
            residual_norm = float(np.linalg.norm(residual, ord=np.inf))
            if residual_norm <= tolerance:
                # Populate BPM readings with the converged orbit after any
                # finite-difference probe tracks performed in prior iterations.
                self._track_coordinates(z, momentum)
                names, positions = self._bpm_positions(bpm_names)
                return ClosedOrbitResult(
                    initial_coordinates=z.copy(),
                    bpm_names=names,
                    bpm_positions=positions,
                    momentum_mev_c=momentum,
                    iterations=iterations,
                    residual_norm=residual_norm,
                )
            if iterations == max_iterations:
                break

            jacobian = np.empty((4, 4), dtype=float)
            for column in range(4):
                offset = np.zeros(4, dtype=float)
                offset[column] = finite_difference_step
                plus = self._track_coordinates(z + offset, momentum)
                minus = self._track_coordinates(z - offset, momentum)
                jacobian[:, column] = (
                    (plus - minus) / (2.0 * finite_difference_step)
                    - np.eye(4)[:, column]
                )
            try:
                step = np.linalg.solve(jacobian, residual)
            except np.linalg.LinAlgError:
                step = np.linalg.lstsq(jacobian, residual, rcond=None)[0]

            accepted = False
            for damping in (1.0, 0.5, 0.25, 0.125, 0.0625):
                candidate = z - damping * step
                try:
                    candidate_residual = (
                        self._track_coordinates(candidate, momentum) - candidate
                    )
                except ParticleLostError:
                    continue
                if np.linalg.norm(candidate_residual, ord=np.inf) < residual_norm:
                    z = candidate
                    accepted = True
                    break
            if not accepted:
                raise ClosedOrbitError(
                    "Closed-orbit Newton iteration did not reduce the residual"
                )

        raise ClosedOrbitError(
            f"Closed orbit did not converge after {max_iterations} iterations; "
            f"residual={residual_norm:.3e}"
        )

    def measure_dispersion(
        self,
        *,
        relative_momentum_step: float = 1e-3,
        bpm_names: Sequence[str] | None = None,
    ) -> DispersionResult:
        """Return periodic dispersion in mm per unit relative momentum error."""
        if relative_momentum_step <= 0:
            raise ValueError("relative_momentum_step must be positive")
        reference = self.find_closed_orbit(bpm_names=bpm_names)
        plus = self.find_closed_orbit(
            momentum_mev_c=self.momentum_mev_c * (1.0 + relative_momentum_step),
            initial_coordinates=reference.initial_coordinates,
            bpm_names=bpm_names,
        )
        minus = self.find_closed_orbit(
            momentum_mev_c=self.momentum_mev_c * (1.0 - relative_momentum_step),
            initial_coordinates=reference.initial_coordinates,
            bpm_names=bpm_names,
        )
        self.find_closed_orbit(
            initial_coordinates=reference.initial_coordinates,
            bpm_names=bpm_names,
        )
        return DispersionResult(
            bpm_names=reference.bpm_names,
            values=(plus.bpm_positions - minus.bpm_positions)
            / (2.0 * relative_momentum_step),
            relative_momentum_step=relative_momentum_step,
        )

    def get_corrector_names(self, plane: str) -> tuple[str, ...]:
        """Return the SAD ZH or ZV correctors in lattice order."""
        prefix = "ZH" if self._plane_index(plane) == 0 else "ZV"
        return tuple(
            corrector.get_name()
            for corrector in self.lattice.get_correctors()
            if corrector.get_name().upper().startswith(prefix)
        )

    @property
    def _p_over_q(self) -> float:
        return self.momentum_mev_c / self.charge

    def get_skew_corrector_names(self) -> tuple[str, ...]:
        """Return the thin SD1R/SF1R skew-K1L actuators in lattice order."""
        return tuple(
            element.get_name()
            for element in self.lattice["*"]
            if element.get_name().endswith("$SKEW")
            and element.get_name().split(".", 1)[0] in {"SD1R", "SF1R"}
        )

    def _get_skew_strength(self, name: str) -> float:
        element = self._single_element(name)
        strengths = np.asarray(element.get_KnL(self._p_over_q), dtype=complex)
        if strengths.size < 2:
            raise ValueError(f"{name} does not expose a quadrupole K1L component")
        return float(np.imag(strengths.reshape(-1)[1]))

    def _set_skew_strength(self, name: str, value: float) -> None:
        element = self._single_element(name)
        strengths = np.asarray(
            element.get_KnL(self._p_over_q), dtype=complex
        ).copy()
        if strengths.size < 2:
            raise ValueError(f"{name} does not expose a quadrupole K1L component")
        flat_strengths = strengths.reshape(-1)
        flat_strengths[1] = flat_strengths[1].real + 1j * float(value)
        element.set_KnL(self._p_over_q, strengths)

    def get_skew_strength(self, name: str) -> float:
        """Return one SD1R/SF1R skew actuator's normalized integrated K1L."""
        if name not in self.get_skew_corrector_names():
            raise ValueError(f"Not an ATF DR skew corrector: {name}")
        return self._get_skew_strength(name)

    def set_skew_strength(self, name: str, value: float) -> None:
        """Set one SD1R/SF1R skew actuator's normalized integrated K1L."""
        if name not in self.get_skew_corrector_names():
            raise ValueError(f"Not an ATF DR skew corrector: {name}")
        self._set_skew_strength(name, value)

    def compute_orbit_response(
        self,
        plane: str,
        *,
        corrector_names: Sequence[str] | None = None,
        bpm_names: Sequence[str] | None = None,
        perturbation: float = 1e-5,
    ) -> OrbitResponseResult:
        """Build a periodic COD response matrix by central differences."""
        plane_index = self._plane_index(plane)
        names = (
            self.get_corrector_names(plane)
            if corrector_names is None
            else tuple(corrector_names)
        )
        if perturbation <= 0:
            raise ValueError("perturbation must be positive")
        if not names:
            raise ValueError("No correctors selected")

        baseline = self.find_closed_orbit(bpm_names=bpm_names)
        original_strengths: dict[str, np.ndarray] = {}
        for name in names:
            element = self._single_element(name)
            if name not in self.get_corrector_names(plane):
                raise ValueError(f"{name} is not a {plane}-plane ATF DR corrector")
            original_strengths[name] = np.asarray(
                element.get_strength(), dtype=float
            ).copy()

        matrix = np.empty((len(baseline.bpm_names), len(names)), dtype=float)
        try:
            for column, name in enumerate(names):
                element = self._single_element(name)
                strength = original_strengths[name]
                plus_strength = strength.copy()
                minus_strength = strength.copy()
                plus_strength[plane_index] += self.actuator_scale * perturbation
                minus_strength[plane_index] -= self.actuator_scale * perturbation

                element.set_strength(*plus_strength)
                plus = self.find_closed_orbit(
                    initial_coordinates=baseline.initial_coordinates,
                    bpm_names=baseline.bpm_names,
                )
                element.set_strength(*minus_strength)
                minus = self.find_closed_orbit(
                    initial_coordinates=baseline.initial_coordinates,
                    bpm_names=baseline.bpm_names,
                )
                element.set_strength(*strength)
                matrix[:, column] = (
                    plus.bpm_positions[:, plane_index]
                    - minus.bpm_positions[:, plane_index]
                ) / (2.0 * perturbation)
        finally:
            for name, strength in original_strengths.items():
                self._single_element(name).set_strength(*strength)
            self.find_closed_orbit(
                initial_coordinates=baseline.initial_coordinates,
                bpm_names=baseline.bpm_names,
            )

        return OrbitResponseResult(
            plane=plane.lower(),
            bpm_names=baseline.bpm_names,
            corrector_names=names,
            matrix=matrix,
            baseline_orbit=baseline.bpm_positions[:, plane_index].copy(),
            perturbation=perturbation,
            actuator_scale=self.actuator_scale,
        )

    def compute_dispersion_response(
        self,
        plane: str,
        *,
        corrector_names: Sequence[str] | None = None,
        bpm_names: Sequence[str] | None = None,
        relative_momentum_step: float = 1e-3,
        perturbation: float = 1e-5,
    ) -> DispersionResponseResult:
        """Build steerer responses for both periodic dispersion and COD."""
        plane_index = self._plane_index(plane)
        names = (
            self.get_corrector_names(plane)
            if corrector_names is None
            else tuple(corrector_names)
        )
        if perturbation <= 0:
            raise ValueError("perturbation must be positive")
        if relative_momentum_step <= 0:
            raise ValueError("relative_momentum_step must be positive")
        if not names:
            raise ValueError("No correctors selected")

        baseline_orbit = self.find_closed_orbit(bpm_names=bpm_names)
        baseline_dispersion = self.measure_dispersion(
            relative_momentum_step=relative_momentum_step,
            bpm_names=baseline_orbit.bpm_names,
        )
        allowed_names = self.get_corrector_names(plane)
        original_strengths: dict[str, np.ndarray] = {}
        for name in names:
            element = self._single_element(name)
            if name not in allowed_names:
                raise ValueError(f"{name} is not a {plane}-plane ATF DR corrector")
            original_strengths[name] = np.asarray(
                element.get_strength(), dtype=float
            ).copy()

        dispersion_matrix = np.empty(
            (len(baseline_orbit.bpm_names), len(names)), dtype=float
        )
        orbit_matrix = np.empty_like(dispersion_matrix)
        try:
            for column, name in enumerate(names):
                element = self._single_element(name)
                strength = original_strengths[name]
                plus_strength = strength.copy()
                minus_strength = strength.copy()
                plus_strength[plane_index] += self.actuator_scale * perturbation
                minus_strength[plane_index] -= self.actuator_scale * perturbation

                element.set_strength(*plus_strength)
                plus_orbit = self.find_closed_orbit(
                    initial_coordinates=baseline_orbit.initial_coordinates,
                    bpm_names=baseline_orbit.bpm_names,
                )
                plus_dispersion = self.measure_dispersion(
                    relative_momentum_step=relative_momentum_step,
                    bpm_names=baseline_orbit.bpm_names,
                )

                element.set_strength(*minus_strength)
                minus_orbit = self.find_closed_orbit(
                    initial_coordinates=baseline_orbit.initial_coordinates,
                    bpm_names=baseline_orbit.bpm_names,
                )
                minus_dispersion = self.measure_dispersion(
                    relative_momentum_step=relative_momentum_step,
                    bpm_names=baseline_orbit.bpm_names,
                )
                element.set_strength(*strength)

                dispersion_matrix[:, column] = (
                    plus_dispersion.values[:, plane_index]
                    - minus_dispersion.values[:, plane_index]
                ) / (2.0 * perturbation)
                orbit_matrix[:, column] = (
                    plus_orbit.bpm_positions[:, plane_index]
                    - minus_orbit.bpm_positions[:, plane_index]
                ) / (2.0 * perturbation)
        finally:
            for name, strength in original_strengths.items():
                self._single_element(name).set_strength(*strength)
            self.find_closed_orbit(
                initial_coordinates=baseline_orbit.initial_coordinates,
                bpm_names=baseline_orbit.bpm_names,
            )

        return DispersionResponseResult(
            plane=plane.lower(),
            bpm_names=baseline_orbit.bpm_names,
            corrector_names=names,
            dispersion_matrix=dispersion_matrix,
            orbit_matrix=orbit_matrix,
            baseline_dispersion=baseline_dispersion.values[:, plane_index].copy(),
            baseline_orbit=baseline_orbit.bpm_positions[:, plane_index].copy(),
            relative_momentum_step=relative_momentum_step,
            perturbation=perturbation,
            actuator_scale=self.actuator_scale,
        )

    def _measure_vertical_response_to_horizontal_probes(
        self,
        probe_corrector_names: Sequence[str],
        *,
        bpm_names: Sequence[str] | None,
        probe_perturbation: float,
        initial_coordinates: Sequence[float] | None = None,
    ) -> tuple[tuple[str, ...], np.ndarray]:
        """Measure ``d y_BPM / d ZH`` for the requested horizontal probes."""
        if probe_perturbation <= 0:
            raise ValueError("probe_perturbation must be positive")
        names = tuple(probe_corrector_names)
        allowed = self.get_corrector_names("x")
        if not names:
            raise ValueError("At least one horizontal probe corrector is required")
        unknown = [name for name in names if name not in allowed]
        if unknown:
            raise ValueError(f"Not horizontal ATF DR correctors: {unknown}")

        baseline = self.find_closed_orbit(
            initial_coordinates=initial_coordinates,
            bpm_names=bpm_names,
        )
        original_strengths = {
            name: np.asarray(self._single_element(name).get_strength(), dtype=float)
            .copy()
            for name in names
        }
        signal = np.empty((len(names), len(baseline.bpm_names)), dtype=float)
        try:
            for index, name in enumerate(names):
                element = self._single_element(name)
                strength = original_strengths[name]
                plus_strength = strength.copy()
                minus_strength = strength.copy()
                plus_strength[0] += self.actuator_scale * probe_perturbation
                minus_strength[0] -= self.actuator_scale * probe_perturbation

                element.set_strength(*plus_strength)
                plus = self.find_closed_orbit(
                    initial_coordinates=baseline.initial_coordinates,
                    bpm_names=baseline.bpm_names,
                )
                element.set_strength(*minus_strength)
                minus = self.find_closed_orbit(
                    initial_coordinates=baseline.initial_coordinates,
                    bpm_names=baseline.bpm_names,
                )
                signal[index] = (plus.y - minus.y) / (2.0 * probe_perturbation)
                element.set_strength(*strength)
        finally:
            for name, strength in original_strengths.items():
                self._single_element(name).set_strength(*strength)
            self.find_closed_orbit(
                initial_coordinates=baseline.initial_coordinates,
                bpm_names=baseline.bpm_names,
            )
        return baseline.bpm_names, signal

    def compute_coupling_response(
        self,
        *,
        probe_corrector_names: Sequence[str] | None = None,
        skew_corrector_names: Sequence[str] | None = None,
        bpm_names: Sequence[str] | None = None,
        probe_perturbation: float = 1e-4,
        skew_perturbation: float = 1e-5,
    ) -> CouplingResponseResult:
        """Build the numerical counterpart of SAD ``CorrectSkewCoupling``.

        A horizontal steerer is changed by ``probe_perturbation`` and the
        resulting vertical periodic COD at every BPM is measured.  Each column
        is the finite-difference derivative of that observable with respect
        to one thin skew-quadrupole K1L actuator at an SD1R or SF1R location.
        """
        probes = (
            self.get_corrector_names("x")[:2]
            if probe_corrector_names is None
            else tuple(probe_corrector_names)
        )
        skews = (
            self.get_skew_corrector_names()
            if skew_corrector_names is None
            else tuple(skew_corrector_names)
        )
        available_skews = self.get_skew_corrector_names()
        unknown_skews = [name for name in skews if name not in available_skews]
        if unknown_skews:
            raise ValueError(f"Not ATF DR skew correctors: {unknown_skews}")
        if not skews:
            raise ValueError("No skew correctors selected")
        if skew_perturbation <= 0:
            raise ValueError("skew_perturbation must be positive")

        baseline_orbit = self.find_closed_orbit(bpm_names=bpm_names)
        bpm_names_out, baseline_signal = (
            self._measure_vertical_response_to_horizontal_probes(
                probes,
                bpm_names=baseline_orbit.bpm_names,
                probe_perturbation=probe_perturbation,
                initial_coordinates=baseline_orbit.initial_coordinates,
            )
        )
        original_strengths = {
            name: self._get_skew_strength(name) for name in skews
        }
        matrix = np.empty((baseline_signal.size, len(skews)), dtype=float)
        try:
            for column, name in enumerate(skews):
                original = original_strengths[name]
                self._set_skew_strength(
                    name, original + self.skew_actuator_scale * skew_perturbation
                )
                _, plus = self._measure_vertical_response_to_horizontal_probes(
                    probes,
                    bpm_names=bpm_names_out,
                    probe_perturbation=probe_perturbation,
                    initial_coordinates=baseline_orbit.initial_coordinates,
                )
                self._set_skew_strength(
                    name, original - self.skew_actuator_scale * skew_perturbation
                )
                _, minus = self._measure_vertical_response_to_horizontal_probes(
                    probes,
                    bpm_names=bpm_names_out,
                    probe_perturbation=probe_perturbation,
                    initial_coordinates=baseline_orbit.initial_coordinates,
                )
                self._set_skew_strength(name, original)
                matrix[:, column] = (plus - minus).reshape(-1) / (
                    2.0 * skew_perturbation
                )
        finally:
            for name, strength in original_strengths.items():
                self._set_skew_strength(name, strength)
            self.find_closed_orbit(
                initial_coordinates=baseline_orbit.initial_coordinates,
                bpm_names=bpm_names_out,
            )

        return CouplingResponseResult(
            bpm_names=bpm_names_out,
            probe_corrector_names=probes,
            skew_corrector_names=skews,
            matrix=matrix,
            baseline_signal=baseline_signal.reshape(-1).copy(),
            probe_perturbation=probe_perturbation,
            skew_perturbation=skew_perturbation,
            actuator_scale=self.actuator_scale,
            skew_actuator_scale=self.skew_actuator_scale,
        )

    @staticmethod
    def _svd_solve(
        matrix: np.ndarray,
        target: np.ndarray,
        rcond: float,
    ) -> tuple[np.ndarray, np.ndarray, int]:
        u, singular_values, vt = np.linalg.svd(matrix, full_matrices=False)
        if singular_values.size == 0:
            return np.zeros(matrix.shape[1]), singular_values, 0
        cutoff = rcond * singular_values[0]
        keep = singular_values > cutoff
        solution = np.zeros(matrix.shape[1], dtype=float)
        if np.any(keep):
            solution = vt[keep].T @ ((u[:, keep].T @ target) / singular_values[keep])
        return solution, singular_values, int(np.count_nonzero(keep))

    @classmethod
    def _greedy_columns(
        cls,
        matrix: np.ndarray,
        target: np.ndarray,
        max_columns: int,
        rcond: float,
    ) -> list[int]:
        """SAD-like forward selection of the steerers reducing residual most."""
        selected: list[int] = []
        remaining = list(range(matrix.shape[1]))
        previous_norm = float(np.linalg.norm(target))
        while remaining and len(selected) < max_columns:
            best_column = None
            best_norm = previous_norm
            for candidate in remaining:
                trial = selected + [candidate]
                solution, _, _ = cls._svd_solve(matrix[:, trial], target, rcond)
                norm = float(np.linalg.norm(target - matrix[:, trial] @ solution))
                if norm < best_norm:
                    best_norm = norm
                    best_column = candidate
            if best_column is None:
                break
            selected.append(best_column)
            remaining.remove(best_column)
            previous_norm = best_norm
        return selected

    def propose_orbit_correction(
        self,
        response: OrbitResponseResult,
        *,
        target_orbit: Sequence[float] | None = None,
        bpm_weights: Sequence[float] | None = None,
        rcond: float = 1e-4,
        gain: float = 1.0,
        max_correctors: int | None = None,
        max_abs_delta: float | None = None,
    ) -> OrbitCorrectionResult:
        """Solve for a COD correction without changing lattice strengths."""
        plane_index = self._plane_index(response.plane)
        measured = self.find_closed_orbit(
            bpm_names=response.bpm_names
        ).bpm_positions[:, plane_index]
        target = (
            np.zeros_like(measured)
            if target_orbit is None
            else np.asarray(target_orbit, dtype=float)
        )
        if target.shape != measured.shape:
            raise ValueError("target_orbit length must match response BPMs")
        weights = (
            np.ones_like(measured)
            if bpm_weights is None
            else np.asarray(bpm_weights, dtype=float)
        )
        if weights.shape != measured.shape or np.any(weights < 0):
            raise ValueError("bpm_weights must be non-negative and match response BPMs")
        if rcond < 0:
            raise ValueError("rcond must be non-negative")

        weighted_matrix = response.matrix * weights[:, None]
        weighted_target = (target - measured) * weights
        if max_correctors is None:
            selected_indices = list(range(response.matrix.shape[1]))
        else:
            if max_correctors <= 0:
                raise ValueError("max_correctors must be positive")
            selected_indices = self._greedy_columns(
                weighted_matrix,
                weighted_target,
                min(max_correctors, response.matrix.shape[1]),
                rcond,
            )
        if not selected_indices:
            raise ClosedOrbitError("No useful corrector columns were selected")

        selected_matrix = weighted_matrix[:, selected_indices]
        selected_delta, singular_values, rank = self._svd_solve(
            selected_matrix, weighted_target, rcond
        )
        selected_delta *= float(gain)
        if max_abs_delta is not None:
            if max_abs_delta <= 0:
                raise ValueError("max_abs_delta must be positive")
            selected_delta = np.clip(
                selected_delta, -float(max_abs_delta), float(max_abs_delta)
            )

        delta = np.zeros(response.matrix.shape[1], dtype=float)
        delta[selected_indices] = selected_delta
        predicted = measured + response.matrix @ delta
        return OrbitCorrectionResult(
            plane=response.plane,
            bpm_names=response.bpm_names,
            corrector_names=response.corrector_names,
            delta_correctors=delta,
            measured_orbit=measured,
            target_orbit=target,
            predicted_orbit=predicted,
            selected_correctors=tuple(
                response.corrector_names[index] for index in selected_indices
            ),
            singular_values=singular_values,
            rank=rank,
        )

    def propose_dispersion_correction(
        self,
        response: DispersionResponseResult,
        *,
        target_dispersion: Sequence[float],
        target_orbit: Sequence[float] | None = None,
        bpm_weights: Sequence[float] | None = None,
        dispersion_weight: float = 0.05,
        orbit_weight: float = 1.0,
        rcond: float = 1e-4,
        gain: float = 1.0,
        max_correctors: int | None = None,
        max_abs_delta: float | None = None,
    ) -> DispersionCorrectionResult:
        """Jointly fit dispersion and COD, matching SAD ``EtaCorrect``."""
        plane_index = self._plane_index(response.plane)
        orbit = self.find_closed_orbit(bpm_names=response.bpm_names)
        dispersion = self.measure_dispersion(
            relative_momentum_step=response.relative_momentum_step,
            bpm_names=response.bpm_names,
        )
        measured_orbit = orbit.bpm_positions[:, plane_index]
        measured_dispersion = dispersion.values[:, plane_index]
        target_dispersion_array = np.asarray(target_dispersion, dtype=float)
        target_orbit_array = (
            np.zeros_like(measured_orbit)
            if target_orbit is None
            else np.asarray(target_orbit, dtype=float)
        )
        if target_dispersion_array.shape != measured_dispersion.shape:
            raise ValueError("target_dispersion length must match response BPMs")
        if target_orbit_array.shape != measured_orbit.shape:
            raise ValueError("target_orbit length must match response BPMs")
        weights = (
            np.ones_like(measured_orbit)
            if bpm_weights is None
            else np.asarray(bpm_weights, dtype=float)
        )
        if weights.shape != measured_orbit.shape or np.any(weights < 0):
            raise ValueError("bpm_weights must be non-negative and match response BPMs")
        if dispersion_weight < 0 or orbit_weight < 0:
            raise ValueError("dispersion_weight and orbit_weight must be non-negative")
        if dispersion_weight == 0 and orbit_weight == 0:
            raise ValueError("At least one correction weight must be positive")
        if rcond < 0:
            raise ValueError("rcond must be non-negative")

        matrix = np.vstack(
            (
                response.dispersion_matrix
                * (dispersion_weight * weights[:, None]),
                response.orbit_matrix * (orbit_weight * weights[:, None]),
            )
        )
        target = np.concatenate(
            (
                (target_dispersion_array - measured_dispersion)
                * (dispersion_weight * weights),
                (target_orbit_array - measured_orbit)
                * (orbit_weight * weights),
            )
        )
        if max_correctors is None:
            selected_indices = list(range(matrix.shape[1]))
        else:
            if max_correctors <= 0:
                raise ValueError("max_correctors must be positive")
            selected_indices = self._greedy_columns(
                matrix,
                target,
                min(max_correctors, matrix.shape[1]),
                rcond,
            )
        if not selected_indices:
            raise ClosedOrbitError("No useful corrector columns were selected")

        selected_delta, singular_values, rank = self._svd_solve(
            matrix[:, selected_indices], target, rcond
        )
        selected_delta *= float(gain)
        if max_abs_delta is not None:
            if max_abs_delta <= 0:
                raise ValueError("max_abs_delta must be positive")
            selected_delta = np.clip(
                selected_delta, -float(max_abs_delta), float(max_abs_delta)
            )

        delta = np.zeros(response.dispersion_matrix.shape[1], dtype=float)
        delta[selected_indices] = selected_delta
        predicted_dispersion = (
            measured_dispersion + response.dispersion_matrix @ delta
        )
        predicted_orbit = measured_orbit + response.orbit_matrix @ delta
        return DispersionCorrectionResult(
            plane=response.plane,
            bpm_names=response.bpm_names,
            corrector_names=response.corrector_names,
            delta_correctors=delta,
            measured_dispersion=measured_dispersion,
            target_dispersion=target_dispersion_array,
            predicted_dispersion=predicted_dispersion,
            measured_orbit=measured_orbit,
            target_orbit=target_orbit_array,
            predicted_orbit=predicted_orbit,
            selected_correctors=tuple(
                response.corrector_names[index] for index in selected_indices
            ),
            singular_values=singular_values,
            rank=rank,
            relative_momentum_step=response.relative_momentum_step,
        )

    def propose_coupling_correction(
        self,
        response: CouplingResponseResult,
        *,
        target_signal: Sequence[float] | None = None,
        measurement_weights: Sequence[float] | None = None,
        rcond: float = 0.4,
        gain: float = 1.0,
        max_skew_correctors: int | None = None,
        max_abs_delta: float | None = None,
    ) -> CouplingCorrectionResult:
        """Solve a skew correction from horizontal-to-vertical COD response.

        The default SVD cutoff (0.4) follows ``skewcor.n``.  SAD's absolute
        SD/SF current limits depend on machine calibration records, which are
        not present in this model, so an explicit ``max_abs_delta`` is used
        when such a limit is needed.
        """
        _, measured_2d = self._measure_vertical_response_to_horizontal_probes(
            response.probe_corrector_names,
            bpm_names=response.bpm_names,
            probe_perturbation=response.probe_perturbation,
        )
        measured = measured_2d.reshape(-1)
        target = (
            np.zeros_like(measured)
            if target_signal is None
            else np.asarray(target_signal, dtype=float)
        )
        if target.shape != measured.shape:
            raise ValueError("target_signal length must match coupling response")
        weights = (
            np.ones_like(measured)
            if measurement_weights is None
            else np.asarray(measurement_weights, dtype=float)
        )
        if weights.shape != measured.shape or np.any(weights < 0):
            raise ValueError(
                "measurement_weights must be non-negative and match coupling response"
            )
        if rcond < 0:
            raise ValueError("rcond must be non-negative")

        weighted_matrix = response.matrix * weights[:, None]
        weighted_target = (target - measured) * weights
        if max_skew_correctors is None:
            selected_indices = list(range(response.matrix.shape[1]))
        else:
            if max_skew_correctors <= 0:
                raise ValueError("max_skew_correctors must be positive")
            selected_indices = self._greedy_columns(
                weighted_matrix,
                weighted_target,
                min(max_skew_correctors, response.matrix.shape[1]),
                rcond,
            )
        if not selected_indices:
            raise ClosedOrbitError("No useful skew-corrector columns were selected")

        selected_delta, singular_values, rank = self._svd_solve(
            weighted_matrix[:, selected_indices], weighted_target, rcond
        )
        selected_delta *= float(gain)
        if max_abs_delta is not None:
            if max_abs_delta <= 0:
                raise ValueError("max_abs_delta must be positive")
            selected_delta = np.clip(
                selected_delta, -float(max_abs_delta), float(max_abs_delta)
            )

        delta = np.zeros(response.matrix.shape[1], dtype=float)
        delta[selected_indices] = selected_delta
        predicted = measured + response.matrix @ delta
        return CouplingCorrectionResult(
            bpm_names=response.bpm_names,
            probe_corrector_names=response.probe_corrector_names,
            skew_corrector_names=response.skew_corrector_names,
            probe_perturbation=response.probe_perturbation,
            delta_skew_correctors=delta,
            measured_signal=measured,
            target_signal=target,
            predicted_signal=predicted,
            selected_skew_correctors=tuple(
                response.skew_corrector_names[index] for index in selected_indices
            ),
            singular_values=singular_values,
            rank=rank,
        )

    def apply_orbit_correction(
        self,
        correction: OrbitCorrectionResult,
    ) -> ClosedOrbitResult:
        """Apply a proposed correction and return the resulting closed orbit."""
        plane_index = self._plane_index(correction.plane)
        for name, delta in zip(
            correction.corrector_names, correction.delta_correctors
        ):
            if delta == 0.0:
                continue
            element = self._single_element(name)
            strength = np.asarray(element.get_strength(), dtype=float)
            strength[plane_index] += self.actuator_scale * float(delta)
            element.set_strength(*strength)
        return self.find_closed_orbit(bpm_names=correction.bpm_names)

    def apply_coupling_correction(
        self,
        correction: CouplingCorrectionResult,
    ) -> np.ndarray:
        """Apply a skew correction and return the remeasured coupling signal."""
        for name, delta in zip(
            correction.skew_corrector_names, correction.delta_skew_correctors
        ):
            if delta == 0.0:
                continue
            self._set_skew_strength(
                name,
                self._get_skew_strength(name)
                + self.skew_actuator_scale * float(delta),
            )
        _, signal = self._measure_vertical_response_to_horizontal_probes(
            correction.probe_corrector_names,
            bpm_names=correction.bpm_names,
            probe_perturbation=correction.probe_perturbation,
        )
        return signal.reshape(-1)

    def apply_dispersion_correction(
        self,
        correction: DispersionCorrectionResult,
    ) -> tuple[ClosedOrbitResult, DispersionResult]:
        """Apply a joint dispersion/COD suggestion and remeasure both."""
        plane_index = self._plane_index(correction.plane)
        for name, delta in zip(
            correction.corrector_names, correction.delta_correctors
        ):
            if delta == 0.0:
                continue
            element = self._single_element(name)
            strength = np.asarray(element.get_strength(), dtype=float)
            strength[plane_index] += self.actuator_scale * float(delta)
            element.set_strength(*strength)
        orbit = self.find_closed_orbit(bpm_names=correction.bpm_names)
        dispersion = self.measure_dispersion(
            relative_momentum_step=correction.relative_momentum_step,
            bpm_names=correction.bpm_names,
        )
        return orbit, dispersion


__all__ = [
    "ATFDRRingCorrection",
    "ClosedOrbitError",
    "ClosedOrbitResult",
    "CouplingCorrectionResult",
    "CouplingResponseResult",
    "DispersionCorrectionResult",
    "DispersionResult",
    "DispersionResponseResult",
    "OrbitCorrectionResult",
    "OrbitResponseResult",
    "ParticleLostError",
]
