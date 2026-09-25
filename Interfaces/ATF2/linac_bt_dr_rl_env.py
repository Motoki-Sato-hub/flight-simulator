"""Dependency-free RL environment for offline Linac--BT--DR orbit studies.

The environment deliberately uses the real SAD-derived RF-Track pipeline at
each ``step``.  It is therefore slow enough that endpoint orbit matching is
the training reward and multi-turn survival is a sparse validation metric.
No control-system interface is imported.  The API mirrors the useful subset of
Gymnasium: ``reset() -> (observation, info)`` and
``step(action) -> (observation, reward, terminated, truncated, info)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from Interfaces.ATF2.ATF2_LinacBTDR_RFTrack import ATF2LinacBTDRRFTrack


MODEL_MATCHED_CAVITY_VOLTAGE_MV = 75.85997836492652
DEFAULT_ACTUATORS = ("ZH1L", "ZH2L", "ZV1L", "ZV2L")
DEFAULT_SYNTHETIC_ERROR = {"ZH5L": 3.0e-4, "ZV5L": -2.0e-4}


@dataclass(frozen=True)
class EpisodeInfo:
    step: int
    orbit_error_mm_mrad: np.ndarray
    correctors_rftrack_strength: dict[str, float]


class LinacBTDRInjectionEnv:
    """Offline continuous-action environment for DR injection-orbit correction.

    Observation is the four-vector ``[x_mm, xp_mrad, y_mm, yp_mrad]`` relative
    to the DR model closed orbit.  Actions are normalised increments in
    ``[-1, 1]`` for each actuator, multiplied by ``action_scale`` in native
    RF-Track corrector-strength units and bounded cumulatively by
    ``max_abs_delta_strength`` around the reset state.  This bound is a
    numerical safety envelope, not a physical supply calibration.  The reward
    is negative normalised orbit norm, with a small action penalty.  A learned
    policy must be compared with the response-matrix solver under identical
    errors, steps, actuator limits, and final multi-turn validation.
    """

    def __init__(
        self,
        *,
        actuators: Sequence[str] = DEFAULT_ACTUATORS,
        synthetic_error: Mapping[str, float] = DEFAULT_SYNTHETIC_ERROR,
        action_scale: float = 1.0e-4,
        orbit_scale_mm_mrad: Sequence[float] = (0.1, 0.1, 0.1, 0.1),
        success_tolerance: float = 1.0e-4,
        max_steps: int = 8,
        max_abs_delta_strength: float = 1.0e-3,
    ):
        self.actuators = tuple(actuators)
        self.synthetic_error = dict(synthetic_error)
        self.action_scale = float(action_scale)
        self.orbit_scale = np.asarray(orbit_scale_mm_mrad, dtype=float)
        self.success_tolerance = float(success_tolerance)
        self.max_steps = int(max_steps)
        self.max_abs_delta_strength = float(max_abs_delta_strength)
        if not self.actuators or not np.isfinite(self.action_scale) or self.action_scale <= 0.0:
            raise ValueError("at least one actuator and a positive finite action_scale are required")
        if self.orbit_scale.shape != (4,) or np.any(~np.isfinite(self.orbit_scale)) or np.any(self.orbit_scale <= 0.0):
            raise ValueError("orbit_scale_mm_mrad must contain four positive finite values")
        if self.max_steps < 1 or self.success_tolerance <= 0.0:
            raise ValueError("max_steps and success_tolerance must be positive")
        if not np.isfinite(self.max_abs_delta_strength) or self.max_abs_delta_strength <= 0.0:
            raise ValueError("max_abs_delta_strength must be positive and finite")

        self.machine = ATF2LinacBTDRRFTrack(
            cavity_voltage_mv=MODEL_MATCHED_CAVITY_VOLTAGE_MV,
            handoff_mode="sad_optics_matched",
            dr_rf_mode="equilibrium",
        )
        all_names = tuple(dict.fromkeys((*self.actuators, *self.synthetic_error)))
        self._initial_strengths = self.machine.get_linac_bt_correctors(all_names)
        self._episode_actuator_centres = {
            name: self._initial_strengths[name] for name in self.actuators
        }
        self._step = 0

    def _orbit_error(self) -> np.ndarray:
        result = self.machine.track(self.machine.make_reference_bunch(), dr_turns=0)
        observed = np.array((
            result.dr_injection.mean_x_mm,
            result.dr_injection.mean_xp_mrad,
            result.dr_injection.mean_y_mm,
            result.dr_injection.mean_yp_mrad,
        ))
        return observed - self.machine.dr_closed_orbit

    def _info(self, orbit_error: np.ndarray) -> dict[str, object]:
        values = self.machine.get_linac_bt_correctors(self.actuators)
        return {
            "step": self._step,
            "orbit_error_mm_mrad": orbit_error.copy(),
            "correctors_rftrack_strength": values,
            "corrector_delta_from_episode_start_rftrack_strength": {
                name: values[name] - self._episode_actuator_centres[name]
                for name in self.actuators
            },
            "max_abs_delta_strength": self.max_abs_delta_strength,
        }

    def reset(self) -> tuple[np.ndarray, dict[str, object]]:
        """Restore synthetic baseline/error settings and return first observation."""
        self.machine.set_linac_bt_correctors(self._initial_strengths)
        self.machine.set_linac_bt_correctors(self.synthetic_error)
        self._episode_actuator_centres = self.machine.get_linac_bt_correctors(self.actuators)
        self._step = 0
        orbit_error = self._orbit_error()
        return orbit_error.copy(), self._info(orbit_error)

    def step(self, action: Sequence[float]) -> tuple[np.ndarray, float, bool, bool, dict[str, object]]:
        """Apply one bounded corrector increment and calculate endpoint reward."""
        action = np.asarray(action, dtype=float)
        if action.shape != (len(self.actuators),):
            raise ValueError(f"action must have shape ({len(self.actuators)},)")
        if np.any(~np.isfinite(action)):
            raise ValueError("action must be finite")
        clipped = np.clip(action, -1.0, 1.0)
        current = self.machine.get_linac_bt_correctors(self.actuators)
        proposed = {
            name: current[name] + self.action_scale * value
            for name, value in zip(self.actuators, clipped)
        }
        # This is a numerical safety envelope for offline RF-Track scans, not
        # a physical ATF supply/current constraint.  In particular, very large
        # native strengths can abort RF-Track's C++ process instead of yielding
        # a recoverable lost-particle result.
        bounded = {
            name: float(np.clip(
                value,
                self._episode_actuator_centres[name] - self.max_abs_delta_strength,
                self._episode_actuator_centres[name] + self.max_abs_delta_strength,
            ))
            for name, value in proposed.items()
        }
        self.machine.set_linac_bt_correctors(bounded)
        self._step += 1
        orbit_error = self._orbit_error()
        norm = float(np.linalg.norm(orbit_error / self.orbit_scale))
        reward = -norm - 1.0e-3 * float(np.dot(clipped, clipped))
        terminated = bool(np.linalg.norm(orbit_error) <= self.success_tolerance)
        truncated = bool(self._step >= self.max_steps and not terminated)
        return orbit_error.copy(), reward, terminated, truncated, self._info(orbit_error)

    def validate_multi_turn(
        self, *, turns: int = 100, turn_history_sample_every: int = 1
    ) -> dict[str, object]:
        """Run sparse long-horizon validation for the current policy result."""
        if turns < 1:
            raise ValueError("turns must be positive")
        if turn_history_sample_every < 1:
            raise ValueError("turn_history_sample_every must be positive")
        result = self.machine.track(
            self.machine.make_reference_bunch(),
            dr_turns=turns,
            record_turn_history=True,
            turn_history_sample_every=turn_history_sample_every,
        )
        history = list(zip(
            result.dr_turn_history_turns,
            (item.survival_fraction_from_input for item in result.dr_turn_history),
        ))
        return {
            "turns_requested": turns,
            "turns_completed": result.completed_dr_turns,
            "ring_survival_after_turns": result.dr_after_turns.survival_fraction_from_input,
            "turn_history_sample_every": result.turn_history_sample_every,
            "first_loss_turn_is_exact": turn_history_sample_every == 1,
            "first_loss_turn": next(
                (turn for turn, value in history if value < 1.0),
                None,
            ),
        }


__all__ = [
    "DEFAULT_ACTUATORS",
    "DEFAULT_SYNTHETIC_ERROR",
    "LinacBTDRInjectionEnv",
]
