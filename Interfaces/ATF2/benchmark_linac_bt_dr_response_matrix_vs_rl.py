"""Compare response-matrix correction with a minimal model-free RL baseline.

Both methods see the same fixed synthetic Linac+BT corrector error, four
actuators, RF-Track endpoint observation, 8-step action horizon, action
scale, and cumulative native-strength safety envelope.  The response-matrix
method gets a finite-difference model; the other method is deliberately a
dependency-free cross-entropy (CEM) policy search over a *constant* bounded
action.  The latter is a small, reproducible model-free RL-style baseline,
not a claim that CEM is a production RL controller or that it has learned a
policy transferable to ATF.

The report includes endpoint interaction counts.  That is essential: a
response matrix should be expected to win this local linear, noiseless model
on sample efficiency.  A future trained feedback policy should only be called
an improvement if it wins a pre-declared noisy/nonlinear robustness test under
the same constraints and independently validates multi-turn survival.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass

import numpy as np

from Interfaces.ATF2.linac_bt_dr_rl_env import LinacBTDRInjectionEnv


@dataclass(frozen=True)
class MethodResult:
    name: str
    final_orbit_error_mm_mrad: np.ndarray
    final_normalised_orbit_norm: float
    actions_applied: int
    endpoint_interactions: int
    final_action: np.ndarray
    corrector_delta_rftrack_strength: dict[str, float]
    multi_turn: dict[str, object]
    detail: dict[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "final_orbit_error_mm_mrad": self.final_orbit_error_mm_mrad.tolist(),
            "final_normalised_orbit_norm": self.final_normalised_orbit_norm,
            "actions_applied": self.actions_applied,
            "endpoint_interactions": self.endpoint_interactions,
            "final_action": self.final_action.tolist(),
            "corrector_delta_rftrack_strength": self.corrector_delta_rftrack_strength,
            "multi_turn": self.multi_turn,
            "detail": self.detail,
        }


def _normalised_norm(env: LinacBTDRInjectionEnv, orbit: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(orbit, dtype=float) / env.orbit_scale))


def _episode_with_constant_action(
    env: LinacBTDRInjectionEnv, action: np.ndarray
) -> tuple[np.ndarray, int]:
    """Reset and apply one stationary bounded policy for the full horizon."""
    orbit, _ = env.reset()
    interactions = 1
    for _ in range(env.max_steps):
        orbit, _, terminated, truncated, _ = env.step(action)
        interactions += 1
        if terminated or truncated:
            break
    return orbit, interactions


def run_response_matrix(
    env: LinacBTDRInjectionEnv, *, response_step: float, svd_rcond: float,
    validation_turns: int,
) -> MethodResult:
    """Use the finite-difference response, while obeying RL action limits."""
    before, reset_info = env.reset()
    response, names = env.machine.linac_bt_dr_orbit_response_matrix(
        env.actuators, step=response_step
    )
    delta, _, rank, singular_values = np.linalg.lstsq(response, -before, rcond=svd_rcond)
    # A correction may not exceed one normalised action per call.  Spread it
    # over the same finite horizon used by the model-free policy.
    required_steps = max(1, int(np.ceil(np.max(np.abs(delta)) / env.action_scale)))
    steps_to_apply = min(required_steps, env.max_steps)
    action = delta / (steps_to_apply * env.action_scale)
    action = np.clip(action, -1.0, 1.0)
    orbit = before
    action_steps = 0
    for _ in range(steps_to_apply):
        orbit, _, terminated, truncated, _ = env.step(action)
        action_steps += 1
        if terminated or truncated:
            break
    final_info = env._info(orbit)
    return MethodResult(
        name="finite_difference_response_matrix",
        final_orbit_error_mm_mrad=orbit,
        final_normalised_orbit_norm=_normalised_norm(env, orbit),
        actions_applied=action_steps,
        # reset endpoint + baseline and four shifted response points + steps
        endpoint_interactions=1 + 1 + len(names) + action_steps,
        final_action=action,
        corrector_delta_rftrack_strength=final_info[
            "corrector_delta_from_episode_start_rftrack_strength"
        ],
        multi_turn=env.validate_multi_turn(turns=validation_turns),
        detail={
            "response_step_rftrack_strength": response_step,
            "response_matrix_mm_mrad_per_strength": response.tolist(),
            "singular_values": singular_values.tolist(),
            "rank": int(rank),
            "unconstrained_total_delta_rftrack_strength": delta.tolist(),
            "steps_required_by_action_limit": required_steps,
            "horizon_limited": required_steps > env.max_steps,
            "reset_orbit_error_mm_mrad": before.tolist(),
            "reset_corrector_strengths": reset_info["correctors_rftrack_strength"],
        },
    )


def run_cem_policy_search(
    env: LinacBTDRInjectionEnv, *, population: int, generations: int,
    elite_fraction: float, seed: int, validation_turns: int,
) -> MethodResult:
    """Train a stationary bounded action policy with black-box CEM returns."""
    if population < 2 or generations < 1 or not 0.0 < elite_fraction < 1.0:
        raise ValueError("population>=2, generations>=1, and 0<elite_fraction<1 are required")
    rng = np.random.default_rng(seed)
    dimensions = len(env.actuators)
    elite_count = max(1, int(np.ceil(population * elite_fraction)))
    mean = np.zeros(dimensions, dtype=float)
    standard_deviation = np.full(dimensions, 0.7, dtype=float)
    best_action = mean.copy()
    best_norm = np.inf
    interactions = 0
    history: list[dict[str, float]] = []
    for generation in range(generations):
        candidates = np.clip(
            rng.normal(mean, standard_deviation, size=(population, dimensions)), -1.0, 1.0
        )
        norms = np.empty(population, dtype=float)
        for index, candidate in enumerate(candidates):
            orbit, used = _episode_with_constant_action(env, candidate)
            interactions += used
            norms[index] = _normalised_norm(env, orbit)
        best_index = int(np.argmin(norms))
        if norms[best_index] < best_norm:
            best_norm = float(norms[best_index])
            best_action = candidates[best_index].copy()
        elite = candidates[np.argsort(norms)[:elite_count]]
        # Smooth updates and retain a small exploratory variance, otherwise a
        # deterministic RF-Track reward can collapse CEM too early.
        mean = 0.3 * mean + 0.7 * np.mean(elite, axis=0)
        standard_deviation = np.maximum(
            0.05, 0.3 * standard_deviation + 0.7 * np.std(elite, axis=0)
        )
        history.append({
            "generation": generation + 1,
            "best_normalised_orbit_norm": float(np.min(norms)),
            "mean_normalised_orbit_norm": float(np.mean(norms)),
        })

    final_orbit, final_used = _episode_with_constant_action(env, best_action)
    interactions += final_used
    final_info = env._info(final_orbit)
    return MethodResult(
        name="model_free_cem_stationary_policy",
        final_orbit_error_mm_mrad=final_orbit,
        final_normalised_orbit_norm=_normalised_norm(env, final_orbit),
        actions_applied=final_used - 1,
        endpoint_interactions=interactions,
        final_action=best_action,
        corrector_delta_rftrack_strength=final_info[
            "corrector_delta_from_episode_start_rftrack_strength"
        ],
        multi_turn=env.validate_multi_turn(turns=validation_turns),
        detail={
            "policy": "constant bounded action repeated over the finite episode",
            "optimizer": "cross-entropy method (black-box episodic return)",
            "population": population,
            "generations": generations,
            "elite_fraction": elite_fraction,
            "random_seed": seed,
            "generation_history": history,
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-turns", type=int, default=10)
    parser.add_argument("--response-step", type=float, default=1.0e-4)
    parser.add_argument("--svd-rcond", type=float, default=1.0e-8)
    parser.add_argument("--population", type=int, default=8)
    parser.add_argument("--generations", type=int, default=4)
    parser.add_argument("--elite-fraction", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()
    if args.validation_turns < 1 or not np.isfinite(args.response_step) or args.response_step == 0.0:
        raise ValueError("validation-turns must be positive and response-step finite/non-zero")

    # Distinct environments ensure one method's final correctors cannot leak
    # into the other's reset, while configuration remains exactly the same.
    response_env = LinacBTDRInjectionEnv()
    cem_env = LinacBTDRInjectionEnv()
    response = run_response_matrix(
        response_env, response_step=args.response_step, svd_rcond=args.svd_rcond,
        validation_turns=args.validation_turns,
    )
    cem = run_cem_policy_search(
        cem_env, population=args.population, generations=args.generations,
        elite_fraction=args.elite_fraction, seed=args.seed,
        validation_turns=args.validation_turns,
    )
    print(json.dumps({
        "simulation_only": True,
        "comparison_contract": {
            "synthetic_error_rftrack_strength": response_env.synthetic_error,
            "actuators": list(response_env.actuators),
            "max_steps": response_env.max_steps,
            "action_scale_rftrack_strength": response_env.action_scale,
            "max_abs_delta_strength": response_env.max_abs_delta_strength,
            "observation": "DR injection [x_mm, xp_mrad, y_mm, yp_mrad] relative to closed orbit",
            "objective": "minimise final normalised endpoint orbit norm; then validate the selected result over DR turns",
        },
        "response_matrix": response.as_dict(),
        "model_free_rl_baseline": cem.as_dict(),
        "interpretation": {
            "expected_local_result": "The response matrix has explicit local derivatives, so it should be more sample-efficient in this deterministic near-linear synthetic task.",
            "not_yet_tested": "Neither method is validated against BPM noise, magnet hysteresis, changing optics, real aperture/loss data, or a surveyed IPZT-to-RING0 map.",
            "decision_rule": "Do not claim RL superiority from this fixture.  Pre-register noisy/nonlinear held-out errors, equal actuator constraints and interaction budget, then require multi-turn and physical-transmission validation.",
        },
    }, indent=2))


if __name__ == "__main__":
    main()
