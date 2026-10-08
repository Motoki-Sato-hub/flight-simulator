"""Fair ORM-versus-RL orbit benchmark on the synthetic ATF BT digital twin.

All controllers see calibrated BPM coordinates and have the same one-step
corrector bound.  They differ only in information and interaction cost:

* nominal ORM: calculated from the uncalibrated model (zero beam probes);
* measured ORM: 44 central-difference pseudo-machine probes (22 correctors);
* REINFORCE: model-free samples from the pseudo machine during training.

The RL policy is deliberately a NumPy-only linear Gaussian REINFORCE agent.
It validates the problem definition and sample accounting, rather than being a
claim that it is the eventual ATF controller (TD3/SAC are natural next agents).
No real-machine connection is used.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass

import numpy as np

from Interfaces.ATF2.InterfaceATF2_BT_RFTrack import InterfaceATF2_BT_RFTrack
from Interfaces.ATF2.synthetic_bt_calibration_audit import (
    BTPseudoMachine,
    MAX_CORRECTION_T_MM,
    SVD_RCOND,
    _measure_response,
    _physical_orbit,
    _rms,
)


@dataclass(frozen=True)
class BenchmarkConfig:
    initial_command_span_t_mm: float = 0.010
    max_action_t_mm: float = MAX_CORRECTION_T_MM
    exploration_std_t_mm: float = 0.004
    state_scale_mm: float = 0.010
    policy_learning_rate: float = 1.0e-5
    batch_size: int = 64
    action_penalty: float = 0.002


class CalibratedPseudoMachineEnvironment:
    """One-step, safe-action environment with hidden BT pseudo-machine errors."""

    def __init__(self, pseudo_machine: BTPseudoMachine, config: BenchmarkConfig):
        self.pseudo_machine = pseudo_machine
        self.config = config
        self.n_actions = len(pseudo_machine.corrector_names)

    def reset(self, commands: np.ndarray) -> np.ndarray:
        commands = np.clip(
            np.asarray(commands, dtype=float),
            -self.config.initial_command_span_t_mm,
            self.config.initial_command_span_t_mm,
        )
        self.pseudo_machine.set_commands(commands)
        return self.pseudo_machine.observe(calibrated=True)

    def step(self, delta: np.ndarray) -> np.ndarray:
        delta = np.clip(
            np.asarray(delta, dtype=float),
            -self.config.max_action_t_mm,
            self.config.max_action_t_mm,
        )
        commands = np.clip(
            self.pseudo_machine.commands + delta,
            -self.config.max_action_t_mm,
            self.config.max_action_t_mm,
        )
        self.pseudo_machine.set_commands(commands)
        return self.pseudo_machine.observe(calibrated=True)

    def restore(self) -> None:
        self.pseudo_machine.set_commands(np.zeros(self.n_actions))


def _nominal_response() -> np.ndarray:
    model = InterfaceATF2_BT_RFTrack(reference_particle=True)
    return _measure_response(
        lambda command: model.set_correctors(model.corrs, command),
        lambda: _physical_orbit(model),
        len(model.corrs),
    )


def _measured_response(env: CalibratedPseudoMachineEnvironment) -> np.ndarray:
    return _measure_response(
        env.pseudo_machine.set_commands,
        lambda: env.pseudo_machine.observe(calibrated=True),
        env.n_actions,
    )


def train_reinforce(
    env: CalibratedPseudoMachineEnvironment,
    episodes: int,
    seed: int,
) -> np.ndarray:
    """Train a continuous-action model-free policy using one plant action/episode."""
    rng = np.random.default_rng(seed)
    n_state, n_action = len(env.reset(np.zeros(env.n_actions))), env.n_actions
    policy = np.zeros((n_action, n_state), dtype=float)
    gradients: list[np.ndarray] = []
    rewards: list[float] = []
    try:
        for episode in range(episodes):
            initial = rng.uniform(
                -env.config.initial_command_span_t_mm,
                env.config.initial_command_span_t_mm,
                n_action,
            )
            state = env.reset(initial)
            scaled_state = state / env.config.state_scale_mm
            noise = rng.normal(0.0, env.config.exploration_std_t_mm, n_action)
            action = policy @ scaled_state + noise
            action = np.clip(action, -env.config.max_action_t_mm, env.config.max_action_t_mm)
            next_state = env.step(action)
            reward = -(
                _rms(next_state) / env.config.state_scale_mm
            ) ** 2 - env.config.action_penalty * float(
                np.mean((action / env.config.max_action_t_mm) ** 2)
            )
            # Score-function gradient of a Gaussian policy.  Very rare clipped
            # actions use the unclipped noise approximation, conservatively.
            gradients.append(
                np.outer(noise / env.config.exploration_std_t_mm**2, scaled_state)
            )
            rewards.append(reward)
            if len(rewards) == env.config.batch_size or episode == episodes - 1:
                advantage = np.asarray(rewards)
                advantage = (advantage - advantage.mean()) / (advantage.std() + 1e-12)
                gradient = np.mean(
                    [a * g for a, g in zip(advantage, gradients)], axis=0
                )
                policy += env.config.policy_learning_rate * gradient
                gradients.clear()
                rewards.clear()
    finally:
        env.restore()
    return policy


def _evaluate_controller(
    env: CalibratedPseudoMachineEnvironment,
    controller,
    trials: int,
    seed: int,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    before, after, action_max = [], [], []
    for _ in range(trials):
        initial = rng.uniform(
            -env.config.initial_command_span_t_mm,
            env.config.initial_command_span_t_mm,
            env.n_actions,
        )
        state = env.reset(initial)
        action = np.clip(
            controller(state),
            -env.config.max_action_t_mm,
            env.config.max_action_t_mm,
        )
        next_state = env.step(action)
        before.append(_rms(state))
        after.append(_rms(next_state))
        action_max.append(float(np.max(np.abs(action))))
    env.restore()
    after_array = np.asarray(after)
    return {
        "before_rms_mm_mean": float(np.mean(before)),
        "after_rms_mm_mean": float(np.mean(after)),
        "success_fraction_0p1mm": float(np.mean(after_array < 0.1)),
        "max_abs_action_t_mm": float(np.max(action_max)),
    }


def run_benchmark(*, episodes: int = 2048, trials: int = 128, seed: int = 20260902) -> dict[str, object]:
    config = BenchmarkConfig()
    pseudo = BTPseudoMachine(seed)
    env = CalibratedPseudoMachineEnvironment(pseudo, config)
    nominal_response = _nominal_response()
    measured_response = _measured_response(env)
    nominal_policy = -np.linalg.pinv(nominal_response, rcond=SVD_RCOND)
    measured_policy = -np.linalg.pinv(measured_response, rcond=SVD_RCOND)
    rl_policy = train_reinforce(env, episodes, seed + 1)
    try:
        results = {
            "nominal_orm": _evaluate_controller(
                env, lambda state: nominal_policy @ state, trials, seed + 2
            ),
            "measured_orm": _evaluate_controller(
                env, lambda state: measured_policy @ state, trials, seed + 2
            ),
            "reinforce": _evaluate_controller(
                env,
                lambda state: rl_policy @ (state / config.state_scale_mm),
                trials,
                seed + 2,
            ),
        }
    finally:
        env.restore()
    return {
        "simulation_only": True,
        "task": "one-step calibrated-BPM corrector-only BT orbit recovery",
        "same_initial_distribution": True,
        "same_action_limit_t_mm": config.max_action_t_mm,
        "results": results,
        "response_matrix": {
            "shape": list(nominal_response.shape),
            "nominal_rank": int(np.linalg.matrix_rank(nominal_response)),
            "measured_rank": int(np.linalg.matrix_rank(measured_response)),
            "measured_orm_plant_perturbations": 2 * env.n_actions,
        },
        "reinforce": {
            "training_plant_actions": episodes,
            "algorithm": "linear Gaussian REINFORCE",
            "not_a_final_controller": True,
        },
        "interpretation": (
            "Compare both residual orbit and interaction count.  A low-dimensional linear "
            "orbit task is expected to favour a measured ORM; RL is being tested for its "
            "future value under nonlinearity, loss, constraints, and time variation."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--episodes", type=int, default=2048)
    parser.add_argument("--trials", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20260902)
    args = parser.parse_args()
    print(json.dumps(run_benchmark(episodes=args.episodes, trials=args.trials, seed=args.seed), indent=2))


if __name__ == "__main__":
    main()
