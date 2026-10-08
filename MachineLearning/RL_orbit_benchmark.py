"""Small, reproducible corrector-only RL benchmark for the ATF Linac RFTrack model.

This intentionally starts with a one-step orbit-recovery task.  It permits a
fair comparison between a response-matrix correction and a genuinely
model-free, continuous-action policy under exactly the same initial-orbit
distribution.  The policy is a Gaussian linear REINFORCE agent; it is a
smoke-test of the environment and safety conventions, not a claim that this
algorithm is the final controller.  TD3/NAF/SAC can use the same environment
once an RL framework is selected.

Example
-------
MPLCONFIGDIR=/tmp/mpl-flight-simulator PYTHONPATH=. \\
  /home/motokisato/rftrack-env/bin/python MachineLearning/RL_orbit_benchmark.py
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass

import numpy as np

from Interfaces.ATF2.InterfaceATF2_Linac_RFTrack import InterfaceATF2_Linac_RFTrack


@dataclass
class BenchmarkConfig:
    finite_difference_step: float = 1e-4
    initial_current_span: float = 5e-3
    max_current_step: float = 5e-3
    exploration_std: float = 1e-3
    policy_learning_rate: float = 5e-6
    policy_batch_size: int = 32
    action_penalty: float = 1e-2


class LinacOrbitEnvironment:
    """Corrector-only one-step environment backed directly by RFTrack."""

    def __init__(self, machine: InterfaceATF2_Linac_RFTrack, config: BenchmarkConfig):
        self.machine = machine
        self.config = config
        self.correctors = list(machine.corrs)
        self.origin = machine.get_correctors(self.correctors)["bdes"].astype(float)
        self.reference_orbit = self.observe()

    def observe(self) -> np.ndarray:
        bpms = self.machine.get_bpms()
        return np.concatenate((bpms["x"].mean(axis=0), bpms["y"].mean(axis=0)))

    def reset(self, initial_delta: np.ndarray) -> np.ndarray:
        self.machine.set_correctors(self.correctors, self.origin + initial_delta)
        return self.observe() - self.reference_orbit

    def step(self, delta: np.ndarray) -> np.ndarray:
        delta = np.clip(delta, -self.config.max_current_step, self.config.max_current_step)
        current = self.machine.get_correctors(self.correctors)["bdes"].astype(float)
        self.machine.set_correctors(self.correctors, current + delta)
        return self.observe() - self.reference_orbit

    def restore(self) -> None:
        self.machine.set_correctors(self.correctors, self.origin)


def measure_response_matrix(env: LinacOrbitEnvironment) -> np.ndarray:
    """Central-difference response matrix in mm per interface-current unit."""
    h = env.config.finite_difference_step
    response = np.empty((env.reference_orbit.size, len(env.correctors)))
    for column in range(len(env.correctors)):
        plus = env.origin.copy()
        plus[column] += h
        env.machine.set_correctors(env.correctors, plus)
        y_plus = env.observe()

        minus = env.origin.copy()
        minus[column] -= h
        env.machine.set_correctors(env.correctors, minus)
        y_minus = env.observe()
        response[:, column] = (y_plus - y_minus) / (2.0 * h)
    env.restore()
    return response


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def train_reinforce(env: LinacOrbitEnvironment, episodes: int, seed: int) -> np.ndarray:
    """Train a Gaussian linear policy directly from RFTrack rewards.

    One reset/action/observation is one episode.  REINFORCE is used here
    because it requires only NumPy and makes the first benchmark runnable in
    the existing RFTrack environment.
    """
    rng = np.random.default_rng(seed)
    n_action, n_state = len(env.correctors), env.reference_orbit.size
    theta = np.zeros((n_action, n_state), dtype=float)
    gradients: list[np.ndarray] = []
    rewards: list[float] = []

    for episode in range(episodes):
        initial = rng.uniform(-env.config.initial_current_span, env.config.initial_current_span, n_action)
        state = env.reset(initial)
        noise = rng.normal(0.0, env.config.exploration_std, n_action)
        action = theta @ state + noise
        action = np.clip(action, -env.config.max_current_step, env.config.max_current_step)
        next_state = env.step(action)
        reward = -rms(next_state) ** 2 - env.config.action_penalty * float(np.mean(action**2))

        # Score-function gradient of log pi(a|s), ignoring rare clipped samples.
        gradients.append(np.outer(noise / env.config.exploration_std**2, state))
        rewards.append(reward)

        if len(rewards) == env.config.policy_batch_size or episode == episodes - 1:
            advantage = np.asarray(rewards, dtype=float)
            advantage = (advantage - advantage.mean()) / (advantage.std() + 1e-12)
            gradient = np.mean([a * g for a, g in zip(advantage, gradients)], axis=0)
            theta += env.config.policy_learning_rate * gradient
            gradients.clear()
            rewards.clear()

    env.restore()
    return theta


def evaluate(
    env: LinacOrbitEnvironment,
    policy: np.ndarray,
    response: np.ndarray,
    trials: int,
    seed: int,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    initial_rms, rm_rms, rl_rms = [], [], []
    for _ in range(trials):
        initial = rng.uniform(-env.config.initial_current_span, env.config.initial_current_span, len(env.correctors))
        state = env.reset(initial)
        initial_rms.append(rms(state))

        rm_action = -np.linalg.pinv(response, rcond=1e-8) @ state
        rm_rms.append(rms(env.step(rm_action)))

        state = env.reset(initial)
        rl_action = policy @ state
        rl_rms.append(rms(env.step(rl_action)))

    env.restore()
    return {
        "initial_rms_mm_mean": float(np.mean(initial_rms)),
        "response_matrix_rms_mm_mean": float(np.mean(rm_rms)),
        "reinforce_rms_mm_mean": float(np.mean(rl_rms)),
        "response_matrix_success_fraction_0p1mm": float(np.mean(np.asarray(rm_rms) < 0.1)),
        "reinforce_success_fraction_0p1mm": float(np.mean(np.asarray(rl_rms) < 0.1)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--episodes", type=int, default=256)
    parser.add_argument("--trials", type=int, default=64)
    parser.add_argument("--nparticles", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260901)
    args = parser.parse_args()

    config = BenchmarkConfig()
    machine = InterfaceATF2_Linac_RFTrack(
        nparticles=args.nparticles,
        bpm_resolution=0.0,
        jitter=0.0,
    )
    machine.log = lambda *unused_args, **unused_kwargs: None
    env = LinacOrbitEnvironment(machine, config)
    try:
        response = measure_response_matrix(env)
        policy = train_reinforce(env, args.episodes, args.seed)
        result = evaluate(env, policy, response, args.trials, args.seed + 1)
        result.update(
            response_shape=list(response.shape),
            response_rank=int(np.linalg.matrix_rank(response, tol=1e-10)),
            rl_algorithm="linear Gaussian REINFORCE (one-step, model-free)",
            rl_training_episodes=args.episodes,
            response_measurement_perturbations=2 * len(env.correctors),
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    finally:
        env.restore()


if __name__ == "__main__":
    main()
