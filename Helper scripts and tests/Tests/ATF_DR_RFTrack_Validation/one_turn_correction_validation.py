"""Deterministic, step-by-step validation of the Kubo ATF DR correction flow.

This is the primary reproducible demonstration of the RF-Track port.  Each
case has one model-unregistered, deterministic injected fault in a virtual
machine; the nominal RF-Track lattice supplies the response matrix and never
receives that fault.
It intentionally excludes random-error robustness studies.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import ATFDRRingCorrection
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import build_atf_dr_lattice


DELTA = 3.5e-3
# Keep the deterministic demonstration inside the closed-orbit solver's
# locally linear range; production scans use the Kubo-style procedure below.
COD_KICK = 1.0e-5
DISPERSION_KICK = 1.0e-5
SKEW_ERROR_K1L = 1.0e-3
PROBES = ("ZH1R", "ZH2R")
COUPLING_PROBE_KICK = 1.0e-5


def _rms(values):
    return float(np.sqrt(np.mean(np.asarray(values, dtype=float) ** 2)))


def _one(correction, name):
    return correction._single_element(name)


def _add_corrector(correction, name, plane, value):
    element = _one(correction, name)
    kick = correction._get_corrector_kick(element)
    kick[plane] += value
    correction._set_corrector_kick(element, kick)


def _apply_correctors(correction, names, plane, values):
    for name, value in zip(names, values):
        if value:
            _add_corrector(correction, name, plane, float(value))


def _solve(matrix, residual, maximum, rcond=1e-4, count=6):
    selected = ATFDRRingCorrection._greedy_columns(
        matrix, -residual, count, rcond
    )
    values, _, rank = ATFDRRingCorrection._svd_solve(
        matrix[:, selected], -residual, rcond
    )
    commands = np.zeros(matrix.shape[1])
    commands[selected] = np.clip(values, -maximum, maximum)
    return commands, selected, rank


def nominal_case():
    correction = ATFDRRingCorrection(build_atf_dr_lattice())
    orbit = correction.find_closed_orbit()
    dispersion = correction.measure_dispersion(relative_momentum_step=DELTA)
    return {
        "bpm_count": len(orbit.bpm_names),
        "cod_x_rms_mm": _rms(orbit.x),
        "cod_y_rms_mm": _rms(orbit.y),
        "dispersion_x_rms_mm": _rms(dispersion.x),
        "dispersion_y_rms_mm": _rms(dispersion.y),
    }


def cod_case(plane, fault_name):
    index = 0 if plane == "x" else 1
    model = ATFDRRingCorrection(build_atf_dr_lattice())
    machine = ATFDRRingCorrection(build_atf_dr_lattice())
    _add_corrector(machine, fault_name, index, COD_KICK)
    before = machine.find_closed_orbit().bpm_positions[:, index]

    candidates = model.get_corrector_names(plane)[1:13]
    response = model.compute_orbit_response(plane, corrector_names=candidates)
    commands, selected, rank = _solve(response.matrix, before, 5e-5)
    _apply_correctors(machine, candidates, index, commands)
    after = machine.find_closed_orbit().bpm_positions[:, index]
    return {
        "fault": f"{fault_name}: +{COD_KICK:g} rad physical kick",
        "selected": [candidates[item] for item in selected],
        "rank": rank,
        "before": before,
        "predicted": before + response.matrix @ commands,
        "after": after,
    }


def vertical_dispersion_case():
    """A hidden vertical steerer setting produces COD and periodic Dy."""
    model = ATFDRRingCorrection(build_atf_dr_lattice())
    machine = ATFDRRingCorrection(build_atf_dr_lattice())
    target = model.measure_dispersion(relative_momentum_step=DELTA)
    _add_corrector(machine, "ZV1R", 1, DISPERSION_KICK)
    orbit_before = machine.find_closed_orbit()
    dispersion_before = machine.measure_dispersion(relative_momentum_step=DELTA)

    candidates = model.get_corrector_names("y")[1:9]
    response = model.compute_dispersion_response(
        "y", corrector_names=candidates, relative_momentum_step=DELTA
    )
    residual = np.concatenate((
        0.05 * (dispersion_before.y - target.y),
        orbit_before.y,
    ))
    matrix = np.vstack((0.05 * response.dispersion_matrix, response.orbit_matrix))
    commands, selected, rank = _solve(matrix, residual, 5e-5)
    _apply_correctors(machine, candidates, 1, commands)
    orbit_after = machine.find_closed_orbit()
    dispersion_after = machine.measure_dispersion(relative_momentum_step=DELTA)
    return {
        "fault": f"ZV1R: +{DISPERSION_KICK:g} rad physical kick",
        "selected": [candidates[item] for item in selected],
        "rank": rank,
        "target": target.y,
        "dispersion_before": dispersion_before.y,
        "dispersion_after": dispersion_after.y,
        "orbit_before": orbit_before.y,
        "orbit_after": orbit_after.y,
    }


def coupling_case():
    model = ATFDRRingCorrection(build_atf_dr_lattice())
    machine = ATFDRRingCorrection(build_atf_dr_lattice())
    available_skews = model.get_skew_corrector_names()
    # Eight phase-distributed candidates retain a fast validation while
    # representing the ring-wide SD1R/SF1R correction family.
    skews = tuple(
        available_skews[index]
        for index in np.linspace(0, len(available_skews) - 2, 8, dtype=int)
    )
    # Deliberately keep the injected skew outside the correction candidate
    # set: this validates compensation by other skew actuators, rather than
    # the trivial reset of the injected actuator.
    injected_skew = available_skews[-1]
    bpms = model.bpm_names
    response = model.compute_coupling_response(
        probe_corrector_names=PROBES,
        skew_corrector_names=skews,
        bpm_names=bpms,
        probe_perturbation=COUPLING_PROBE_KICK,
        skew_perturbation=1e-5,
    )
    machine.set_skew_strength(injected_skew, SKEW_ERROR_K1L)
    _, before_2d = machine._measure_vertical_response_to_horizontal_probes(
        PROBES, bpm_names=bpms, probe_perturbation=COUPLING_PROBE_KICK
    )
    before = before_2d.reshape(-1)
    commands, selected, rank = _solve(response.matrix, before, 2e-3, rcond=0.4)
    for name, command in zip(skews, commands):
        if command:
            machine.set_skew_strength(name, machine.get_skew_strength(name) + command)
    _, after_2d = machine._measure_vertical_response_to_horizontal_probes(
        PROBES, bpm_names=bpms, probe_perturbation=COUPLING_PROBE_KICK
    )
    return {
        "fault": f"{injected_skew}: +{SKEW_ERROR_K1L:g} K1L",
        "selected": [skews[item] for item in selected],
        "rank": rank,
        "before": before,
        "predicted": before + response.matrix @ commands,
        "after": after_2d.reshape(-1),
    }


def _summary(result):
    rows = []
    for plane in ("x", "y"):
        item = result["cod"][plane]
        rows.append((f"COD {plane}", _rms(item["before"]), _rms(item["after"]), "mm"))
    item = result["vertical_dispersion"]
    rows.append(("vertical dispersion", _rms(item["dispersion_before"] - item["target"]), _rms(item["dispersion_after"] - item["target"]), "mm/delta"))
    rows.append(("vertical COD with dispersion", _rms(item["orbit_before"]), _rms(item["orbit_after"]), "mm"))
    item = result["coupling"]
    rows.append(("coupling probe response", _rms(item["before"]), _rms(item["after"]), "mm/kick"))
    return rows


def _serializable(result):
    serial = {"nominal": result["nominal"], "cases": {}}
    for name in ("cod", "vertical_dispersion", "coupling"):
        value = result[name]
        if name == "cod":
            serial["cases"][name] = {
                plane: {key: (item.tolist() if isinstance(item, np.ndarray) else item)
                        for key, item in data.items()}
                for plane, data in value.items()
            }
        else:
            serial["cases"][name] = {
                key: (item.tolist() if isinstance(item, np.ndarray) else item)
                for key, item in value.items()
            }
    serial["summary"] = [list(row) for row in _summary(result)]
    return serial


def main():
    result = {
        "nominal": nominal_case(),
        "cod": {"x": cod_case("x", "ZH1R"), "y": cod_case("y", "ZV1R")},
        "vertical_dispersion": vertical_dispersion_case(),
        "coupling": coupling_case(),
    }
    rows = _summary(result)
    for name, before, after, unit in rows:
        print(f"{name:30s} {before:.4e} -> {after:.4e} {unit}")

    # Model and virtual machine are independent objects, but differ by one
    # deliberately injected fault that is not registered in the model.
    # COD uses a restricted set of twelve model steerers; the other two tests
    # are expected to converge more tightly in this deterministic setting.
    assert all(after < 0.20 * before for _, before, after, _ in rows[:2])
    assert rows[2][2] < 0.02 * rows[2][1]  # vertical dispersion
    assert rows[3][2] < 0.20 * rows[3][1]  # accompanying vertical COD
    assert rows[4][2] < 0.10 * rows[4][1]  # coupling probe response

    analysis = Path(__file__).resolve().parents[4] / "analysis" / "DR-RFTrack"
    (analysis / "kubo_step_by_step_result.json").write_text(
        json.dumps(_serializable(result), indent=2), encoding="utf-8"
    )

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.0), constrained_layout=True)
    for axis, plane in zip(axes[0], ("x", "y")):
        item = result["cod"][plane]
        axis.plot(item["before"], label="injected fault")
        axis.plot(item["predicted"], "--", label="model prediction")
        axis.plot(item["after"], ".", ms=3, label="after correction")
        axis.set(title=f"Step {1 if plane == 'x' else 2}: {plane}-COD", xlabel="BPM index", ylabel="orbit [mm]")
        axis.legend(fontsize=8)
    item = result["vertical_dispersion"]
    axes[1, 0].plot(item["dispersion_before"] - item["target"], label="before")
    axes[1, 0].plot(item["dispersion_after"] - item["target"], label="after")
    axes[1, 0].set(title="Step 3: vertical dispersion residual", xlabel="BPM index", ylabel="Dy residual [mm/delta]")
    axes[1, 0].legend(fontsize=8)
    item = result["coupling"]
    count = len(item["before"]) // 2
    axes[1, 1].plot(item["before"][:count], label="ZH1R before")
    axes[1, 1].plot(item["after"][:count], ".", ms=3, label="ZH1R after")
    axes[1, 1].plot(item["before"][count:], label="ZH2R before")
    axes[1, 1].plot(item["after"][count:], ".", ms=3, label="ZH2R after")
    axes[1, 1].set(title="Step 4: coupling probe response", xlabel="BPM index", ylabel="y / kick [mm]")
    axes[1, 1].legend(fontsize=8)
    fig.suptitle("RF-Track implementation of Kubo-style ATF DR correction", fontsize=15)
    fig.savefig(analysis / "kubo_step_by_step_result.png", dpi=180)


if __name__ == "__main__":
    main()
