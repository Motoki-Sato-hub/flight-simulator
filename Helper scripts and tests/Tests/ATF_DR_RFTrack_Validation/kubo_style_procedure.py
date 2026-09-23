"""Kubo (2003)-style ATF DR correction simulation using RF-Track.

This validation deliberately follows the simulation topology in K. Kubo,
Phys. Rev. ST Accel. Beams 6, 092801 (2003): nominal responses are calculated
once without alignment/skew errors, then applied to a separate lattice with
magnet and BPM errors through COD -> vertical COD/dispersion -> coupling.

The radiation-equilibrium evaluation is performed separately by
``ATF_DR_RFTrack_emittance`` after this transverse correction flow is stable.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import RF_Track as rft

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import ATFDRRingCorrection
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import (
    build_atf_dr_lattice,
    get_bpm_nearest_magnet_names,
    get_lattice_metadata,
    set_quadrupole_roll_error,
    set_sextupole_roll_error,
)


# Kubo 2003, Table I.  These are *additional* Gaussian errors; for the 2011
# compatible lattice the published Fig. 1/2 alignment is loaded separately.
SEED = int(os.environ.get("ATF_DR_KUBO_SEED", "2003"))
# ``1`` is the published Table-I strength.  The preliminary orbit search below
# is deliberately run before the three correction stages, as in Sec. III C of
# Kubo (2003); it keeps its steerer settings while the random errors are ramped.
MAGNET_ERROR_SCALE = float(os.environ.get("ATF_DR_KUBO_MAGNET_ERROR_SCALE", "1.0"))
if not 0.0 < MAGNET_ERROR_SCALE <= 1.0:
    raise ValueError("ATF_DR_KUBO_MAGNET_ERROR_SCALE must be in (0, 1]")
# RF-Track reports particle transverse coordinates in mm, but its element
# ``set_offsets`` interface is SI: Kubo's 30 micrometre value is 30e-6 m.
SIGMA_MAGNET_OFFSET_M = 30e-6
SIGMA_MAGNET_ROLL_RAD = 300e-6
# Diagnostic multiplier: zero isolates BPM/readback effects from magnet-offset
# feed-down.  The published Table-I baseline is one.
MAGNET_OFFSET_MULTIPLIER = float(
    os.environ.get("ATF_DR_KUBO_MAGNET_OFFSET_MULTIPLIER", "1.0")
)
if MAGNET_OFFSET_MULTIPLIER < 0.0:
    raise ValueError("ATF_DR_KUBO_MAGNET_OFFSET_MULTIPLIER must be non-negative")
# Diagnostic multiplier only.  ``1`` is Table I; values other than one are
# explicitly recorded so offset- and roll-driven loss of a periodic orbit can
# be separated without changing the published standard condition.
MAGNET_ROLL_MULTIPLIER = float(
    os.environ.get("ATF_DR_KUBO_MAGNET_ROLL_MULTIPLIER", "1.0")
)
if MAGNET_ROLL_MULTIPLIER < 0.0:
    raise ValueError("ATF_DR_KUBO_MAGNET_ROLL_MULTIPLIER must be non-negative")
SIGMA_BPM_OFFSET_MM = 0.3
SIGMA_BPM_ROLL_RAD = 20e-3
DELTA = 3.5e-3
# Match the 0.01-mrad 2011-lattice linear probe used for the nominal map.
PROBE = 1e-5
KUBO_GAIN_FIRST = 0.7
DISPERSION_WEIGHT = 0.05
# The full inverse retains numerical modes that demand multi-mrad kicks for a
# sub-mm BPM-offset residual.  Kubo specifies SVD response inversion but not a
# published singular-value cutoff; retain a well-conditioned subset and report
# this implementation choice rather than claiming it is a Kubo parameter.
RCOND_COD = float(os.environ.get("ATF_DR_KUBO_COD_RCOND", "1e-2"))
RCOND_COUPLING = float(os.environ.get("ATF_DR_KUBO_COUPLING_RCOND", "0.4"))
if not 0.0 <= RCOND_COUPLING < 1.0:
    raise ValueError("ATF_DR_KUBO_COUPLING_RCOND must be in [0, 1)")
# ``sad/operation/lib/cod.n`` and ``coddispersion.n`` use a constrained,
# greedy corrector selection rather than a single unconstrained inverse.  It
# is now the default for the Kubo reproduction path.  ``svd`` is retained as
# an explicitly labelled diagnostic comparison with the earlier RFTrack work.
COD_DISPERSION_SOLVER = os.environ.get(
    "ATF_DR_KUBO_COD_DISPERSION_SOLVER", "sad_greedy"
)
if COD_DISPERSION_SOLVER not in {"sad_greedy", "svd"}:
    raise ValueError("ATF_DR_KUBO_COD_DISPERSION_SOLVER must be sad_greedy or svd")
# ``nominal`` is the Kubo/digital-twin baseline.  ``local_orm`` is a
# diagnostic upper bound: it measures the virtual-machine ORM by the same
# finite corrector changes used on a real machine.  It is intentionally never
# confused with the nominal-model result.
RESPONSE_SOURCE = os.environ.get("ATF_DR_KUBO_RESPONSE_SOURCE", "nominal")
if RESPONSE_SOURCE not in {"nominal", "local_orm"}:
    raise ValueError("ATF_DR_KUBO_RESPONSE_SOURCE must be nominal or local_orm")
LOCAL_ORM_DIAGNOSTIC_ONLY = os.environ.get(
    "ATF_DR_KUBO_LOCAL_ORM_DIAGNOSTIC_ONLY", "0"
) == "1"
if LOCAL_ORM_DIAGNOSTIC_ONLY and RESPONSE_SOURCE != "local_orm":
    raise ValueError("local_orm diagnostic-only mode requires RESPONSE_SOURCE=local_orm")
# SAD cod.n / coddispersion.n: maxk0=0.16138E-02.  This is the allowed
# excursion from the starting setting of each horizontal/vertical steerer.
SAD_MAX_STEERER_KICK_RAD = 0.16138e-2
LATTICE_DATA_PATH = os.environ.get("ATF_DR_KUBO_LATTICE_DATA", "")
LATTICE_LABEL = (
    Path(LATTICE_DATA_PATH).stem.replace("ATF_DR_", "").replace("_RFTrack_lattice", "")
    if LATTICE_DATA_PATH
    else "20111111b"
)
DEFAULT_ALIGNMENT_MAP = (
    Path(__file__).resolve().parents[3]
    / "Interfaces/ATF2/DR_ATF2/ATF_DR_20111111b_Kubo2003_digitized_alignment.json"
)
# The digitised map is deliberately opt-in: it is associated to the 2011
# compatible lattice by longitudinal order, not read from an original 2003
# SAD file.  Set this variable to its path (or ``default``) for the explicit
# Fig. 1/2 comparison experiment.
_alignment_map_value = os.environ.get("ATF_DR_KUBO_ALIGNMENT_MAP", "")
ALIGNMENT_MAP_PATH = (
    DEFAULT_ALIGNMENT_MAP if _alignment_map_value == "default"
    else Path(_alignment_map_value) if _alignment_map_value else None
)
# Kubo used 0.1 mrad.  The 2011 reference lattice loses the horizontal
# periodic orbit at that amplitude because its sextupoles/working point differ
# from the 2003 lattice.  0.01 mrad remains in the linear response regime and
# is the largest symmetric finite-difference probe with a periodic orbit.
# The response is divided by this kick, so it retains Kubo's unit-kick map.
STEERER_RESPONSE_KICK_RAD = float(
    os.environ.get("ATF_DR_KUBO_RESPONSE_KICK_RAD", "1e-5")
)
if STEERER_RESPONSE_KICK_RAD <= 0.0:
    raise ValueError("ATF_DR_KUBO_RESPONSE_KICK_RAD must be positive")
# The full-error lattice has a much smaller locally periodic perturbation
# range than the nominal lattice.  This is only for the optional local-ORM
# diagnostic; it does not alter the nominal Kubo response probe above.
LOCAL_ORM_PROBE_RAD = float(
    os.environ.get("ATF_DR_KUBO_LOCAL_ORM_PROBE_RAD", "1e-6")
)
if LOCAL_ORM_PROBE_RAD <= 0.0:
    raise ValueError("ATF_DR_KUBO_LOCAL_ORM_PROBE_RAD must be positive")
RESPONSE_CACHE_VERSION = 8
# SAD's skew-correction configuration selects this pair; it also satisfies the
# phase-separation criterion stated in Kubo (2003), Sec. IV E.
PROBES = ("ZH28R", "ZH30R")
SKEW_FAMILY = os.environ.get("ATF_DR_KUBO_SKEW_FAMILY", "SF1R")
if SKEW_FAMILY not in {"SD1R", "SF1R"}:
    raise ValueError("ATF_DR_KUBO_SKEW_FAMILY must be 'SD1R' or 'SF1R'")
LATTICE_SOURCE_SHA256 = get_lattice_metadata(LATTICE_DATA_PATH or None)["source_sha256"]
PRELIMINARY_X_LIMIT_MM = 2.0
PRELIMINARY_Y_LIMIT_MM = 1.0
# The published rough-COD acceptance limits are deliberately loose.  On the
# current 2011 lattice, waiting for those limits reaches a lost-orbit branch
# before feedback starts.  These tighter continuation thresholds preserve the
# periodic branch without weakening Table-I errors or changing final limits.
PRELIMINARY_CONTINUATION_X_MM = float(
    os.environ.get("ATF_DR_KUBO_CONTINUATION_X_MM", "0.5")
)
PRELIMINARY_CONTINUATION_Y_MM = float(
    os.environ.get("ATF_DR_KUBO_CONTINUATION_Y_MM", "0.5")
)
if PRELIMINARY_CONTINUATION_X_MM <= 0 or PRELIMINARY_CONTINUATION_Y_MM <= 0:
    raise ValueError("Kubo continuation limits must be positive")
PRELIMINARY_MAX_FEEDBACK_STEPS = 20
ROUGH_RESPONSE_KICK_RAD = 1e-6
# RF-Track coordinates are mm/mrad.  The default corresponds to 0.1 nm or
# 0.1 nrad.  A larger value is only an explicit numerical-solver study; it is
# recorded in the JSON and never confused with a BPM measurement error.
CLOSED_ORBIT_TOLERANCE = float(
    os.environ.get("ATF_DR_KUBO_CLOSED_ORBIT_TOLERANCE", "1e-7")
)
if CLOSED_ORBIT_TOLERANCE <= 0.0:
    raise ValueError("ATF_DR_KUBO_CLOSED_ORBIT_TOLERANCE must be positive")


def _load_published_alignment_mm():
    """Load the vector-digitised Kubo Fig. 1/2 main-magnet offsets, if selected."""
    if ALIGNMENT_MAP_PATH is None:
        return {}
    with ALIGNMENT_MAP_PATH.open(encoding="utf-8") as stream:
        payload = json.load(stream)
    rows = payload.get("magnets", [])
    if len(rows) != 204:
        raise ValueError(f"Expected 204 published main-magnet offsets in {ALIGNMENT_MAP_PATH}")
    return {
        row["name"]: np.asarray((row["dx_mm"], row["dy_mm"]), dtype=float) * 1e-3
        for row in rows
    }


PUBLISHED_ALIGNMENT_M = _load_published_alignment_mm()


def _build_lattice():
    """Build the selected nominal lattice without changing the 2011 default."""
    return build_atf_dr_lattice(
        lattice_data_path=LATTICE_DATA_PATH or None,
    )


def _bpm_reference_offsets_mm(lattice, magnet_offsets_mm, bpm_local_offsets_mm):
    """Return each BPM's offset from the design coordinates in Kubo Eq. (3)."""
    nearest = get_bpm_nearest_magnet_names(LATTICE_DATA_PATH or None)
    bpm_names = _paper_bpm_names(ATFDRRingCorrection(lattice))
    if not set(bpm_names).issubset(nearest):
        raise RuntimeError("SAD-derived BPM-to-magnet mapping does not match RF-Track lattice")
    return np.asarray([
        magnet_offsets_mm[name] for name in (nearest[bpm] for bpm in bpm_names)
    ]) + bpm_local_offsets_mm


def _paper_bpm_names(machine):
    """Select the 96 arc-cell BPMs used by Kubo, excluding R/CAV monitors."""
    names = tuple(name for name in machine.bpm_names if name.startswith("MB"))
    if len(names) != 96:
        raise RuntimeError(f"Expected the paper's 96 BPMs, got {len(names)}")
    return names


def _paper_corrector_names(machine, plane):
    """Select the 48 H / 50 V operational steerers in Kubo's lattice."""
    names = tuple(machine.get_corrector_names(plane))
    if plane == "x":
        names = tuple(name for name in names if name not in {"ZH100R", "ZH101R"})
        expected = 48
    elif plane == "y":
        names = tuple(name for name in names if name != "ZV100R")
        expected = 50
    else:
        raise ValueError("plane must be 'x' or 'y'")
    if len(names) != expected:
        raise RuntimeError(f"Expected {expected} Kubo {plane}-plane steerers, got {len(names)}")
    return names


def _rms(values):
    return float(np.sqrt(np.mean(np.asarray(values, dtype=float) ** 2)))


def _one(correction, name):
    return correction._single_element(name)


def _add_correctors(correction, names, plane, values, *, max_abs_kick_rad=None):
    for name, value in zip(names, values):
        if value == 0.0:
            continue
        element = _one(correction, name)
        kick = correction._get_corrector_kick(element)
        kick[plane] += float(value)
        if max_abs_kick_rad is not None:
            kick[plane] = np.clip(
                kick[plane], -float(max_abs_kick_rad), float(max_abs_kick_rad)
            )
        correction._set_corrector_kick(element, kick)


def _apply_magnet_errors(lattice, scale):
    """Apply alignment errors and explicit multipole-field roll errors.

    RF-Track placement rotations do not create the skew field of a rolled
    quadrupole, hence its field-roll equivalent is set separately.  The
    sextupole field-roll equivalent is also included explicitly.  Combined-
    function bend rolls remain a separately tracked limitation because their
    dipole, quadrupole and edge-field rotations cannot be replaced by one
    thin multipole without a map-level validation.
    """
    rng = np.random.default_rng(SEED)
    counts = {"quadrupoles": 0, "sextupoles": 0, "bends": 0, "correctors": 0}
    magnet_offsets_mm = {}
    for element in lattice["*"]:
        if isinstance(element, (rft.Quadrupole, rft.Sextupole, rft.SBend)):
            random_dx, random_dy = rng.normal(
                0.0, SIGMA_MAGNET_OFFSET_M * scale * MAGNET_OFFSET_MULTIPLIER, 2
            )
            roll = rng.normal(
                0.0, SIGMA_MAGNET_ROLL_RAD * scale * MAGNET_ROLL_MULTIPLIER
            )
            base_dx, base_dy = PUBLISHED_ALIGNMENT_M.get(
                element.get_name(), np.zeros(2)
            )
            dx, dy = base_dx + random_dx, base_dy + random_dy
            element.set_offsets(dx, dy, 0.0, 0.0, 0.0, 0.0, "center")
            magnet_offsets_mm[element.get_name()] = np.array((dx, dy)) * 1e3
            if isinstance(element, rft.Quadrupole):
                set_quadrupole_roll_error(lattice, element.get_name(), roll)
                counts["quadrupoles"] += 1
            elif isinstance(element, rft.Sextupole):
                set_sextupole_roll_error(lattice, element.get_name(), roll)
                counts["sextupoles"] += 1
            elif isinstance(element, rft.SBend):
                counts["bends"] += 1
    if PUBLISHED_ALIGNMENT_M and len(magnet_offsets_mm) != len(PUBLISHED_ALIGNMENT_M):
        raise RuntimeError("Published alignment map does not match the selected RF-Track lattice")
    return counts, magnet_offsets_mm


def _restore_correctors(machine, names, plane, values):
    """Restore a trial steerer change after an unhelpful feedback step."""
    _add_correctors(machine, names, plane, -np.asarray(values, dtype=float))


def _set_corrector_plane(machine, names, plane, values):
    """Set physical corrector kicks exactly, preserving the other plane."""
    for name, value in zip(names, values):
        element = _one(machine, name)
        kick = machine._get_corrector_kick(element)
        kick[plane] = float(value)
        machine._set_corrector_kick(element, kick)


def _apply_dispersion_step_with_backtracking(
    machine, names, proposal, gain, guess, offsets, rolls, target_dy
):
    """Apply a nonlinear-safe vertical COD/Dy correction step.

    The nominal response is deliberately reused, as in the Kubo procedure.
    On the 2011 lattice a full linear step can leave the nearby periodic-orbit
    branch at Table-I strength.  This line search changes only the applied
    fraction of that same response solution and accepts a step only when the
    measured stacked COD/Dy objective decreases.
    """
    _, _, before_orbit = _orbit(machine, guess, offsets, rolls)
    _, _, before_dispersion = _dispersion(machine, guess, offsets, rolls)
    before = np.concatenate((
        before_orbit[:, 1],
        DISPERSION_WEIGHT * (before_dispersion[:, 1] - target_dy),
    ))
    before_norm = float(np.dot(before, before))
    original = np.asarray([
        machine._get_corrector_kick(_one(machine, name))[1]
        for name in names
    ])
    # First retain Kubo's requested 0.7/full step; subsequent values are a
    # numerical safeguard for the current (not historical-2003) lattice.
    fractions = tuple(dict.fromkeys(
        float(gain) * value for value in (1.0, 0.7, 0.5, 0.3, 0.1, 0.03)
    ))
    for fraction in fractions:
        candidate = np.clip(
            original + fraction * np.asarray(proposal, dtype=float),
            -SAD_MAX_STEERER_KICK_RAD,
            SAD_MAX_STEERER_KICK_RAD,
        )
        _set_corrector_plane(machine, names, 1, candidate)
        try:
            candidate_guess, _, candidate_orbit = _orbit(
                machine, guess, offsets, rolls
            )
            candidate_guess, _, candidate_dispersion = _dispersion(
                machine, candidate_guess, offsets, rolls
            )
            residual = np.concatenate((
                candidate_orbit[:, 1],
                DISPERSION_WEIGHT * (candidate_dispersion[:, 1] - target_dy),
            ))
            if float(np.dot(residual, residual)) < before_norm * (1.0 - 1e-10):
                return candidate_guess, fraction
        except Exception:
            pass
        _set_corrector_plane(machine, names, 1, original)
    _set_corrector_plane(machine, names, 1, original)
    return guess, 0.0


def _local_orbit_response(machine, guess, offsets, rolls, names, plane):
    """Measure an error-lattice ORM in the BPM readback coordinate system."""
    plane_index = 0 if plane == "x" else 1
    response = np.empty((len(offsets), len(names)))
    for column, name in enumerate(names):
        element = _one(machine, name)
        original = machine._get_corrector_kick(element)
        plus = original.copy()
        minus = original.copy()
        plus[plane_index] += LOCAL_ORM_PROBE_RAD
        minus[plane_index] -= LOCAL_ORM_PROBE_RAD
        try:
            machine._set_corrector_kick(element, plus)
            _, _, positive = _orbit(machine, guess, offsets, rolls)
            machine._set_corrector_kick(element, minus)
            _, _, negative = _orbit(machine, guess, offsets, rolls)
            response[:, column] = (
                positive[:, plane_index] - negative[:, plane_index]
            ) / (2.0 * LOCAL_ORM_PROBE_RAD)
        finally:
            machine._set_corrector_kick(element, original)
    return response


def _local_vertical_dispersion_response(machine, guess, offsets, rolls, names):
    r"""Measure \(d(D_y^{BPM})/d\theta_y\) on the virtual error lattice."""
    response = np.empty((len(offsets), len(names)))
    for column, name in enumerate(names):
        element = _one(machine, name)
        original = machine._get_corrector_kick(element)
        plus = original.copy()
        minus = original.copy()
        plus[1] += LOCAL_ORM_PROBE_RAD
        minus[1] -= LOCAL_ORM_PROBE_RAD
        try:
            machine._set_corrector_kick(element, plus)
            _, _, positive = _dispersion(machine, guess, offsets, rolls)
            machine._set_corrector_kick(element, minus)
            _, _, negative = _dispersion(machine, guess, offsets, rolls)
            response[:, column] = (
                positive[:, 1] - negative[:, 1]
            ) / (2.0 * LOCAL_ORM_PROBE_RAD)
        finally:
            machine._set_corrector_kick(element, original)
    return response


def _local_orm_response_set(machine, guess, offsets, rolls, nominal):
    """Return a diagnostic response set measured on the virtual machine.

    This deliberately keeps the Kubo nominal response untouched.  The returned
    set is only used when ``ATF_DR_KUBO_RESPONSE_SOURCE=local_orm`` to quantify
    the benefit of an experimentally measured ORM over the digital-twin ORM.
    """
    response = dict(nominal)
    response["x"] = _local_orbit_response(
        machine, guess, offsets, rolls, response["x_names"], "x"
    )
    response["y_orbit"] = _local_orbit_response(
        machine, guess, offsets, rolls, response["y_names"], "y"
    )
    response["y_dispersion"] = _local_vertical_dispersion_response(
        machine, guess, offsets, rolls, response["y_names"]
    )
    diagnostics = {}
    for key in ("x", "y_orbit", "y_dispersion"):
        model = np.asarray(nominal[key], dtype=float)
        measured = np.asarray(response[key], dtype=float)
        model_norm = float(np.linalg.norm(model))
        measured_norm = float(np.linalg.norm(measured))
        overlap = float(np.vdot(model, measured).real)
        diagnostics[key] = {
            "shape": list(model.shape),
            "relative_frobenius_difference": float(
                np.linalg.norm(measured - model) / model_norm
            ) if model_norm else float("nan"),
            "matrix_cosine_similarity": overlap / (model_norm * measured_norm),
            "nominal_rank_rcond_1e-2": _svd_rank(model, RCOND_COD),
            "local_rank_rcond_1e-2": _svd_rank(measured, RCOND_COD),
        }
    return response, diagnostics


def _preliminary_orbit_search(response):
    """Reproduce the preliminary all-steerer orbit search of Kubo Sec. III C.

    The alignment errors are increased adiabatically.  At every error level a
    small response-matrix feedback is accepted only if it reduces the *maximum*
    BPM readback and a periodic orbit still exists.  Crucially, the accumulated
    steerer strengths are applied when the lattice is rebuilt at the next level.
    This is the part that was absent from the earlier pilot.
    """
    rng = np.random.default_rng(SEED + 1)
    n_bpms = len(_paper_bpm_names(ATFDRRingCorrection(_build_lattice())))
    offsets = rng.normal(0.0, SIGMA_BPM_OFFSET_MM, (n_bpms, 2))
    rolls = rng.normal(0.0, SIGMA_BPM_ROLL_RAD, n_bpms)
    x_commands = np.zeros(len(response["x_names"]))
    y_commands = np.zeros(len(response["y_names"]))
    guess = np.zeros(4)
    fraction = 0.0
    increment = 0.05
    records = []
    counts = None

    while fraction < 1.0 - 1e-12:
        trial_fraction = min(1.0, fraction + increment)
        lattice = _build_lattice()
        counts, magnet_offsets_mm = _apply_magnet_errors(
            lattice, MAGNET_ERROR_SCALE * trial_fraction
        )
        machine = ATFDRRingCorrection(lattice)
        if len(_paper_bpm_names(machine)) != len(offsets):
            raise RuntimeError("Unexpected RF-Track BPM count in Kubo procedure")
        _add_correctors(machine, response["x_names"], 0, x_commands)
        _add_correctors(machine, response["y_names"], 1, y_commands)
        reference_offsets = _bpm_reference_offsets_mm(
            lattice, magnet_offsets_mm, offsets
        )

        try:
            guess, _, measured = _orbit(machine, guess, reference_offsets, rolls)
        except Exception as error:
            # A smaller error increment is the numerical analogue of bringing
            # the machine up gradually.  It does not silently discard a fault.
            increment *= 0.5
            if increment < 1e-4:
                raise RuntimeError(
                    "Preliminary orbit search cannot continue at "
                    f"{trial_fraction:.5f} of the requested Table-I errors"
                ) from error
            continue

        accepted_steps = 0
        for _ in range(PRELIMINARY_MAX_FEEDBACK_STEPS):
            maximum = np.max(np.abs(measured), axis=0)
            if (
                maximum[0] < PRELIMINARY_CONTINUATION_X_MM
                and maximum[1] < PRELIMINARY_CONTINUATION_Y_MM
            ):
                break
            dx = _least_squares(response["x"], measured[:, 0], RCOND_COD)
            dy = _least_squares(response["y_orbit"], measured[:, 1], RCOND_COD)
            accepted = False
            for gain in (0.30, 0.15, 0.075, 0.03, 0.01):
                _add_correctors(machine, response["x_names"], 0, gain * dx)
                _add_correctors(machine, response["y_names"], 1, gain * dy)
                try:
                    trial_guess, _, trial_measured = _orbit(
                        machine, guess, reference_offsets, rolls
                    )
                    trial_maximum = np.max(np.abs(trial_measured), axis=0)
                except Exception:
                    trial_maximum = np.array((np.inf, np.inf))
                if np.max(
                    trial_maximum / np.array((
                        PRELIMINARY_CONTINUATION_X_MM,
                        PRELIMINARY_CONTINUATION_Y_MM,
                    ))
                ) < np.max(
                    maximum / np.array((
                        PRELIMINARY_CONTINUATION_X_MM,
                        PRELIMINARY_CONTINUATION_Y_MM,
                    ))
                ):
                    x_commands += gain * dx
                    y_commands += gain * dy
                    guess, measured = trial_guess, trial_measured
                    accepted_steps += 1
                    accepted = True
                    break
                _restore_correctors(machine, response["x_names"], 0, gain * dx)
                _restore_correctors(machine, response["y_names"], 1, gain * dy)
            if not accepted:
                break

        maximum = np.max(np.abs(measured), axis=0)
        records.append({
            "error_fraction": trial_fraction,
            "max_measured_x_mm": float(maximum[0]),
            "max_measured_y_mm": float(maximum[1]),
            "accepted_feedback_steps": accepted_steps,
        })
        print(
            f"preliminary fraction={trial_fraction:.3f}, "
            f"max BPM=({maximum[0]:.3f}, {maximum[1]:.3f}) mm, "
            f"feedback={accepted_steps}",
            flush=True,
        )
        fraction = trial_fraction
        increment = min(0.10, increment * 1.35)

    # Kubo's rough-COD stage uses SAD's free-parameter search on the *error
    # lattice*, not the nominal response matrix used by the following stages.
    # Finish the adiabatic ramp with the corresponding local virtual-machine
    # response.  The very small 1-urad difference is only an optimizer probe;
    # it is not the 0.1-mrad nominal-response convention of Sec. IV C.
    # Disabled by default until the 2011 lattice's local closed-orbit branch
    # selection is validated.  The original adiabatic search above remains the
    # reproducible baseline; setting this environment flag enables the more
    # expensive actual-response refinement for solver development.
    for iteration in range(
        PRELIMINARY_MAX_FEEDBACK_STEPS
        if os.environ.get("ATF_DR_KUBO_LOCAL_ROUGH", "0") == "1"
        else 0
    ):
        guess, _, measured = _orbit(machine, guess, reference_offsets, rolls)
        maximum = np.max(np.abs(measured), axis=0)
        if (
            maximum[0] < PRELIMINARY_X_LIMIT_MM
            and maximum[1] < PRELIMINARY_Y_LIMIT_MM
        ):
            break
        local_x = machine.compute_orbit_response(
            "x", corrector_names=response["x_names"],
            perturbation=ROUGH_RESPONSE_KICK_RAD,
            initial_coordinates=guess,
            bpm_names=_paper_bpm_names(machine),
            central_difference=False,
        ).matrix
        local_y = machine.compute_orbit_response(
            "y", corrector_names=response["y_names"],
            perturbation=ROUGH_RESPONSE_KICK_RAD,
            initial_coordinates=guess,
            bpm_names=_paper_bpm_names(machine),
            central_difference=False,
        ).matrix
        dx = _least_squares(local_x, measured[:, 0], RCOND_COD)
        dy = _least_squares(local_y, measured[:, 1], RCOND_COD)
        accepted = False
        for gain in (1.0, 0.7, 0.35, 0.15, 0.07, 0.03):
            _add_correctors(machine, response["x_names"], 0, gain * dx)
            _add_correctors(machine, response["y_names"], 1, gain * dy)
            try:
                candidate_guess, _, candidate = _orbit(machine, guess, reference_offsets, rolls)
                candidate_maximum = np.max(np.abs(candidate), axis=0)
            except Exception:
                candidate_maximum = np.array((np.inf, np.inf))
            if np.max(candidate_maximum / np.array((PRELIMINARY_X_LIMIT_MM, PRELIMINARY_Y_LIMIT_MM))) < np.max(maximum / np.array((PRELIMINARY_X_LIMIT_MM, PRELIMINARY_Y_LIMIT_MM))):
                x_commands += gain * dx
                y_commands += gain * dy
                guess, measured = candidate_guess, candidate
                accepted = True
                records.append({
                    "error_fraction": 1.0,
                    "phase": "local_virtual_machine_search",
                    "max_measured_x_mm": float(candidate_maximum[0]),
                    "max_measured_y_mm": float(candidate_maximum[1]),
                    "accepted_feedback_steps": iteration + 1,
                })
                break
            _restore_correctors(machine, response["x_names"], 0, gain * dx)
            _restore_correctors(machine, response["y_names"], 1, gain * dy)
        if not accepted:
            break

    return machine, guess, reference_offsets, rolls, counts, records, x_commands, y_commands


def _make_machine(response):
    """Create the full-error virtual machine after Kubo preliminary search."""
    (
        machine, guess, offsets, rolls, counts, records, x_commands, y_commands
    ) = _preliminary_orbit_search(response)
    return machine, guess, offsets, rolls, counts, records, x_commands, y_commands


def _bpm_measurement(positions, offsets, rolls):
    """BPM offset and roll measurement model; inputs/outputs are mm."""
    values = np.asarray(positions, dtype=float) - np.asarray(offsets, dtype=float)
    cosine = np.cos(rolls)[:, None]
    sine = np.sin(rolls)[:, None]
    rotation = np.concatenate((
        cosine * values[:, :1] + sine * values[:, 1:2],
        -sine * values[:, :1] + cosine * values[:, 1:2],
    ), axis=1)
    return rotation


def _orbit(machine, guess, offsets, rolls):
    closed = machine.find_closed_orbit(
        initial_coordinates=guess,
        bpm_names=_paper_bpm_names(machine),
        max_iterations=30,
        tolerance=CLOSED_ORBIT_TOLERANCE,
    )
    return closed.initial_coordinates, closed.bpm_positions, _bpm_measurement(
        closed.bpm_positions, offsets, rolls
    )


def _dispersion(machine, guess, offsets, rolls):
    bpm_names = _paper_bpm_names(machine)
    reference = machine.find_closed_orbit(
        initial_coordinates=guess, bpm_names=bpm_names, max_iterations=30,
        tolerance=CLOSED_ORBIT_TOLERANCE,
    )
    plus = machine.find_closed_orbit(
        momentum_mev_c=machine.momentum_mev_c * (1.0 + DELTA),
        initial_coordinates=reference.initial_coordinates,
        bpm_names=bpm_names,
        max_iterations=30,
        tolerance=CLOSED_ORBIT_TOLERANCE,
    )
    minus = machine.find_closed_orbit(
        momentum_mev_c=machine.momentum_mev_c * (1.0 - DELTA),
        initial_coordinates=reference.initial_coordinates,
        bpm_names=bpm_names,
        max_iterations=30,
        tolerance=CLOSED_ORBIT_TOLERANCE,
    )
    true = (plus.bpm_positions - minus.bpm_positions) / (2.0 * DELTA)
    measured = (
        _bpm_measurement(plus.bpm_positions, offsets, rolls)
        - _bpm_measurement(minus.bpm_positions, offsets, rolls)
    ) / (2.0 * DELTA)
    return reference.initial_coordinates, true, measured


def _coupling(machine, guess, offsets, rolls):
    """Two-probe Kubo Cxy measurement, including BPM rotation but not offset."""
    bpm_names = _paper_bpm_names(machine)
    baseline = machine.find_closed_orbit(
        initial_coordinates=guess, bpm_names=bpm_names, max_iterations=30,
        tolerance=CLOSED_ORBIT_TOLERANCE,
    )
    horizontal, vertical = [], []
    for name in PROBES:
        element = _one(machine, name)
        initial = machine._get_corrector_kick(element)
        plus, minus = initial.copy(), initial.copy()
        plus[0] += PROBE
        minus[0] -= PROBE
        machine._set_corrector_kick(element, plus)
        positive = machine.find_closed_orbit(
            initial_coordinates=baseline.initial_coordinates,
            bpm_names=bpm_names,
            max_iterations=30,
            tolerance=CLOSED_ORBIT_TOLERANCE,
        )
        machine._set_corrector_kick(element, minus)
        negative = machine.find_closed_orbit(
            initial_coordinates=baseline.initial_coordinates,
            bpm_names=bpm_names,
            max_iterations=30,
            tolerance=CLOSED_ORBIT_TOLERANCE,
        )
        machine._set_corrector_kick(element, initial)
        true = (positive.bpm_positions - negative.bpm_positions) / (2.0 * PROBE)
        measured = (
            _bpm_measurement(positive.bpm_positions, offsets, rolls)
            - _bpm_measurement(negative.bpm_positions, offsets, rolls)
        ) / (2.0 * PROBE)
        horizontal.append(measured[:, 0])
        vertical.append(measured[:, 1])
    horizontal, vertical = np.asarray(horizontal), np.asarray(vertical)
    cxy = float(np.sqrt(np.mean(vertical**2) / np.mean(horizontal**2)))
    return baseline.initial_coordinates, vertical.reshape(-1), cxy


def _load_or_build_responses(cache):
    """Kubo computes one error-free response set shared by all seeds."""
    if cache.exists():
        data = np.load(cache, allow_pickle=True)
        if (
            int(data.get("cache_version", -1)) == RESPONSE_CACHE_VERSION
            and str(data.get("skew_family", "")) == SKEW_FAMILY
            and str(data.get("lattice_source_sha256", "")) == LATTICE_SOURCE_SHA256
        ):
            return {
                "x_names": tuple(data["x_names"].tolist()),
                "x": data["x"],
                "y_names": tuple(data["y_names"].tolist()),
                "y_orbit": data["y_orbit"],
                "y_dispersion": data["y_dispersion"],
                "skew_names": tuple(data["skew_names"].tolist()),
                "coupling": data["coupling"],
            }

    # Do not convert older caches: RF-Track Corrector.set_kick uses mrad while
    # the correction API is deliberately in rad.  A response made before that
    # boundary conversion has a factor-of-1000 actuator error and is invalid.

    model = ATFDRRingCorrection(_build_lattice())
    bpm_names = _paper_bpm_names(model)
    x_names = _paper_corrector_names(model, "x")
    y_names = _paper_corrector_names(model, "y")
    x = model.compute_orbit_response(
        "x", corrector_names=x_names, bpm_names=bpm_names,
        perturbation=STEERER_RESPONSE_KICK_RAD
    )
    y = model.compute_dispersion_response(
        "y", corrector_names=y_names, bpm_names=bpm_names,
        relative_momentum_step=DELTA,
        perturbation=STEERER_RESPONSE_KICK_RAD,
    )
    # Kubo 2003 used all 34 correctors of one sextupole family.  The archived
    # SAD configuration selects SF1R; this can be made explicit in reruns.
    skew_names = tuple(
        name for name in model.get_skew_corrector_names()
        if name.startswith(SKEW_FAMILY)
    )
    if len(skew_names) != 34:
        raise RuntimeError(f"Expected 34 {SKEW_FAMILY} skew correctors")
    coupling = model.compute_coupling_response(
        probe_corrector_names=PROBES,
        skew_corrector_names=skew_names,
        bpm_names=bpm_names,
        probe_perturbation=STEERER_RESPONSE_KICK_RAD,
        skew_perturbation=1e-5,
    )
    np.savez_compressed(
        cache,
        cache_version=RESPONSE_CACHE_VERSION,
        x_names=np.asarray(x_names), x=x.matrix,
        y_names=np.asarray(y_names), y_orbit=y.orbit_matrix, y_dispersion=y.dispersion_matrix,
        skew_names=np.asarray(skew_names), coupling=coupling.matrix,
        skew_family=np.asarray(SKEW_FAMILY),
        lattice_source_sha256=np.asarray(LATTICE_SOURCE_SHA256),
    )
    return {
        "x_names": x_names, "x": x.matrix, "y_names": y_names,
        "y_orbit": y.orbit_matrix, "y_dispersion": y.dispersion_matrix,
        "skew_names": skew_names, "coupling": coupling.matrix,
    }


def _least_squares(matrix, residual, rcond):
    return ATFDRRingCorrection._svd_solve(matrix, -residual, rcond)[0]


def _sad_greedy_least_squares(matrix, residual, max_abs_delta):
    """SAD ``NewCorrectCOD``-style bounded, greedy response inversion.

    ``sad/operation/lib/cod.n`` and ``coddispersion.n`` append one candidate
    steerer at a time, retain the candidate with the smallest residual, and
    enforce ``maxk0=0.16138E-02``.  The SAD code treats the live BPM data as
    its right-hand side; here ``residual`` is the same measurement minus the
    desired value, hence the required correction has the opposite sign.

    The implementation is intentionally deterministic: ties retain the first
    lattice-order corrector.  It exposes the selected count to the result JSON
    so that this is auditable against the SAD configuration.
    """
    matrix = np.asarray(matrix, dtype=float)
    target = -np.asarray(residual, dtype=float)
    if matrix.ndim != 2 or target.shape != (matrix.shape[0],):
        raise ValueError("response matrix and residual have incompatible shapes")
    if max_abs_delta <= 0.0:
        raise ValueError("max_abs_delta must be positive")

    selected = []
    remaining = list(range(matrix.shape[1]))
    best_full = np.zeros(matrix.shape[1])
    best_norm = float(np.dot(target, target))
    # The SAD caller can request up to all available steerers.  A candidate
    # that has already saturated is removed, just as SAD removes it from its
    # active set before solving the remaining system.
    while remaining:
        chosen = None
        for candidate in remaining:
            trial = selected + [candidate]
            solution, *_ = np.linalg.lstsq(matrix[:, trial], target, rcond=None)
            solution = np.clip(solution, -max_abs_delta, max_abs_delta)
            trial_residual = target - matrix[:, trial] @ solution
            norm = float(np.dot(trial_residual, trial_residual))
            if chosen is None or norm < chosen[0]:
                chosen = (norm, candidate, solution)
        assert chosen is not None
        norm, candidate, solution = chosen
        # Further columns cannot improve a numerically exact fit.  In the
        # non-exact case keep going through the full actuator set, matching
        # the all-steerer setting used for the Table-I study.
        if norm > best_norm * (1.0 - 1e-13):
            break
        selected.append(candidate)
        remaining.remove(candidate)
        best_norm = norm
        best_full = np.zeros(matrix.shape[1])
        best_full[selected] = solution
        if norm <= 1e-24:
            break
    return best_full, tuple(selected), best_norm


def _cod_dispersion_solution(matrix, residual):
    if COD_DISPERSION_SOLVER == "svd":
        solution = _least_squares(matrix, residual, RCOND_COD)
        return solution, tuple(np.flatnonzero(np.abs(solution) > 0.0)), None
    return _sad_greedy_least_squares(
        matrix, residual, SAD_MAX_STEERER_KICK_RAD
    )


def _svd_rank(matrix, rcond):
    return ATFDRRingCorrection._svd_solve(
        matrix, np.zeros(matrix.shape[0]), rcond
    )[2]


def _stage_metrics(machine, guess, offsets, rolls):
    guess, true_orbit, measured_orbit = _orbit(machine, guess, offsets, rolls)
    guess, true_dispersion, measured_dispersion = _dispersion(machine, guess, offsets, rolls)
    guess, coupling_signal, cxy = _coupling(machine, guess, offsets, rolls)
    return guess, {
        "closed_orbit_initial_coordinates_mm_mrad": guess.tolist(),
        "true_cod_x_mm": _rms(true_orbit[:, 0]),
        "true_cod_y_mm": _rms(true_orbit[:, 1]),
        "measured_cod_x_mm": _rms(measured_orbit[:, 0]),
        "measured_cod_y_mm": _rms(measured_orbit[:, 1]),
        # The 96 arc-cell BPMs used by Kubo are selected in _dispersion.
        "true_dy_all_bpm_mm": _rms(true_dispersion[:, 1]),
        "measured_dy_all_bpm_mm": _rms(measured_dispersion[:, 1]),
        "measured_cxy": cxy,
    }


def _actuator_snapshot(machine, response):
    """Return physical rad/K1L settings for radiation-model replay."""
    return {
        "x_kick_rad": [
            float(machine._get_corrector_kick(_one(machine, name))[0])
            for name in response["x_names"]
        ],
        "y_kick_rad": [
            float(machine._get_corrector_kick(_one(machine, name))[1])
            for name in response["y_names"]
        ],
        "skew_k1l": [
            float(machine.get_skew_strength(name))
            for name in response["skew_names"]
        ],
    }


def run_case(response):
    (
        machine, guess, offsets, rolls, counts, preliminary, x_commands, y_commands
    ) = _make_machine(response)
    local_orm_diagnostics = None
    if RESPONSE_SOURCE == "local_orm":
        print("Diagnostic: measuring local error-lattice ORM", flush=True)
        response, local_orm_diagnostics = _local_orm_response_set(
            machine, guess, offsets, rolls, response
        )
    if LOCAL_ORM_DIAGNOSTIC_ONLY:
        return {
            "status": "local ORM mismatch diagnostic; no correction applied",
            "seed": SEED,
            "nominal_lattice_label": LATTICE_LABEL,
            "magnet_error_scale_of_kubo_table_i": MAGNET_ERROR_SCALE,
            "magnet_offset_multiplier": MAGNET_OFFSET_MULTIPLIER,
            "magnet_roll_multiplier": MAGNET_ROLL_MULTIPLIER,
            "response_settings": {
                "source": RESPONSE_SOURCE,
                "local_orm_probe_mrad": LOCAL_ORM_PROBE_RAD * 1e3,
            },
            "preliminary_final_bpm_max_mm": preliminary[-1],
            "local_orm_mismatch": local_orm_diagnostics,
        }
    stages = {}
    actuator_settings = {}

    print("Kubo C: COD", flush=True)
    # Kubo C: all steerers, 0.7 of calculated correction then full correction.
    cod_selected_counts = []
    for gain in (KUBO_GAIN_FIRST, 1.0):
        guess, _, measured = _orbit(machine, guess, offsets, rolls)
        dx, x_selected, _ = _cod_dispersion_solution(response["x"], measured[:, 0])
        dy, y_selected, _ = _cod_dispersion_solution(response["y_orbit"], measured[:, 1])
        cod_selected_counts.append({"horizontal": len(x_selected), "vertical": len(y_selected)})
        _add_correctors(
            machine, response["x_names"], 0, gain * dx,
            max_abs_kick_rad=SAD_MAX_STEERER_KICK_RAD,
        )
        _add_correctors(
            machine, response["y_names"], 1, gain * dy,
            max_abs_kick_rad=SAD_MAX_STEERER_KICK_RAD,
        )
    guess, stages["after_cod"] = _stage_metrics(machine, guess, offsets, rolls)
    actuator_settings["after_cod"] = _actuator_snapshot(machine, response)

    print("Kubo D: vertical COD and dispersion", flush=True)
    # Kubo D: all vertical steerers; vertical COD and Dy are fitted together.
    model_target_dy = np.zeros(len(_paper_bpm_names(machine)))
    matrix = np.vstack((response["y_orbit"], DISPERSION_WEIGHT * response["y_dispersion"]))
    dispersion_selected_counts = []
    dispersion_applied_gains = []
    for gain in (KUBO_GAIN_FIRST, 1.0):
        guess, _, measured_orbit = _orbit(machine, guess, offsets, rolls)
        guess, _, measured_dispersion = _dispersion(machine, guess, offsets, rolls)
        residual = np.concatenate((
            measured_orbit[:, 1],
            DISPERSION_WEIGHT * (measured_dispersion[:, 1] - model_target_dy),
        ))
        dy, selected, _ = _cod_dispersion_solution(matrix, residual)
        dispersion_selected_counts.append(len(selected))
        guess, applied_gain = _apply_dispersion_step_with_backtracking(
            machine, response["y_names"], dy, gain, guess, offsets, rolls,
            model_target_dy,
        )
        dispersion_applied_gains.append(applied_gain)
    guess, stages["after_cod_dispersion"] = _stage_metrics(machine, guess, offsets, rolls)
    actuator_settings["after_cod_dispersion"] = _actuator_snapshot(machine, response)

    print("Kubo E: coupling", flush=True)
    # Kubo E: two horizontal probes and all SD1R-family skew correctors.
    guess, measured_signal, _ = _coupling(machine, guess, offsets, rolls)
    skew_delta = _least_squares(response["coupling"], measured_signal, RCOND_COUPLING)
    for name, value in zip(response["skew_names"], skew_delta):
        machine.set_skew_strength(name, machine.get_skew_strength(name) + value)
    guess, stages["after_coupling"] = _stage_metrics(machine, guess, offsets, rolls)
    actuator_settings["after_coupling"] = _actuator_snapshot(machine, response)

    return {
        "seed": SEED,
        "nominal_lattice_label": LATTICE_LABEL,
        "nominal_lattice_data_path": LATTICE_DATA_PATH or "ATF_DR_RFTrack_lattice.json (default)",
        "magnet_error_scale_of_kubo_table_i": MAGNET_ERROR_SCALE,
        "errors": {
            "magnet_offset_um": (
                SIGMA_MAGNET_OFFSET_M * MAGNET_ERROR_SCALE
                * MAGNET_OFFSET_MULTIPLIER * 1e6
            ),
            "magnet_offset_multiplier": MAGNET_OFFSET_MULTIPLIER,
            "magnet_roll_urad": (
                SIGMA_MAGNET_ROLL_RAD * MAGNET_ERROR_SCALE
                * MAGNET_ROLL_MULTIPLIER * 1e6
            ),
            "magnet_roll_multiplier": MAGNET_ROLL_MULTIPLIER,
            "bpm_offset_um": SIGMA_BPM_OFFSET_MM * 1e3,
            "bpm_roll_mrad": SIGMA_BPM_ROLL_RAD * 1e3,
        },
        "components_with_magnet_errors": counts,
        "published_alignment": {
            "enabled": bool(PUBLISHED_ALIGNMENT_M),
            "path": str(ALIGNMENT_MAP_PATH) if ALIGNMENT_MAP_PATH else None,
            "status": (
                "vector-digitised Kubo Fig. 1/2 points associated by longitudinal order"
                if PUBLISHED_ALIGNMENT_M else "not applied"
            ),
        },
        "preliminary_orbit_search": {
            "reference": "Kubo 2003, Sec. III C: all steerers; target |x| < 2 mm, |y| < 1 mm",
            "continuation_limits_mm": {
                "x": PRELIMINARY_CONTINUATION_X_MM,
                "y": PRELIMINARY_CONTINUATION_Y_MM,
            },
            "continuation_note": (
                "Tighter than the published final acceptance limits; used only "
                "to retain the nearby periodic-orbit branch while ramping the "
                "unchanged Table-I random errors on the 2011 lattice."
            ),
            "records": preliminary,
            "final_x_command_rms_mrad": _rms(x_commands) * 1e3,
            "final_y_command_rms_mrad": _rms(y_commands) * 1e3,
            "final_x_command_max_mrad": float(np.max(np.abs(x_commands)) * 1e3),
            "final_y_command_max_mrad": float(np.max(np.abs(y_commands)) * 1e3),
        },
        "cod_dispersion_solver": {
            "name": COD_DISPERSION_SOLVER,
            "max_steerer_excursion_mrad": SAD_MAX_STEERER_KICK_RAD * 1e3,
            "cod_selected_counts": cod_selected_counts,
            "dispersion_selected_counts": dispersion_selected_counts,
            "dispersion_requested_gains": [KUBO_GAIN_FIRST, 1.0],
            "dispersion_applied_gains": dispersion_applied_gains,
            "dispersion_line_search": "same nominal response; backtracking on periodic orbit and stacked COD/Dy residual",
            "source": "sad/operation/lib/cod.n and coddispersion.n",
        },
        "response_dimensions": {
            "x_cod": list(response["x"].shape),
            "y_cod_dispersion": list(response["y_orbit"].shape),
            "coupling": list(response["coupling"].shape),
        },
        "skew_corrector_family": SKEW_FAMILY,
        "response_settings": {
            "source": RESPONSE_SOURCE,
            "steerer_response_kick_mrad": STEERER_RESPONSE_KICK_RAD * 1e3,
            "local_orm_probe_mrad": (
                LOCAL_ORM_PROBE_RAD * 1e3
                if RESPONSE_SOURCE == "local_orm" else None
            ),
            "paper_steerer_response_kick_mrad": 0.1,
            "response_probe_note": (
                "Matches Kubo's 0.1-mrad nominal-response probe."
                if np.isclose(STEERER_RESPONSE_KICK_RAD, 1e-4)
                else "A non-paper probe is used because the selected lattice "
                "does not retain a periodic orbit at 0.1 mrad."
            ),
            "cod_svd_rcond": RCOND_COD,
            "coupling_svd_rcond": RCOND_COUPLING,
            "closed_orbit_tolerance_mm_mrad": CLOSED_ORBIT_TOLERANCE,
            "x_cod_svd_rank": _svd_rank(response["x"], RCOND_COD),
            "y_cod_svd_rank": _svd_rank(response["y_orbit"], RCOND_COD),
        },
        "local_orm_mismatch": local_orm_diagnostics,
        "stages": stages,
        "actuator_settings": actuator_settings,
        "scope": (
            f"Kubo-2003 correction sequence on {LATTICE_LABEL} for one seed at "
            f"{MAGNET_ERROR_SCALE:g} times the Table-I random magnet errors; "
            "validates COD, vertical-dispersion and coupling intermediate "
            "observables. The response SVD cutoff, selected-lattice-to-2003 "
            "actuator calibration, radiation-equilibrium model, and 500-seed "
            "Table-II statistic remain to be validated."
        ),
    }


def main():
    analysis = Path(__file__).resolve().parents[4] / "analysis" / "DR-RFTrack"
    probe_label = f"{STEERER_RESPONSE_KICK_RAD * 1e3:g}".replace(".", "p")
    rcond_label = f"{RCOND_COD:.0e}".replace("-", "m")
    cache = analysis / (
        f"kubo2003_nominal_kick_response_cache_{LATTICE_LABEL}_"
        f"{SKEW_FAMILY}_probe_{probe_label}mrad.npz"
    )
    response = _load_or_build_responses(cache)
    result = run_case(response)
    suffix = (
        f"_{LATTICE_LABEL}_{SKEW_FAMILY}_probe_{probe_label}mrad_"
        f"solver_{COD_DISPERSION_SOLVER}_response_{RESPONSE_SOURCE}_"
        f"offset_{MAGNET_OFFSET_MULTIPLIER:g}_roll_{MAGNET_ROLL_MULTIPLIER:g}_"
        f"rcond_{rcond_label}_seed_{SEED}"
    ) + (
        "" if MAGNET_ERROR_SCALE == 1.0
        else f"_scale_{MAGNET_ERROR_SCALE:g}".replace(".", "p")
    )
    result_path = analysis / f"kubo2003_rftrack_procedure_result{suffix}.json"
    result_path.write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    if "stages" not in result:
        print(json.dumps(result, indent=2))
        return

    labels = ("after COD", "after COD + Dy", "after coupling")
    stages = [result["stages"][key] for key in ("after_cod", "after_cod_dispersion", "after_coupling")]
    figure, axes = plt.subplots(1, 3, figsize=(11.2, 3.4), constrained_layout=True)
    plots = (
        ("true_cod_y_mm", "vertical COD RMS [mm]"),
        ("true_dy_all_bpm_mm", "vertical dispersion RMS [mm/delta]"),
        ("measured_cxy", "measured Cxy"),
    )
    for axis, (key, ylabel) in zip(axes, plots):
        values = np.array([stage[key] for stage in stages])
        bars = axis.bar(labels, values, color=("0.55", "tab:blue", "tab:green"))
        axis.bar_label(bars, labels=[f"{value:.2e}" for value in values], padding=3, fontsize=8)
        axis.set(title=key.replace("_", " "), ylabel=ylabel)
        axis.tick_params(axis="x", rotation=15, labelsize=8)
        axis.grid(True, axis="y", alpha=0.25)
    figure.suptitle(
        "RF-Track Kubo-2003-style sequential correction "
        f"({LATTICE_LABEL}, {MAGNET_ERROR_SCALE:g} Table-I magnet errors)",
        fontsize=14,
    )
    figure.savefig(
        analysis / f"kubo2003_rftrack_procedure_result{suffix}.png", dpi=180
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
