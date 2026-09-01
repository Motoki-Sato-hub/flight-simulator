"""Generate and load linear Twiss data for the SAD-derived ATF DR lattice.

The values are calculated from the RF-Track lattice itself, linearized about
its periodic closed orbit.  They are not imported from MAD-X or SAD.  The
generated TFS file contains one row per RF-Track lattice element and is intended as the
nominal optics reference for the DR correction-validation notebook.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import ATFDRRingCorrection
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import (
    NOMINAL_MOMENTUM_MEV_C,
    REFERENCE_SAD_RELATIVE_PATH,
    build_atf_dr_lattice,
)


TWISS_FILENAME = "ATF_DR_20111111b_RFTrack_twiss.tfs"
DEFAULT_FINITE_DIFFERENCE_STEP = 1e-5


def twiss_path() -> Path:
    """Return the checked-in RF-Track-generated DR Twiss file path."""
    return Path(__file__).with_name(TWISS_FILENAME)


def _track_coordinates(
    lattice,
    coordinates: np.ndarray,
    momentum_mev_c: float = NOMINAL_MOMENTUM_MEV_C,
) -> np.ndarray:
    import RF_Track as rft

    phase_space = np.concatenate(
        (np.asarray(coordinates, dtype=float), [0.0, float(momentum_mev_c)])
    )
    bunch = rft.Bunch6d(rft.electronmass, 0.0, -1.0, phase_space)
    tracked = lattice.track(bunch)
    if tracked.size() != 1:
        raise RuntimeError("RF-Track lost a finite-difference Twiss probe")
    return np.asarray(tracked.get_phase_space()[0, :4], dtype=float)


def _linear_map(lattice, closed_orbit: np.ndarray, step: float) -> np.ndarray:
    """Return the 4x4 map linearized about ``closed_orbit``.

    RF-Track coordinates use mm and mrad, so the usual Twiss matrix formulas
    retain their numerical metre/radian interpretation for each 2x2 plane.
    """
    matrix = np.empty((4, 4), dtype=float)
    for column in range(4):
        delta = np.zeros(4, dtype=float)
        delta[column] = step
        plus = _track_coordinates(lattice, closed_orbit + delta)
        minus = _track_coordinates(lattice, closed_orbit - delta)
        matrix[:, column] = (plus - minus) / (2.0 * step)
    return matrix


def _periodic_twiss(matrix: np.ndarray) -> tuple[float, float, float]:
    """Return beta, alpha and tune from a stable uncoupled 2x2 map."""
    cosine = float(np.clip(np.trace(matrix) / 2.0, -1.0, 1.0))
    sine_magnitude = float(np.sqrt(max(0.0, 1.0 - cosine**2)))
    sign_source = matrix[0, 1] if abs(matrix[0, 1]) > 1e-14 else -matrix[1, 0]
    sine = float(np.copysign(sine_magnitude, sign_source))
    if abs(sine) < 1e-12:
        raise RuntimeError("Cannot determine periodic Twiss parameters at integer/half-integer tune")
    beta = float(matrix[0, 1] / sine)
    alpha = float((matrix[0, 0] - matrix[1, 1]) / (2.0 * sine))
    if beta <= 0.0:
        raise RuntimeError("The linearized one-turn map is not stable")
    tune = float(np.arctan2(sine, cosine) / (2.0 * np.pi))
    if tune < 0.0:
        tune += 1.0
    return beta, alpha, tune


def _propagate_twiss(map_to_point: np.ndarray, beta0: float, alpha0: float) -> tuple[float, float]:
    gamma0 = (1.0 + alpha0**2) / beta0
    sigma0 = np.array(((beta0, -alpha0), (-alpha0, gamma0)), dtype=float)
    sigma = map_to_point @ sigma0 @ map_to_point.T
    return float(sigma[0, 0]), float(-sigma[0, 1])


def _keyword(element) -> str:
    """Map RF-Track element classes to familiar MAD-X TFS keywords."""
    element_type = type(element).__name__
    if element_type == "Bpm":
        return "MONITOR"
    if element_type == "Quadrupole":
        return "QUADRUPOLE"
    if element_type == "Sextupole":
        return "SEXTUPOLE"
    if element_type == "SBend":
        return "SBEND"
    if element_type == "Corrector":
        return "HKICKER" if element.get_name().startswith("ZH") else "VKICKER"
    if element_type == "Multipole":
        return "MULTIPOLE"
    return "DRIFT"


def _magnetic_columns(element) -> tuple[float, float, float, float, float, float, float, float, float, float, float]:
    """Return the TFS fields required to reconstruct RF-Track elements."""
    p_over_q = -NOMINAL_MOMENTUM_MEV_C
    k0l = k0sl = k1l = k1sl = k2l = k2sl = angle = e1 = e2 = 0.0
    element_type = type(element).__name__
    if element_type == "Quadrupole":
        k1l = float(element.get_K1L(p_over_q))
    elif element_type == "Sextupole":
        k2l = float(element.get_K2L(p_over_q))
    elif element_type == "SBend":
        angle = float(element.get_angle())
        k0l = float(element.get_K0())
        k1l = float(element.get_K1L())
        e1 = e2 = angle / 2.0
    elif element_type == "Multipole":
        strengths = np.asarray(element.get_KnL(p_over_q), dtype=complex).reshape(-1)
        if strengths.size > 0:
            k0l, k0sl = float(strengths[0].real), float(strengths[0].imag)
        if strengths.size > 1:
            k1l, k1sl = float(strengths[1].real), float(strengths[1].imag)
        if strengths.size > 2:
            k2l, k2sl = float(strengths[2].real), float(strengths[2].imag)
    return k0l, k0sl, k1l, k1sl, k2l, k2sl, angle, e1, e2, 0.0, 0.0


def generate_atf_dr_twiss(
    output_path: str | Path | None = None,
    *,
    finite_difference_step: float = DEFAULT_FINITE_DIFFERENCE_STEP,
) -> Path:
    """Generate an all-element TFS Twiss table from the nominal RF-Track model."""
    if finite_difference_step <= 0.0:
        raise ValueError("finite_difference_step must be positive")

    output = twiss_path() if output_path is None else Path(output_path)
    import RF_Track as rft

    lattice = build_atf_dr_lattice()
    correction = ATFDRRingCorrection(lattice)
    closed_orbit = correction.find_closed_orbit()
    relative_momentum_step = 1e-3
    orbit_plus = correction.find_closed_orbit(
        momentum_mev_c=NOMINAL_MOMENTUM_MEV_C * (1.0 + relative_momentum_step),
        initial_coordinates=closed_orbit.initial_coordinates,
    )
    orbit_minus = correction.find_closed_orbit(
        momentum_mev_c=NOMINAL_MOMENTUM_MEV_C * (1.0 - relative_momentum_step),
        initial_coordinates=closed_orbit.initial_coordinates,
    )

    one_turn_map = _linear_map(
        lattice, closed_orbit.initial_coordinates, finite_difference_step
    )
    coupling = max(
        float(np.max(np.abs(one_turn_map[:2, 2:]))),
        float(np.max(np.abs(one_turn_map[2:, :2]))),
    )
    if coupling > 1e-8:
        raise RuntimeError(
            f"The nominal one-turn map is coupled ({coupling:.3e}); uncoupled Twiss export is invalid"
        )
    beta_x0, alpha_x0, tune_x = _periodic_twiss(one_turn_map[:2, :2])
    beta_y0, alpha_y0, tune_y = _periodic_twiss(one_turn_map[2:, 2:])

    rows = []
    map_to_element = np.eye(4)
    reference = closed_orbit.initial_coordinates.copy()
    plus = orbit_plus.initial_coordinates.copy()
    minus = orbit_minus.initial_coordinates.copy()
    for element in lattice["*"]:
        # RF-Track exposes tracking on Lattice, not directly on Element.
        # Linearize one element at a time and compose its map.  This gives an
        # all-element table without repeatedly tracking every growing prefix.
        element_lattice = rft.Lattice()
        element_lattice.append(element)
        element_map = _linear_map(
            element_lattice, reference, finite_difference_step
        )
        map_to_element = element_map @ map_to_element
        reference = _track_coordinates(element_lattice, reference)
        plus = _track_coordinates(
            element_lattice,
            plus,
            NOMINAL_MOMENTUM_MEV_C * (1.0 + relative_momentum_step),
        )
        minus = _track_coordinates(
            element_lattice,
            minus,
            NOMINAL_MOMENTUM_MEV_C * (1.0 - relative_momentum_step),
        )
        name = element.get_name()
        beta_x, alpha_x = _propagate_twiss(
            map_to_element[:2, :2], beta_x0, alpha_x0
        )
        beta_y, alpha_y = _propagate_twiss(
            map_to_element[2:, 2:], beta_y0, alpha_y0
        )
        local_dispersion = (plus - minus) / (2.0 * relative_momentum_step)
        magnetic_columns = _magnetic_columns(element)
        rows.append(
            (
                name,
                _keyword(element),
                element.get_S("exit"),
                element.get_length(),
                beta_x,
                alpha_x,
                beta_y,
                alpha_y,
                *reference,
                *local_dispersion,
                *magnetic_columns,
            )
        )

    if len(rows) != len(lattice["*"]):
        raise RuntimeError("Generated Twiss table does not contain every RF-Track element")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as stream:
        stream.write('@ NAME             %05s "ATF_DR_20111111b_RFTRACK"\n')
        stream.write('@ TYPE             %05s "TWISS"\n')
        stream.write('@ SEQUENCE         %05s "RING0"\n')
        stream.write('@ PARTICLE         %08s "ELECTRON"\n')
        stream.write('@ MASS             %le       0.00051099895\n')
        stream.write('@ CHARGE           %le                  -1\n')
        stream.write('@ ENERGY           %le                 1.3\n')
        stream.write(f'@ REFERENCE_SAD    %30s "{REFERENCE_SAD_RELATIVE_PATH}"\n')
        stream.write('@ ORIGIN           %30s "RF-Track periodic-orbit linearization"\n')
        stream.write(f"@ PC               %le {NOMINAL_MOMENTUM_MEV_C / 1000.0:18.10g}\n")
        stream.write(f"@ Q1               %le {tune_x:18.12g}\n")
        stream.write(f"@ Q2               %le {tune_y:18.12g}\n")
        stream.write(f"@ FD_STEP          %le {finite_difference_step:18.12g}\n")
        stream.write(
            "* NAME               KEYWORD                             S                  L               BETX               ALFX"
            "               BETY               ALFY                  X                 PX                  Y                 PY"
            "                 DX                DPX                 DY                DPY"
            "                K0L               K0SL                K1L               K1SL                K2L               K2SL"
            "              ANGLE                 E1                 E2              HKICK              VKICK\n"
        )
        stream.write("$ " + " ".join(["%s", "%s"] + ["%le"] * 25) + "\n")
        for row in rows:
            (
                name, keyword, s, length, beta_x, alpha_x, beta_y, alpha_y,
                x, px, y, py, dx, dpx, dy, dpy,
                k0l, k0sl, k1l, k1sl, k2l, k2sl, angle, e1, e2, hkick, vkick,
            ) = row
            stream.write(f' "{name}"'.ljust(21))
            stream.write(f' "{keyword}"'.ljust(37))
            stream.write(
                f"{s:18.12g}{length:19.12g}{beta_x:19.12g}{alpha_x:19.12g}"
                f"{beta_y:19.12g}{alpha_y:19.12g}{x:19.12g}{px:19.12g}"
                f"{y:19.12g}{py:19.12g}{dx:19.12g}{dpx:19.12g}"
                f"{dy:19.12g}{dpy:19.12g}"
            )
            stream.write(
                f"{k0l:19.12g}{k0sl:19.12g}{k1l:19.12g}{k1sl:19.12g}"
                f"{k2l:19.12g}{k2sl:19.12g}{angle:19.12g}{e1:19.12g}{e2:19.12g}"
                f"{hkick:19.12g}{vkick:19.12g}\n"
            )
    return output


def load_atf_dr_twiss(path: str | Path | None = None) -> np.ndarray:
    """Load the generated all-element Twiss table as a structured array."""
    source = twiss_path() if path is None else Path(path)
    lines = source.read_text(encoding="utf-8").splitlines()
    header_index = next(
        (index for index, line in enumerate(lines) if line.startswith("* ")),
        None,
    )
    if header_index is None or header_index + 2 > len(lines):
        raise ValueError(f"Invalid TFS Twiss table: {source}")
    names = lines[header_index].split()[1:]
    data_lines = [line for line in lines[header_index + 2 :] if line.strip()]
    table = np.atleast_1d(np.genfromtxt(
        data_lines,
        names=names,
        dtype=None,
        encoding="utf-8",
    ))
    for field in ("NAME", "KEYWORD"):
        if field in table.dtype.names:
            table[field] = np.char.strip(table[field], '"')
    return table


if __name__ == "__main__":
    generated = generate_atf_dr_twiss()
    print(f"Generated {generated}")
