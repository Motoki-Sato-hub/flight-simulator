"""RF-Track lattice builder for the ATF damping ring.

Reference SAD daihon: ``atfdr-design-20111111b.sad``

The checked-in JSON is generated from the SAD element definitions and
``LINE RING0`` by ``generate_atf_dr_rftrack_lattice.py``.  SAD itself is not
required to build or track this lattice.

The default model is deliberately transverse-only: synchrotron radiation and
the thin SAD RF cavity are disabled.  ``rf_mode='equilibrium'`` is an explicit
opt-in pilot which adds incoherent synchrotron radiation in every SBend and a
short pillbox representation of SAD's 714-MHz CAV.  It is for turn-by-turn
equilibrium-emittance studies, not for the closed-orbit response routines.
"""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
from typing import Any

import numpy as np


REFERENCE_SAD_DAIHON = "atfdr-design-20111111b.sad"
REFERENCE_SAD_RELATIVE_PATH = "operation/daihon/atfdr-design-20111111b.sad"
LATTICE_DATA_FILENAME = "ATF_DR_RFTrack_lattice.json"
NOMINAL_MOMENTUM_MEV_C = 1299.9999
EQUILIBRIUM_RF_CAVITY_LENGTH_M = 1e-3
EQUILIBRIUM_RF_PHASE_DEG = 180.0


def load_lattice_data(
    lattice_data_path: str | Path | None = None,
) -> dict[str, Any]:
    """Load nominal data, or an explicitly selected compatible SAD export."""
    path = (
        Path(lattice_data_path)
        if lattice_data_path is not None
        else Path(__file__).with_name(LATTICE_DATA_FILENAME)
    )
    with path.open(encoding="utf-8") as stream:
        data = json.load(stream)
    metadata = data.get("metadata", {})
    if metadata.get("ring_line") != "RING0":
        raise ValueError(
            f"Expected a RING0 SAD export in {path}, got {metadata!r}"
        )
    return data


def get_lattice_metadata(lattice_data_path: str | Path | None = None) -> dict[str, Any]:
    """Return a copy of the source and geometry metadata."""
    return dict(load_lattice_data(lattice_data_path)["metadata"])


def get_bpm_nearest_magnet_names(
    lattice_data_path: str | Path | None = None,
) -> dict[str, str]:
    """Map each BPM to its nearest quadrupole or sextupole centre.

    Kubo's BPM-error model defines the BPM offset with respect to the field
    centre of the nearest quadrupole or sextupole.  The mapping is derived
    from the generated SAD sequence, not from a hard-coded BPM numbering.
    """
    data = load_lattice_data(lattice_data_path)
    definitions = data["definitions"]
    circumference = float(data["metadata"]["circumference_m"])
    occurrences: Counter[str] = Counter()
    bpm_index = 0
    drift_index = 0
    s = 0.0
    bpm_positions: list[tuple[str, float]] = []
    magnet_positions: list[tuple[str, float]] = []

    for source_name in data["sequence"]:
        definition = definitions[source_name]
        element_type = str(definition["type"])
        attributes = definition["attributes"]
        length = float(attributes.get("L", 0.0))
        occurrences[source_name] += 1
        if element_type == "DRIFT" or (
            element_type == "BEND"
            and (abs(float(attributes.get("ANGLE", 0.0))) == 0.0)
        ):
            drift_index += 1
        if source_name == "M":
            bpm_index += 1
        instance_name = _instance_name(
            source_name, element_type, occurrences, bpm_index, drift_index
        )
        centre = s + 0.5 * length
        if source_name == "M":
            bpm_positions.append((instance_name, centre))
        if element_type in {"QUAD", "SEXT"} and length > 0.0:
            magnet_positions.append((instance_name, centre))
        s += length

    if not bpm_positions or not magnet_positions:
        raise ValueError("The SAD export did not contain BPMs and quadrupole/sextupole magnets")
    return {
        bpm_name: min(
            magnet_positions,
            key=lambda item: min(
                abs(position - item[1]), circumference - abs(position - item[1])
            ),
        )[0]
        for bpm_name, position in bpm_positions
    }


def get_magnet_centres(
    lattice_data_path: str | Path | None = None,
) -> list[dict[str, Any]]:
    """Return SAD-order centres of the 204 main ring magnets.

    This deliberately excludes correctors and monitors.  It is used to attach
    Kubo's published Fig. 1/2 alignment points to the corresponding main
    quadrupole, sextupole, and finite-angle bend in a compatible SAD export.
    """
    data = load_lattice_data(lattice_data_path)
    definitions = data["definitions"]
    occurrences: Counter[str] = Counter()
    bpm_index = 0
    drift_index = 0
    s = 0.0
    centres: list[dict[str, Any]] = []

    for source_name in data["sequence"]:
        definition = definitions[source_name]
        element_type = str(definition["type"])
        attributes = definition["attributes"]
        length = float(attributes.get("L", 0.0))
        angle = float(attributes.get("ANGLE", 0.0))
        occurrences[source_name] += 1
        if element_type == "DRIFT" or (element_type == "BEND" and angle == 0.0):
            drift_index += 1
        if source_name == "M":
            bpm_index += 1
        instance_name = _instance_name(
            source_name, element_type, occurrences, bpm_index, drift_index
        )
        if (
            (element_type in {"QUAD", "SEXT"} and length > 0.0)
            or (element_type == "BEND" and abs(angle) > 0.0)
        ):
            centres.append(
                {
                    "name": instance_name,
                    "type": element_type,
                    "s_m": s + 0.5 * length,
                }
            )
        s += length

    return centres


# Historical ATF design values quoted as *half apertures* in the NLC/ATF
# damping-ring design report, Chapter 3.3.3:
# https://cds.cern.ch/record/450524/files/slac-r-559.pdf
#
# * 12 mm standard aperture in the arcs;
# *  6 mm mask in the south straight;
# *  5 mm photon masks in wiggler sections;
# *  5 mm extraction-kicker vacuum chamber.
#
# The 2011 SAD daihon does not encode a chamber/survey table, so these values
# must not be expanded into a complete loss model by assumption.  Only the
# extraction-kicker restriction can be located unambiguously in this export:
# its line sequence is ``IEX KIX KIX`` and yields KIX.1/KIX.2 in RF-Track.
HISTORICAL_ATF_DR_DESIGN_HALF_APERTURES_MM = {
    "arc_standard": 12.0,
    "south_straight_mask": 6.0,
    "wiggler_photon_mask": 5.0,
    "extraction_kicker_chamber": 5.0,
}
HISTORICAL_ATF_DR_APERTURE_SOURCE = (
    "NLC/ATF damping-ring design report, Chapter 3.3.3, "
    "https://cds.cern.ch/record/450524/files/slac-r-559.pdf"
)


def get_historical_extraction_kicker_apertures(
    lattice_data_path: str | Path | None = None,
    *,
    shape: str = "circular",
) -> dict[str, tuple[float, float, str]]:
    """Return the safely located historical 5-mm KIX half-aperture screen.

    This intentionally returns only ``KIX`` instances.  It is a partial
    historical-design loss screen, not an ATF transmission aperture table:
    the source values for arcs, wiggler photon masks, and the south-straight
    mask cannot be mapped to all corresponding 2011 SAD elements without a
    survey/engineering chamber map.  The source quotes a half aperture but
    does not establish the cross-section geometry; ``shape='circular'`` is a
    stated RF-Track study assumption, not an engineering claim.  Replace it
    when chamber drawings/survey data are available.
    """
    data = load_lattice_data(lattice_data_path)
    definitions = data["definitions"]
    occurrences: Counter[str] = Counter()
    bpm_index = 0
    drift_index = 0
    apertures: dict[str, tuple[float, float]] = {}
    for source_name in data["sequence"]:
        definition = definitions[source_name]
        element_type = str(definition["type"])
        attributes = definition["attributes"]
        angle = float(attributes.get("ANGLE", 0.0))
        occurrences[source_name] += 1
        if element_type == "DRIFT" or (element_type == "BEND" and angle == 0.0):
            drift_index += 1
        if source_name == "M":
            bpm_index += 1
        instance_name = _instance_name(
            source_name, element_type, occurrences, bpm_index, drift_index
        )
        if source_name == "KIX":
            half_aperture = HISTORICAL_ATF_DR_DESIGN_HALF_APERTURES_MM[
                "extraction_kicker_chamber"
            ]
            apertures[instance_name] = (half_aperture, half_aperture, shape)
    if not apertures:
        raise ValueError("The selected SAD export has no KIX extraction-kicker element")
    return apertures


def _instance_name(
    source_name: str,
    element_type: str,
    occurrences: Counter[str],
    bpm_index: int,
    drift_index: int,
) -> str:
    if source_name == "M":
        return f"MB{bpm_index}R"
    if element_type == "DRIFT":
        return f"DRIFT_{drift_index}"
    if source_name.upper().startswith(("ZH", "ZV")):
        return source_name.upper()
    if source_name == "IDUM":
        return f"IDUMB{occurrences[source_name]}R"
    return f"{source_name.upper()}.{occurrences[source_name]}"


def build_atf_dr_lattice(
    momentum_mev_c: float = NOMINAL_MOMENTUM_MEV_C,
    charge: float = -1.0,
    *,
    rf_mode: str = "disabled",
    radiation_quantum: bool = True,
    radiation_steps: int = 10,
    rf_phase_deg: float | None = None,
    rf_voltage_scale: float = 1.0,
    lattice_data_path: str | Path | None = None,
):
    """Build the ATF DR ``RING0`` lattice using RF-Track elements.

    Parameters
    ----------
    momentum_mev_c:
        Reference momentum used to convert normalized multipole strengths.
    charge:
        Particle charge in elementary-charge units.  The default is electron.
    rf_mode:
        ``"disabled"`` keeps the transverse response model.  ``"equilibrium"``
        adds radiation with quantum excitation to each SBend and represents
        the zero-length SAD CAV as a 1-mm 714-MHz pillbox.  Its phase must be
        synchronized by the emittance workflow before multi-turn tracking.
    radiation_quantum:
        Enable quantum excitation when ``rf_mode='equilibrium'``.  Set false
        only for deterministic RF-phase synchronization and damping checks.
    radiation_steps:
        Collective-effect integration steps per SBend in equilibrium mode.
    rf_phase_deg:
        Pillbox phase used only by ``rf_mode='equilibrium'``.  ``None`` uses
        the 180-degree 2011-reference setting; historical optics may require
        an explicit phase scan before radiation equilibrium is evaluated.
    rf_voltage_scale:
        Multiplicative diagnostic factor for the SAD cavity voltage.  The
        default is one; any other value is an explicit RF-model scan.
    lattice_data_path:
        Optional generated SAD export.  The default is the checked-in 2011
        lattice; an explicit path is intended for historical comparisons and
        never changes that default.
    """
    if rf_mode not in {"disabled", "equilibrium"}:
        raise ValueError("rf_mode must be 'disabled' or 'equilibrium'")
    if radiation_steps < 1:
        raise ValueError("radiation_steps must be positive")
    if charge == 0:
        raise ValueError("charge must be non-zero")
    if rf_voltage_scale <= 0.0:
        raise ValueError("rf_voltage_scale must be positive")

    import RF_Track as rft

    data = load_lattice_data(lattice_data_path)
    definitions = data["definitions"]
    p_over_q = float(momentum_mev_c) / float(charge)
    lattice = rft.Lattice()

    start = rft.Drift(0.0)
    start.set_name("RING0$START")
    lattice.append(start)

    occurrences: Counter[str] = Counter()
    bpm_index = 0
    drift_index = 0

    for source_name in data["sequence"]:
        definition = definitions[source_name]
        element_type = definition["type"]
        attributes = definition["attributes"]
        length = float(attributes.get("L", 0.0))
        occurrences[source_name] += 1

        if element_type == "DRIFT":
            element = rft.Drift(length)
            drift_index += 1
        elif element_type == "BEND":
            angle = float(attributes.get("ANGLE", 0.0))
            if source_name.upper().startswith(("ZH", "ZV")):
                element = rft.Corrector(length)
            elif abs(angle) > 0.0:
                # SAD BEND is represented as a rectangular bend.  RF-Track's
                # face rotations use +ANGLE/2 for the equivalent map.
                element = rft.SBend(
                    length,
                    angle,
                    p_over_q,
                    angle / 2.0,
                    angle / 2.0,
                )
                element.set_K1(float(attributes.get("K1", 0.0)))
                if rf_mode == "equilibrium":
                    element.set_cfx_nsteps(int(radiation_steps))
                    element.add_collective_effect(
                        rft.IncoherentSynchrotronRadiation(
                            quantum=bool(radiation_quantum)
                        )
                    )
            else:
                # Zero-angle injection/extraction bends and BHE are passive in
                # this model; operational DR steerers are the ZH/ZV elements.
                element = rft.Drift(length)
                drift_index += 1
        elif element_type == "QUAD":
            # SAD's K1 in this daihon is the integrated strength K1*L.
            # RF-Track's Quadrupole constructor expects K1 [m^-2].
            k1 = float(attributes.get("K1", 0.0)) / length if length else 0.0
            element = rft.Quadrupole(
                length,
                p_over_q,
                k1,
            )
        elif element_type == "SEXT":
            # Likewise SAD K2 is integrated K2*L, while RF-Track expects
            # normalized K2 [m^-3].
            k2 = float(attributes.get("K2", 0.0)) / length if length else 0.0
            if length:
                element = rft.Sextupole(length, p_over_q, k2)
            else:
                # Preserve SAD's zero-length sextupole-compensation markers
                # without producing NaN from Sextupole(K2L/L) at L=0.
                element = rft.Multipole()
                element.set_strengths(np.zeros(3, dtype=complex))
        elif element_type == "MONI":
            element = rft.Bpm(length)
            if source_name == "M":
                bpm_index += 1
        elif element_type == "MARK":
            element = rft.Drift(length)
            drift_index += 1
        elif element_type == "CAVI":
            if rf_mode == "disabled":
                element = rft.Drift(0.0)
                drift_index += 1
            else:
                voltage_v = float(attributes["VOLT"]) * float(rf_voltage_scale)
                frequency_hz = float(attributes["FREQ"])
                # SAD's CAVI is thin.  A short pillbox preserves the location
                # while making a longitudinal RF kick available to RF-Track.
                element = rft.Pillbox_Cavity(
                    np.array([[voltage_v / EQUILIBRIUM_RF_CAVITY_LENGTH_M]]),
                    frequency_hz,
                    EQUILIBRIUM_RF_CAVITY_LENGTH_M,
                    1,
                )
                element.set_t0(0.0)
                element.set_phid(
                    EQUILIBRIUM_RF_PHASE_DEG
                    if rf_phase_deg is None else float(rf_phase_deg)
                )
        else:  # pragma: no cover - protected by the generator
            raise ValueError(f"Unsupported SAD element type: {element_type}")

        instance_name = _instance_name(
            source_name,
            element_type,
            occurrences,
            bpm_index,
            drift_index,
        )
        element.set_name(instance_name)
        lattice.append(element)

        if element_type == "QUAD":
            # A physical quadrupole roll is represented explicitly as a thin
            # skew-K1L companion.  RF-Track Element.set_offsets changes the
            # placement frame but does not provide the field-roll error needed
            # for ATF coupling/emittance studies.  The companion is zero in
            # the nominal lattice and is driven by set_quadrupole_roll_error.
            roll = rft.Multipole(0.0)
            roll.set_KnL(p_over_q, np.zeros(2, dtype=complex))
            roll.set_name(f"{instance_name}$ROLL")
            lattice.append(roll)

        if element_type == "SEXT" and length:
            # A sextupole roll produces a skew-sextupole component.  It is
            # kept separate from the normal Sextupole so that a Table-I
            # magnet-roll error can be changed without modifying the SAD
            # reference field.
            roll = rft.Multipole(0.0)
            roll.set_KnL(p_over_q, np.zeros(3, dtype=complex))
            roll.set_name(f"{instance_name}$ROLL")
            lattice.append(roll)

        if source_name in {"SD1R", "SF1R"}:
            # SAD's skew-coupling correction drives the skew-quadrupole
            # component of these sextupole locations (skewcor.n).  Keep the
            # physical normal sextupole as above and add a thin, initially
            # zero skew K1L component.  This preserves the daihon geometry
            # while making the same 34 SD1R and 34 SF1R actuators available
            # to the RF-Track correction model.
            skew = rft.Multipole(0.0)
            skew.set_KnL(p_over_q, np.zeros(2, dtype=complex))
            skew.set_name(f"{instance_name}$SKEW")
            lattice.append(skew)

    end = rft.Drift(0.0)
    end.set_name("RING0$END")
    lattice.append(end)
    return lattice


def set_quadrupole_roll_error(
    lattice,
    quadrupole_name: str,
    roll_rad: float,
    *,
    momentum_mev_c: float = NOMINAL_MOMENTUM_MEV_C,
    charge: float = -1.0,
) -> None:
    """Set the thin skew-K1L equivalent of a rolled quadrupole.

    To first order, a normal quadrupole ``K1L`` rolled by ``theta`` contains a
    skew component ``2 theta K1L``.  The normal quadrupole and zero-length
    companion are adjacent, so this is the appropriate linear error model for
    the small (300-urad) Kubo Table-I rotations.
    """
    quadrupole = lattice[quadrupole_name]
    if isinstance(quadrupole, list):
        if len(quadrupole) != 1:
            raise ValueError(f"Expected one quadrupole named {quadrupole_name}")
        quadrupole = quadrupole[0]
    companion_name = f"{quadrupole_name}$ROLL"
    companion = lattice[companion_name]
    if isinstance(companion, list):
        if len(companion) != 1:
            raise ValueError(f"Expected one roll companion named {companion_name}")
        companion = companion[0]
    p_over_q = float(momentum_mev_c) / float(charge)
    strengths = np.asarray(companion.get_KnL(p_over_q), dtype=complex).copy()
    strengths.reshape(-1)[1] = 2j * float(roll_rad) * quadrupole.get_K1L(p_over_q)
    companion.set_KnL(p_over_q, strengths)


def set_sextupole_roll_error(
    lattice,
    sextupole_name: str,
    roll_rad: float,
    *,
    momentum_mev_c: float = NOMINAL_MOMENTUM_MEV_C,
    charge: float = -1.0,
) -> None:
    """Set the thin skew-K2L equivalent of a rolled sextupole.

    To first order, rolling a normal sextupole by ``theta`` produces the
    skew component ``3 theta K2L``.  RF-Track placement-frame rotations do
    not rotate the multipole field itself, so a zero-length companion is used
    just as for quadrupoles.
    """
    sextupole = lattice[sextupole_name]
    if isinstance(sextupole, list):
        if len(sextupole) != 1:
            raise ValueError(f"Expected one sextupole named {sextupole_name}")
        sextupole = sextupole[0]
    companion_name = f"{sextupole_name}$ROLL"
    companion = lattice[companion_name]
    if isinstance(companion, list):
        if len(companion) != 1:
            raise ValueError(f"Expected one roll companion named {companion_name}")
        companion = companion[0]
    p_over_q = float(momentum_mev_c) / float(charge)
    strengths = np.asarray(companion.get_KnL(p_over_q), dtype=complex).copy()
    strengths.reshape(-1)[2] = 3j * float(roll_rad) * sextupole.get_K2L(p_over_q)
    companion.set_KnL(p_over_q, strengths)


__all__ = [
    "REFERENCE_SAD_DAIHON",
    "REFERENCE_SAD_RELATIVE_PATH",
    "NOMINAL_MOMENTUM_MEV_C",
    "build_atf_dr_lattice",
    "set_quadrupole_roll_error",
    "set_sextupole_roll_error",
    "get_lattice_metadata",
    "get_bpm_nearest_magnet_names",
    "get_magnet_centres",
    "load_lattice_data",
]
