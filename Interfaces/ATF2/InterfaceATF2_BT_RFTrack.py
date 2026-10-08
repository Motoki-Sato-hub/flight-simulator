"""ATF Linac-to-DR transport (BT) RF-Track interface from the SAD daihon.

The historical BT optics are stored at 1.542282 GeV/c, whereas present ATF
operation uses 1.3 GeV/c.  The lattice geometry and normalized strengths are
read from the SAD-derived TFS, then all main quadrupoles and sector bends are
retuned at the requested reference momentum.  This preserves the design
linear optics and bend angles under the usual magnetic-rigidity scaling.

The converted TFS is presently maintained in the sibling
``ATF2_LinacBT_RFTrack`` project.  Its path can be overridden with
``ATF2_BT_RFT_TFS``.  Keeping that source explicit is intentional: the
historical SAD-to-RF-Track optics equivalence still needs further validation
before the generated TFS is vendored as a frozen Flight Simulator asset.
"""

from __future__ import annotations

import os
import shlex
from pathlib import Path

import numpy as np
import RF_Track as rft

from Interfaces.AbstractMachineInterface import AbstractMachineInterface


BT_SOURCE_MOMENTUM_MEV_C = 1542.282
BT_DEFAULT_MOMENTUM_MEV_C = 1300.0
_BPM_ALIASES = {"ML10T": "MB10T", "ML11T": "MB11T"}


def _default_tfs_path() -> Path:
    override = os.environ.get("ATF2_BT_RFT_TFS")
    if override:
        return Path(override).expanduser()
    projects_root = Path(__file__).resolve().parents[3]
    return projects_root / "ATF2_LinacBT_RFTrack" / "generated" / "atf2_bt_1p542282GeV_RFTrack.twiss"


def _read_tfs_rows(path: Path) -> list[dict[str, float | str]]:
    """Read just the TFS columns required to rigidity-scale BT magnets."""
    columns: list[str] = []
    types: list[str] = []
    rows: list[dict[str, float | str]] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("*"):
            columns = shlex.split(line)[1:]
        elif line.startswith("$"):
            types = shlex.split(line)[1:]
        elif columns and types and not line.startswith("@"):
            values = shlex.split(line)
            if len(values) != len(columns):
                raise ValueError(f"{path}: malformed TFS row")
            row: dict[str, float | str] = {}
            for name, type_code, value in zip(columns, types, values):
                row[name] = value if type_code.endswith("s") else float(value)
            rows.append(row)
    if not rows:
        raise ValueError(f"{path}: no TFS table rows found")
    return rows


class InterfaceATF2_BT_RFTrack(AbstractMachineInterface):
    """Corrector-capable 1.3 GeV/c BT model built from the SAD lattice."""

    def get_name(self):
        return "ATF2_BT_RFT"

    def __init__(
        self,
        population: float = 1e10,
        jitter: float = 0.0,
        bpm_resolution: float = 0.0,
        nsamples: int = 1,
        nparticles: int = 1000,
        momentum_mev_c: float = BT_DEFAULT_MOMENTUM_MEV_C,
        tfs_path: str | os.PathLike[str] | None = None,
        reference_particle: bool = False,
        apertures_mm: dict[str, tuple[float, float] | tuple[float, float, str]] | None = None,
    ):
        super().__init__()
        self.log = print
        self.population = float(population)
        self.jitter = float(jitter)
        self.nsamples = int(nsamples)
        self.nparticles = int(nparticles)
        # Use this for linear optics/ORM identification only.  Finite bunches
        # are needed separately for transmission and aperture studies.
        self.reference_particle = bool(reference_particle)
        self.Pref0 = float(momentum_mev_c)
        self.Pref = self.Pref0
        self.Q = -1.0
        self.electronmass = rft.electronmass
        self.dfs_test_energy = 0.995
        self.wfs_test_charge = 0.90
        self._beam_mode = "nominal"
        self.tfs_path = Path(tfs_path) if tfs_path is not None else _default_tfs_path()
        if not self.tfs_path.is_file():
            raise FileNotFoundError(
                "BT SAD-derived TFS not found: "
                f"{self.tfs_path}. Set ATF2_BT_RFT_TFS to its location."
            )

        self._tfs_rows = _read_tfs_rows(self.tfs_path)
        self.lattice = rft.Lattice(str(self.tfs_path))
        self._retune_main_magnets_for_momentum(self.Pref0)
        self._configured_apertures = self._apply_apertures(apertures_mm or {})
        self.lattice.set_bpm_resolution(float(bpm_resolution))
        self._model_bpms = [element.get_name() for element in self.lattice.get_bpms()]
        self.bpms = [_BPM_ALIASES.get(name, name) for name in self._model_bpms]
        self._interface_to_model_bpm = {
            interface: model for interface, model in zip(self.bpms, self._model_bpms)
        }
        self.corrs = [element.get_name() for element in self.lattice.get_correctors()]
        self.quadrupoles = list(dict.fromkeys(element.get_name() for element in self.lattice.get_quadrupoles()))
        self.sequence = [
            _BPM_ALIASES.get(element.get_name(), element.get_name())
            for element in self.lattice["*"]
        ]
        self._setup_beam(population=self.population, momentum_mev_c=self.Pref0)
        self._track_bunch()
        self.machine_name = "ATF2"
        self.lattice.align_elements()

    def _retune_main_magnets_for_momentum(self, momentum_mev_c: float) -> None:
        """Set physical fields for a new reference rigidity.

        The source TFS stores SAD normalized K1L at 1.542282 GeV/c.  Passing
        the same normalized K1 with a lower design momentum to RF-Track scales
        the physical quadrupole gradient by p/p_source; likewise changing the
        bend P/Q preserves its nominal bend angle.
        """
        grouped: dict[tuple[str, str], list[dict[str, float | str]]] = {}
        for row in self._tfs_rows:
            keyword = str(row["KEYWORD"])
            if keyword in {"QUADRUPOLE", "SBEND"}:
                grouped.setdefault((keyword, str(row["NAME"])), []).append(row)
        p_over_q = float(momentum_mev_c) / self.Q
        for (keyword, name), rows in grouped.items():
            elements = self.lattice.get_elements_by_name(name)
            if not isinstance(elements, list):
                elements = [elements]
            if len(elements) != len(rows):
                raise ValueError(f"{name}: TFS/RF-Track occurrence count mismatch")
            for element, row in zip(elements, rows):
                if keyword == "QUADRUPOLE":
                    length = float(row["L"])
                    k1 = float(row["K1L"]) / length if length else 0.0
                    element.set_K1(-k1, momentum_mev_c)
                else:
                    element.set_P_over_Q(p_over_q)

    def _apply_apertures(
        self,
        apertures_mm: dict[str, tuple[float, float] | tuple[float, float, str]],
    ) -> dict[str, tuple[float, float, str]]:
        """Apply supplied half-apertures in RF-Track millimetres.

        The historical BT TFS contains no aperture information.  This method
        intentionally accepts only explicit caller-provided values, normally
        from a survey or an approved engineering aperture table.  A two-value
        tuple is circular; a third value selects an RF-Track shape such as
        ``"rectangular"`` or ``"elliptical"``.
        """
        configured: dict[str, tuple[float, float, str]] = {}
        for name, specification in apertures_mm.items():
            if len(specification) not in (2, 3):
                raise ValueError(f"{name}: aperture must be (x_mm, y_mm[, shape])")
            x_mm, y_mm = float(specification[0]), float(specification[1])
            shape = str(specification[2]) if len(specification) == 3 else "circular"
            if not np.isfinite(x_mm) or not np.isfinite(y_mm) or x_mm <= 0.0 or y_mm <= 0.0:
                raise ValueError(f"{name}: aperture dimensions must be positive finite mm values")
            elements = self.lattice.get_elements_by_name(name)
            if not isinstance(elements, list):
                elements = [elements]
            if not elements:
                raise ValueError(f"{name}: no BT lattice element exists for aperture assignment")
            for element in elements:
                if len(specification) == 3:
                    element.set_aperture(x_mm, y_mm, shape)
                else:
                    element.set_aperture(x_mm, y_mm)
            configured[name] = (x_mm, y_mm, shape)
        return configured

    def get_model_coverage(self) -> dict[str, object]:
        """State which prerequisites for physical transmission are present."""
        return {
            "aperture_model": bool(self._configured_apertures),
            "configured_apertures_mm": dict(self._configured_apertures),
            "transmission_interpretation": (
                "aperture-aware survival at explicitly configured elements; "
                "not necessarily a complete BT transmission model" if self._configured_apertures
                else "geometric tracking survival only; no BT aperture data loaded"
            ),
        }

    def _setup_beam(self, *, population: float, momentum_mev_c: float) -> None:
        if self.reference_particle:
            self.B0 = rft.Bunch6d(
                self.electronmass,
                population,
                self.Q,
                np.array([[0.0, 0.0, 0.0, 0.0, 0.0, momentum_mev_c]]),
            )
            return
        twiss = rft.Bunch6d_twiss()
        # A conservative ATF DR-like injected bunch.  Actual injection Twiss
        # and charge should replace these in digital-twin calibration.
        twiss.emitt_x = 5.2
        twiss.emitt_y = 0.03
        twiss.beta_x = 4.0
        twiss.beta_y = 4.0
        twiss.alpha_x = 0.0
        twiss.alpha_y = 0.0
        twiss.sigma_t = 8.0
        twiss.sigma_pt = 0.8
        self.B0 = rft.Bunch6d_QR(
            self.electronmass,
            population,
            self.Q,
            momentum_mev_c,
            twiss,
            self.nparticles,
        )

    def _track_bunch(self) -> None:
        if self.reference_particle:
            self.B1 = self.lattice.track(self.B0)
            return
        info = self.B0.get_info()
        bunch = self.B0.displaced(
            self.jitter * info.sigma_x,
            self.jitter * info.sigma_y,
            0.0,
            0.0,
            0.0,
            self.jitter * info.sigma_py,
            self.jitter * info.sigma_px,
        )
        self.B1 = self.lattice.track(bunch)

    def get_sequence(self):
        return list(self.sequence)

    def get_hcorrectors_names(self):
        return [name for name in self.corrs if name.startswith(("ZH", "ZX"))]

    def get_vcorrectors_names(self):
        return [name for name in self.corrs if name.startswith(("ZV", "ZY"))]

    def get_bpms(self, names=None):
        x = np.zeros((self.nsamples, len(self.bpms)))
        y = np.zeros_like(x)
        tmit = np.zeros_like(x)
        for sample in range(self.nsamples):
            for index, model_name in enumerate(self._model_bpms):
                reading = self.lattice[model_name].get_reading()
                x[sample, index] = reading[0]
                y[sample, index] = reading[1]
                tmit[sample, index] = self.lattice[model_name].get_total_charge()
        result = {"names": list(self.bpms), "x": x, "y": y, "tmit": tmit}
        return self._select_bpms(result, names)

    def _select_bpms(self, bpms, names):
        if isinstance(names, str):
            names = [names]
        if names is None:
            return bpms
        indices = [index for index, name in enumerate(bpms["names"]) if name in names]
        return {
            "names": [bpms["names"][index] for index in indices],
            "x": bpms["x"][:, indices],
            "y": bpms["y"][:, indices],
            "tmit": bpms["tmit"][:, indices],
        }

    def get_icts(self, names=None):
        bpms = self.get_bpms(names)
        return {"names": bpms["names"], "charge": bpms["tmit"].mean(axis=0)}

    def get_transmission(self) -> float:
        """Model survival fraction at the BT exit.

        The current SAD-derived TFS has no physical aperture/chamber data, so
        this is a tracking-survival diagnostic, *not* a prediction of measured
        ATF transmission until surveyed apertures and loss monitors are added.
        """
        input_charge = abs(float(self.B0.get_total_charge()))
        if input_charge == 0.0:
            return 0.0
        return abs(float(self.B1.get_total_charge())) / input_charge

    def get_correctors(self, names=None):
        values = np.zeros(len(self.corrs), dtype=float)
        for index, name in enumerate(self.corrs):
            strength = self.lattice[name].get_strength()
            values[index] = 10.0 * (strength[0] if name.startswith(("ZH", "ZX")) else strength[1])
        result = {"names": list(self.corrs), "bdes": values, "bact": values.copy()}
        return self._select_magnets(result, names)

    def set_correctors(self, names, values):
        if isinstance(names, str):
            names = [names]
        if np.isscalar(values):
            values = [values]
        if len(names) != len(values):
            raise ValueError("names and values must have the same length")
        for name, value in zip(names, values):
            if name not in self.corrs:
                raise ValueError(f"Unknown BT corrector: {name}")
            if name.startswith(("ZH", "ZX")):
                self.lattice[name].set_strength(float(value) / 10.0, 0.0)
            else:
                self.lattice[name].set_strength(0.0, float(value) / 10.0)
        self._track_bunch()

    def vary_correctors(self, names, values):
        if isinstance(names, str):
            names = [names]
        if np.isscalar(values):
            values = [values]
        current = self.get_correctors(names)["bdes"]
        self.set_correctors(names, current + np.asarray(values, dtype=float))

    def get_quadrupoles(self, names=None):
        values = np.zeros(len(self.quadrupoles), dtype=float)
        for index, name in enumerate(self.quadrupoles):
            elements = self.lattice.get_elements_by_name(name)
            if not isinstance(elements, list):
                elements = [elements]
            strength = elements[0].get_K1(self.Pref0 / self.Q)
            values[index] = float(strength[0] if isinstance(strength, (list, tuple, np.ndarray)) else strength)
        result = {"names": list(self.quadrupoles), "bdes": values, "bact": values.copy()}
        return self._select_magnets(result, names)

    def set_quadrupoles(self, names, values):
        if isinstance(names, str):
            names = [names]
        if np.isscalar(values):
            values = [values]
        if len(names) != len(values):
            raise ValueError("names and values must have the same length")
        for name, value in zip(names, values):
            elements = self.lattice.get_elements_by_name(name)
            if not isinstance(elements, list):
                elements = [elements]
            for element in elements:
                element.set_K1(self.Pref0 / self.Q, float(value))
        self._track_bunch()

    @staticmethod
    def _select_magnets(values, names):
        if isinstance(names, str):
            names = [names]
        if names is None:
            return values
        indices = [index for index, name in enumerate(values["names"]) if name in names]
        return {
            "names": [values["names"][index] for index in indices],
            "bdes": values["bdes"][indices],
            "bact": values["bact"][indices],
        }

    def _reference_orbit_at_momentum(self, momentum_mev_c: float) -> tuple[np.ndarray, np.ndarray]:
        """Track one probe particle and return its BPM orbit in millimetres."""
        probe = rft.Bunch6d(
            self.electronmass,
            1.0,
            self.Q,
            np.array([[0.0, 0.0, 0.0, 0.0, 0.0, momentum_mev_c]]),
        )
        tracked = self.lattice.track(probe)
        if tracked.size() != 1:
            raise RuntimeError("BT reference probe was lost while measuring dispersion")
        bpms = self.get_bpms()
        return np.asarray(bpms["x"][0], dtype=float), np.asarray(bpms["y"][0], dtype=float)

    def get_model_dispersion(
        self,
        names=None,
        relative_momentum_step: float = 1.0e-3,
    ) -> dict[str, np.ndarray | list[str] | float]:
        """Return model BPM dispersion in mm per relative momentum offset.

        Magnets remain at the nominal BT rigidity.  The method tracks symmetric
        reference probes at ``p0 * (1 +/- delta)`` and restores the interface's
        bunch/BPM state before returning.
        """
        if relative_momentum_step <= 0.0 or relative_momentum_step >= 1.0:
            raise ValueError("relative_momentum_step must be between 0 and 1")
        try:
            x_plus, y_plus = self._reference_orbit_at_momentum(
                self.Pref0 * (1.0 + relative_momentum_step)
            )
            x_minus, y_minus = self._reference_orbit_at_momentum(
                self.Pref0 * (1.0 - relative_momentum_step)
            )
        finally:
            # Reference-probe tracks update BPM readbacks; restore the public
            # interface state even if a probe fails.
            self._track_bunch()
        result = {
            "names": list(self.bpms),
            "x": (x_plus - x_minus) / (2.0 * relative_momentum_step),
            "y": (y_plus - y_minus) / (2.0 * relative_momentum_step),
            "relative_momentum_step": float(relative_momentum_step),
        }
        if isinstance(names, str):
            names = [names]
        if names is None:
            return result
        indices = [index for index, name in enumerate(self.bpms) if name in names]
        return {
            "names": [self.bpms[index] for index in indices],
            "x": result["x"][indices],
            "y": result["y"][indices],
            "relative_momentum_step": float(relative_momentum_step),
        }

    def get_target_dispersion(self, names=None):
        dispersion = self.get_model_dispersion(names)
        return dispersion["x"], dispersion["y"]

    def change_energy(self):
        self._setup_beam(population=self.population, momentum_mev_c=self.dfs_test_energy * self.Pref0)
        self._track_bunch()
        self._beam_mode = "energy_changed"
        return self.dfs_test_energy - 1.0

    def reset_energy(self):
        self._setup_beam(population=self.population, momentum_mev_c=self.Pref0)
        self._track_bunch()
        self._beam_mode = "nominal"

    def change_intensity(self):
        self._setup_beam(population=self.wfs_test_charge * self.population, momentum_mev_c=self.Pref0)
        self._track_bunch()
        self._beam_mode = "intensity_changed"

    def reset_intensity(self):
        self._setup_beam(population=self.population, momentum_mev_c=self.Pref0)
        self._track_bunch()
        self._beam_mode = "nominal"
