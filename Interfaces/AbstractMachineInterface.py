from abc import ABC, abstractmethod
from Backend.State import State
import time
import numpy as np

class AbstractMachineInterface(ABC):

    def _wait_for_readback(self, read_value, target, *, description, tolerance=1e-4, timeout=10.0, poll_interval=0.05):
        t0 = time.perf_counter()
        last_value = np.nan
        while time.perf_counter() - t0 < timeout:
            try:
                value = np.asarray(read_value())
                last_value = float(value.flat[0]) if value.size else np.nan
            except Exception:
                last_value = np.nan
            if np.isfinite(last_value) and abs(last_value - float(target)) <= tolerance:
                return True
            time.sleep(poll_interval)
        log = getattr(self, "log", print)
        log(
            f"Warning: {description} did not reach target {float(target):.6g} "
            f"within {timeout:.2f}s. Last readback = {last_value:.6g}"
        )
        return False

    @abstractmethod
    def get_name(self):
        pass

    @abstractmethod
    def get_sequence(self):
        pass

    @abstractmethod
    def get_correctors(self, names=None):
        pass

    @abstractmethod
    def get_bpms(self, names=None):
        pass

    @abstractmethod
    def get_icts(self, names=None):
        pass

    def get_quadrupoles(self, names=None):
        return {
            "names": [],
            "bdes": np.array([]),
            "bact": np.array([]),
        }

    def get_sextupoles(self, names=None):
        return {
            "names": [],
            "bdes": np.array([]),
            "bact": np.array([]),
        }

    def set_sextupoles(self, names, values):
        raise NotImplementedError(f"{self.get_name()} does not implement set_sextupoles")

    def get_screens(self, names=None):
        return {"names": [], "hpixel": np.array([]), "vpixel": np.array([]), "x":np.array([]),"y":np.array([]), "sigx":np.array([]), "sigy":np.array([]),"sum":np.array([]),"hedges":[],"vedges":[],"images":[],"S":np.array([])}

    def get_target_dispersion(self, names=None):
        if names is None:
            names = self.bpms
        if isinstance(names, str):
            names = [names]
        return np.zeros(len(names)), np.zeros(len(names))

    @abstractmethod
    def set_correctors(self, names, values):
        pass

    @abstractmethod
    def change_energy(self):
        pass

    @abstractmethod
    def reset_energy(self):
        pass

    @abstractmethod
    def change_intensity(self):
        pass

    @abstractmethod
    def reset_intensity(self):
        pass

    def set_quadrupoles(self, names, values):
        raise NotImplementedError(f"{self.get_name()} does not implement set_quadrupoles")

    def get_beam_settings(self):
        if hasattr(self, "_beam_mode"):
            return {
                "simulation": {
                    "beam_mode": self._beam_mode,
                    "dfs_energy_scale": float(self.dfs_test_energy),
                    "wfs_charge_scale": float(self.wfs_test_charge),
                }
            }
        return {}

    def restore_beam_settings(self, settings):
        simulation = (settings or {}).get("simulation")
        if simulation is None:
            return False
        self.dfs_test_energy = float(simulation.get("dfs_energy_scale", self.dfs_test_energy))
        self.wfs_test_charge = float(simulation.get("wfs_charge_scale", self.wfs_test_charge))
        mode = simulation.get("beam_mode", "nominal")
        self.reset_energy()
        if mode == "energy_changed":
            self.change_energy()
        elif mode == "intensity_changed":
            self.change_intensity()
        return True

    def get_state(self):
        return State(
            correctors=self.get_correctors(),
            bpms=self.get_bpms(),
            icts=self.get_icts(),
            sequence=self.get_sequence(),
            hcorrectors_names=self.get_hcorrectors_names(),
            vcorrectors_names=self.get_vcorrectors_names(),
            quadrupoles=self.get_quadrupoles(),
            sextupoles=self.get_sextupoles(),
            beam_settings=self.get_beam_settings(),
            interface_id=f"{type(self).__module__}.{type(self).__name__}",
        )

    def restore_correctors_state(self, state):
        correctors = state.get_correctors()
        self.set_correctors(correctors["names"], correctors["bdes"])

    def restore_quadrupoles_state(self, state):
        quadrupoles = state.get_quadrupoles()
        if len(quadrupoles["names"]) > 0:
            try:
                self.set_quadrupoles(quadrupoles["names"], quadrupoles["bdes"])
            except NotImplementedError:
                pass

        if not hasattr(self, "apply_qmag_xyroll"):
            return
        xdes = np.asarray(quadrupoles.get("xdes", []), dtype=float)
        ydes = np.asarray(quadrupoles.get("ydes", []), dtype=float)
        rolldes = np.asarray(quadrupoles.get("rolldes", []), dtype=float)
        valid = np.isfinite(xdes) & np.isfinite(ydes) & np.isfinite(rolldes)
        if np.any(valid):
            names = np.asarray(quadrupoles["names"], dtype=object)[valid].tolist()
            self.apply_qmag_xyroll(names, xdes[valid], ydes[valid], rolldes[valid], wait=True)

    def restore_sextupoles_one_by_one(self, state, callback=None):
        sextupoles = state.get_sextupoles()
        names = list(sextupoles["names"])
        values = np.asarray(sextupoles["bdes"], dtype=float)
        for name, value in zip(names, values):
            self.set_sextupoles([name], [value])
            if callback is not None:
                callback(name)

    def restore_sextupoles_state(self, state):
        sextupoles = state.get_sextupoles()
        if len(sextupoles["names"]) > 0:
            try:
                self.set_sextupoles(sextupoles["names"], sextupoles["bdes"])
            except NotImplementedError:
                pass
