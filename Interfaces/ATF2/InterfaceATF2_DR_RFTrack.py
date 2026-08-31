import RF_Track as rft
import numpy as np
import time,os
from Backend.LogConsole import LogConsole
from datetime import datetime
from Interfaces.AbstractMachineInterface import AbstractMachineInterface
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import (
    REFERENCE_SAD_DAIHON,
    build_atf_dr_lattice,
)
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import (
    ATFDRRingCorrection,
)

class InterfaceATF2_DR_RFTrack(AbstractMachineInterface):

    def get_name(self):
        return 'ATF2_DR_RFT'

    def __init__(self, population=2e10, jitter=0.0, bpm_resolution=0.0, nsamples=1):
        self.log = print
        self.twiss_path=os.path.join(os.path.dirname(__file__),'DR_ATF2','ATF_DR_twiss_file.tws')
        self.lattice = build_atf_dr_lattice(rf_mode="disabled")
        self.lattice_reference = REFERENCE_SAD_DAIHON
        self.lattice.set_bpm_resolution(bpm_resolution)
        self.sequence = [ e.get_name() for e in self.lattice['*']]
        self.bpms = [ e.get_name() for e in self.lattice.get_bpms()]
        self.corrs = [ e.get_name() for e in self.lattice.get_correctors()]
        self.screens = []
        self.quadrupoles = list(dict.fromkeys(e.get_name() for e in self.lattice.get_quadrupoles()))
        self.Pref = 1.2999999e3 # 1.3 GeV/c
        self.population = population
        self.jitter = jitter
        self.nsamples = nsamples
        # The SAD DR dispersion correction uses an RF-frequency perturbation
        # corresponding to a sub-percent dp/p.  A 2% step is outside the
        # validated momentum acceptance of this transverse-only ring model.
        self.dfs_test_energy = 0.995
        self.wfs_test_charge = 0.90
        self.Q=-1
        self.electronmass = rft.electronmass
        # Interface corrector values are converted with strength = value / 10.
        # Using the same scale here makes response/correction deltas compatible
        # with get/set/vary_correctors.
        self.ring_correction = ATFDRRingCorrection(
            self.lattice,
            self.Pref,
            self.Q,
            population=self.population,
            actuator_scale=0.1,
        )
        self._nominal_model_dispersion = self.ring_correction.measure_dispersion(
            relative_momentum_step=1e-3,
        )
        self._closed_orbit = None
        self.__setup_beam0()
        self.__track_bunch()
        self._saved_sextupoles_state = None
        self.sextupoles = self._get_element_names_from_twiss_types({"SEXTUPOLE"})

    def _get_element_names_from_twiss_types(self, allowed_types): # because rf track doesn't have get sextupoles
        with open(self.twiss_path, "r") as file:
            lines = [line.strip() for line in file if line.strip()]

        header_idx = None
        format_idx = None
        name_column = None
        type_column = None

        for i, line in enumerate(lines):
            if not line.startswith("*"):
                continue

            columns = line.lstrip("*").split()
            if "NAME" not in columns:
                continue

            type_col_name = None
            if "KEYWORD" in columns:
                type_col_name = "KEYWORD"
            elif "TYPE" in columns:
                type_col_name = "TYPE"
            else:
                continue

            header_idx = i
            name_column = columns.index("NAME")
            type_column = columns.index(type_col_name)

            for j in range(i + 1, len(lines)):
                if lines[j].startswith("$"):
                    format_idx = j
                    break
            break

        if header_idx is None or format_idx is None:
            return []

        names = []
        seen = set()

        for line in lines[format_idx + 1:]:
            if line.startswith("@") or line.startswith("*") or line.startswith("$"):
                continue

            data = line.split()
            if len(data) <= max(name_column, type_column):
                continue

            elem_name = data[name_column].strip('"')
            elem_type = data[type_column].strip('"').upper()

            if elem_type in allowed_types and elem_name not in seen:
                names.append(elem_name)
                seen.add(elem_name)

        return names

    def log_messages(self,console):
        self.log=console or print

    def __setup_beam0(self):
        self._tracking_momentum = self.Pref
        self.ring_correction.population = self.population
        T = rft.Bunch6d_twiss()
        T.emitt_x = 5.2 # mm.mrad normalised emittance
        T.emitt_y = 0.03 # mm.mrad
        T.alpha_x = -4.46977512
        T.alpha_y = 1.286939438
        T.beta_x = 5.329581905 # m
        T.beta_y = 1.710527297 # m
        T.disp_x = 0.1327935005 # m
        T.disp_px = 0.1112759177
        T.sigma_t = 8 # mm/c
        T.sigma_pt = 0.8 # permille
        Nparticles = 10000 # number of macroparticles
        self.B0 = rft.Bunch6d_QR(rft.electronmass, self.population, -1, self.Pref, T, Nparticles)
        
    def __setup_beam1(self):
        # Beam for DFS - Reduced energy
        Pref = self.dfs_test_energy * self.Pref
        self._tracking_momentum = Pref
        self.ring_correction.population = self.population
        T = rft.Bunch6d_twiss()
        T.emitt_x = 5.2 # mm.mrad normalised emittance
        T.emitt_y = 0.03 # mm.mrad
        T.alpha_x = -4.46977512
        T.alpha_y = 1.286939438
        T.beta_x = 5.329581905 # m
        T.beta_y = 1.710527297 # m
        T.disp_x = 0.1327935005 # m
        T.disp_px = 0.1112759177
        T.sigma_t = 8 # mm/c
        T.sigma_pt = 0.8 # permille
        Nparticles = 10000 # number of macroparticles
        self.B0 = rft.Bunch6d_QR(rft.electronmass, self.population, -1, Pref, T, Nparticles)

    def __setup_beam2(self):
        # Beam for WFS - Reduced bunch charge
        population = self.wfs_test_charge * self.population
        self._tracking_momentum = self.Pref
        self.ring_correction.population = population
        #population = 0.90 * self.population # 90% of nominal charge
        T = rft.Bunch6d_twiss()
        T.emitt_x = 5.2 # mm.mrad normalised emittance
        T.emitt_y = 0.03 # mm.mrad
        T.alpha_x = -4.46977512
        T.alpha_y = 1.286939438
        T.beta_x = 5.329581905 # m
        T.beta_y = 1.710527297 # m
        T.disp_x = 0.1327935005 # m
        T.disp_px = 0.1112759177
        T.sigma_t = 8 # mm/c
        T.sigma_pt = 0.8 # permille
        Nparticles = 10000 # number of macroparticles
        self.B0 = rft.Bunch6d_QR(rft.electronmass, population, -1, self.Pref, T, Nparticles)

    def __track_bunch(self):
        I0 = self.B0.get_info()
        dx = self.jitter*I0.sigma_x
        dy = self.jitter*I0.sigma_y
        dz, roll, dt = 0.0, 0.0, 0.0
        pitch = self.jitter*I0.sigma_py
        yaw   = self.jitter*I0.sigma_px
        B0_offset = self.B0.displaced(dx, dy, dz, dt, roll, pitch, yaw)
        tracked_bunch = self.lattice.track(B0_offset)
        # A storage ring BPM observes the periodic closed orbit, not a one-pass
        # trajectory launched at the lattice origin.  This final probe track
        # also leaves the RF-Track BPM objects with the correct COD readings.
        self._closed_orbit = self.ring_correction.find_closed_orbit(
            momentum_mev_c=self._tracking_momentum,
        )
        I=tracked_bunch.get_info()
        self.log("Emittance after tracking:")
        self.log(f"εx = {I.emitt_x}[mm.rad]")
        self.log(f"εy = {I.emitt_y}[mm.rad]")
        self.log(f"εz = {I.emitt_z}[mm.permille]")

    def get_beam_factors(self):
        gamma_rel = np.sqrt((self.Pref / self.electronmass) ** 2 + 1.0)
        beta_rel = np.sqrt(1.0 - 1.0 / gamma_rel ** 2)
        return gamma_rel, beta_rel

    def change_energy(self):
        self.__setup_beam1()
        self.__track_bunch()
        dP_P = self.dfs_test_energy - 1.0
        return dP_P

    def reset_energy(self):
        self.__setup_beam0()
        self.__track_bunch()

    def change_intensity(self): #reduced charge
        self.__setup_beam2()
        self.__track_bunch()

    def reset_intensity(self):
        self.__setup_beam0()
        self.__track_bunch()

    def get_sequence(self):
        return self.sequence

    def get_closed_orbit(self, names=None, momentum_mev_c=None):
        """Return the periodic RF-Track COD in mm at the selected BPMs."""
        if isinstance(names, str):
            names = [names]
        result = self.ring_correction.find_closed_orbit(
            momentum_mev_c=momentum_mev_c,
            bpm_names=names,
        )
        self._closed_orbit = result
        return {
            "names": list(result.bpm_names),
            "x": result.x.copy(),
            "y": result.y.copy(),
            "initial_coordinates": result.initial_coordinates.copy(),
            "residual_norm": result.residual_norm,
            "iterations": result.iterations,
            "momentum_mev_c": result.momentum_mev_c,
        }

    def get_model_dispersion(self, names=None, relative_momentum_step=1e-3):
        """Return periodic model dispersion in mm per unit dp/p."""
        if isinstance(names, str):
            names = [names]
        result = self.ring_correction.measure_dispersion(
            relative_momentum_step=relative_momentum_step,
            bpm_names=names,
        )
        return {
            "names": list(result.bpm_names),
            "x": result.x.copy(),
            "y": result.y.copy(),
            "relative_momentum_step": result.relative_momentum_step,
        }

    def get_nominal_model_dispersion(self, names=None):
        """Return the nominal SAD-derived RF-Track dispersion target in mm."""
        if isinstance(names, str):
            names = [names]
        target = self._nominal_model_dispersion
        if names is None:
            indices = np.arange(len(target.bpm_names))
        else:
            name_to_index = {
                name: index for index, name in enumerate(target.bpm_names)
            }
            unknown = [name for name in names if name not in name_to_index]
            if unknown:
                raise ValueError(f"Unknown BPMs: {unknown}")
            indices = np.array([name_to_index[name] for name in names])
        return {
            "names": [target.bpm_names[index] for index in indices],
            "x": target.x[indices].copy(),
            "y": target.y[indices].copy(),
            "relative_momentum_step": target.relative_momentum_step,
        }

    def compute_periodic_orbit_response(
        self,
        plane,
        corrector_names=None,
        bpm_names=None,
        perturbation=1e-4,
    ):
        """Build a periodic COD response in mm per interface corrector unit."""
        return self.ring_correction.compute_orbit_response(
            plane,
            corrector_names=corrector_names,
            bpm_names=bpm_names,
            perturbation=perturbation,
        )

    def suggest_periodic_orbit_correction(self, response, **kwargs):
        """Calculate, but do not apply, a periodic COD correction."""
        return self.ring_correction.propose_orbit_correction(response, **kwargs)

    def apply_periodic_orbit_correction(self, correction):
        """Apply a suggestion and return the resulting periodic COD."""
        result = self.ring_correction.apply_orbit_correction(correction)
        self._closed_orbit = result
        return result

    def compute_periodic_dispersion_response(
        self,
        plane,
        corrector_names=None,
        bpm_names=None,
        relative_momentum_step=1e-3,
        perturbation=1e-4,
    ):
        """Build joint periodic dispersion and COD steerer responses."""
        return self.ring_correction.compute_dispersion_response(
            plane,
            corrector_names=corrector_names,
            bpm_names=bpm_names,
            relative_momentum_step=relative_momentum_step,
            perturbation=perturbation,
        )

    def suggest_periodic_dispersion_correction(
        self,
        response,
        target_dispersion=None,
        **kwargs,
    ):
        """Calculate, but do not apply, a joint dispersion/COD correction."""
        if target_dispersion is None:
            target = self.get_nominal_model_dispersion(response.bpm_names)
            target_dispersion = target[response.plane]
        return self.ring_correction.propose_dispersion_correction(
            response,
            target_dispersion=target_dispersion,
            **kwargs,
        )

    def apply_periodic_dispersion_correction(self, correction):
        """Apply a joint dispersion/COD suggestion and remeasure it."""
        orbit, dispersion = self.ring_correction.apply_dispersion_correction(
            correction
        )
        self._closed_orbit = orbit
        return orbit, dispersion

    def get_skew_correctors_names(self):
        """Return the SD1R/SF1R skew-K1L actuators used for coupling control."""
        return list(self.ring_correction.get_skew_corrector_names())

    def get_skew_correctors(self, names=None):
        """Read skew-quadrupole K1L values (normalized integrated strength)."""
        available = self.get_skew_correctors_names()
        if isinstance(names, str):
            names = [names]
        selected = available if names is None else list(names)
        unknown = [name for name in selected if name not in available]
        if unknown:
            raise ValueError(f"Unknown skew correctors: {unknown}")
        values = np.array(
            [self.ring_correction.get_skew_strength(name) for name in selected]
        )
        return {"names": selected, "bdes": values, "bact": values.copy()}

    def set_skew_correctors(self, names, values, track=True):
        """Set SD1R/SF1R skew-K1L values and optionally refresh BPM readings."""
        if isinstance(names, str):
            names = [names]
        if not isinstance(values, (list, tuple, np.ndarray)):
            values = [values]
        for name, value in zip(names, values):
            self.ring_correction.set_skew_strength(name, float(value))
        if track:
            self.__track_bunch()

    def compute_periodic_coupling_response(
        self,
        probe_corrector_names=None,
        skew_corrector_names=None,
        bpm_names=None,
        probe_perturbation=1e-4,
        skew_perturbation=1e-5,
    ):
        """Build the SAD-style horizontal-kick to vertical-COD skew response."""
        return self.ring_correction.compute_coupling_response(
            probe_corrector_names=probe_corrector_names,
            skew_corrector_names=skew_corrector_names,
            bpm_names=bpm_names,
            probe_perturbation=probe_perturbation,
            skew_perturbation=skew_perturbation,
        )

    def suggest_periodic_coupling_correction(self, response, **kwargs):
        """Calculate, but do not apply, a skew-coupling correction."""
        return self.ring_correction.propose_coupling_correction(response, **kwargs)

    def apply_periodic_coupling_correction(self, correction):
        """Apply a skew-coupling suggestion and return its measured signal."""
        return self.ring_correction.apply_coupling_correction(correction)

    def get_hcorrectors_names(self):
        return [string for string in self.corrs if (string.lower().startswith('zh')) or (string.lower().startswith('zx'))]

    def get_vcorrectors_names(self):
        return [string for string in self.corrs if string.lower().startswith('zv')]

    def get_elements_indices(self, names):
        if isinstance(names, str):
            names = [names]
        name_to_index = {string: index for index, string in enumerate(self.sequence)}
        return [name_to_index.get(name, np.nan) for name in names]

    def get_target_dispersion(self, names=None):
        # The generic correction GUIs expect target dispersion in metres.
        # Use the nominal lattice paired with this RF-Track model instead of
        # the legacy TWS file from a different SAD revision.
        target = self.get_nominal_model_dispersion(names)
        return target["x"] / 1e3, target["y"] / 1e3

    def get_icts(self, names=None):
        self.log("Reading ict's...")
        icts = {
            "names": self.bpms,
            "charge": np.array([bpm.get_total_charge() for bpm in self.lattice.get_bpms()])
        }
        if isinstance(names, str):
            names = [names]
        if names is not None:
            idx = np.array([i for i, s in enumerate(icts["names"]) if s in names])
            icts = {
                "names": np.array(icts["names"])[idx],
                "charge": np.array(icts["charge"])[idx],
            }

        return icts

    def get_correctors(self, names=None):
        self.log("Reading correctors' strengths...")
        bdes = np.zeros(len(self.corrs))
        for i, corrector in enumerate(self.corrs):
            if corrector[:2] == "ZH" or corrector[:2] == "ZX":
                bdes[i] = self.lattice[corrector].get_strength()[0] * 10
            elif corrector[:2] == "ZV":
                bdes[i] = self.lattice[corrector].get_strength()[1] * 10

        correctors = {"names": self.corrs, "bdes": bdes, "bact": bdes.copy()}

        if isinstance(names, str):
            names = [names]
        if names is not None:
            idx = np.array([i for i, s in enumerate(correctors["names"]) if s in names])
            correctors = {
                "names": np.array(correctors["names"])[idx],
                "bdes": np.array(correctors["bdes"])[idx],
                "bact": np.array(correctors["bact"])[idx],
            }

        return correctors

    def get_bpms(self, names=None):
        self.log('Reading bpms...')
        x = np.zeros((self.nsamples, len(self.bpms)))
        y = np.zeros(x.shape)
        tmit = np.zeros(x.shape)

        for i in range(self.nsamples):
            for j, bpm in enumerate(self.bpms):
                b = self.lattice[bpm]
                reading = b.get_reading()
                x[i, j] = reading[0]
                y[i, j] = reading[1]
                tmit[i, j] = b.get_total_charge()

        bpms = {"names": self.bpms, "x": x, "y": y, "tmit": tmit}

        if isinstance(names, str):
            names = [names]
        if names is not None:
            idx = np.array([i for i, s in enumerate(bpms["names"]) if s in names])
            bpms = {
                "names": np.array(bpms["names"])[idx],
                "x": np.array(bpms["x"])[:, idx],
                "y": np.array(bpms["y"])[:, idx],
                "tmit": np.array(bpms["tmit"])[:, idx],
            }

        return bpms

    def get_quadrupoles(self, names=None):
        self.log("Reading quadrupoles' strengths...")
        bdes = np.zeros(len(self.quadrupoles), dtype=float)

        for i, quadrupole_name in enumerate(self.quadrupoles):
            elements = self.lattice[quadrupole_name]
            if not isinstance(elements, list):
                elements = [elements]

            k1_values = []
            for element in elements:
                try:
                    strength = element.get_K1(self.Pref / self.Q)
                except Exception:
                    continue
                if isinstance(strength, (list, tuple, np.ndarray)):
                    if len(strength) > 0:
                        k1_values.append(float(strength[0]))
                else:
                    k1_values.append(float(strength))

            bdes[i] = k1_values[0] if k1_values else 0.0

        quadrupoles = {"names": self.quadrupoles, "bdes": bdes, "bact": bdes.copy()}

        if isinstance(names, str):
            names = [names]
        if names is not None:
            idx = np.array([i for i, s in enumerate(quadrupoles["names"]) if s in names])
            quadrupoles = {
                "names": np.array(quadrupoles["names"])[idx],
                "bdes": np.array(quadrupoles["bdes"])[idx],
                "bact": np.array(quadrupoles["bact"])[idx],
            }

        return quadrupoles

    def set_quadrupoles(self, names, values_range, track = True):
        if isinstance(names, str):
            names = [names]
        if not (isinstance(values_range, (list, tuple, np.ndarray))):
            values_range = [values_range]
        for quadrupole_name, value in zip(names, values_range):
            elements = self.lattice[quadrupole_name]
            if not isinstance(elements, (list)): elements = [elements]
            for element in elements:
                element.set_K1(self.Pref / self.Q,float(value))
        if track:
            self.__track_bunch()

    def set_correctors(self, names, corr_vals):
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        for corr, val in zip(names, corr_vals):
            if corr[:2] == "ZH" or corr[:2] == "ZX":
                self.lattice[corr].set_strength(val / 10, 0.0)
            elif corr[:2] == "ZV":
                self.lattice[corr].set_strength(0.0, val / 10)
        self.__track_bunch()

    def vary_correctors(self, names, corr_vals):
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        for corr, val in zip(names, corr_vals):
            if corr[:2] == "ZH" or corr[:2] == "ZX":
                self.lattice[corr].vary_strength(val / 10, 0.0)
            elif corr[:2] == "ZV":
                self.lattice[corr].vary_strength(0.0, val / 10)
        self.__track_bunch()

    def vary_quadrupoles(self, names, delta_values):
        if not isinstance(names, list):
            names = [names]
        if not isinstance(delta_values, (list, tuple, np.ndarray)):
            delta_values = [delta_values]
        for quadrupole_name, val in zip(names, delta_values):
            elements = self.lattice[quadrupole_name]
            if not isinstance(elements, list):
                elements = [elements]
            current_values=[]
            for element in elements:
                current=element.get_K1(self.Pref / self.Q)
                current=float(current[0]) if isinstance(current, (list, tuple,np.ndarray)) else float(current)
                current_values.append(current)
            if len(current_values)>1 and not np.allclose(current_values, current_values[0], rtol=0.0, atol=1e-12):
                self.log(f"Parts of quadrupole {quadrupole_name} have different values")
            target_value=(current_values[0] if len(current_values)>0 else 0.0) +float(val)
            for element in elements:
                element.set_K1(self.Pref / self.Q,target_value)

        self.__track_bunch()

    def align_everything(self):
        self.lattice.align_elements()
        self.__track_bunch()

    def misalign_quadrupoles(self,sigma_x=0.02,sigma_y=0.02):
        self.lattice.scatter_elements('quadrupole', sigma_x, sigma_y, 0, 0, 0, 0, 'center')
        self.__track_bunch()

    def misalign_bpms(self,sigma_x=0.100,sigma_y=0.100):
        self.lattice.scatter_elements('bpm', sigma_x, sigma_y, 0, 0, 0, 0, 'center')
        self.__track_bunch()

    def get_sextupoles(self, names = None):
        self.log("Reading sextupoles' strengths...")
        bdes = np.zeros(len(self.sextupoles), dtype=float)

        for i, sextupole_name in enumerate(self.sextupoles):
            elements = self.lattice[sextupole_name]
            if not isinstance(elements, list):
                elements = [elements]

            k2_values = []
            for element in elements:
                try:
                    k2_values.append(float(element.get_K2(self.Pref / self.Q)))
                    continue
                except Exception:
                    pass
                try:
                    strengths = np.asarray(
                        element.get_strengths(), dtype=complex
                    ).ravel()
                except Exception:
                    continue
                k2_values.append(
                    float(np.real(strengths[2])) if strengths.size >= 3 else 0.0
                )
            if len(k2_values) > 1 and not np.allclose(k2_values, k2_values[0], rtol = 0.0, atol = 1e-12):
                self.log(f"Parts of sextupole {sextupole_name} are not consistent.")

            bdes[i] = k2_values[0] if k2_values else 0.0

        sextupoles = {"names": self.sextupoles, "bdes": bdes, "bact": bdes.copy()}

        if isinstance(names, str):
            names = [names]
        if names is not None:
            idx = np.array([i for i, s in enumerate(sextupoles["names"]) if s in names])
            sextupoles = {
                "names": np.array(sextupoles["names"])[idx],
                "bdes": np.array(sextupoles["bdes"])[idx],
                "bact": np.array(sextupoles["bact"])[idx],
            }
        return sextupoles

    def set_sextupoles(self, names, values_range):
        if isinstance(names, str):
            names = [names]
        if not (isinstance(values_range, (list, tuple, np.ndarray))):
            values_range = [values_range]
        for sextupole_name, value in zip(names, values_range):
            elements = self.lattice[sextupole_name]
            if not isinstance(elements, (list)): elements = [elements]
            for element in elements:
                try:
                    element.set_K2(self.Pref / self.Q, float(value))
                    continue
                except Exception:
                    pass
                try:
                    strengths = element.get_strengths()
                except Exception:
                    continue
                strengths = np.asarray(strengths, dtype = complex).ravel()

                if strengths.size < 3:
                    padded = np.zeros(3, dtype = complex)
                    padded[:strengths.size] = strengths
                    strengths = padded
                strengths[2] = complex(float(value), 0.0)
                element.set_strengths(strengths)
        self.__track_bunch()

    def _get_optics_from_twiss_file(self,names=None):
        with open(self.twiss_path, "r") as file:
            lines = [line.strip() for line in file if line.strip()]

        star_symbol = next(i for i, line in enumerate(lines) if line.startswith("*"))
        dollar_sign = next(i for i, line in enumerate(lines) if line.startswith("$") and i > star_symbol)
        columns = lines[star_symbol].lstrip("*").split()

        cols=["NAME","S","BETX","ALFX","BETY","ALFY","MUX","MUY","L"]
        index={}
        for col in cols:
            try:
                index[col]=columns.index(col)
            except ValueError:
                raise RuntimeError(f"Column {col} not found in twiss file")
        result={k: [] for k in cols}

        duplicated={}
        for line in lines[dollar_sign + 1:]:
            data = line.split()
            if len(data) <= max(index.values()):  # if a line has less column than needed, it is omitted
                continue
            name = data[index["NAME"]].strip('"')
            result["NAME"].append(name)

            if name not in duplicated:
                duplicated[name]={col:[] for col in cols[1:]}
            for col in cols[1:]:
                try:
                    duplicated[name][col].append(float(data[index[col]]))
                except ValueError:
                    duplicated[name][col].append(float("nan"))
        result_names=list(duplicated.keys())
        result={col:[] for col in cols[1:]}

        for name in result_names:
            vals=duplicated[name]
            result["S"].append(vals["S"][0])
            result["BETX"].append(vals["BETX"][0])
            result["ALFX"].append(vals["ALFX"][0])
            result["MUX"].append(vals["MUX"][0])
            result["L"].append(sum(vals["L"]))
            result["BETY"].append(vals["BETY"][0])
            result["ALFY"].append(vals["ALFY"][0])
            result["MUY"].append(vals["MUY"][0])
        return {
            "names":result_names,
            "S": np.array(result["S"]),
            "betx": np.array(result["BETX"]),
            "alfx": np.array(result["ALFX"]),
            "bety": np.array(result["BETY"]),
            "alfy": np.array(result["ALFY"]),
            "mux": np.array(result["MUX"]),
            "muy": np.array(result["MUY"]),
            "L": np.array(result["L"]),
        }

    def get_twiss_at_element(self,name):
        optics=self._get_optics_from_twiss_file()
        names=list(optics["names"])
        if name not in names:
            raise ValueError(f"Element {name} not found in twiss file")
        i=names.index(name)
        return {
            "name": names[i],
            "S": float(optics["S"][i]),
            "betx": float(optics["betx"][i]),
            "alfx": float(optics["alfx"][i]),
            "bety": float(optics["bety"][i]),
            "alfy": float(optics["alfy"][i]),

        }
