import numpy as np
import time, os
from datetime import datetime
from Interfaces.AbstractMachineInterface import AbstractMachineInterface
import RF_Track as rft
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "FACET2"))
import FACET2

class InterfaceFACET2_Linac_RFTrack(AbstractMachineInterface):
    def get_name(self):
        return 'FACET2_Linac_RFT'

    def __init__(self, population=rft.nC, jitter=0.0, bpm_resolution=0.0, nsamples=1, nparticles=1000):
        super().__init__()
        self.log = print
        self.lattice = FACET2.load_FACET()
        self.lattice.set_bpm_resolution(bpm_resolution)
        self.sequence = [ e.get_name() for e in self.lattice['*']]
        self.bpms = [ e.get_name() for e in self.lattice.get_bpms()]
        self.corrs = [ e.get_name() for e in self.lattice.get_correctors()]
        self.screens = [ e.get_name() for e in self.lattice.get_screens() if e.get_name().startswith(("IM", "WS"))]
        self.sextupoles = []
        self.quadrupoles = list(dict.fromkeys(e.get_name() for e in self.lattice.get_quadrupoles()))
        self.Pref = np.sqrt(125.0**2 - rft.electronmass**2)
        self.bunch_length_ps = 3
        self.nparticles = nparticles
        self.population = population
        self.jitter = jitter
        self.nsamples = nsamples
        self.dfs_test_energy = 0.98
        self.wfs_test_charge = 0.90
        self._beam_mode = "nominal"
        self.Q = -1
        self.__setup_beam0()
        self.__track_bunch()
        self._saved_sextupoles_state = None
        self.electronmass = rft.electronmass
        self.machine_name = "FACET"
        self.lattice.align_elements()
        '''Uncomment lines below to scatter elements in the lattice.'''
        # self.lattice.scatter_elements('bpm', 0.100, 0.100, 0, 0, 0, 0, 'center')
        # self.lattice.scatter_elements('quadrupole', 0.100, 0.100, 0, 0, 0, 0, 'center')


    def log_messages(self,console):
        self.log=console or print

    def __setup_beam0(self):
        T = rft.Bunch6d_twiss()
        T.beta_x = 5.6
        T.alpha_x = -2.11
        T.beta_y = 2.9
        T.alpha_y = 0.0
        T.sigma_t = self.bunch_length_ps * rft.ps
        T.sigma_pt = 0.0
        T.emitt_x = 3.2
        T.emitt_y = 3.2
        T.sigma_pt = 0.8 # permille
        self.B0 = rft.Bunch6d_QR(rft.electronmass, self.population, -1, self.Pref, T, self.nparticles)

    def __setup_beam1(self):
        # Beam for DFS - Reduced energy
        Pref = self.dfs_test_energy * self.Pref
        T = rft.Bunch6d_twiss()
        T.beta_x = 5.6
        T.alpha_x = -2.11
        T.beta_y = 2.9
        T.alpha_y = 0.0
        T.sigma_t = self.bunch_length_ps * rft.ps
        T.sigma_pt = 0.0
        T.emitt_x = 3.2
        T.emitt_y = 3.2
        T.sigma_pt = 0.8 # permille
        self.B0 = rft.Bunch6d_QR(rft.electronmass, self.population, -1, Pref, T, self.nparticles)

    def __setup_beam2(self):
        # Beam for WFS - Reduced bunch charge
        population = self.wfs_test_charge * self.population
        T = rft.Bunch6d_twiss()
        T.beta_x = 5.6
        T.alpha_x = -2.11
        T.beta_y = 2.9
        T.alpha_y = 0.0
        T.sigma_t = self.bunch_length_ps * rft.ps
        T.sigma_pt = 0.0
        T.emitt_x = 3.2
        T.emitt_y = 3.2
        T.sigma_pt = 0.8 # permille
        self.B0 = rft.Bunch6d_QR(rft.electronmass, population, -1, self.Pref, T, self.nparticles)

    def __track_bunch(self):
        I0 = self.B0.get_info()
        dx = self.jitter*I0.sigma_x
        dy = self.jitter*I0.sigma_y
        dz, dt, roll = 0.0, 0.0, 0.0
        pitch = self.jitter*I0.sigma_py
        yaw   = self.jitter*I0.sigma_px
        B0_offset = self.B0.displaced(dx, dy, dz, dt, roll, pitch, yaw)
        B1 = self.lattice.track(B0_offset)
        I = B1.get_info()
        # self.log("Emittance after tracking:")
        # self.log(f"εx = {I.emitt_x}[mm.mrad]")
        # self.log(f"εy = {I.emitt_y}[mm.mrad]")
        # self.log(f"εz = {I.emitt_z}[mm.permille]")

    def get_beam_factors(self):
        gamma_rel = np.sqrt((self.Pref / self.electronmass) ** 2 + 1.0)
        beta_rel = np.sqrt(1.0 - 1.0 / gamma_rel ** 2)
        return gamma_rel, beta_rel

    def change_energy(self):
        self.__setup_beam1()
        self.__track_bunch()
        self._beam_mode = "energy_changed"
        dP_P = self.dfs_test_energy - 1.0
        return dP_P

    def reset_energy(self):
        self.__setup_beam0()
        self.__track_bunch()
        self._beam_mode = "nominal"

    def change_intensity(self): #reduced charge
        self.__setup_beam2()
        self.__track_bunch()
        self._beam_mode = "intensity_changed"

    def reset_intensity(self):
        self.__setup_beam0()
        self.__track_bunch()
        self._beam_mode = "nominal"

    def get_sequence(self):
        return self.sequence

    def get_hcorrectors_names(self):
        return [string for string in self.corrs if string.lower().startswith('x')]

    def get_vcorrectors_names(self):
        return [string for string in self.corrs if string.lower().startswith('y')]

    def get_elements_position(self,names):
        return [index for index, string in enumerate(self.sequence) if string in names]

    def get_icts(self, names=None):
        self.log("Reading ict's...")
        icts = {
            "names": self.bpms,
            "charge": np.array([bpm.get_total_charge() for bpm in self.lattice.get_bpms()])
        }
        if isinstance(names, str):
            names = [names]
        if names is not None:
            idx = [i for i, s in enumerate(icts["names"]) if s in names]
            icts = {
                "names": [icts["names"][i] for i in idx],
                "charge": np.asarray(icts["charge"])[idx],
            }
        return icts

    def get_quadrupoles(self, names=None):
        #self.log("Reading quadrupoles' strengths...")
        bdes = np.zeros(len(self.quadrupoles), dtype=float)

        for i, quadrupole_name in enumerate(self.quadrupoles):
            elements = self.lattice[quadrupole_name]
            if not isinstance(elements, list):
                elements = [elements]

            k1l_values = []
            for element in elements:
                try:
                    strength = element.get_K1(self.Pref / self.Q)
                except Exception:
                    continue
                if isinstance(strength, (list, tuple, np.ndarray)):
                    if len(strength) > 0:
                        k1l_values.append(float(strength[0]))
                else:
                    k1l_values.append(float(strength))

            bdes[i] = k1l_values[0] if k1l_values else 0.0

        quadrupoles = {"names": self.quadrupoles, "bdes": bdes, "bact": bdes.copy()}

        if isinstance(names, str):
            names = [names]
        if names is not None:
            idx = [i for i, s in enumerate(quadrupoles["names"]) if s in names]
            quadrupoles = {
                "names": [quadrupoles["names"][i] for i in idx],
                "bdes": np.asarray(quadrupoles["bdes"])[idx],
                "bact": np.asarray(quadrupoles["bact"])[idx],
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
        # AS A TEST!
            self.__track_bunch()



    def get_correctors(self,names=None):
        self.log("Reading correctors' strengths...")
        bdes = np.zeros(len(self.corrs))
        for i,corrector in enumerate(self.corrs):
            if corrector[:2] == "XC":
                bdes[i] = (self.lattice[corrector].get_strength()[0]*10)  # gauss*m
            elif corrector[:2] == "YC":
                bdes[i] = (self.lattice[corrector].get_strength()[1]*10)  # gauss*m
        correctors = { "names": self.corrs, "bdes": bdes, "bact": bdes }

        if isinstance(names, str):
            names = [names]
        if names is not None:
            idx = [i for i, s in enumerate(correctors["names"]) if s in names]
            correctors = {
                "names": [correctors["names"][i] for i in idx],
                "bdes": np.asarray(correctors["bdes"])[idx],
                "bact": np.asarray(correctors["bact"])[idx],
            }

        return correctors

    def get_bpms(self,names=None):
        self.log('Reading bpms...')
        x = np.zeros((self.nsamples, len(self.bpms)))
        y = np.zeros(x.shape)
        tmit = np.zeros(x.shape)
        for i in range(self.nsamples):
            for j,bpm in enumerate(self.bpms):
                b = self.lattice[bpm]
                reading = b.get_reading()
                x[i,j] = reading[0]
                y[i,j] = reading[1]
                tmit[i,j] = b.get_total_charge()
        bpms = { "names": self.bpms, "x": x, "y": y, "tmit": tmit }
        if isinstance(names, str):
            names = [names]
        if names is not None:
            idx = [i for i, s in enumerate(bpms["names"]) if s in names]
            bpms = {
                "names": [bpms["names"][i] for i in idx],
                "x": np.asarray(bpms["x"])[:, idx],
                "y": np.asarray(bpms["y"])[:, idx],
                "tmit": np.asarray(bpms["tmit"])[:, idx],
            }

        return bpms

    def _give_elements_to_show_beamline(self, quad_selected):
        start_quad_element_name = quad_selected
        return start_quad_element_name

    def _get_elements_positions_show_beamline(self, names=None):
        if isinstance(names, str):
            names = [names]
        all_names = []
        all_s = []
        for element in self.lattice['*']:
            element_name = element.get_name()
            if names is None or element_name in names:
                all_names.append(element_name)
                all_s.append(float(element.get_S("entrance")))
        return {
            "names": all_names,
            "S": np.array(all_s, dtype=float),
        }

    def get_screens(self, names=None):
        if isinstance(names, str):
            names = [names]
        hpixel = 0.001  # mm
        vpixel = 0.001  # mm

        selected_screens = [
            screen for screen in self.screens
            if names is None or screen in names
        ]

        hpixel_list = []
        vpixel_list = []
        xb_list = []
        yb_list = []
        sigx_list = []
        sigy_list = []
        sum_list = []
        images = []
        hedges_all = []
        vedges_all = []
        screen_names = []
        s_list = []

        for screen_name in selected_screens:
            screen = self.lattice[screen_name]
            screen_names.append(screen_name)
            s_list.append(float(screen.get_S("entrance")))
            hpixel_list.append(hpixel)
            vpixel_list.append(vpixel)
            bunch = screen.get_bunch()
            if bunch is None:
                xb_list.append(np.nan)
                yb_list.append(np.nan)
                sigx_list.append(np.nan)
                sigy_list.append(np.nan)
                sum_list.append(0.0)
                images.append(np.zeros((1, 1)))
                hedges_all.append(np.array([0.0, hpixel]))
                vedges_all.append(np.array([0.0, vpixel]))
                continue

            m = bunch.get_phase_space('%x %y')

            if m is None or len(m) == 0:
                xb_list.append(np.nan)
                yb_list.append(np.nan)
                sigx_list.append(np.nan)
                sigy_list.append(np.nan)
                sum_list.append(0.0)
                images.append(np.zeros((1, 1)))
                hedges_all.append(np.array([0.0, hpixel]))
                vedges_all.append(np.array([0.0, vpixel]))
                continue

            sumw = len(m[:, 0])

            xb_list.append(np.mean(m[:, 0]))
            yb_list.append(np.mean(m[:, 1]))
            sigx_list.append(np.std(m[:, 0]))
            sigy_list.append(np.std(m[:, 1]))
            sum_list.append(sumw)

            nx = int(np.ceil(np.ptp(m[:, 0]) / hpixel)) if np.ptp(m[:, 0]) > 0 else 1
            ny = int(np.ceil(np.ptp(m[:, 1]) / vpixel)) if np.ptp(m[:, 1]) > 0 else 1

            nx = int(np.clip(nx, 10, 400))
            ny = int(np.clip(ny, 10, 400))

            image, hedges, vedges = np.histogram2d(
                m[:, 0],
                m[:, 1],
                bins=(nx, ny),
            )

            images.append(image)
            hedges_all.append(hedges)
            vedges_all.append(vedges)

        screens = {
            "names": screen_names,
            "hpixel": np.array(hpixel_list, dtype=float),
            "vpixel": np.array(vpixel_list, dtype=float),
            "x": np.array(xb_list, dtype=float),
            "y": np.array(yb_list, dtype=float),
            "sigx": np.array(sigx_list, dtype=float),
            "sigy": np.array(sigy_list, dtype=float),
            "sum": np.array(sum_list, dtype=float),
            "hedges": hedges_all,
            "vedges": vedges_all,
            "images": images,
            "S": np.array(s_list, dtype=float),
        }

        return screens

    def set_correctors(self, names, corr_vals):
        if not isinstance(names, list):
            names = [ names ] # makes it a list
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        for corr, val in zip(names, corr_vals):
            if corr[:2] == "XC":
                self.lattice[corr].set_strength(val/10, 0.0)  # T*mm
            elif corr[:2] == "YC":
                self.lattice[corr].set_strength(0.0, val/10)  # T*mm
        self.__track_bunch()

    def vary_correctors(self, names, corr_vals):
        if not isinstance(names, list):
            names = [ names ] # makes it a list
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        for corr, val in zip(names, corr_vals):
            if corr[:2] == "XC":
                self.lattice[corr].vary_strength(val/10, 0.0)  # T*mm
            elif corr[:2] == "YC":
                self.lattice[corr].vary_strength(0.0, val/10)  # T*mm
        self.__track_bunch()

    def _get_elements_positions(self, names=None):
        if isinstance(names, str):
            names = [names]
        all_names = []
        all_s = []
        all_l = []

        for element in self.lattice['*']:
            element_name = element.get_name()
            if names is None or element_name in names:
                all_names.append(element_name)
                all_s.append(float(element.get_S("entrance")))
                all_l.append(float(element.get_length()))

        return {
            "names": all_names,
            "S": np.array(all_s, dtype=float),
            "L": np.array(all_l, dtype=float),
        }

    def _build_bunch_from_guesses(self, emit_x, emit_y, beta_x0, beta_y0, alpha_x0, alpha_y0):
        T = rft.Bunch6d_twiss()
        T.emitt_x = float(emit_x)
        T.emitt_y = float(emit_y)
        T.beta_x = float(beta_x0)
        T.beta_y = float(beta_y0)
        T.alpha_x = float(alpha_x0)
        T.alpha_y = float(alpha_y0)
        T.sigma_t = self.bunch_length_ps * rft.ps
        T.sigma_pt = 0.8
        return rft.Bunch6d_QR(rft.electronmass, self.population, self.Q, self.Pref, T, self.nparticles)

    def _read_tracked_bunch_screen_sigmas(self, screens):
        screen_data = self.get_screens(names=screens)
        name_to_index = {name: i for i, name in enumerate(screen_data["names"])}
        sigx = np.full(len(screens), np.nan, dtype=float)
        sigy = np.full(len(screens), np.nan, dtype=float)

        for i, screen in enumerate(screens):
            idx = name_to_index.get(screen)
            if idx is None:
                continue
            sigx[i] = float(screen_data["sigx"][idx])
            sigy[i] = float(screen_data["sigy"][idx])

        return sigx, sigy

    def get_elements_indices(self, names):
        if isinstance(names, str):
            names = [names]
        name_to_index = {string: index for index, string in enumerate(self.sequence)}
        return [name_to_index.get(name, np.nan) for name in names]

    def predict_emittance_scan_response(self, quad_name, screens, k1l_values, emit_x, emit_y, beta_x0, beta_y0, alpha_x0, alpha_y0, stop_checker=None, reference_screen=None):
        screens = list(screens)
        k1l_values = np.asarray(k1l_values, dtype=float)

        if len(screens) == 0:
            raise RuntimeError("No screens provided for emittance scan prediction.")
        if quad_name not in self.quadrupoles:
            raise ValueError(f"Quadrupole {quad_name} not found in quadrupoles")
        missing_screens = [screen for screen in screens if screen not in self.screens]
        if missing_screens:
            raise ValueError(f"Screens not found: {missing_screens}")

        original_quad = self.get_quadrupoles(names=[quad_name])
        if len(original_quad["bdes"]) == 0:
            raise RuntimeError(f"Could not find original strength for quad {quad_name}")
        K1L_original = float(original_quad["bdes"][0])
        B0_original = self.B0

        output_x = np.full((len(k1l_values), len(screens)), np.nan, dtype=float)
        output_y = np.full((len(k1l_values), len(screens)), np.nan, dtype=float)

        try:
            for k, K1L in enumerate(k1l_values):
                if callable(stop_checker) and stop_checker():
                    raise RuntimeError("__OPTIMIZATION_STOP__")

                self.set_quadrupoles([quad_name], [float(K1L)], track=False)
                self.B0 = self._build_bunch_from_guesses(
                    emit_x=emit_x,
                    emit_y=emit_y,
                    beta_x0=beta_x0,
                    beta_y0=beta_y0,
                    alpha_x0=alpha_x0,
                    alpha_y0=alpha_y0,
                )
                self.__track_bunch()
                output_x[k, :], output_y[k, :] = self._read_tracked_bunch_screen_sigmas(screens)
        finally:
            self.B0 = B0_original
            self.set_quadrupoles([quad_name], [K1L_original], track=False)
            self.__track_bunch()

        return output_x, output_y
