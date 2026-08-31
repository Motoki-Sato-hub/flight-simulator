import numpy as np
import matplotlib.pyplot as plt
import RF_Track as rft
from scipy.optimize import minimize
import os
from Backend.State import State
from Interfaces.AbstractMachineInterface import AbstractMachineInterface

class InterfaceCLEAR_RFTrack(AbstractMachineInterface):
    def get_name(self):
        return 'CLEAR_RFT'

    def get_ITF(self, I):
        return 1.29404711e-2 - 2.59458259e-07 * I  # T/A

    def get_grad(self, I, Lquad=0.226):
        G_0 = I * self.get_ITF(I) / Lquad  # T/m
        return G_0

    def get_Quad_K(self, G_0, Pref):
        K = 299.8 * G_0 / Pref  # 1/m^2
        return K

    def get_Quad_K_from_I(self, I, Lquad, Pref):
        G_0 = self.get_grad(I, Lquad)
        K = self.get_Quad_K(G_0, Pref)
        return K # 1/m^2

    @staticmethod
    def _replace_btv_monitors_with_screens(lattice):
        for element in list(lattice['*']):
            if not element.get_name().startswith("CA.BTV"):
                continue
            if hasattr(element, "get_bunch"):
                continue
            screen = rft.Screen()
            screen.set_name(element.get_name())
            screen.set_length(element.get_length())
            element.replace_with(screen)

    def __init__(self, population=300 * rft.pC, jitter=0.0, bpm_resolution=0.0, nsamples=1, nparticles=10000):
        self.sigmaCut = 2.0
        self.Pref = 198 # MeV/c
        self.Q=-1
        self.population = population
        self.jitter = jitter
        self.nsamples = nsamples
        self.nparticles = nparticles
        self.electronmass=rft.electronmass
        self.is_simulation = True
        survey_path = os.path.join(os.path.dirname(__file__),"Setup_files", "twissinit.tfs")
        self.lattice = rft.Lattice(survey_path)
        for element in self.lattice["*ICT*"]:
            drift = rft.Drift(element.get_length())
            drift.set_name(element.get_name())
            element.replace_with(drift)
        for element in self.lattice["*SCT*"]:
            drift = rft.Drift(element.get_length())
            drift.set_name(element.get_name())
            element.replace_with(drift)
        self._replace_btv_monitors_with_screens(self.lattice)
        elements_in_lattice = list(self.lattice['*'])
        self.sequence = [element.get_name() for element in elements_in_lattice]
        self.start = self.sequence[0]
        self.end = self.sequence[-1]
        self.lattice.set_bpm_resolution(bpm_resolution)
        self.lattice.set_tt_nsteps(0)
        self.log = print
        self.bpms = [element.get_name() for element in self.lattice.get_bpms()]
        self.corrs = [element.get_name() for element in self.lattice.get_correctors()]
        self.screens = [element.get_name() for element in self.lattice.get_screens()]
        self.quadrupoles = [element.get_name() for element in self.lattice.get_quadrupoles()]
        self.sextupoles = []
        self.__setup_beam0()
        '''Uncomment lines below to scatter elements in the lattice.'''
        self.lattice.scatter_elements('bpm', 0.100, 0.100, 0, 0, 0, 0, 'center')
        self.lattice.scatter_elements('quadrupole', 0.100, 0.100, 0, 0, 0, 0, 'center')
        self.freq=2.997e9
        self.nr_quad=11
        self.Lquad=0.226 #magnetic length of the quadrupole in [m]
        self.nominal_K=0.7752883624676146 #3.35  # 1/m
        self.machine_name = "CLEAR"
        self.lattice.align_elements()
        self.chosen_ict = "CA.BPM0890"
        # qfd520 = self.lattice["CA.QFD0520"]
        #                     # dx   # dy   #dz  # roll  # pitch # yaw
        # qfd520.set_offsets(0.0, 0.0, 0.0, 0, 0.0, 0.0, "center")

        # for element in self.lattice['*']:
        #     element.set_aperture(2e-3, 2e-3, "circular") # in reality, CLEAR has 40mm, or 30-20mm diameter

        '''test of bpm invertion'''
        self.invert_bpm = True


        self.__track_bunch()

    def __setup_beam0(self):
        T = rft.Bunch6d_twiss()
        T.emitt_x = 12.7000 # mm.mrad normalised emittance
        T.emitt_y = 4.3400   # mm.mrad
        T.beta_x = 46.1766 # m
        T.beta_y = 109.1117  # m; back-propagated from QFD350 entrance
        T.alpha_x = 1.3126
        T.alpha_y = 7.2462

        T.sigma_t = 0#10*rft.ps # mm/c
        T.sigma_pt = 0#10 # permille
        T.mean_xp = 0.0
        T.mean_yp = 0.0
        sigmaCut = 2.0
        self.P0 = rft.Bunch6d_QR(rft.electronmass, self.population, 1, self.Pref, T, self.nparticles, sigmaCut) # reference particle
        self.B0 = rft.Bunch6d_QR(rft.electronmass, self.population, self.Q, self.Pref, T, self.nparticles, sigmaCut) # reference bunch
        self.dfs_test_energy = 0.90 #0.963
        self.wfs_test_charge = 0.90
        self._beam_mode = "nominal"

    def __setup_beam1(self):
        # Beam for DFS - Reduced energy
        Pref = self.dfs_test_energy * self.Pref
        T = rft.Bunch6d_twiss()

        T.emitt_x = 12.7000 # mm.mrad normalised emittance
        T.emitt_y = 4.3400   # mm.mrad
        T.beta_x = 46.244900 # m
        T.beta_y = 109.143468  # m; back-propagated from QFD350 entrance
        T.alpha_x = 1.314308
        T.alpha_y = 7.248155

        T.sigma_t = 0#10*rft.ps # mm/c
        T.sigma_pt = 0#10 # permille
        T.mean_xp=0.0
        T.mean_yp=0.0
        self.B0 = rft.Bunch6d_QR(rft.electronmass, self.population, self.Q , Pref, T, self.nparticles, self.sigmaCut)
        self.P0 = rft.Bunch6d_QR(rft.electronmass, self.population,  1 , Pref, T, self.nparticles, self.sigmaCut)

    def __setup_beam2(self):
        # Beam for WFS - Reduced bunch charge
        population = self.wfs_test_charge * self.population
        T = rft.Bunch6d_twiss()

        T.emitt_x = 12.7000 # mm.mrad normalised emittance
        T.emitt_y = 4.3400   # mm.mrad
        T.beta_x = 46.244900 # m
        T.beta_y = 109.143468  # m; back-propagated from QFD350 entrance
        T.alpha_x = 1.314308
        T.alpha_y = 7.248155

        T.sigma_t = 0#10*rft.ps # mm/c
        T.sigma_pt = 0#10 # permille
        T.mean_xp=0.0
        T.mean_yp=0.0
        self.B0 = rft.Bunch6d_QR(rft.electronmass, population, self.Q, self.Pref, T, self.nparticles, self.sigmaCut)
        self.P0 = rft.Bunch6d_QR(rft.electronmass, population,  1, self.Pref, T, self.nparticles, self.sigmaCut)

    def get_screens(self, names=None):
        if isinstance(names, str): names = [names]
        hpixel = 0.001
        vpixel = 0.001
        selected_screens = [screen for screen in self.screens if names is None or screen in names]
        hpixel_list = []
        vpixel_list = []
        xb_list = []
        yb_list = []
        sigx_list = []
        sigy_list = []
        sigxy_list = []
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
            m = bunch.get_phase_space('%x %y')
            x_mean = float(np.mean(m[:, 0]))
            y_mean = float(np.mean(m[:, 1]))
            x_centered = m[:, 0] - x_mean # deviation from center of the beam; if the beam centroid changes, this will stay the same
            y_centered = m[:, 1] - y_mean # it has information about the tilt of ellipse

            xb_list.append(float(np.mean(m[:, 0])))
            yb_list.append(float(np.mean(m[:, 1])))
            sigx_list.append(float(np.std(m[:, 0])))
            sigy_list.append(float(np.std(m[:, 1])))
            sigxy_list.append(float(np.mean(x_centered * y_centered)))
            sum_list.append(float(len(m[:, 0])))

            nx = int(np.ceil(np.ptp(m[:, 0]) / hpixel)) if np.ptp(m[:, 0]) > 0 else 1
            ny = int(np.ceil(np.ptp(m[:, 1]) / vpixel)) if np.ptp(m[:, 1]) > 0 else 1
            nx = int(np.clip(nx, 10, 400))
            ny = int(np.clip(ny, 10, 400))

            image, hedges, vedges = np.histogram2d(m[:, 0], m[:, 1], bins=(nx, ny))
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
            "sigxy": np.array(sigxy_list, dtype=float),
            "sum": np.array(sum_list, dtype=float),
            "hedges": hedges_all,
            "vedges": vedges_all,
            "images": images,
            "S": np.array(s_list, dtype=float),
        }

        return screens

    def __track_bunch(self):
        I0 = self.B0.get_info()
        dx = self.jitter * I0.sigma_x
        dy = self.jitter * I0.sigma_y
        dz, dt, roll = 0.0, 0.0, 0.0
        pitch = self.jitter * I0.sigma_py
        yaw = self.jitter * I0.sigma_px
        B0_offset = self.B0.displaced(dx, dy, dz, dt, roll, pitch, yaw)
        B1=self.lattice.track(B0_offset)
        I = B0_offset.get_info()

    def get_beam_factors(self):
        gamma_rel = np.sqrt((self.Pref / self.electronmass) ** 2 + 1.0)
        beta_rel = np.sqrt(1.0 - 1.0 / gamma_rel ** 2)
        beta_gamma = gamma_rel * beta_rel
        if not np.isfinite(beta_gamma) or beta_gamma <= 0:
            raise RuntimeError("Invalid beam factors")
        return gamma_rel, beta_rel, beta_gamma

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

    def _get_elements_positions(self, names=None):
        if isinstance(names, str):
            names = [names]
        selected = [name for name in self.sequence if names is None or name in names]
        return {
            "names": selected,
            "S": np.array([self.lattice[name].get_S("entrance") for name in selected], dtype=float),
            "L": np.array([self.lattice[name].get_length() for name in selected], dtype=float),
        }

    def get_sequence(self):
        return self.sequence

    def get_hcorrectors_names(self):
        return [string for string in self.corrs if "DHG" in string]

    def get_vcorrectors_names(self):
        return [string for string in self.corrs if ("DHJ" in string) or ("SDV" in string) ]

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

    def get_correctors(self,names=None):
        #self.log("Reading correctors' strengths...")
        bdes = np.zeros(len(self.corrs))
        for i,corrector in enumerate(self.corrs):
            c=self.lattice[corrector]
            hx,hy=c.get_strength()
            if "DHG" in corrector: #horizontal
                bdes[i] = (hx*10)  # gauss*m
            elif ("SDV" in corrector) or ("DHJ" in corrector): #vertical
                bdes[i] = (hy*10)  # gauss*m
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

    def get_bpms(self, names=None):
        self.log('Reading bpms...')
        x = np.zeros((self.nsamples, len(self.bpms)))
        y = np.zeros_like(x)
        tmit = np.zeros_like(x)

        for i in range(self.nsamples):
            for j, bpm_name in enumerate(self.bpms):
                bpm = self.lattice[bpm_name]
                reading = bpm.get_reading()
                '''test of inverted bpm'''
                sign = -1.0 if (self.invert_bpm and bpm_name == "CA.BPM0890") else 1.0

                x[i, j] = sign * reading[0]
                y[i, j] = sign * reading[1]
                # x[i, j] = reading[0]
                # y[i, j] = reading[1]
                tmit[i, j] = bpm.get_total_charge()

        bpms = {"names": self.bpms, "x": x, "y": y, "tmit": tmit}

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

    def set_correctors(self, names, corr_vals):
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        for corr, val in zip(names, corr_vals):
            c = self.lattice[corr]
            if "DHG" in corr:
                c.set_strength(val / 10, 0.0)
            elif ("DHJ" in corr) or ("SDV" in corr):
                c.set_strength(0.0, val / 10)
        self.__track_bunch()

    def vary_correctors(self, names, corr_vals):
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        for corr, val in zip(names, corr_vals):
            c = self.lattice[corr]
            if "DHG" in corr:
                c.vary_strength(val / 10, 0.0)
            elif ("DHJ" in corr) or ("SDV" in corr):
                c.vary_strength(0.0, val / 10)
        self.__track_bunch()

    def _build_bunch_from_guesses(self, emit_x, emit_y, beta_x0, beta_y0, alpha_x0, alpha_y0):
        T = rft.Bunch6d_twiss()
        T.emitt_x = float(emit_x)
        T.emitt_y = float(emit_y)
        T.beta_x = float(beta_x0)
        T.beta_y = float(beta_y0)
        T.alpha_x = float(alpha_x0)
        T.alpha_y = float(alpha_y0)
        T.sigma_t = 0#10*rft.ps # mm/c
        T.sigma_pt = 0#10 # permille
        T.mean_xp = 0.0
        T.mean_yp = 0.0
        return rft.Bunch6d_QR(rft.electronmass, self.population, self.Q, self.Pref, T, self.nparticles, self.sigmaCut)

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

    def predict_emittance_scan_response(self, quad_name, screens, K1L_values, emit_x, emit_y, beta_x0, beta_y0, alpha_x0, alpha_y0, stop_checker = None, reference_screen = None):
        screens = list(screens)
        K1L_values = np.asarray(K1L_values, dtype=float)
        screens = list(screens)
        if reference_screen is None:
            reference_screen = screens[0]
        start_element_name = str(quad_name)
        end_element_name = str(screens[-1])
        original_quads = self.get_quadrupoles(names=[quad_name])
        if len(original_quads["bdes"]) == 0:
            raise RuntimeError(f"Could not find original strength for quad {quad_name}")
        K1L_original = float(original_quads["bdes"][0])

        B0_original = self.B0
        output_x = np.full((len(K1L_values),len(screens)), np.nan, dtype=float)
        output_y = np.full((len(K1L_values),len(screens)), np.nan, dtype=float)

        try:
            for k,K1L in enumerate(K1L_values):
                if callable(stop_checker) and stop_checker():
                    raise RuntimeError("__OPTIMIZATION_STOP__")
                self.set_quadrupoles([quad_name], [float(K1L)], track=False)
                start_element = self.lattice[start_element_name]
                if isinstance(start_element, list):
                    start_element = start_element[0]

                end_element = self.lattice[end_element_name]
                if isinstance(end_element, list):
                    end_element = end_element[-1]

                temp_bunch = self._build_bunch_from_guesses(
                    emit_x=float(emit_x), emit_y=float(emit_y),
                    beta_x0=float(beta_x0), beta_y0=float(beta_y0),
                    alpha_x0=float(alpha_x0), alpha_y0=float(alpha_y0),
                )

                lattice_view = rft.Lattice_view(self.lattice, start_element, end_element)
                tracked_to_last_screen = lattice_view.track(temp_bunch)

                for si, screen_name in enumerate(screens):
                    screen_elem = self.lattice[screen_name]
                    if isinstance(screen_elem, list):
                        screen_elem = screen_elem[-1]
                    bunch_at_screen = None
                    try:
                        bunch_at_screen = screen_elem.get_bunch()
                    except Exception:
                        bunch_at_screen = None
                    if bunch_at_screen is None and str(screen_name) == end_element_name:
                        bunch_at_screen = tracked_to_last_screen
                    if bunch_at_screen is None:
                        continue
                    m = bunch_at_screen.get_phase_space('%x %y')
                    if m is not None and len(m) > 0:
                        output_x[k, si] = float(np.std(m[:, 0]))
                        output_y[k, si] = float(np.std(m[:, 1]))

        finally:
            self.set_quadrupoles([quad_name], [float(K1L_original)], track=False)
            self.B0 = B0_original
            self.__track_bunch()

        return output_x, output_y

    def _predict_scan_response_full(self, quad_name, screens, K1L_values, emit_x, emit_y, beta_x0, beta_y0, alpha_x0, alpha_y0, quad_dx0=None, quad_dy0=None, quad_roll=None, stop_checker=None, reference_screen=None):
        screens = list(screens)
        K1L_values = np.asarray(K1L_values, dtype=float)
        if reference_screen is None: reference_screen = screens[0]
        end_element_name = str(screens[-1])
        start_element_name = str(quad_name)
        B0_original = self.B0
        lattice_reference = self.lattice
        self.lattice = lattice_reference.clone()
        start_element = self.lattice[start_element_name]
        if isinstance(start_element, list):
            start_element = start_element[0]
        override_offsets = any(value is not None for value in (quad_dx0, quad_dy0, quad_roll))
        dx = 0.0 if quad_dx0 is None else float(quad_dx0)
        dy = 0.0 if quad_dy0 is None else float(quad_dy0)
        roll = 0.0 if quad_roll is None else float(quad_roll)
        nK1L, nscreens = len(K1L_values), len(screens)
        sigma_x = np.full((nK1L, nscreens), np.nan, dtype=float)
        sigma_y = np.full((nK1L, nscreens), np.nan, dtype=float)
        x_mean = np.full((nK1L, nscreens), np.nan, dtype=float)
        y_mean = np.full((nK1L, nscreens), np.nan, dtype=float)
        sigma_xy = np.full((nK1L, nscreens), np.nan, dtype=float)

        try:
            if override_offsets:
                start_element.set_offsets(dx, dy, 0.0, roll, 0.0, 0.0, "center")
            for k,K1L in enumerate(K1L_values):
                if callable(stop_checker) and stop_checker():
                    raise RuntimeError("__OPTIMIZATION_STOP__")
                self.set_quadrupoles([quad_name], [float(K1L)], track=False)

                end_element = self.lattice[end_element_name]
                if isinstance(end_element, list):
                    end_element = end_element[-1]

                temp_bunch = self._build_bunch_from_guesses(emit_x=float(emit_x), emit_y=float(emit_y), beta_x0=float(beta_x0), beta_y0=float(beta_y0), alpha_x0=float(alpha_x0), alpha_y0=float(alpha_y0))
                lattice_view = rft.Lattice_view(self.lattice, start_element, end_element)
                tracked_to_last_screen = lattice_view.track(temp_bunch)

                for si, screen_name in enumerate(screens):
                    screen_elem = self.lattice[screen_name]
                    if isinstance(screen_elem, list):
                        screen_elem = screen_elem[-1]
                    bunch_at_screen = screen_elem.get_bunch()
                    if bunch_at_screen is None and str(screen_name) == end_element_name:
                        bunch_at_screen = tracked_to_last_screen
                    m = bunch_at_screen.get_phase_space('%x %y')
                    if m is not None and len(m) > 0:
                        xs = m[:, 0]
                        ys = m[:, 1]
                        xm = float(np.mean(xs))
                        ym = float(np.mean(ys))
                        sigma_x[k, si] = float(np.std(xs))
                        sigma_y[k, si] = float(np.std(ys))
                        x_mean[k, si] = xm
                        y_mean[k, si] = ym
                        sigma_xy[k, si] = float(np.mean((xs - xm) * (ys - ym))) # covariance

        finally:
            self.lattice = lattice_reference
            self.B0 = B0_original
            self.__track_bunch()

        return {
            "sigma_x": sigma_x, "sigma_y": sigma_y,
            "x_mean": x_mean, "y_mean": y_mean,
            "sigma_xy": sigma_xy,
        }

    def predict_emittance_scan_response(self, quad_name, screens, K1L_values, emit_x, emit_y, beta_x0, beta_y0, alpha_x0, alpha_y0, stop_checker = None, reference_screen = None):
        full = self._predict_scan_response_full(quad_name, screens, K1L_values, emit_x, emit_y, beta_x0, beta_y0, alpha_x0, alpha_y0, stop_checker=stop_checker, reference_screen=reference_screen)
        return full["sigma_x"], full["sigma_y"]

    def predict_emittance_scan_response_full(self, quad_name, screens, K1L_values, emit_x, emit_y, beta_x0, beta_y0, alpha_x0, alpha_y0, quad_dx0=None, quad_dy0=None, quad_roll=None, stop_checker=None, reference_screen=None):
        return self._predict_scan_response_full(quad_name, screens, K1L_values, emit_x, emit_y, beta_x0, beta_y0, alpha_x0, alpha_y0, quad_dx0=quad_dx0, quad_dy0=quad_dy0, quad_roll=quad_roll, stop_checker=stop_checker, reference_screen=reference_screen)


    def get_phase_space_transport_to_screens(self, reference_screen=None, screens=None):
        if screens is None:
            screens = list(self.screens)
        if isinstance(screens, str):
            screens = [screens]
        screens = list(screens)

        if reference_screen is None:
            reference_screen = screens[0]

        original_bunch = self.B0

        result = {
            "reference_screen": str(reference_screen),
            "screens": [str(s) for s in screens],
            "x": {"R11": [], "R12": [], "R21": [], "R22": []},
            "y": {"R33": [], "R34": [], "R43": [], "R44": []},
        }

        try:
            start_element = self.lattice[reference_screen]
            if isinstance(start_element, list):
                start_element = start_element[-1]

            end_element = self.lattice[screens[-1]]
            if isinstance(end_element, list):
                end_element = end_element[-1]

            bx = np.array([
                [1.0, 0.0, 0.0, 0.0, 0.0, self.Pref],
                [0.0, 1.0, 0.0, 0.0, 0.0, self.Pref],
            ], dtype=float)

            bunch_x = rft.Bunch6d(rft.electronmass, 0.0, self.Q, bx)
            lattice_view = rft.Lattice_view(self.lattice, start_element, end_element)
            tracked_x = lattice_view.track(bunch_x)

            for screen in screens:
                screen_element = self.lattice[screen]
                if isinstance(screen_element, list):
                    screen_element = screen_element[-1]

                b = screen_element.get_bunch()
                if b is None and str(screen) == str(screens[-1]):
                    b = tracked_x

                ps = np.asarray(b.get_phase_space("%x %xp"), dtype=float)

                result["x"]["R11"].append(float(ps[0, 0]))
                result["x"]["R12"].append(float(ps[1, 0]))
                result["x"]["R21"].append(float(ps[0, 1]))
                result["x"]["R22"].append(float(ps[1, 1]))

            by = np.array([
                [0.0, 0.0, 1.0, 0.0, 0.0, self.Pref],
                [0.0, 0.0, 0.0, 1.0, 0.0, self.Pref],
            ], dtype=float)

            bunch_y = rft.Bunch6d(rft.electronmass, 0.0, self.Q, by)
            lattice_view = rft.Lattice_view(self.lattice, start_element, end_element)
            tracked_y = lattice_view.track(bunch_y)

            for screen in screens:
                screen_element = self.lattice[screen]
                if isinstance(screen_element, list):
                    screen_element = screen_element[-1]

                b = screen_element.get_bunch()
                if b is None and str(screen) == str(screens[-1]):
                    b = tracked_y

                ps = np.asarray(b.get_phase_space("%y %yp"), dtype=float)

                # [
                #     [x_of_particle_0, xp_of_particle_0],
                #     [x_of_particle_1, xp_of_particle_1],
                # ]

                result["y"]["R33"].append(float(ps[0, 0]))
                result["y"]["R34"].append(float(ps[1, 0]))
                result["y"]["R43"].append(float(ps[0, 1]))
                result["y"]["R44"].append(float(ps[1, 1]))

        finally:
            self.B0 = original_bunch
            self.__track_bunch()

        return result

    def get_R_matrix_scan(self, quad_name, screens, K1L_values):
        screens = list(screens)
        K1L_values = np.asarray(K1L_values, dtype=float)
        original_quads = self.get_quadrupoles(names=[quad_name])
        K1L_original = float(original_quads["bdes"][0])
        end_element_name = str(screens[-1])

        n_k1l = len(K1L_values)
        n_screens = len(screens)
        R11 = np.full((n_k1l, n_screens), np.nan, dtype=float)
        R12 = np.full((n_k1l, n_screens), np.nan, dtype=float)
        R33 = np.full((n_k1l, n_screens), np.nan, dtype=float)
        R34 = np.full((n_k1l, n_screens), np.nan, dtype=float)

        original_bunch = self.B0
        try:
            for k, K1L in enumerate(K1L_values):
                self.set_quadrupoles([quad_name], [float(K1L)], track=False)
                start_element = self.lattice[quad_name]
                if isinstance(start_element, list): start_element = start_element[0]
                end_element = self.lattice[end_element_name]
                if isinstance(end_element, list): end_element = end_element[-1]

                bx = np.array([
                    [1.0, 0.0, 0.0, 0.0, 0.0, self.Pref],  # x = 1, x' = 0
                    [0.0, 1.0, 0.0, 0.0, 0.0, self.Pref],  # x = 0, x' = 1
                ], dtype=float)

                bunch_x = rft.Bunch6d(rft.electronmass, 0.0, self.Q, bx)
                tracked_x = rft.Lattice_view(self.lattice, start_element, end_element).track(bunch_x)

                for si, screen_name in enumerate(screens):
                    screen_elem = self.lattice[screen_name]
                    if isinstance(screen_elem, list):
                        screen_elem = screen_elem[-1]
                    b_x = None
                    try:
                        b_x = screen_elem.get_bunch()  # Read the test particles cached at this screen.
                    except Exception:
                        b_x = None
                    if b_x is None and str(screen_name) == end_element_name:
                        b_x = tracked_x
                    if b_x is None:
                        continue
                    phase_space_x = np.asarray(b_x.get_phase_space("%x %xp"), dtype=float)
                    R11[k, si] = phase_space_x[0, 0]  # k: K1L index, si: screen index.
                    R12[k, si] = phase_space_x[1, 0]

                by = np.array([
                    [0.0, 0.0, 1.0, 0.0, 0.0, self.Pref],  # y = 1, y' = 0
                    [0.0, 0.0, 0.0, 1.0, 0.0, self.Pref],  # y = 0, y' = 1
                ], dtype=float)
                bunch_y = rft.Bunch6d(rft.electronmass, 0.0, self.Q, by)
                tracked_y = rft.Lattice_view(self.lattice, start_element, end_element).track(bunch_y)

                for si, screen_name in enumerate(screens):
                    screen_elem = self.lattice[screen_name]
                    if isinstance(screen_elem, list):
                        screen_elem = screen_elem[-1]
                    b_y = None
                    try:
                        b_y = screen_elem.get_bunch()  # Read the test particles cached at this screen.
                    except Exception:
                        b_y = None
                    if b_y is None and str(screen_name) == end_element_name:
                        b_y = tracked_y
                    if b_y is None:
                        continue
                    phase_space_y = np.asarray(b_y.get_phase_space("%y %yp"), dtype=float)
                    R33[k, si] = phase_space_y[0, 0]  # k: K1L index, si: screen index.
                    R34[k, si] = phase_space_y[1, 0]
        finally:
            self.set_quadrupoles([quad_name], [float(K1L_original)], track=False)
            self.B0 = original_bunch
            self.__track_bunch()

        return {"R11": R11, "R12": R12, "R33": R33, "R34": R34, "screens": screens, "K1L_values": K1L_values}

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
                    strength = element.get_K1L(self.Pref / self.Q)
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

    def _get_elements_positions_show_beamline(self, names=None):
        if isinstance(names, str):
            names = [names]
        all_names = [name for name in self.sequence if names is None or name in names]

        return {
            "names": all_names,
            "S": np.array([self.lattice[name].get_S("entrance") for name in all_names], dtype=float),
        }

    def _give_elements_to_show_beamline(self, quad_selected):
        start_quad_element_name = quad_selected
        return start_quad_element_name

    def set_quadrupoles(self, names, values_range, track=True):
        if isinstance(names, str):
            names = [names]
        if not (isinstance(values_range, (list, tuple, np.ndarray))):
            values_range = [values_range]
        for quadrupole_name, value in zip(names, values_range):
            elements = self.lattice[quadrupole_name]
            if not isinstance(elements, (list)): elements = [elements]
            for element in elements:
                element.set_K1L(self.Pref / self.Q,float(value))
        if track:
            self.__track_bunch()

    def get_elements_indices(self, names):
        if isinstance(names, str):
            names = [names]
        name_to_index = {string: index for index, string in enumerate(self.sequence)}
        return [name_to_index.get(name, np.nan) for name in names]
