import RF_Track as rft
import numpy as np
import time, os, re, copy
from Backend.LogConsole import LogConsole
from datetime import datetime
from Interfaces.AbstractMachineInterface import AbstractMachineInterface

class InterfaceATF2_Linac_RFTrack(AbstractMachineInterface):

    def get_name(self):
        return 'ATF2_Linac_RFT'

    def __init__(self, population=1e+10, jitter=0.0, bpm_resolution=0.0, nsamples=1, nparticles=1000):
        super().__init__()
        self.log = print
        self.twiss_path = os.path.join(os.path.dirname(__file__), 'Linac_ATF2', 'linacend.tws')
        self.madx_path = os.path.join(os.path.dirname(__file__), 'Linac_ATF2', 'End_linac.madx')
        self.Pref0 = 80
        self.Pref = self.Pref0
        self.nparticles = nparticles
        self.electronmass = rft.electronmass
        self.population = population
        self.jitter = jitter
        self.nsamples = nsamples
        self.Q = -1
        self.dfs_test_energy = 0.98
        self.wfs_test_charge = 0.90
        self._beam_mode = "nominal"
        self.freq = 2855.9822615999997e+6 # Hz
        self.lattice = self._build_lattice()
        self.lattice.set_bpm_resolution(bpm_resolution)
        self.sequence = [e.get_name() for e in self.lattice['*']]
        self.bpms = [e.get_name() for e in self.lattice.get_bpms()]
        self.corrs = [e.get_name() for e in self.lattice.get_correctors()]
        self.quadrupoles = list(dict.fromkeys(e.get_name() for e in self.lattice.get_quadrupoles()))
        self.__setup_beam0()
        self.__track_bunch()
        self.machine_name = "ATF2"
        self.lattice.align_elements()
        '''Uncomment lines below to scatter elements in the lattice.'''
        # self.lattice.scatter_elements('bpm', 0.100, 0.100, 0, 0, 0, 0, 'center')
        # self.lattice.scatter_elements('quadrupole', 0.100, 0.100, 0, 0, 0, 0, 'center')


    def get_beam_factors(self):
        gamma_rel = np.sqrt((self.Pref0 / self.electronmass) ** 2 + 1.0)
        beta_rel = np.sqrt(1.0 - 1.0 / gamma_rel ** 2)
        return gamma_rel, beta_rel

    def log_messages(self, console):
        self.log = console or print

    def __setup_beam0(self):
        T = rft.Bunch6d_twiss()
        T.emitt_x = 24.7e-4  # mm.mrad normalised emittance
        T.emitt_y = 24.7e-4  # mm.mrad
        T.beta_x = 0.451  # m
        T.beta_y = 2.236  # m
        T.alpha_x = -3.599
        T.alpha_y = -10.419
        T.sigma_t = 8  # mm/c
        T.sigma_pt = 0.8  # permille
        self.B0 = rft.Bunch6d_QR(rft.electronmass, self.population, self.Q, self.Pref, T, self.nparticles)

    def __setup_beam1(self):
        # Beam for DFS - Reduced energy
        Pref = self.dfs_test_energy * self.Pref0
        T = rft.Bunch6d_twiss()
        T.emitt_x = 24.7e-4  # mm.mrad normalised emittance
        T.emitt_y = 24.7e-4  # mm.mrad
        T.beta_x = 0.451  # m
        T.beta_y = 2.236  # m
        T.alpha_x = -3.599
        T.alpha_y = -10.419
        T.sigma_t = 8  # mm/c
        T.sigma_pt = 0.8  # permille
        Nparticles = 1000  # number of macroparticles
        self.B0 = rft.Bunch6d_QR(rft.electronmass, self.population, self.Q, Pref, T, self.nparticles)

    def __setup_beam2(self):
        # Beam for WFS - Reduced bunch charge
        population = self.wfs_test_charge * self.population
        T = rft.Bunch6d_twiss()
        T.emitt_x = 24.7e-4  # mm.mrad normalised emittance
        T.emitt_y = 24.7e-4  # mm.mrad
        T.beta_x = 0.451  # m
        T.beta_y = 2.236  # m
        T.alpha_x = -3.599
        T.alpha_y = -10.419
        T.sigma_t = 8  # mm/c
        T.sigma_pt = 0.8  # permille
        Nparticles = 1000  # number of macroparticles
        self.B0 = rft.Bunch6d_QR(rft.electronmass, population, self.Q, self.Pref, T, self.nparticles)

    def __track_bunch(self):
        I0 = self.B0.get_info()
        dx = self.jitter * I0.sigma_x
        dy = self.jitter * I0.sigma_y
        dz, dt, roll = 0.0, 0.0, 0.0
        pitch = self.jitter * I0.sigma_py
        yaw = self.jitter * I0.sigma_px
        B0_offset = self.B0.displaced(dx, dy, dz, dt, roll, pitch, yaw)
        self.B1 = self.lattice.track(B0_offset)
        I = self.B1.get_info()

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

    def change_intensity(self):  # reduced charge
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
        return [string for string in self.corrs if
                (string.lower().startswith('zh')) or (string.lower().startswith('zx'))]

    def get_vcorrectors_names(self):
        return [string for string in self.corrs if string.lower().startswith('zv')]

    def _set_name(self, element, name):
        try:
            element.set_name(name)
        except Exception:
            pass
        return element

    def _make_tw(self, name):
        a0 = 25.5e6 # V/m, principal Fourier coefficient
        ph_adv = 2 * np.pi / 3  # radian, phase advance per cell
        n_cells = -85  # number of cells, negative sign indicates a start from the beginning of the cell
        TW = rft.TW_Structure(a0, 0, self.freq, ph_adv, n_cells)  # 0=el primer elemento del armonico
        phid = -91.67 + 90  # LAG -.2546388888888886 [2*pi]
        TW.set_phid(phid)
        self._set_name(TW, name)
        return TW

    def _make_twb1(self, name):
        BCph_adv = 2.0908741;
        BCBn_cells = -85;
        TWB1 = rft.TW_Structure(0, 0, self.freq, BCph_adv, BCBn_cells)  # TW Bunching
        BCphid = -90 + 90
        TWB1.set_phid(BCphid)
        self._set_name(TWB1, name)
        return TWB1

    def _make_twb2(self, name):
        BC2ph_adv = 2.0958035;
        BC2Bn_cells = -85;
        TWB2 = rft.TW_Structure(0, 0, self.freq, BC2ph_adv, BC2Bn_cells)  # TW Bunching 2
        BCphid = -90 + 90
        TWB2.set_phid(BCphid)
        self._set_name(TWB2, name)
        return TWB2

    def _append_drift(self, lattice, length, name = None):
        if length <= 1e-12:
            return
        Drift = rft.Drift(float(length))
        if name is not None:
            self._set_name(Drift, name)
        lattice.append(Drift)

    def _append_main_linac(self, lattice):
        charge = -1
        P_Q = self.Pref0 / charge
        pos_actual = 0.0

        quad_specifics = [
            # name, LQ, pos_quad, strength
            ('QA1L', 0.132, 1.356, 1.1835109401030304),
            ('QA2L', 0.182, 1.703, -9.264115425347802),
            ('QA3L', 0.182, 2.128, 7.601904693021979),
            ('QA4L', 0.182, 4.774, -8.309547432909891),
            ('QA5L', 0.132, 5.171, 10.368563413319697),
            ('QD1L', 0.132, 9.548499999999997, -6.394111798396212),
            ('QF1L', 0.182, 9.948499999999996, 8.468997394013186),
            ('QD2L', 0.132, 10.348499999999994, -6.415397177518939),
            ('QF2L', 0.125, 14.143499999999994, 4.6362775786152),
            ('QD3L', 0.125, 14.894999999999994, -4.7315774972976),
            ('QF3L', 0.125, 18.986999999999995, 5.8107864654864),
            ('QD4L', 0.125, 19.73699999999999, -5.1676165334704),
            ('QF4L', 0.125, 23.84099999999999, 7.6235715827536),
            ('QD5L', 0.125, 24.491999999999987, -6.7358412875152),
            ('QF5L', 0.125, 28.486499999999985, 7.6898162668192),
            ('QD6L', 0.125, 29.13649999999998, -6.8687757100048),
            ('QF6L', 0.125, 36.74699999999999, 4.4820119116136),
            ('QD7L', 0.125, 37.397, -4.4944090948816),
            ('QS5L', 0.125, 45.298760000000016, 0.9884761807816),
            ('QS6L', 0.125, 53.137260000000026, -1.1188273572312),
            ('QS7L', 0.125, 60.66226000000003, 0.8013307329232),
            ('QS8L', 0.125, 71.47405000000003, -0.522816766032),
            ('QM1L', 0.175, 78.59605000000002, 1.4909714285714286),
            ('QM2L', 0.325, 79.61655000000002, -1.894584615384615),
            ('QM3L', 0.175, 81.59355000000002, 3.646114285714286),
        ]

        cavity_specifics = [
            ('CA1L', 7.3835, 'TW'),
            ('CA2L', 12.142714999999995, 'TW'),
            ('CA3L', 17.006999999999994, 'TW'),
            ('CA4L', 21.87099999999999, 'TW'),
            ('CB1L', 26.520499999999988, 'TWB1'),
            ('CA5L', 31.180499999999984, 'TW'),
            ('CA6L', 34.78399999999999, 'TW'),
            ('CA7L', 39.43000000000001, 'TW'),
            ('CA8L', 43.03397000000001, 'TW'),
            ('CA9L', 47.26726000000002, 'TW'),
            ('CA10L', 50.872260000000026, 'TW'),
            ('CA11L', 55.10826000000003, 'TW'),
            ('CA12L', 58.397760000000034, 'TW'),
            ('CB2L', 62.63276000000003, 'TWB2'),
            ('CA13L', 65.91976000000003, 'TW'),
            ('CA14L', 69.20926000000003, 'TW'),
            ('CA15L', 73.44305000000003, 'TW'),
            ('CA16L', 76.73305000000002, 'TW'),
        ]

        corrector_specifics = [
            ("ZH1L", 1.356), ("ZV1L", 1.356),
            ("ZH2L", 1.703), ("ZV2L", 1.703),
            ("ZH3L", 2.128), ("ZV3L", 2.128),
            ("ZH4L", 4.774), ("ZV4L", 4.774),
            ("ZH5L", 5.171), ("ZV5L", 5.171),
        ]

        actions = []

        for name, length, pos, k1 in quad_specifics:
            actions.append(
                {
                    'element_type': 'quadrupole',
                    'position': pos,
                    'name': name,
                    'length': length,
                    'k1': k1,
                }
            )
        tw = self._make_tw('tw')
        twb1 = self._make_twb1('twb1')
        twb2 = self._make_twb2('twb2')

        cavity_length = {
            'TW' : tw.get_length(),
            'TWB1' : twb1.get_length(),
            'TWB2' : twb2.get_length(),
        }

        for name, pos, ctype in cavity_specifics:
            actions.append(
                {
                    'element_type': 'cavity',
                    'position': pos,
                    'name': name,
                    'ctype': ctype,
                    'length': cavity_length[ctype],
                }
            )

        for name, pos in corrector_specifics:
            actions.append({
                "element_type": "corrector",
                "position": pos,
                "name": name,
                "length": 0.0,
            })

        actions.sort(key = lambda item: item['position'])

        drift_id = 1

        for item in actions:
            drift = item['position'] - item['length'] / 2.0 - pos_actual
            self._append_drift(lattice, drift, f'DRIFT_LINAC_{drift_id:03d}')
            drift_id += 1

            if item['element_type'] == 'quadrupole':
                quadrupole = rft.Quadrupole(item['length'], P_Q, item['k1'])
                self._set_name(quadrupole, item['name'])
                lattice.append(quadrupole)

            elif item['element_type'] == 'corrector':
                corrector = rft.Corrector()
                self._set_name(corrector, item['name'])
                lattice.append(corrector)

            else:
                if item['ctype'] == 'TW':
                    cavity = self._make_tw(item['name'])
                elif item['ctype'] == 'TWB1':
                    cavity = self._make_twb1(item['name'])
                else:
                    cavity = self._make_twb2(item['name'])
                lattice.append(cavity)
                try:
                    lattice.unset_t0()
                except Exception:
                    pass

            pos_actual = item['position'] + item['length'] / 2.0
        pos_IPZL = 83.40465000000002
        self._append_drift(lattice, pos_IPZL - pos_actual, 'DRIFT_TO_IPZL')
        ipzl = rft.Drift(0.0)
        self._set_name(ipzl, 'IPZL')
        lattice.append(ipzl)

    def _get_endlinac_twiss(self):
        try:
            from cpymad.madx import Madx
            madx = Madx(stdout=True)
            madx.call(self.madx_path);
        except Exception as e:
            self.log("Could not get endlinac twiss data")

    def _append_endlinac(self, lattice):
        self._get_endlinac_twiss()
        if not os.path.isfile(self.twiss_path):
            raise FileNotFoundError("Could not find endlinac twiss file")
        lattice.append(rft.Lattice(self.twiss_path))

    def _adjust_strength_with_autophase(self, lattice):
        mass = rft.electronmass
        charge = self.Q
        population = self.population
        Pi = self.Pref0
        P_Q = Pi / charge
        P0 = rft.Bunch6d(mass, population, charge, np.array([0, 0, 0, 0, 0, Pi]).T)

        TWL = rft.Lattice()
        TWL.append(self._make_tw('AUTOPHASE_TW'))
        try:
            TWL.unset_t0()
        except Exception:
            pass

        try:
            Pmax = TWL.autophase(P0)
            p_gain = Pmax - Pi
        except Exception as e:
            self.log(f"Autophase failed, keeping nominal main linac quadrupoles: {e}")
            return

        groups = [
            ['QA1L', 'QA2L', 'QA3L', 'QA4L', 'QA5L'],
            ['QD1L', 'QF1L', 'QD2L'],
            ['QF2L', 'QD3L'],
            ['QF3L', 'QD4L'],
            ['QF4L', 'QD5L', 'QF5L', 'QD6L'],
            ['QF6L', 'QD7L'],
            ['QS5L'],
            ['QS6L'],
            ['QS7L'],
            ['QS8L'],
            ['QM1L', 'QM2L', 'QM3L'],
        ]
        increments_after_group = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0]

        p = Pi
        for names, inc in zip(groups, increments_after_group):
            for name in names:
                try:
                    element = lattice[name]
                    current_K1 = element.get_K1(P_Q)
                    if isinstance(current_K1, (list, tuple, np.ndarray)):
                        current_K1 = float(current_K1[0])
                    else:
                        current_K1 = float(current_K1)
                    element.set_K1(p / charge, current_K1)
                except Exception as e:
                    self.log(f"Could not set K1 for {name}: {e}")
            p += inc * p_gain

    def _build_lattice(self):
        lattice = rft.Lattice()
        self._append_main_linac(lattice)
        self._append_endlinac(lattice)
        self._adjust_strength_with_autophase(lattice)
        return lattice

    def get_target_dispersion(self, names=None): # for DR too
        if names is None:
            names = self.bpms
        if isinstance(names, str):
            names = [names]
        twiss_path = os.path.join(os.path.dirname(__file__), 'Linac_ATF2', 'linacend.tws')
        with open(twiss_path, "r") as file:
            lines = [line.strip() for line in file if line.strip()]

        star_symbol = next(i for i, line in enumerate(lines) if line.startswith("*"))
        dollar_sign = next(i for i, line in enumerate(lines) if line.startswith("$") and i > star_symbol)
        columns = lines[star_symbol].lstrip("*").split()
        try:
            DX_column = columns.index("DX")
            DY_column = columns.index("DY")
            elements_names = columns.index("NAME")
        except ValueError as e:
            raise RuntimeError("There are no such columns in the twiss file")
        disp_values = {}
        target_disp_x, target_disp_y = [], []
        for line in lines[dollar_sign + 1:]:
            data = line.split()
            if len(data) <= max(DX_column, DY_column,
                                elements_names):  # if a line has less column than needed, it is omitted
                continue
            bpms_name = data[elements_names].strip('"')
            try:
                disp_values[bpms_name] = (float(data[DX_column]), float(data[DY_column]))
            except ValueError:
                continue
        for bpm in names:
            if bpm in disp_values:
                dx, dy = disp_values[bpm]
            else:
                dx, dy = float("nan"), float("nan")
            target_disp_x.append(dx)
            target_disp_y.append(dy)
        return target_disp_x, target_disp_y

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
            idx = [i for i, s in enumerate(bpms["names"]) if s in names]
            bpms = {
                "names": [bpms["names"][i] for i in idx],
                "x": np.asarray(bpms["x"])[:, idx],
                "y": np.asarray(bpms["y"])[:, idx],
                "tmit": np.asarray(bpms["tmit"])[:, idx],
            }

        return bpms

    def get_quadrupoles(self, names=None):
        self.log("Reading quadrupoles' strengths...")
        bdes = np.zeros(len(self.quadrupoles), dtype=float)

        for i, quadrupole_name in enumerate(self.quadrupoles):
            elements = self.lattice[quadrupole_name]
            if not isinstance(elements, list):
                elements = [elements]

            k1l_values = []
            for element in elements:
                try:
                    strength = element.get_K1(self.Pref0 / self.Q)
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
                element.set_K1(self.Pref0 / self.Q,float(value))
        if track:
            self.__track_bunch()

    def set_correctors(self, names, corr_vals):
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        for corr, val in zip(names, corr_vals):
            if corr[:2] == "ZH" or corr[:2] == "ZX":
                self.lattice[corr].set_strength(val / 10, 0.0)  # T*mm
            elif corr[:2] == "ZV":
                self.lattice[corr].set_strength(0.0, val / 10)  # T*mm
        self.__track_bunch()

    def vary_correctors(self, names, corr_vals):
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        for corr, val in zip(names, corr_vals):
            if corr[:2] == "ZH" or corr[:2] == "ZX":
                self.lattice[corr].vary_strength(val / 10, 0.0)  # T*mm
            elif corr[:2] == "ZV":
                self.lattice[corr].vary_strength(0.0, val / 10)  # T*mm
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

    def _get_optics_from_twiss_file(self,names=None):

        quad_specifics = [
            # name, LQ, pos_quad, strength
            ('QA1L', 0.132, 1.356, 1.1835109401030304),
            ('QA2L', 0.182, 1.703, -9.264115425347802),
            ('QA3L', 0.182, 2.128, 7.601904693021979),
            ('QA4L', 0.182, 4.774, -8.309547432909891),
            ('QA5L', 0.132, 5.171, 10.368563413319697),
            ('QD1L', 0.132, 9.548499999999997, -6.394111798396212),
            ('QF1L', 0.182, 9.948499999999996, 8.468997394013186),
            ('QD2L', 0.132, 10.348499999999994, -6.415397177518939),
            ('QF2L', 0.125, 14.143499999999994, 4.6362775786152),
            ('QD3L', 0.125, 14.894999999999994, -4.7315774972976),
            ('QF3L', 0.125, 18.986999999999995, 5.8107864654864),
            ('QD4L', 0.125, 19.73699999999999, -5.1676165334704),
            ('QF4L', 0.125, 23.84099999999999, 7.6235715827536),
            ('QD5L', 0.125, 24.491999999999987, -6.7358412875152),
            ('QF5L', 0.125, 28.486499999999985, 7.6898162668192),
            ('QD6L', 0.125, 29.13649999999998, -6.8687757100048),
            ('QF6L', 0.125, 36.74699999999999, 4.4820119116136),
            ('QD7L', 0.125, 37.397, -4.4944090948816),
            ('QS5L', 0.125, 45.298760000000016, 0.9884761807816),
            ('QS6L', 0.125, 53.137260000000026, -1.1188273572312),
            ('QS7L', 0.125, 60.66226000000003, 0.8013307329232),
            ('QS8L', 0.125, 71.47405000000003, -0.522816766032),
            ('QM1L', 0.175, 78.59605000000002, 1.4909714285714286),
            ('QM2L', 0.325, 79.61655000000002, -1.894584615384615),
            ('QM3L', 0.175, 81.59355000000002, 3.646114285714286),
        ]

        cavity_specifics = [
            ('CA1L', 7.3835, 'TW'),
            ('CA2L', 12.142714999999995, 'TW'),
            ('CA3L', 17.006999999999994, 'TW'),
            ('CA4L', 21.87099999999999, 'TW'),
            ('CB1L', 26.520499999999988, 'TWB1'),
            ('CA5L', 31.180499999999984, 'TW'),
            ('CA6L', 34.78399999999999, 'TW'),
            ('CA7L', 39.43000000000001, 'TW'),
            ('CA8L', 43.03397000000001, 'TW'),
            ('CA9L', 47.26726000000002, 'TW'),
            ('CA10L', 50.872260000000026, 'TW'),
            ('CA11L', 55.10826000000003, 'TW'),
            ('CA12L', 58.397760000000034, 'TW'),
            ('CB2L', 62.63276000000003, 'TWB2'),
            ('CA13L', 65.91976000000003, 'TW'),
            ('CA14L', 69.20926000000003, 'TW'),
            ('CA15L', 73.44305000000003, 'TW'),
            ('CA16L', 76.73305000000002, 'TW'),
        ]

        tw = self._make_tw('tw')
        twb1 = self._make_tw('twb1')
        twb2 = self._make_tw('twb2')
        cavity_length = {
            'TW': float(tw.get_length()),
            'TWB1': float(twb1.get_length()),
            'TWB2': float(twb2.get_length()),
        }
        combined = {}

        for name, length, pos, k1 in quad_specifics:
            combined[name] = {
                'S': float(pos),
                'betx': np.nan,
                'alfx': np.nan,
                'bety': np.nan,
                'alfy': np.nan,
                'mux': np.nan,
                'muy': np.nan,
                'L': float(length),
            }

        for name, pos, ctype in cavity_specifics:
            combined[name] = {
                'S': float(pos),
                'betx': np.nan,
                'alfx': np.nan,
                'bety': np.nan,
                'alfy': np.nan,
                'mux': np.nan,
                'muy': np.nan,
                'L': float(cavity_length[ctype]),
            }

        combined['IPZL'] = {
            'S': 83.40465000000002,
            'betx': np.nan,
            'alfx': np.nan,
            'bety': np.nan,
            'alfy': np.nan,
            'mux': np.nan,
            'muy': np.nan,
            'L': 0.0,
        }

        with open(self.twiss_path, "r") as file:
            lines = [line.strip() for line in file if line.strip()]

        star_symbol = next(i for i, line in enumerate(lines) if line.startswith("*"))
        dollar_sign = next(i for i, line in enumerate(lines) if line.startswith("$") and i > star_symbol)
        columns = lines[star_symbol].lstrip("*").split()

        cols=["NAME", "S" , "BETX" , "ALFX" , "BETY" , "ALFY" , "MUX" , "MUY" , "L"]
        index={}
        for col in cols:
            try:
                index[col]=columns.index(col)
            except ValueError:
                raise RuntimeError(f"Column {col} not found in twiss file")
        duplicated={}
        for line in lines[dollar_sign + 1:]:
            data = line.split()
            if len(data) <= max(index.values()):  # if a line has less column than needed, it is omitted
                continue
            name = data[index["NAME"]].strip('"')

            if name not in duplicated:
                duplicated[name]={col:[] for col in cols[1:]}
            for col in cols[1:]:
                try:
                    duplicated[name][col].append(float(data[index[col]]))
                except ValueError:
                    duplicated[name][col].append(float("nan"))
        endlinac_s_offset = 83.40465000000002

        for name, vals in duplicated.items():
            local_s = float(vals["S"][0]) if len(vals["S"]) > 0 else np.nan
            global_s = local_s + endlinac_s_offset if np.isfinite(local_s) else np.nan

            combined[name] = {
                'S': global_s,
                'betx': float(vals["BETX"][0]) if len(vals["BETX"]) > 0 else np.nan,
                'alfx': float(vals["ALFX"][0]) if len(vals["ALFX"]) > 0 else np.nan,
                'bety': float(vals["BETY"][0]) if len(vals["BETY"]) > 0 else np.nan,
                'alfy': float(vals["ALFY"][0]) if len(vals["ALFY"]) > 0 else np.nan,
                'mux': float(vals["MUX"][0]) if len(vals["MUX"]) > 0 else np.nan,
                'muy': float(vals["MUY"][0]) if len(vals["MUY"]) > 0 else np.nan,
                'L': float(np.nansum(vals["L"])) if len(vals["L"]) > 0 else np.nan,
            }

        out_names = list(combined.keys())
        if isinstance(names, str):
            names = [names]
        if names is not None:
            requested = set(names)
            out_names = [name for name in out_names if name in requested]

        return {
            "names": out_names,
            "S": np.array([combined[name]["S"] for name in out_names], dtype=float),
            "betx": np.array([combined[name]["betx"] for name in out_names], dtype=float),
            "alfx": np.array([combined[name]["alfx"] for name in out_names], dtype=float),
            "bety": np.array([combined[name]["bety"] for name in out_names], dtype=float),
            "alfy": np.array([combined[name]["alfy"] for name in out_names], dtype=float),
            "mux": np.array([combined[name]["mux"] for name in out_names], dtype=float),
            "muy": np.array([combined[name]["muy"] for name in out_names], dtype=float),
            "L": np.array([combined[name]["L"] for name in out_names], dtype=float),
        }

    def get_twiss_at_element(self,name):
        optics=self._get_optics_from_twiss_file(name)
        names=list(optics["names"])
        if name not in names:
            raise ValueError(f"Element {name} not found in optics map")
        i=names.index(name)
        return {
            "name": names[i],
            "S": float(optics["S"][i]),
            "betx": float(optics["betx"][i]),
            "alfx": float(optics["alfx"][i]),
            "bety": float(optics["bety"][i]),
            "alfy": float(optics["alfy"][i]),

        }

    def get_elements_indices(self, names):
        if isinstance(names, str):
            names = [names]
        name_to_index = {string: index for index, string in enumerate(self.sequence)}
        return [name_to_index.get(name, np.nan) for name in names]
