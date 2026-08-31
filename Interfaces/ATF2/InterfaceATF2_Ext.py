import RF_Track as rft
import sys, time, math, os, threading, struct, ctypes, subprocess
import numpy as np
from epics import PV, ca, caget
from Interfaces.AbstractMachineInterface import AbstractMachineInterface
from collections import defaultdict
from Interfaces.ATF2.InterfaceATF2_Ext_RFTrack import InterfaceATF2_Ext_RFTrack

class CurrentDropToZeroError(RuntimeError):
    def __init__(self, message, *, target=None, readback=None, magnets=None):
        super().__init__(message)
        self.target = dict(target or {})
        self.readback = dict(readback or {})
        self.magnets = list(magnets or [])

class MagKiWrapper:
    MODE_K_TO_I = 1
    MODE_I_TO_K = 2

    def __init__(self, library_path):
        self.library_path = os.path.abspath(str(library_path))
        self.lib = ctypes.CDLL(self.library_path)
        self._func = self.lib.mag_ki_main
        self._func.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_float,
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
        ]
        self._func.restype = ctypes.c_int

    def _call(self, mode, name, energy_GeV, k_main=0.0, current_main=0.0):
        kvalue = (ctypes.c_float * 2)(float(k_main), 0.0)
        current = (ctypes.c_float * 2)(float(current_main), 0.0)
        field = (ctypes.c_float * 2)(0.0, 0.0)
        efflen = ctypes.c_float(0.0)
        status = int(
            self._func(int(mode), str(name).encode("ascii"), ctypes.c_float(float(energy_GeV)), kvalue, current, ctypes.byref(efflen), field))
        if status != 1:
            raise RuntimeError(f"mag_ki_main failed for {name}, mode={mode}, status={status}")
        return {
            "k": float(kvalue[0]),
            "current": float(current[0]),
            "efflen": float(efflen.value),
            "field": float(field[0]),
        }

    def current_to_k1l(self, name, current_A, energy_GeV):
        return self._call(self.MODE_I_TO_K, name, energy_GeV, current_main=current_A)["k"]

    def k1l_to_current(self, name, k1, energy_GeV):
        return self._call(self.MODE_K_TO_I, name, energy_GeV, k_main=k1)["current"]

class InterfaceATF2_Ext(AbstractMachineInterface):
    def get_name(self):
        return 'ATF2_Ext'

    def __init__(self, bg_shots=10.0, nsamples=10, nominal_intensity=0.15, wfs_intensity=0.1):
        self.screen_backgrounds = {}
        self.nsamples = nsamples
        self.bpm_sample_interval_s = 0.5
        self.twiss_path = os.path.join(os.path.dirname(__file__), 'Ext_ATF2', 'ATF2_EXT_FF_v5.2.twiss')
        self.lattice = rft.Lattice(self.twiss_path)
        self.lattice.set_tt_nsteps(0)
        self.electronmass = 0.51099895 # MeV/c^2
        self.Pref = 1.2999999e3 # MeV/c, until a PV is specified
        self.screen_names = ['OTR0X', 'OTR1X', 'OTR2X', 'OTR3X']
        self.screen_pv_names = {
            'OTR0X': 'mOTR0',
            'OTR1X': 'mOTR1',
            'OTR2X': 'mOTR2',
            'OTR3X': 'mOTR3'
        }
        self.is_simulation = False
        self.bpm_sample_interval_s = 0.5
        self.screen_image_shape = (960, 1280) # image size = 1280 x 960
        self.tracking_interface = InterfaceATF2_Ext_RFTrack()
        # Bpms and correctors in beamline order
        sequence = [
            "MB2X", "ZV1X", "QF1X", "ZV2X", "QD2X", "QF3X", "ZH1X", "ZV3X", "QF4X",
            "ZH2X", "QD5X", "ZV4X", "ZV5X", "QF6X", "QF7X", "ZH3X", "QD8X", "ZV6X",
            "QF9X", "ZH4X", "FONTK1", "ZV7X", "FONTP1", "QD10X", "ZH5X", "QF11X",
            "FONTK2", "ZV8X", "FONTP2", "QD12X", "ZH6X", "QF13X", "QD14X", "FONTP3",
            "ZH7X", "QF15X", "ZV9X", "QD16X", "ZH8X", "QF17X", "ZV10X", "QD18X","OTR0X",
            "ZH9X", "QF19X","OTR1X", "ZV11X", "QD20X","OTR2X" , "ZH10X", "QF21X", "OTR3X","IPT1", "IPT2",
            "IPT3", "IPT4", "QM16FF", "ZH1FF", "ZV1FF", "QM15FF", "QM14FF", "FB2FF",
            "QM13FF", "QM12FF", "QM11FF", "QD10BFF", "QD10AFF", "QF9BFF",
            "MSF6FF", "QF9AFF", "QD8FF", "QF7FF", "QD6FF", "QF5BFF", "MSF5FF",
            "QF5AFF", "QD4BFF", "MSD4FF", "QD4AFF", "QF3FF", "QD2BFF", "QD2AFF",
            "MSF1FF", "QF1FF", "MSD0FF", "QD0FF", "PREIP", "IPA", "IPB", "IPC", "M-PIP"
        ]

        quadrupoles = [
            "QF1X", "QD2X", "QF3X","QF4X", "QD5X", "QF6X", "QF7X","QD8X","QF9X", "QD10X",
            "QF11X", "QD12X", "QF13X", "QD14X", "QF15X","QD16X", "QF17X","QD18X","QF19X",
            "QD20X", "QF21X", "QM16FF", "QM15FF", "QM14FF", "QM13FF", "QM12FF", "QM11FF", "QD10BFF", "QD10AFF", "QF9BFF",
            "QF9AFF", "QD8FF", "QF7FF", "QD6FF", "QF5BFF", "QF5AFF", "QD4BFF", "QD4AFF", "QF3FF", "QD2BFF", "QD2AFF",
            "QF1FF",  "QD0FF"
        ]
        self.quadrupoles = list(quadrupoles)

        screens = ['OTR0X','OTR1X','OTR2X','OTR3X']
        self.screens = list(screens)
        # ATF2' BPMs Epics names
        # https://atf.kek.jp/atfbin/view/ATF/EPICS_DATABASE
        monitors = [
            "MB1X", "MB2X", "MQF1X", "MQD2X", "MQF3X", "MQF4X", "MQD5X", "MQF6X",
            "MQF7X", "MQD8X", "MQF9X", "MQD10X", "MQF11X", "MQD12X", "MQF13X",
            "MQD14X", "MQF15X", "MQD16X", "MQF17X", "MQD18X", "MQF19X",
            "MQF21X", "MQM16FF",
            "MQM15FF", "MQM14FF",
            "MQM11FF", "MQD10BFF", "MQD10AFF", "MQF9AFF",
            "MQD8FF", "MQF7FF", "MQF5BFF", "MQD4BFF", "MQF3FF", "MQD2BFF", "MQD2AFF",
            "MSF1FF", "MPREIP", "MIPB"
            ]

        # keep monitor indexes from MONITOR_INDEX_TO_NAME as the source of truth.
        self.MONITOR_INDEX_TO_NAME = {
            0: "MB1X", 1: "MB2X", 2: "MQF1X", 3: "MQD2X", 4: "MQF3X", 5: "MQF4X", 6: "MQD5X", 7: "MQF6X", 8: "MQF7X",
            9: "MQD8X", 10: "MQF9X", 11: "MQD10X", 12: "MQF11X", 13: "MQD12X", 14: "MQF13X", 15: "MQD14X", 16: "MQF15X", 17: "MQD16X",
            18: "MQF17X", 19: "MQD18X", 20: "MQF19X", 21: "MQD20X", 22: "MQF21X", 23: "IPBPM1", 24: "IPBPM2", 25: "nBPM1",
            26: "nBPM2", 27: "nBPM3", 28: "MQM16FF", 29: "MQM15FF", 30: "MQM14FF", 31: "MFB2FF", 32: "MQM13FF", 33: "MQM12FF",
            34: "MFB1FF", 35: "MQM11FF", 36: "MQD10BFF", 37: "MQD10AFF", 38: "MQF9BFF", 39: "MSF6FF", 40: "MQF9AFF", 41: "MQD8FF",
            42: "MQF7FF", 43: "MQD6FF", 44: "MQF5BFF", 45: "MSF5FF", 46: "MQF5AFF", 47: "MQD4BFF", 48: "MSD4FF", 49: "MQD4AFF",
            50: "MQF3FF", 51: "MQD2BFF", 52: "MQD2AFF", 53: "MSF1FF", 54: "MQF1FF", 55: "MSD0FF", 56: "MQD0FF", 57: "M1&2IP", 58: "MPIP",
            59: "MDUMP", 60: "ICT1X", 61: "ICTDUMP", 62: "MW1X", 63: "MW1IP", 64: "MPREIP", 65: "MIPA", 66: "MIPB"}

        self.sextupoles = ["SF6FF", "SK4FF", "SK3FF", "SF5FF", "SD4FF", "SK2FF", "SK1FF", "SF1FF", "SD0FF"]
        self.screens = ['OTR0X', 'OTR1X','OTR2X','OTR3X']

        # Use list comprehension to filter out strings starting with 'Z' or 'z'
        monitors_from_sequence = [string for string in sequence if not string.lower().startswith('z')]

        # Check if the bpms in the config files are known to Epics
        bpm_ok = all(bpm in monitors for bpm in monitors_from_sequence)
        if not bpm_ok:
            bpms_unknown = [bpm for bpm in monitors_from_sequence if bpm not in monitors]
            print(f'Unknown bpms {bpms_unknown} removed from list')

        # Only retain BPMs in config file which are known by Epics
        sequence_filtered = [element for element in sequence if (element in monitors) or element.lower().startswith('z')]

        # Subset of BPMs and correctors from the config file
        self.sequence = sequence_filtered
        self.sequence_raw = list(sequence)
        self.movable_magnets = [string for string in self.sequence_raw if string.upper().startswith(('MQ', 'MS'))]
        #self.bpms = [string for string in self.sequence if not string.lower().startswith('z')]
        self.corrs = [string for string in self.sequence if string.lower().startswith('z')]
        self.qmag_alias_to_canonical = {}
        for name in self.movable_magnets:
            self.qmag_alias_to_canonical[name] = name
        q_aliases = [string for string in self.sequence_raw if string.upper().startswith('Q')]
        for qname in q_aliases:
            candidate = f"M{qname}"
            if candidate in self.movable_magnets:
                self.qmag_alias_to_canonical[qname] = candidate
            else:
                # Keep alias as-is if no canonical counterpart exists.
                self.qmag_alias_to_canonical[qname] = qname

        # # Index of the selected BPMs in the Epics PV ATF2:monitors
        # self.bpm_indexes = [index for index, string in enumerate(monitors) if string in self.bpms]

        '''
           Sato-san's way:
           BPM order must follow MONITOR_INDEX_TO_NAME, independent from sequence.
           
        '''
        name_to_monitor_index = {name: index for index, name in self.MONITOR_INDEX_TO_NAME.items()}
        #self.bpms = [element for element in self.sequence if not element.lower().startswith('z') and element in name_to_monitor_index]
        #self.bpm_indexes = np.array([name_to_monitor_index[name] for name in self.bpms], dtype=int)
        self.bpm_indexes = np.array(sorted(self.MONITOR_INDEX_TO_NAME.keys()), dtype=int)
        self.bpms = [self.MONITOR_INDEX_TO_NAME[i] for i in self.bpm_indexes]

        # Bunch current monitors
        self.ict_names = [
            'gun:GUNcharge', 'l0:L0charge', 'linacbt:LNEcharge', 'linacbt:BTMcharge',
            'ext:EXTcharge', 'linacbt:BTEcharge', 'BIM:DR:nparticles', 'BIM:IP:nparticles'
        ]
        self.nominal_laser_intensity = nominal_intensity
        self.laser_intensity = PV('RFGun:LaserIntensity1:Read').get()
        self.test_laser_intensity = wfs_intensity
        #PV('RFGun:LaserIntensity1:Read').get()

        # k_T_per_A : integrated-gradient slope GL/I [T/A]
        # L_m       : magnetic length
        self.mag_ki = None
        mag_ki_library_candidates = []

        env_mag_ki_library_path = os.environ.get("ATF2_MAG_KI_LIB", "")
        if env_mag_ki_library_path:
            mag_ki_library_candidates.append(env_mag_ki_library_path)

        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
        libmagnet_dir = os.path.join(repo_root, "Machine specifics, user implementations", "ATF2", "libmagnet")
        mag_ki_library_candidates.extend([
            os.path.join(libmagnet_dir, "libmagnet.dylib"),
            os.path.join(libmagnet_dir, "libmagnet.so"),
        ])

        for mag_ki_library_path in mag_ki_library_candidates:
            if not mag_ki_library_path or not os.path.exists(mag_ki_library_path):
                continue
            try:
                self.mag_ki = MagKiWrapper(mag_ki_library_path)
                print(f"Loaded ATF2 mag_ki library: {mag_ki_library_path}")
                break
            except Exception as exc:
                print(f"ATF2 mag_ki library '{mag_ki_library_path}': {exc} not loaded")

        if self.mag_ki is None:
            print("ATF2 mag_ki library not loaded.")

        # IPBSM hooks and FF knob definitions (real machine)
        self.error_history = []
        self.pv_trigger = PV("IPBSM:FringeScan:RemoteScanStart")
        self.pv_end = PV("IPBSM:FringeScan:RemoteScanFinish")
        self._ipbsm_lock = threading.Lock()
        self.datafile = "/atf/data/ipbsm/knob/knob_fringe_result_v2.dat"

        self.linear_matrix = {
            "Ax": {"SD0FF": (-116.5, 0), "SF1FF": (-35.2, 0), "SD4FF": (-37.8, 0), "SF5FF": (0, 0), "SF6FF": (-623.4, 0)},
            "Ex": {"SD0FF": (-813.0, 0), "SF1FF": (793.5, 0), "SD4FF": (-137.9, 0), "SF5FF": (0, 0), "SF6FF": (1492.9, 0)},
            "Ay": {"SD0FF": (-98.9, 0), "SF1FF": (-9.6, 0), "SD4FF": (-246.7, 0), "SF5FF": (0, 0), "SF6FF": (120.0, 0)},
            "Ey": {"SD0FF": (0, 374.1), "SF1FF": (0, -120.0), "SD4FF": (0, -451.3), "SF5FF": (0, 0), "SF6FF": (0, -64.3)},
            "Coup1": {"SD0FF": (0, -100), "SF1FF": (0, 100), "SD4FF": (0, 0), "SF5FF": (0, 0), "SF6FF": (0, 0)},
            "Coup2": {"SD0FF": (0, 107.8), "SF1FF": (0, 4.2), "SD4FF": (0, 152.6), "SF5FF": (0, 0), "SF6FF": (0, 90.4)},
            "Spare1": {"SD0FF": (-676.0, 0), "SF1FF": (-316.0, 0), "SD4FF": (252.0, 0), "SF5FF": (587.0, 0), "SF6FF": (593.0, 0)},
            "Spare2": {"SD0FF": (0, 78.0), "SF1FF": (0, 188.0), "SD4FF": (0, -55.0), "SF5FF": (-179.0, 0), "SF6FF": (0, 38.0)},
            "Spare3": {"SD0FF": (0, 0), "SF1FF": (0, 0), "SD4FF": (0, 1), "SF5FF": (0, 0), "SF6FF": (0, 0)},
        }
        self.nonlinear_matrix = {
            "Y24": {"SK1FF": 0.0, "SK2FF": 0.0, "SK3FF": 0.0, "SK4FF": 0.0, "SD0FF": 0.119, "SF1FF": -0.013, "SD4FF": -0.554, "SF5FF": -0.083, "SF6FF": -0.175},
            "Y46": {"SK1FF": 0.0, "SK2FF": 0.0, "SK3FF": 0.0, "SK4FF": 0.0, "SD0FF": 0.259, "SF1FF": -0.057, "SD4FF": 1.049, "SF5FF": -0.106, "SF6FF": -0.056},
            "Y22": {"SK1FF": -1.629, "SK2FF": 0.174, "SK3FF": 1.024, "SK4FF": 2.435, "SD0FF": 0.0, "SF1FF": 0.0, "SD4FF": 0.0, "SF5FF": 0.0, "SF6FF": 0.0},
            "Y26": {"SK1FF": 1.763, "SK2FF": -0.126, "SK3FF": 0.463, "SK4FF": -0.701, "SD0FF": 0.0, "SF1FF": 0.0, "SD4FF": 0.0, "SF5FF": 0.0, "SF6FF": 0.0},
            "Y66": {"SK1FF": 5.571, "SK2FF": -0.207, "SK3FF": -4.668, "SK4FF": -6.673, "SD0FF": 0.0, "SF1FF": 0.0, "SD4FF": 0.0, "SF5FF": 0.0, "SF6FF": 0.0},
            "Y44": {"SK1FF": 0.037, "SK2FF": 1.614, "SK3FF": -0.458, "SK4FF": -0.186, "SD0FF": 0.0, "SF1FF": 0.0, "SD4FF": 0.0, "SF5FF": 0.0, "SF6FF": 0.0},
            "Spare": {"SK1FF": 0.0, "SK2FF": 0.0, "SK3FF": 0.0, "SK4FF": 0.0, "SD0FF": 0.0, "SF1FF": 0.0, "SD4FF": 0.0, "SF5FF": 0.0, "SF6FF": 0.0},
        }
        self.corrector_knob_pvs = {
            "corrector 1": {
                "pv_set": "IPBPM:PMC:D:setPosWrite",
                "pv_move": "IPBPM:PMC:D:Move.PROC",
                "pv_busy": "IPBPM:PMC:D:busyStatus",
            },
            "Abe chamber": {
                "pv_set": "WAKECHAMBER:VER:SETVAL",
                "pv_move": "WAKECHAMBER:VER:CHANGE_ABS_AND_MOVE:PROC",
                "pv_busy": "WAKECHAMBER:VER:BUSY",
            },
        }
        self.zscan_knob_name = "Z scan knob"
        self._zscan_knob_values = {self.zscan_knob_name: 0.0}
        self._linear_knob_values = {k: 0.0 for k in self.linear_matrix.keys()}
        self._nonlinear_knob_values = {k: 0.0 for k in self.nonlinear_matrix.keys()}
        self._corrector_knob_values = {k: float("nan") for k in self.corrector_knob_pvs.keys()}
        self.zscan_mode_axes = {
            "2-8": "M8LY",
            "30": "M30LY",
            "174": "M174LY",
        }
        for name in self.movable_magnets:
            self.qmag_alias_to_canonical[name] = name
        q_aliases = [string for string in self.sequence_raw if string.upper().startswith('Q')]
        for qname in q_aliases:
            candidate = f"M{qname}"
            if candidate in self.movable_magnets:
                self.qmag_alias_to_canonical[qname] = candidate
            else:
                # Keep alias as-is if no canonical counterpart exists.
                self.qmag_alias_to_canonical[qname] = qname
        # Prefer Q* aliases in UIs; fallback to canonical names when no alias exists.
        qm_set = set(q_aliases)
        if not qm_set:
            qm_set = set(self.movable_magnets)
        self.qmags = sorted(qm_set, key=lambda n: self.sequence_raw.index(n) if n in self.sequence_raw else 10 ** 9)
        self.qmag_pv = {}
        for alias in self.qmag_alias_to_canonical.keys():
            # Use the alias itself for mover PVs (e.g. QM16FF:MAG:DES:X).
            self.qmag_pv[alias] = self._build_qmag_pv_names(alias)
        self.bg_shots = int(bg_shots)

    def current_to_k1l(self, name, current_A):
        current_A = float(current_A)
        if not np.isfinite(current_A):
            return np.nan
        canonical = self.qmag_alias_to_canonical.get(name, name)
        mag_name_for_library = canonical[1:] if canonical.startswith("M") else canonical
        if self.mag_ki is not None:
            return self.mag_ki.current_to_k1l(mag_name_for_library, current_A, self.Pref / 1e3)
        else: raise KeyError(f"No A-K1 calibration for '{name}' ")

    def k1l_to_current(self, name, k1):
        k1 = float(k1)
        if not np.isfinite(k1):
            raise ValueError(f"Cannot convert K1 for quadrupole '{name}': {k1}")
        canonical = self.qmag_alias_to_canonical.get(name, name)
        mag_name_for_library = canonical[1:] if canonical.startswith("M") else canonical
        if self.mag_ki is not None:
            return self.mag_ki.k1l_to_current(mag_name_for_library, k1, self.Pref / 1e3)
        else:
            raise KeyError(f"No A-K1 calibration for '{name}' ")

    def update_tracking_model_with_clibmagnet(self, nominal_bdes_value = None, nominal_bact_value = None, names=None):
        if isinstance(names, str):
            names = [names]
        bact = np.array(nominal_bact_value, dtype=float)
        bdes = np.array(nominal_bdes_value, dtype=float)
        self.log(f"Updating tracking model for {names}..."
                 f"Current nominal bdes value from the machine [1/m]: {bdes}, "
                 f"Current nominal bact value from the machine [1/m]: {bact}")
        self.tracking_interface.set_quadrupoles(names, bact)
        self.log(f"Reading new values set in the lattice:")
        print(self.tracking_interface.get_quadrupoles(names)["bdes"])
        self.log(f"Changed model values to match the lattice!")

    def insert_screen(self, screen_name):
        screen_pv_name = self.screen_pv_names.get(screen_name)
        if screen_pv_name is None:
            raise ValueError(f"Unknown screen: {screen_name}")
        self.log(f"Getting the state of screen {screen_name}...")
        status = PV(f'{screen_pv_name}:Target:READ:INOUT').get()
        self.log(f"Current status of screen {screen_name} is {status}...")
        PV(f"{screen_pv_name}:Target:WRITE:IN").put(1)
        if not self._wait_for_screen_target_position(screen_name, 1):
            raise RuntimeError(f"Screen {screen_name} was not inserted within time.")

    def extract_screen(self, screen_name):
        screen_pv_name = self.screen_pv_names.get(screen_name)
        if screen_pv_name is None:
            raise ValueError(f"Unknown screen: {screen_name}")
        self.log(f"Getting the state of screen {screen_name}...")
        status = PV(f'{screen_pv_name}:Target:READ:INOUT').get()
        self.log(f"Current status of screen {screen_name} is {status}...")
        PV(f"{screen_pv_name}:Target:WRITE:OUT").put(1)
        if not self._wait_for_screen_target_position(screen_name, 0):
            raise RuntimeError(f"Screen {screen_name} was not extracted within time.")

    def _map_quadrupoles_names_from_lattice(self, name):
        name = str(name)
        elements =  []
        for quad_name in [name, f"{name}_1", f"{name}_2"]:
            try:
                quad_parts = self.lattice[quad_name]
            except Exception as e:
                continue
            if not isinstance(quad_parts, list):
                quad_parts = [quad_parts]
            elements.extend(quad_parts) # adding each element of the list, not the whole list as an element
        return elements

    def _get_elements_positions_show_beamline(self, names=None):
        if isinstance(names, str):
            names = [names]

        all_names = []
        all_s = []
        all_l = []
        s_pos = 0.0

        for element in self.lattice['*']:
            element_name = element.get_name()
            try:
                element_length = float(element.get_length())
            except Exception:
                element_length = 0.0

            if names is None or element_name in names:
                all_names.append(element_name)
                all_s.append(s_pos)
                all_l.append(element_length)

            s_pos += element_length

        return {
            "names": all_names,
            "S": np.array(all_s, dtype=float),
        }

    def _give_elements_to_show_beamline(self, quad_selected):
        mapped_quad_elements = self._map_quadrupoles_names_from_lattice(quad_selected)
        if not isinstance(mapped_quad_elements, list):
            mapped_quad_elements = [mapped_quad_elements]
        start_quad_element_name = mapped_quad_elements[0]
        return start_quad_element_name

    def get_beam_factors(self):
        # TO BE REPLACED WITH A PV OF REAL BEAM ENERGY
        gamma_rel = np.sqrt((self.Pref / self.electronmass) ** 2 + 1.0)
        beta_rel = np.sqrt(1.0 - 1.0 / gamma_rel ** 2)
        beta_gamma = gamma_rel * beta_rel
        if not np.isfinite(beta_gamma) or beta_gamma <= 0:
            raise RuntimeError("Invalid beam factors")
        return gamma_rel, beta_rel, beta_gamma

    def _read_twiss_file(self):
        with open(self.twiss_path, "r") as file:
            lines = [line.strip() for line in file if line.strip()]
        star_symbol = next(i for i, line in enumerate(lines) if line.startswith("*"))
        dollar_sign = next(i for i, line in enumerate(lines) if line.startswith("$") and i > star_symbol)
        columns = lines[star_symbol].lstrip("*").split()
        return lines, columns, dollar_sign

    def  _get_twiss_s_positions(self, names):
        names = list(names)
        lines, columns, dollar_sign = self._read_twiss_file()
        try:
            name_column = columns.index("NAME")
            s_column = columns.index("S")
        except ValueError:
            return [np.nan] * len(names)
        s_pos = {}

        for line in lines[dollar_sign + 1:]:
            data = line.split()
            if len(data) <= max(name_column, s_column):
                continue
            elem_name = data[name_column].strip('"')
            try:
                s_pos[elem_name] = float(data[s_column])
            except ValueError:
                continue
        return [s_pos.get(name, np.nan) for name in names]

    @staticmethod
    def make_safe_float(value, default = np.nan): #so even if pv returns none, empty array or whatever, interface still works
        try:
            if value is None:
                return float(default)
            arr = np.asarray(value)
            if arr.size == 0:
                return float(default)
            return float(arr.flat[0])
        except Exception:
            return float(default)

    def _wait_for_pv_readback(self, pv_name, target, tolerance=1e-3, timeout=10.0):
        return self._wait_for_readback(lambda: self.make_safe_float(PV(pv_name).get(), default=np.nan), target, description=pv_name, tolerance=tolerance, timeout=timeout)

    def _valid_pv_value(self, pv_names, default = np.nan):
        for pv_name in pv_names:
            try:
                value = caget(pv_name)
            except Exception:
                continue
            value = self.make_safe_float(value, default = np.nan)
            if np.isfinite(value):
                return value
        return float(default)

    def get_movable_magnets_names(self):
        return self.movable_magnets

    def get_quadrupole_movers_names(self):
        return list(self.qmags)

    def predict_emittance_scan_response(self, *args, **kwargs): # from optimizer, delegates all variables to the simulation
        return self.tracking_interface.predict_emittance_scan_response(*args, **kwargs)

    def get_R_matrix_scan(self, *args, **kwargs):
        return self.tracking_interface.get_R_matrix_scan(*args, **kwargs)

    def get_phase_space_transport_to_screens(self, *args, **kwargs):
        return self.tracking_interface.get_phase_space_transport_to_screens(*args, **kwargs)

    def get_twiss_evolution(self, *args, **kwargs):
        return self.tracking_interface.get_twiss_evolution(*args, **kwargs)

    def _quadrupole_current_pv_name(self,name):
        if name.startswith("M") and name[1:].startswith(("QF", "QD", "QM")):
            return name[1:]
        return name

    def _quad_mover_pv_name(self,name):
        return self.qmag_alias_to_canonical.get(name, name)

    def get_quadrupoles(self, names=None, include_pv_names=False):
        print(" 'get_quadrupoles' running...")
        if names is None:
            names = self.quadrupoles # quadrupoles names
        if type(names) == str:
            names = [names]
        names = [name for name in names if name in self.quadrupoles]
        ides, iact = [], []
        xdes, ydes, rolldes = [], [], []
        xact, yact, rollact = [], [], []

        for name in names:
            print(" 'Getting quadrupoles' PVs ...")
            current_name = self._quadrupole_current_pv_name(name)
            mover_name = self._quad_mover_pv_name(name)
            mover_pv = self.qmag_pv.get(mover_name)
            if mover_pv is None:
                mover_pv = self._build_qmag_pv_names(mover_name)
                self.qmag_pv[mover_name] = mover_pv
            desired_current = self._pv_get(f"{current_name}:currentWrite", default=np.nan, timeout=0.7)
            actual_current = self._pv_get(f"{current_name}:current", default=np.nan, timeout=0.7)
            if not np.isfinite(actual_current):
                actual_current = self._pv_get(f"{current_name}:currentRead", default=np.nan, timeout=0.7)
            if not np.isfinite(desired_current):
                desired_current = actual_current

            ides.append(desired_current)
            iact.append(actual_current)
            xdes.append(self._pv_get(mover_pv.get("pv_set_x"), default=np.nan))
            ydes.append(self._pv_get(mover_pv.get("pv_set_y"), default=np.nan))
            rolldes.append(self._pv_get(mover_pv.get("pv_set_roll"), default=np.nan))
            xact.append(self._pv_get(mover_pv.get("pv_read_enc_x"), default=np.nan))
            yact.append(self._pv_get(mover_pv.get("pv_read_enc_y"), default=np.nan))
            rollact.append(self._pv_get(mover_pv.get("pv_read_enc_roll"), default=np.nan))

        ides = np.array(ides, dtype=float)
        iact = np.array(iact, dtype=float)

        try:
            bdes = np.array([self.current_to_k1l(n, i) for n, i in zip(names, ides)], dtype=float)
            bact = np.array([self.current_to_k1l(n, i) for n, i in zip(names, iact)], dtype=float)
            self.update_tracking_model_with_clibmagnet(names=names, nominal_bdes_value = bdes, nominal_bact_value = bact)
        except Exception as exc:
            print(f"current to K1 conversion failed for quadrupoles {names}: {exc}")
            bdes = np.array([np.nan for _ in names], dtype=float)
            bact = np.array([np.nan for _ in names], dtype=float)

        data = {
            "names": names,
            "bdes": bdes, # 1/m,
            "bact": bact,  # 1/m
            "ides": np.array(ides, dtype=float),
            "iact": np.array(iact, dtype=float),
            "xdes": np.array(xdes, dtype=float),
            "ydes": np.array(ydes, dtype=float),
            "rolldes": np.array(rolldes, dtype=float),
            "xact": np.array(xact, dtype=float),
            "yact": np.array(yact, dtype=float),
            "rollact": np.array(rollact, dtype=float),
        }
        if include_pv_names:
            data["pvs"] = {name: dict(self.qmag_pv[name]) for name in names}
        return data

    def set_quadrupoles(self, names, k1l_values, track=True):
        if type(names) == str:
            names = [names]
        if not isinstance(k1l_values, (list, tuple, np.ndarray)):
            k1l_values = [k1l_values]
        if len(names) != len(k1l_values):
            raise ValueError(f"len(names)={len(names)} != len(k1l_values)={len(k1l_values)} in set_quadrupoles")

        for name, k1 in zip(names, k1l_values):
            if name not in self.quadrupoles:
                raise ValueError(f"Quadrupole '{name}' is not magnet list.")
            current_name = self._quadrupole_current_pv_name(name)
            print(name)
            target_current = self.k1l_to_current(current_name, float(k1))  # A
            self._pv_put(f"{current_name}:currentWrite", float(target_current))
            self._wait_for_magnet_readback(current_name, float(target_current))

    """Methods for OTRs from mOTRs_measurement.py"""

    # --- Gaussian Fit ---
    def gaussian(self, x, amplitude, mean, stddev, offset):
        """1D Gaussian function for curve fitting."""
        return amplitude * np.exp(-((x - mean) / (2 * stddev)) ** 2) + offset

    def plot_otr_analysis_inset(self, ax_main, img_data, title_str, h_factor_um_px=1.0, v_factor_um_px=1.0):

        im = ax_main.imshow(img_data, cmap='gray', origin='lower')
        ax_main.set_title(title_str)
        ax_main.axis('off')

        plt.colorbar(im, ax=ax_main, fraction=0.046, pad=0.04, label='Pixel Intensity (Counts)')

        h, w = img_data.shape
        proj_x = np.sum(img_data, axis=0)
        proj_y = np.sum(img_data, axis=1)
        x_coords = np.arange(w)
        y_coords = np.arange(h)
        p0_x = [np.max(proj_x), np.argmax(proj_x), w / 10, np.min(proj_x)]
        p0_y = [np.max(proj_y), np.argmax(proj_y), h / 10, np.min(proj_y)]

        sigma_h_um = None
        sigma_v_um = None

        try:
            popt_x, pcov_x = curve_fit(self.gaussian, x_coords, proj_x, p0=p0_x)
            amp_x, mean_x, stddev_x, offset_x = popt_x
            fit_x = self.gaussian(x_coords, *popt_x)
            beam_size_x_px = abs(stddev_x) * 2
            sigma_h_um = beam_size_x_px * h_factor_um_px

            popt_y, pcov_y = curve_fit(self.gaussian, y_coords, proj_y, p0=p0_y)
            amp_y, mean_y, stddev_y, offset_y = popt_y
            fit_y = self.gaussian(y_coords, *popt_y)
            beam_size_y_px = abs(stddev_y) * 2
            sigma_v_um = beam_size_y_px * v_factor_um_px

            # PLOTS
            # Horizontal histogram (top margin)
            ax_hist_x = ax_main.inset_axes([0.0, 1.05, 1.0, 0.2], transform=ax_main.transAxes)
            ax_hist_x.plot(x_coords, proj_x, label='H Projection')
            ax_hist_x.plot(x_coords, fit_x, 'r--', label=f'Size: {sigma_h_um:.2f} $\\mu$m')
            ax_hist_x.legend(fontsize='x-small', loc='upper right')
            ax_hist_x.axis('off')

            # Vertical histogram (right margin)
            ax_hist_y = ax_main.inset_axes([1.05, 0.0, 0.2, 1.0], transform=ax_main.transAxes)
            ax_hist_y.plot(proj_y, y_coords, label='V Projection')
            ax_hist_y.plot(fit_y, y_coords, 'r--', label=f'Size: {sigma_v_um:.2f} $\\mu$m')
            ax_hist_y.legend(fontsize='x-small', loc='lower right')
            ax_hist_y.axis('off')

        except RuntimeError:
            print(f"Could not fit Gaussian for {title_str}")
            ax_hist_x = ax_main.inset_axes([0.0, 1.05, 1.0, 0.2], transform=ax_main.transAxes)
            ax_hist_y = ax_main.inset_axes([1.05, 0.0, 0.2, 1.0], transform=ax_main.transAxes)
            ax_hist_x.axis('off')
            ax_hist_y.axis('off')

        return proj_x, proj_y, sigma_h_um, sigma_v_um

    # --- Data Acquisition Functions ---
    def get_pixel_calibrations(self, screen_pv_name):
        h_pv_name = f'{screen_pv_name}:H:x1:Calibration:Factor'
        v_pv_name = f'{screen_pv_name}:V:x1:Calibration:Factor'
        h_factor = PV(h_pv_name).get()
        v_factor = PV(v_pv_name).get()
        if h_factor is None or v_factor is None:
            h_factor = 1.0
            v_factor = 1.0
        return h_factor, v_factor

    def acquire_screen_image(self, screen_name, min_total_intensity=135000, max_retries=3):
        """
        Acquires an OTR image with background subtraction.
        Takes the background once, then loops the beam acquisition if the
        total integrated intensity is below min_total_intensity.
        """
        screen_pv_name = self.screen_pv_names.get(screen_name)
        pv_in_name = f'{screen_pv_name}:Target:WRITE:IN'
        pv_out_name = f'{screen_pv_name}:Target:WRITE:OUT'
        pv_img_data_name = f'{screen_pv_name}:IMAGE:ArrayData'
        pv_acquire_name = f'{screen_pv_name}:CAMERA:Acquire'
        otr_in_pv = PV(pv_in_name)
        otr_out_pv = PV(pv_out_name)
        image_data_pv = PV(pv_img_data_name)
        image_acquire_pv = PV(pv_acquire_name)
        bg_img = self.screen_backgrounds[screen_name]
        print(f"Inserting screen {screen_name} for beam...")
        self.insert_screen(screen_name)
        time.sleep(5)
        beam_frames = []
        for attempt in range(1, max_retries + 1):
            beam_frames = []
            for i in range(5):
                image_acquire_pv.put(1)
                time.sleep(3)
                beam_data = image_data_pv.get()
                image_acquire_pv.put(0)
                beam_frames.append(beam_data.reshape(960, 1280).astype(np.float64))
            beam_img = np.median(beam_frames, axis=0)
            subtracted_img = beam_img - bg_img
            subtracted_img[subtracted_img < 0] = 0
            tot_intensity = np.sum(subtracted_img)
            if tot_intensity >= min_total_intensity:
                break
            else:
                print(f" -> WARNING: Total intensity too low ({tot_intensity:,.0f} < threshold of {min_total_intensity:,.0f}).")
                if attempt < max_retries:
                    print(" -> Retaking beam image (keeping screen inserted)...\n")
                else:
                    print(" -> Max retries reached. Proceeding with the empty/low-intensity frame.")
        print(" -> Retracting screen...")
        otr_out_pv.put(1)
        return subtracted_img, bg_img, beam_img

    def _screen_data_from_image(self, image, hpixel, vpixel): # better be subtracted!
        if image is None or np.asarray(image, dtype=float).ndim!=2 or np.asarray(image, dtype=float).size==0:
            return np.nan, np.nan, np.nan, np.nan, 0.0, np.zeros((1, 1)), np.array([0.0, 1.0]), np.array([0.0, 1.0])

        # I_xi_yi, I = I_beam - I_background
        img = np.asarray(image, dtype=float).copy()
        img[~np.isfinite(img)] = 0.0
        img[img < 0] = 0.0
        ny, nx = img.shape # e.g. ny = no. of rows, nx = no. of columns
        j = np.arange(nx)
        i = np.arange(ny)
        summed_intensity = np.sum(img)
        if summed_intensity <= 0:
            hedges = (np.arange(nx + 1) - nx / 2) * hpixel
            vedges = (np.arange(ny + 1) - ny / 2) * vpixel
            return np.nan, np.nan, np.nan, np.nan, 0.0, img, hedges, vedges
        proj_x = np.sum(img, axis=0)
        proj_y = np.sum(img, axis=1)

        # center of each pixel can be expressed as: x_j = (j - (N_x - 1)/2)*hpixel, y_i = (i - (N_i - 1)/2)vhpixel

        '''
        x_centers and y_centers are arrays containing the physical coordinates of the centre of each pixel. 
        They allow to convert the image from pixel intensities into beam position and beam size.
        '''
        x_pixels_positions = (j - (nx - 1) / 2) * hpixel # coordinates of centre of each pixel
        y_pixels_positions = (i - (ny - 1) / 2) * vpixel
        x_mean_positions = np.sum((x_pixels_positions * proj_x) / summed_intensity)
        y_mean_positions = np.sum((y_pixels_positions * proj_y) / summed_intensity)

        sigx = np.sqrt(np.sum((x_pixels_positions - x_mean_positions)**2 * proj_x) / summed_intensity)
        sigy = np.sqrt(np.sum((y_pixels_positions - y_mean_positions)**2 * proj_y) / summed_intensity)

        hedges = (np.arange(nx+1) - nx / 2) * hpixel
        vedges = (np.arange(ny+1) - ny / 2) * vpixel


        '''
        Logic commented is for reading sigx and sigy from precomputed PVs, not from image analysis.
        
        # # mOTR:analyzer:dispersion:selectedmotr
        # # hack to avoid background subtraction
        # command = f"caput mOTR:analyzer:dispersion:selectedmotr {screen_pv_name[-1]}"
        # result = subprocess.run(command,shell=True)
        # time.sleep(1)
        # result = subprocess.run(command,shell=True)
        # time.sleep(1)
        # result = subprocess.run(command,shell=True)
        # time.sleep(10)
        # sigx_pv = f"mOTR:analyzer:size:H"
        # sigy_pv = f"mOTR:analyzer:size:V"
        # sigx_pv_value = self.make_safe_float(PV(sigx_pv).get(), default=np.nan)
        # sigy_pv_value = self.make_safe_float(PV(sigy_pv).get(), default=np.nan)
        # 
        # sigx_from_image = float(np.sqrt(max(np.sum(((x_centers - x_mean) ** 2) * proj_x) / total, 0.0))) # mm
        # sigy_from_image = float(np.sqrt(max(np.sum(((y_centers - y_mean) ** 2) * proj_y) / total, 0.0))) # mm

        # # sigx_prev = self.make_safe_float(PV(sigx_pv).get(), default=np.nan)
        # # sigy_prev = self.make_safe_float(PV(sigy_pv).get(), default=np.nan)
        # # max_retries = 5
        # # for attempt in range(max_retries):
        # #     time.sleep(5)
        # #     sigx_new = self.make_safe_float(PV(sigx_pv).get(), default=np.nan)
        # #     sigy_new = self.make_safe_float(PV(sigy_pv).get(), default=np.nan)
        # #     if not np.isfinite(sigx_prev) or not np.isfinite(sigy_prev):
        # #         sigx_prev, sigy_prev = sigx_new, sigy_new
        # #         continue
        # #     if sigx_prev <= 0 or sigy_prev <= 0 or sigx_new <= 0 or sigy_new <= 0:
        # #         sigx_prev, sigy_prev = sigx_new, sigy_new
        # #         continue
        # #     change_x = max(sigx_new / sigx_prev, sigx_prev / sigx_new)
        # #     change_y = max(sigy_new / sigy_prev, sigy_prev / sigy_new)
        # #     if change_x <= 8 and change_y <= 8:
        # #         sigx = sigx_new / 1000.0
        # #         sigy = sigy_new / 1000.0
        # #         break
        # #     print("Screen size changed too much between measurements of sigx and sigy. Remeasuring...")
        # #     sigx_prev, sigy_prev = sigx_new, sigy_new
        # # else:
        # #     sigx = sigx_prev / 1000.0 if np.isfinite(sigx_prev) else sigx
        # #     sigy = sigy_prev / 1000.0 if np.isfinite(sigy_prev) else sigy
        # 
        # print("sigx from precomputed PV: ", sigx_pv_value)
        # print("sigy from precomputed PV: ", sigy_pv_value)
        # print("sigx from image analysis: ", sigx_from_image)
        # print("sigy from image analysis: ", sigy_from_image)
        # # np.average(h[1:], weights=np.sum(i,axis=0)) # andrea's suggestion
        # # mOTR:analyzer:size
        '''
        return x_mean_positions, y_mean_positions, sigx, sigy, summed_intensity, img, hedges, vedges

    def _wait_for_screen_target_position(self, screen_name, target, timeout=10.0, poll_interval=0.05):
        screen_pv_name = self.screen_pv_names.get(screen_name)
        t0 = time.perf_counter()
        status = np.nan
        while time.perf_counter() - t0 < timeout:
            status = self.make_safe_float(PV(f'{screen_pv_name}:Target:READ:INOUT').get(), default=np.nan)
            if status == 0 and target == 0:
                return True
            if status > 0 and target > 0:
                return True
            time.sleep(poll_interval)
        self.log(
            f'Warning: {screen_name} did not reach target state = {target:.6g} '
            f'within {timeout:.2f}s. Last readback = {status:.6g}'
        )
        return False

    def get_screens(self, names=None):
        print('Reading screens...')
        if isinstance(names, str):
            names = [names]
        selected_names = self.screen_names if names is None else [name for name in self.screen_names if name in names]
        s_positions = self._get_twiss_s_positions(selected_names)
        hpixel_list = []
        vpixel_list = []
        xb_list = []
        yb_list = []
        sigx_list = []
        sigy_list = []
        sum_list = []
        images = []
        background_images = []
        beam_images = []
        hedges_all = []
        vedges_all = []
        inout_list = [] # is screen in or out

        for screen_name in selected_names:
            screen_pv_name = self.screen_pv_names.get(screen_name)
            otr_id = screen_name.replace('OTR', '')
            hpixel_um, vpixel_um = self.get_pixel_calibrations(screen_pv_name)
            hpixel = hpixel_um / 1000.0
            vpixel = vpixel_um / 1000.0

            try:
                if screen_name not in self.screen_backgrounds:
                    self.log(f"Acquiring background image for {screen_name}.")
                    self.acquire_screen_background(screen_name, frames=10)
                subtracted_img, bg_img, beam_img = self.acquire_screen_image(screen_name)
                x_mean, y_mean, sigx, sigy, total, img, hedges, vedges = self._screen_data_from_image(subtracted_img, hpixel, vpixel)

            except Exception as e:
                self.log(f"Couldn't acquire screen image for {screen_name}, because: {e}")
                x_mean = np.nan
                y_mean = np.nan
                sigx = np.nan
                sigy = np.nan
                total = 0.0
                subtracted_img = np.zeros((1, 1), dtype=float)
                bg_img = np.zeros((1, 1), dtype=float)
                beam_img = np.zeros((1, 1), dtype=float)
                hedges = np.array([0.0, 1.0], dtype=float)
                vedges = np.array([0.0, 1.0], dtype=float)

            status = PV(f'{screen_pv_name}:Target:READ:INOUT').get()
            hpixel_list.append(hpixel)
            vpixel_list.append(vpixel)
            xb_list.append(x_mean)
            yb_list.append(y_mean)
            sigx_list.append(sigx)
            sigy_list.append(sigy)
            sum_list.append(total)
            images.append(np.asarray(subtracted_img, dtype=float))
            background_images.append(np.asarray(bg_img, dtype=float))
            beam_images.append(np.asarray(beam_img, dtype=float))
            hedges_all.append(np.asarray(hedges, dtype=float))
            vedges_all.append(np.asarray(vedges, dtype=float))
            inout_list.append(status)

        screens = {
            "names": selected_names,
            "hpixel": np.asarray(hpixel_list, dtype=float), # mm/pixel
            "vpixel": np.asarray(vpixel_list, dtype=float), # mm/pixel
            "x": np.asarray(xb_list, dtype=float), # mm
            "y": np.asarray(yb_list, dtype=float), # mm
            "sigx": np.asarray(sigx_list, dtype=float), # mm
            "sigy": np.asarray(sigy_list, dtype=float), # mm
            "sum": np.asarray(sum_list, dtype=float),
            "hedges": hedges_all, # mm
            "vedges": vedges_all, # mm
            "images": images,
            "background_images": background_images,
            "beam_images": beam_images,
            "S": np.asarray(s_positions, dtype=float), # "S": np.full(len(selected_names), np.nan)
            "inout": np.asarray(inout_list, dtype=float),
        }
        return screens

    def acquire_screen_background(self, screen_name, frames = None):
        if frames is None: frames = self.bg_shots
        self.extract_screen(screen_name)
        screen_pv_name = self.screen_pv_names.get(screen_name)
        pv_in_name = f'{screen_pv_name}:Target:WRITE:IN'
        pv_out_name = f'{screen_pv_name}:Target:WRITE:OUT'
        pv_img_data_name = f'{screen_pv_name}:IMAGE:ArrayData'
        pv_acquire_name = f'{screen_pv_name}:CAMERA:Acquire'
        otr_in_pv = PV(pv_in_name)
        otr_out_pv = PV(pv_out_name)
        image_data_pv = PV(pv_img_data_name)
        image_acquire_pv = PV(pv_acquire_name)
        background_frames = []
        if frames == 0:
            bg_img = np.zeros((960, 1280), dtype=float)
        else:
            for frame in range(frames):
                self.log(f"Acquiring {frame}/{frames} background frames...")
                image_acquire_pv.put(1)
                time.sleep(1)
                bg_data = image_data_pv.get()
                image_acquire_pv.put(0)
                background_frames.append(bg_data.reshape(960, 1280).astype(np.float64))
                self.log(f"Acquired {frame}/{frames} background frames...")
            bg_img = np.median(background_frames, axis=0)
            if not background_frames:
                raise RuntimeError(f"No background frames available for {screen_name}")
            self.log(f"Acquired {frames} background frames. Calculating the median...")
            bg_img = np.median(np.stack(background_frames, axis=0), axis=0)
            self.log(f"Median calculated.")
            self.screen_backgrounds[screen_name] = bg_img

        return bg_img

    def change_energy(self, delta_freq=4):
        PV('RAMP:CONTROL_ON_SW').put(1)
        self._wait_for_pv_readback('RAMP:CONTROL_ON_SW', 1)
        # delta_freq in kHz
        ### delta_freq MUST MATCH :MI2: to EPICS --> means "MINUS2"
        # PV('RAMP:MI2:ONOFF_SW').put(1)
        PV('RAMP:PL4:ONOFF_SW').put(1)
        self._wait_for_pv_readback('RAMP:PL4:ONOFF_SW', 1)
        DR_freq = 714e3  # 714 MHz in kHz
        DR_momentum_compaction = 2.1e-3
        dP_P = -float(delta_freq) / DR_freq / DR_momentum_compaction
        return dP_P

    def reset_energy(self):
        PV('RAMP:CONTROL_OFF_SW').put(0)
        self._wait_for_pv_readback('RAMP:CONTROL_ON_SW', 0)

    def change_intensity(self, laserintensity=None):
        if laserintensity is None:
            laserintensity = self.test_laser_intensity
        print(f'Changing laser intensity to {laserintensity}...')
        self.laser_intensity = self.make_safe_float(PV('RFGun:LaserIntensity1:Read').get(), default=np.nan)
        laser_intensity = float(laserintensity) * 100 * 5
        PV('RFGun:LaserIntensity1:Write').put(laser_intensity)
        self._wait_for_pv_readback('RFGun:LaserIntensity1:Read', laser_intensity)
        return self

    def reset_intensity(self):
        new_laser_intensity = self.nominal_laser_intensity
        print(f'Resetting laser intensity to {new_laser_intensity}...')
        laser_intensity = new_laser_intensity * 100 * 5 # Korysko dixit: 100 for percent, 5 convesion factor
        PV('RFGun:LaserIntensity1:Write').put(laser_intensity)
        self._wait_for_pv_readback('RFGun:LaserIntensity1:Read', laser_intensity)
        return self

    def get_beam_settings(self):
        settings = {"energy": {}, "intensity": {}}
        for name, pv_name in (("ramp_control", "RAMP:CONTROL_ON_SW"), ("ramp_pl4", "RAMP:PL4:ONOFF_SW")):
            settings["energy"][name] = float(PV(pv_name).get())
            settings["intensity"]["laser_intensity1"] = float(PV('RFGun:LaserIntensity1:Read').get())
        return settings

    def restore_beam_settings(self, settings):
        settings = settings or {}
        energy = settings.get("energy", {})
        for name, pv_name in (("ramp_control", "RAMP:CONTROL_ON_SW"), ("ramp_pl4", "RAMP:PL4:ONOFF_SW")):
            if name in energy:
                target = float(energy[name])
                PV(pv_name).put(target)
                self._wait_for_pv_readback(pv_name, target)
        intensity = settings.get("intensity", {}).get("laser_intensity1")
        if intensity is not None:
            target = float(intensity)
            PV('RFGun:LaserIntensity1:Write').put(target)
            self._wait_for_pv_readback('RFGun:LaserIntensity1:Read', target)
        return True

    def get_sequence(self):
        return self.sequence

    def get_hcorrectors_names(self):
        return [string for string in self.corrs if (string.lower().startswith('zh')) or (string.lower().startswith('zx'))]

    def get_vcorrectors_names(self):
        return [string for string in self.corrs if string.lower().startswith('zv')]

    def get_elements_indices(self, names):
        if isinstance(names, str):
            names = [names]
        sequence_for_index = list(getattr(self, "sequence_raw", self.sequence))
        name_to_index = {string: index for index, string in enumerate(sequence_for_index)}

        for alias, canonical in getattr(self, "qmag_alias_to_canonical", {}).items():
            if alias in name_to_index and canonical not in name_to_index:
                name_to_index[canonical] = name_to_index[alias]
            if canonical in name_to_index and alias not in name_to_index:
                name_to_index[alias] = name_to_index[canonical]
        return [name_to_index.get(name, np.nan) for name in names]

    def get_target_dispersion(self, names=None): # for DR too
        if names is None:
            names = self.bpms
        if isinstance(names, str):
            names = [names]
        twiss_path = self.twiss_path
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
        print("Reading ict's...")
        charge = []
        for ict in self.ict_names:
            pv = PV(f'{ict}')
            charge.append(pv.get())

        icts = {
            "names": self.ict_names,
            "charge": np.array(charge),
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
        print("Reading correctors' strengths...")

        if isinstance(names, str):
            names = [names]

        corr_names = self.corrs if names is None else [corr for corr in self.corrs if corr in names]
        bdes, bact = [], []
        for corrector in corr_names:
            des_name = f'{corrector}:currentWrite'
            act_name = f'{corrector}:currentRead'
            des_val = self._pv_get(des_name, default=np.nan, timeout=0.7)
            act_val = self._pv_get(act_name, default=np.nan, timeout=0.7)
            if np.isnan(des_val) or np.isnan(act_val):
                print(f"Corrector PV read failed: {corrector} ({des_name}, {act_name})")
            bdes.append(des_val)
            bact.append(act_val)

        return {
            "names": corr_names,
            "bdes": np.array(bdes, dtype=float),
            "bact": np.array(bact, dtype=float),
        }

    def get_bpms(self, names=None):
        print('Reading bpms...')
        if isinstance(names, str):
            names = [names]
        x_list, y_list, tmit_list = [], [], []
        s_ref = None
        bpm_names = None
        p = PV("ATF2:monitors")
        sample = 0
        while sample < self.nsamples:
            try:
                print(f"Reading BPM sample {sample + 1}/{self.nsamples}...")
                raw = np.asarray(p.get(), dtype=float)
                a = raw.reshape((-1, 10))
                bpm_indexes = np.array([i for i in self.bpm_indexes if i < a.shape[0]], dtype=int)

                current_names = [self.MONITOR_INDEX_TO_NAME[i] for i in bpm_indexes]
                bpm = a[bpm_indexes]

                status = bpm[:, 0]
                status = np.where(status == 1, 1.0, 0.0)

                x = bpm[:, 1]
                y = bpm[:, 2]
                tmit = status * bpm[:, 3]
                s = bpm[:, 4]

                if s_ref is None:
                    s_ref = s.copy()
                    bpm_names = current_names

                x_list.append(x)
                y_list.append(y)
                tmit_list.append(tmit)

                sample += 1

                if sample < self.nsamples:
                    time.sleep(self.bpm_sample_interval_s)

            except Exception as e:
                print(f'An error occurred while reading BPMs: {e}')
                time.sleep(self.bpm_sample_interval_s)

        bpms = {
            "names": list(bpm_names),
            "x": np.vstack(x_list) / 1e3,
            "y": np.vstack(y_list) / 1e3,
            "tmit": np.vstack(tmit_list),
            "S": np.asarray(s_ref, dtype=float),
        }

        if names is not None:
            idx = [i for i, s in enumerate(bpms["names"]) if s in names]
            bpms = {
                "names": [bpms["names"][i] for i in idx],
                "x": bpms["x"][:, idx],
                "y": bpms["y"][:, idx],
                "tmit": bpms["tmit"][:, idx],
                "S": bpms["S"][idx],
            }

        return bpms

    def get_sextupoles(self, names=None):
        if names is None:
            names = self.sextupoles
        if isinstance(names, str):
            names = [names]
        bdes, bact = [], []

        for name in names:
            bdes.append(self._pv_get(f"{name}:currentWrite"))
            bact.append(self._pv_get(f"{name}:currentRead"))
        return {
            "names": names,
            "bdes": np.array(bdes, dtype=float),
            "bact": np.array(bact, dtype=float),
        }

    def set_sextupoles(self, names, values):
        if isinstance(names, str):
            names = [names]
        values = np.asarray(values, dtype=float).reshape(-1)
        if len(names) != len(values):
            raise ValueError("len(names) != len(values) in set_sextupoles")
        for name, value in zip(names, values):
            self._pv_put(f"{name}:currentWrite", float(value))
            self._wait_for_magnet_readback(name, value)

    def _wait_for_magnet_readback(self, magnet, target, tolerance=1e-4, timeout=1.0, poll_interval=0.05):
        readback_pv = PV(f'{magnet}:currentRead')
        t0 = time.perf_counter()
        last_value = np.nan
        while time.perf_counter() - t0 < timeout:
            try:
                last_value = self.make_safe_float(readback_pv.get(), default=np.nan)
            except Exception:
                last_value = np.nan
            if np.isfinite(last_value) and abs(last_value - float(target)) <= tolerance:
                return True
            time.sleep(poll_interval)
        print(
            f'Warning: {magnet}:currentRead did not reach target {float(target):.6g} '
            f'within {timeout:.2f}s. Last readback = {last_value:.6g}'
        )
        return False

    def set_correctors(self, names, corr_vals):
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        if len(names) != len(corr_vals):
            print('Error: len(names) != len(corr_vals) in set_correctors(names, corr_vals)')
        for corrector, corr_val in zip(names, corr_vals):
            pv_des = PV(f'{corrector}:currentWrite')
            pv_des.put(corr_val)
            self._wait_for_magnet_readback(corrector, corr_val)

    def vary_correctors(self, names, corr_vals):
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        if len(names) != len(corr_vals):
            print('Error: len(names) != len(corr_vals) in vary_correctors(names, corr_vals)')
        for corrector, corr_val in zip(names, corr_vals):
            pv_des = PV(f'{corrector}:currentWrite')
            curr_val = self.make_safe_float(pv_des.get(), default=np.nan)
            target = curr_val + float(corr_val)
            pv_des.put(target)
            self._wait_for_magnet_readback(corrector, target)


    '''
    METHODS FROM SATO-SAN'S REPO:
    '''

    def get_ipbsm(self, timeout=300, file_wait=330.0, poll=0.1):
        dat = self.get_ipbsm_full(timeout=timeout, file_wait=file_wait, poll=poll)
        return float(dat["modulation"]), abs(float(dat["error"]))

    def get_ipbsm_state(self):
        dat = self.get_ipbsm_full()
        return {
            "modulation": float(dat.get("modulation", float("nan"))),
            "modulation_error": float(dat.get("error", float("nan"))),
            "sigma_y_m": float(dat.get("beamsize", float("nan"))),
            "angle_deg": float(dat.get("phase", float("nan"))),
            "average": float(dat.get("average", float("nan"))),
            "ict_average": float(dat.get("ict_average", float("nan"))),
            "filename": str(dat.get("filename", "")),
        }

    def _zscan_pvs_for_mode(self, scan_mode_label):
        mode = self._normalize_scan_mode_label(scan_mode_label)
        axis = str(self.zscan_mode_axes[mode])
        base = f"IPBSM:{axis}"
        return {
            "mode": mode,
            "axis": axis,
            "pv_position": f"{base}:Position",
            "pv_state": f"{base}:State",
            "pv_busy": f"{base}:Busy",
        }

    def get_zscan_status(self, scan_mode_label="30"):
        pvs = self._zscan_pvs_for_mode(scan_mode_label)
        pos = float(self._pv_get(pvs["pv_position"], default=np.nan, timeout=0.7))
        busy = float(self._pv_get(pvs["pv_busy"], default=np.nan, timeout=0.7))
        return {
            "mode": pvs["mode"],
            "axis": pvs["axis"],
            "position": pos,
            "busy": busy,
            "pv_position": pvs["pv_position"],
            "pv_state": pvs["pv_state"],
            "pv_busy": pvs["pv_busy"],
        }

    def apply_qmag_current(self, names, currents):
        # cannot be used as se_quadrupoles, because it adds up values,
        # wouldn't be suitable for a scan
        if type(currents) is float:
            currents = np.array([currents])
        if type(names) == str:
            names = [names]
        if len(names) != currents.size:
            print('Error: len(names) != len(currents) in apply_qmag_current(names, currents)')
        for qmag, current in zip(names, currents):
            canonical = self.qmag_alias_to_canonical.get(qmag, qmag)
            pv_des = PV(f'{canonical}:currentWrite')
            target = self.make_safe_float(pv_des.get(), default=np.nan) + float(current)
            pv_des.put(target)
            self._wait_for_magnet_readback(canonical, target)

    def apply_qmag_xyroll(self, names, x_um, y_um, roll_m, wait=True, max_attempts=5, attempt_timeout=30.0, settle_dt=0.5, tol_um=15.0):
        if type(names) == str:
            names = [names]

        nmag = len(names)

        def _expand(values):
            if np.isscalar(values):
                return np.full(nmag, float(values), dtype=float)
            arr = np.asarray(values, dtype=float).reshape(-1)
            if arr.size != nmag:
                raise ValueError(f"Length mismatch in apply_qmag_xyroll: expected {nmag}, got {arr.size}")
            return arr

        xs = _expand(x_um)
        ys = _expand(y_um)
        rolls = _expand(roll_m)

        entries = []
        for name, x_target, y_target, roll_target in zip(names, xs, ys, rolls):
            if name not in self.qmag_pv:
                self.qmag_pv[name] = self._build_qmag_pv_names(name)
            pv = self.qmag_pv[name]
            entries.append({
                "name": name,
                "pv": pv,
                "x_target": float(x_target),
                "y_target": float(y_target),
                "roll_target": float(roll_target),
                "last_x": np.nan,
                "last_y": np.nan,
                "last_roll": np.nan,
                "any_timeout": False,
            })

        for entry in entries:
            print(
                f"Moving {entry['name']}: target x={entry['x_target']:.3f}, "
                f"y={entry['y_target']:.3f}, roll={entry['roll_target']:.6g}"
            )
            self._pv_put(entry["pv"]["pv_set_x"], entry["x_target"])
            self._pv_put(entry["pv"]["pv_set_y"], entry["y_target"])
            self._pv_put(entry["pv"]["pv_set_roll"], entry["roll_target"])

        if not wait:
            return

        time.sleep(settle_dt)

        pending = list(entries)
        attempt = 0
        while pending and attempt < max_attempts:
            attempt += 1

            for entry in pending:
                try:
                    self._pv_put(entry["pv"]["pv_dotrim"], 1)
                except Exception:
                    pass
            time.sleep(settle_dt)

            deadline = time.time() + float(attempt_timeout)
            timed_out_names = set()
            while True:
                busy = []
                for entry in pending:
                    state = int(self._pv_get(entry["pv"]["pv_dotrim"], default=0.0) or 0)
                    if state != 0:
                        busy.append(entry)
                if not busy:
                    break
                if time.time() >= deadline:
                    timed_out_names = {entry["name"] for entry in busy}
                    for entry in busy:
                        entry["any_timeout"] = True
                        try:
                            self._pv_put(entry["pv"]["pv_dotrim"], 0)
                        except Exception:
                            pass
                    time.sleep(1.0)
                    break
                time.sleep(settle_dt)

            next_pending = []
            for entry in pending:
                entry["last_x"] = self._pv_get(entry["pv"]["pv_set_mag_x"], default=0.0)
                entry["last_y"] = self._pv_get(entry["pv"]["pv_set_mag_y"], default=0.0)
                entry["last_roll"] = self._pv_get(entry["pv"]["pv_set_mag_roll"], default=0.0)
                settled = (
                        abs(entry["last_x"] - entry["x_target"]) <= tol_um and
                        abs(entry["last_y"] - entry["y_target"]) <= tol_um and
                        abs(entry["last_roll"] - entry["roll_target"]) <= tol_um
                )
                if settled:
                    print(
                        f"Mover {entry['name']} settled: "
                        f"x={entry['last_x']:.3f}, y={entry['last_y']:.3f}, roll={entry['last_roll']:.6g}"
                    )
                elif entry["name"] not in timed_out_names or attempt < max_attempts:
                    next_pending.append(entry)
            pending = next_pending

        for entry in entries:
            settled = (
                    abs(entry["last_x"] - entry["x_target"]) <= tol_um and
                    abs(entry["last_y"] - entry["y_target"]) <= tol_um and
                    abs(entry["last_roll"] - entry["roll_target"]) <= tol_um
            )
            if not settled:
                if entry["any_timeout"]:
                    print(
                        f"[WARN] {entry['name']}: MAG attempt-timeout after {max_attempts} tries "
                        f"(last x={entry['last_x']}, y={entry['last_y']}, roll={entry['last_roll']}, "
                        f"target x={entry['x_target']}, y={entry['y_target']}, "
                        f"roll={entry['roll_target']}, tol={tol_um})"
                    )
                else:
                    print(
                        f"[WARN] {entry['name']}: did not converge within tolerance "
                        f"(last x={entry['last_x']}, y={entry['last_y']}, roll={entry['last_roll']}, "
                        f"target x={entry['x_target']}, y={entry['y_target']}, "
                        f"roll={entry['roll_target']}, tol={tol_um})"
                    )
            else:
                print(
                    f"Mover {entry['name']} ready for measurement: "
                    f"x={entry['last_x']:.3f}, y={entry['last_y']:.3f}, roll={entry['last_roll']:.6g}"
                )

    def apply_qf1ff_qd0ff(self, qf1_dx, qf1_dy, qf1_droll, qf1_dk1, qd0_dx, qd0_dy, qd0_droll, qd0_dk1):
        q = self.get_quadrupoles(["QF1FF", "QD0FF"])
        names = list(q.get("names", []))
        if not names:
            raise RuntimeError("QF1FF/QD0FF were not found in quadrupoles table.")
        by_name = {name: idx for idx, name in enumerate(names)}

        def _target(name, dx, dy, droll):
            idx = by_name[name]
            x0 = float(q["xdes"][idx]) if np.isfinite(q["xdes"][idx]) else 0.0
            y0 = float(q["ydes"][idx]) if np.isfinite(q["ydes"][idx]) else 0.0
            r0 = float(q["rolldes"][idx]) if np.isfinite(q["rolldes"][idx]) else 0.0
            return x0 + float(dx), y0 + float(dy), r0 + float(droll) * 1e-6

        qf1_x, qf1_y, qf1_r = _target("QF1FF", qf1_dx, qf1_dy, qf1_droll)
        qd0_x, qd0_y, qd0_r = _target("QD0FF", qd0_dx, qd0_dy, qd0_droll)
        self.apply_qmag_xyroll(["QF1FF", "QD0FF"], [qf1_x, qd0_x], [qf1_y, qd0_y], [qf1_r, qd0_r], wait=True)
        self.apply_qmag_current(["QF1FF", "QD0FF"], np.array([float(qf1_dk1), float(qd0_dk1)]))

    def apply_zscan_knobs(self, knob_values: dict, scan_mode_label="30", timeout=30.0, poll=0.05, settle_dt=0.05,
                          base_values=None):
        for knob, value in knob_values.items():
            key = str(knob)
            if key != self.zscan_knob_name:
                raise KeyError(f"Unknown Z scan knob: {knob}")
            delta = float(value)
            if isinstance(base_values, dict) and key in base_values and np.isfinite(float(base_values[key])):
                target = float(base_values[key]) + delta
            else:
                pvs = self._zscan_pvs_for_mode(scan_mode_label)
                pos0 = self._pv_get(pvs["pv_position"], default=np.nan, timeout=0.7)
                if not np.isfinite(pos0):
                    raise RuntimeError(f"Z scan position read failed: {pvs['pv_position']}")
                target = float(pos0) + delta
            self._move_zscan_knob_absolute(
                key,
                target,
                scan_mode_label=scan_mode_label,
                timeout=timeout,
                poll=poll,
                settle_dt=settle_dt,
            )

    def apply_sum_knob(self, dA):
        self.apply_qmag_current(["QS1X", "QS2X"], np.array([float(dA), float(dA)]))

    def get_qf1ff_qd0ff_state(self):
        q = self.get_quadrupoles(["QF1FF", "QD0FF"])
        out = {}
        for idx, name in enumerate(q.get("names", [])):
            out[name] = {
                "dx_um": float(q.get("xact", [np.nan])[idx]),
                "dy_um": float(q.get("yact", [np.nan])[idx]),
                "roll_urad": float(q.get("rollact", [np.nan])[idx]) * 1e6,
                "k1": float(q.get("iact", [np.nan])[idx]),
            }
        return out


    def _normalize_scan_mode_label(self, scan_mode_label):
        text = str(scan_mode_label or "30").strip()
        token = text.replace(" ", "")
        if token in ("2-8", "2/8", "28"):
            return "2-8"
        if token in ("30",):
            return "30"
        if token in ("174",):
            return "174"
        raise ValueError(f"Unsupported scan mode label: {scan_mode_label}")

    def _normalize_linear_knob_name(self, knob_name: str) -> str:
        aliases = {
            "A_X": "Ax",
            "E_X": "Ex",
            "A_Y": "Ay",
            "E_Y": "Ey",
        }
        k = str(knob_name)
        if k in self.linear_matrix:
            return k
        key = k.replace(" ", "").upper()
        return aliases.get(key, k)

    def get_linear_knob_names(self):
        return list(self.linear_matrix.keys())

    def get_nonlinear_knob_names(self):
        return list(self.nonlinear_matrix.keys())

    def get_corrector_knob_names(self):
        return list(self.corrector_knob_pvs.keys())

    def get_zscan_knob_names(self):
        return [self.zscan_knob_name]

    def set_linear_knob(self, knob_name: str, value: float):
        key = self._normalize_linear_knob_name(knob_name)
        if key not in self._linear_knob_values:
            raise KeyError(f"Unknown linear knob: {knob_name}")
        self._linear_knob_values[key] = float(value)
        self.apply_linear_knobs({key: float(value)})

    def set_nonlinear_knob(self, knob_name: str, value: float):
        key = str(knob_name)
        if key not in self._nonlinear_knob_values:
            raise KeyError(f"Unknown nonlinear knob: {knob_name}")
        self._nonlinear_knob_values[key] = float(value)
        self.apply_nonlinear_knobs({key: float(value)})

    def set_corrector_knob(self, knob_name: str, value: float):
        key = str(knob_name)
        if key not in self._corrector_knob_values:
            raise KeyError(f"Unknown corrector knob: {knob_name}")
        self._corrector_knob_values[key] = float(value)
        self.apply_corrector_knobs({key: float(value)})

    def set_zscan_knob(self, knob_name: str, value: float, scan_mode_label="30"):
        key = str(knob_name)
        if key not in self._zscan_knob_values:
            raise KeyError(f"Unknown Z scan knob: {knob_name}")
        self._zscan_knob_values[key] = float(value)
        self.apply_zscan_knobs({key: float(value)}, scan_mode_label=scan_mode_label)

    def _build_linear_deltas(self, knob_values: dict):
        dpos = defaultdict(lambda: [0.0, 0.0])
        for knob, kval in knob_values.items():
            key = self._normalize_linear_knob_name(knob)
            if key not in self.linear_matrix:
                raise KeyError(f"Unknown linear knob: {knob}")
            for mag, (ax, ay) in self.linear_matrix[key].items():
                dpos[mag][0] += ax * float(kval)
                dpos[mag][1] += ay * float(kval)
        return {m: (v[0], v[1]) for m, v in dpos.items()}

    def _build_nonlinear_deltas(self, knob_values: dict):
        dcur = defaultdict(float)
        for knob, kval in knob_values.items():
            key = str(knob)
            if key not in self.nonlinear_matrix:
                raise KeyError(f"Unknown nonlinear knob: {knob}")
            for mag, coeff in self.nonlinear_matrix[key].items():
                dcur[mag] += float(coeff) * float(kval)
        return dict(dcur)

    def _read_linear_magnet_state(self, mag: str):
        return {
            "x_des": float(self._pv_get(f"{mag}:MAG:DES:X", default=np.nan)),
            "y_des": float(self._pv_get(f"{mag}:MAG:DES:Y", default=np.nan)),
            "x_read": float(self._pv_get(f"{mag}:MAG:X", default=np.nan)),
            "y_read": float(self._pv_get(f"{mag}:MAG:Y", default=np.nan)),
            "current_write": float(self._pv_get(f"{mag}:currentWrite", default=np.nan)),
            "current_read": float(self._pv_get(f"{mag}:currentRead", default=np.nan)),
        }

    def _read_current_magnet_state(self, mag: str):
        return {
            "current_write": float(self._pv_get(f"{mag}:currentWrite", default=np.nan)),
            "current_read": float(self._pv_get(f"{mag}:currentRead", default=np.nan)),
            "current_rb": float(self._pv_get(f"{mag}:current", default=np.nan)),
        }

    def _linear_magnets_for_knobs(self, knob_names):
        mags = []
        for knob in knob_names:
            key = self._normalize_linear_knob_name(knob)
            if key not in self.linear_matrix:
                continue
            mags.extend(self.linear_matrix[key].keys())
        return sorted(set(mags))

    def _nonlinear_magnets_for_knobs(self, knob_names):
        mags = []
        for knob in knob_names:
            key = str(knob)
            if key not in self.nonlinear_matrix:
                continue
            mags.extend(self.nonlinear_matrix[key].keys())
        return sorted(set(mags))

    def capture_knob_origin(self, knob_names, scan_mode_label="30"):
        knob_names = list(knob_names or [])
        linear_names = set(self.get_linear_knob_names())
        nonlinear_names = set(self.get_nonlinear_knob_names())
        corrector_names = set(self.get_corrector_knob_names())
        zscan_names = set(self.get_zscan_knob_names())
        selected_linear = [k for k in knob_names if k in linear_names]
        selected_nonlinear = [k for k in knob_names if k in nonlinear_names]
        selected_correctors = [k for k in knob_names if k in corrector_names]
        selected_zscan = [k for k in knob_names if k in zscan_names]
        linear_magnets = self._linear_magnets_for_knobs(selected_linear)
        nonlinear_magnets = self._nonlinear_magnets_for_knobs(selected_nonlinear)

        linear_state = {}
        for mag in linear_magnets:
            linear_state[mag] = self._read_linear_magnet_state(mag)

        nonlinear_state = {}
        for mag in nonlinear_magnets:
            nonlinear_state[mag] = self._read_current_magnet_state(mag)

        corrector_state = {}
        for knob in selected_correctors:
            pvs = self.corrector_knob_pvs[knob]
            corrector_state[knob] = {
                "set_value": float(self._pv_get(pvs["pv_set"], default=np.nan, timeout=0.7)),
                "busy": float(self._pv_get(pvs["pv_busy"], default=np.nan, timeout=0.7)),
                "pv_set": pvs["pv_set"],
                "pv_move": pvs["pv_move"],
                "pv_busy": pvs["pv_busy"],
            }

        zscan_state = {}
        zscan_base_values = {}
        if selected_zscan:
            zscan_pvs = self._zscan_pvs_for_mode(scan_mode_label)
            current_pos = float(self._pv_get(zscan_pvs["pv_position"], default=np.nan, timeout=0.7))
            current_busy = float(self._pv_get(zscan_pvs["pv_busy"], default=np.nan, timeout=0.7))
            for knob in selected_zscan:
                zscan_state[knob] = {
                    "position": current_pos,
                    "busy": current_busy,
                    "mode": zscan_pvs["mode"],
                    "axis": zscan_pvs["axis"],
                    "pv_position": zscan_pvs["pv_position"],
                    "pv_state": zscan_pvs["pv_state"],
                    "pv_busy": zscan_pvs["pv_busy"],
                }
                if np.isfinite(current_pos):
                    zscan_base_values[knob] = current_pos

        linear_base_positions = {}
        for mag in linear_magnets:
            st = linear_state[mag]
            x0 = float(st.get("x_read", np.nan))
            y0 = float(st.get("y_read", np.nan))
            if not np.isfinite(x0):
                x0 = float(st.get("x_des", np.nan))
            if not np.isfinite(y0):
                y0 = float(st.get("y_des", np.nan))
            linear_base_positions[mag] = {"x": x0, "y": y0}

        return {
            "captured_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "params": knob_names,
            "linear_knobs": selected_linear,
            "nonlinear_knobs": selected_nonlinear,
            "corrector_knobs": selected_correctors,
            "zscan_knobs": selected_zscan,
            "scan_mode_label": self._normalize_scan_mode_label(scan_mode_label),
            "linear_magnets": linear_magnets,
            "nonlinear_magnets": nonlinear_magnets,
            "linear_state": linear_state,
            "nonlinear_state": nonlinear_state,
            "corrector_state": corrector_state,
            "zscan_state": zscan_state,
            "linear_base_positions": linear_base_positions,
            "nonlinear_base_currents": {
                mag: float(nonlinear_state[mag]["current_write"])
                for mag in nonlinear_magnets
            },
            "corrector_base_values": {
                knob: float(corrector_state[knob]["set_value"])
                for knob in selected_correctors
            },
            "zscan_base_values": {str(knob): float(val) for knob, val in zscan_base_values.items()},
        }

    def restore_knob_origin(self, origin_state, *, pos_tol=15.0, pos_timeout=30.0, current_tol=0.05, current_timeout=15.0, poll=0.05, settle_dt=0.5, use_trim=True, scan_mode_label=None):
        if not isinstance(origin_state, dict):
            raise ValueError("origin_state must be a dict")

        linear_base_positions = dict(origin_state.get("linear_base_positions", {}) or {})
        nonlinear_base_currents = dict(origin_state.get("nonlinear_base_currents", {}) or {})
        corrector_base_values = dict(origin_state.get("corrector_base_values", {}) or {})
        zscan_base_values = dict(origin_state.get("zscan_base_values", {}) or {})
        zscan_mode_label = str(scan_mode_label or origin_state.get("scan_mode_label", "30"))

        if linear_base_positions:
            raw_target_positions = {
                str(mag): (
                    float(item.get("x", np.nan)),
                    float(item.get("y", np.nan)),
                )
                for mag, item in linear_base_positions.items()
            }
            target_positions = {
                mag: (xt, yt)
                for mag, (xt, yt) in raw_target_positions.items()
                if np.isfinite(xt) and np.isfinite(yt)
            }
            if not target_positions and raw_target_positions:
                raise RuntimeError("No valid saved linear origins were available.")
            for mag, (xt, yt) in target_positions.items():
                self._pv_put(f"{mag}:MAG:DES:X", xt)
                self._pv_put(f"{mag}:MAG:DES:Y", yt)
            time.sleep(settle_dt)
            if use_trim:
                for mag in target_positions.keys():
                    self._pv_put(f"{mag}:TRIM", 1)
                time.sleep(settle_dt)

            deadline = time.time() + float(pos_timeout)
            pending = set(target_positions.keys())
            while pending:
                done = []
                for mag in list(pending):
                    if use_trim:
                        st = self._pv_get(f"{mag}:TRIM", default=1)
                        if not np.isfinite(st) or int(st) != 0:
                            continue
                    xm = self._pv_get(f"{mag}:MAG:X", default=np.nan)
                    ym = self._pv_get(f"{mag}:MAG:Y", default=np.nan)
                    xt, yt = target_positions[mag]
                    if np.isfinite(xm) and np.isfinite(ym) and abs(float(xm) - xt) <= pos_tol and abs(
                            float(ym) - yt) <= pos_tol:
                        done.append(mag)
                for mag in done:
                    pending.remove(mag)
                if not pending:
                    break
                if time.time() >= deadline:
                    raise TimeoutError(f"Restore position timeout. pending={sorted(pending)}")
                time.sleep(poll)

        if nonlinear_base_currents:
            target_currents = {str(mag): float(val) for mag, val in nonlinear_base_currents.items()}
            for mag, target_i in target_currents.items():
                if not np.isfinite(target_i):
                    raise RuntimeError(f"Invalid saved current origin for {mag}")
                self._pv_put(f"{mag}:currentWrite", target_i)

            deadline = time.time() + float(current_timeout)
            pending = set(target_currents.keys())
            while pending:
                done = []
                for mag in list(pending):
                    rb = self._pv_get(f"{mag}:current", default=np.nan)
                    if np.isfinite(rb) and abs(float(rb) - target_currents[mag]) <= current_tol:
                        done.append(mag)
                for mag in done:
                    pending.remove(mag)
                if not pending:
                    break
                if time.time() >= deadline:
                    raise TimeoutError(f"Restore current timeout. pending={sorted(pending)}")
                time.sleep(poll)

        if corrector_base_values:
            self.apply_corrector_knobs(
                {str(knob): float(val) for knob, val in corrector_base_values.items()},
                timeout=current_timeout,
                poll=poll,
            )

        if zscan_base_values:
            target_base = {str(knob): float(val) for knob, val in zscan_base_values.items()}
            for knob, val in target_base.items():
                if not np.isfinite(float(val)):
                    raise RuntimeError(f"Invalid saved Z scan origin for {knob}")
            self.apply_zscan_knobs(
                {knob: 0.0 for knob in target_base.keys()},
                scan_mode_label=zscan_mode_label,
                timeout=current_timeout,
                poll=poll,
                settle_dt=max(0.05, float(settle_dt) * 0.2),
                base_values=target_base,
            )

    def _apply_positions_batch(self, name_to_dxy: dict, *, tol, timeout, poll, settle_dt, use_trim,
                               base_positions=None):
        target = {}
        for mag, (dx, dy) in name_to_dxy.items():
            pv_des_x = f"{mag}:MAG:DES:X"
            pv_des_y = f"{mag}:MAG:DES:Y"
            if isinstance(base_positions, dict) and mag in base_positions:
                x0 = float(base_positions[mag].get("x", np.nan))
                y0 = float(base_positions[mag].get("y", np.nan))
            else:
                x0 = self._pv_get(pv_des_x, default=np.nan)
                y0 = self._pv_get(pv_des_y, default=np.nan)
            if not np.isfinite(x0) or not np.isfinite(y0):
                raise RuntimeError(f"DES read failed: {pv_des_x} / {pv_des_y}")
            xt = float(x0) + float(dx)
            yt = float(y0) + float(dy)
            target[mag] = (xt, yt)
            self._pv_put(pv_des_x, xt)
            self._pv_put(pv_des_y, yt)

        if not target:
            return

        time.sleep(settle_dt)
        if use_trim:
            for mag in target.keys():
                self._pv_put(f"{mag}:TRIM", 1)
            time.sleep(settle_dt)

        deadline = time.time() + float(timeout)
        pending = set(target.keys())
        while pending:
            done = []
            for mag in list(pending):
                if use_trim:
                    st = self._pv_get(f"{mag}:TRIM", default=1)
                    if not np.isfinite(st) or int(st) != 0:
                        continue
                xm = self._pv_get(f"{mag}:MAG:X", default=np.nan)
                ym = self._pv_get(f"{mag}:MAG:Y", default=np.nan)
                if not np.isfinite(xm) or not np.isfinite(ym):
                    continue
                xt, yt = target[mag]
                if abs(float(xm) - xt) <= tol and abs(float(ym) - yt) <= tol:
                    done.append(mag)
            for mag in done:
                pending.remove(mag)
            if not pending:
                return
            if time.time() >= deadline:
                raise TimeoutError(f"Position timeout. pending={sorted(pending)}")
            time.sleep(poll)

    def _apply_currents_batch(self, name_to_di: dict, *, tol, timeout, poll, base_currents=None):
        target = {}
        initial = {}
        for mag, dcur in name_to_di.items():
            pv_set = f"{mag}:currentWrite"
            if isinstance(base_currents, dict) and mag in base_currents:
                cur0 = float(base_currents[mag])
            else:
                cur0 = self._pv_get(pv_set, default=np.nan)
            if not np.isfinite(cur0):
                raise RuntimeError(f"Current setpoint read failed: {pv_set}")
            initial[mag] = float(cur0)
            itarget = float(cur0) + float(dcur)
            target[mag] = itarget
            self._pv_put(pv_set, itarget)

        if not target:
            return

        deadline = time.time() + float(timeout)
        pending = set(target.keys())
        drop_zero_limit = 0.02
        drop_target_limit = max(5.0 * float(tol), 0.1)
        last_readback = {}

        def _is_suspicious_zero_drop(mag: str, rb: float) -> bool:
            init_i = float(initial.get(mag, np.nan))
            target_i = float(target[mag])
            near_zero_target = abs(target_i) < drop_target_limit
            sign_flip_move = (
                    np.isfinite(init_i)
                    and abs(init_i) >= drop_target_limit
                    and abs(target_i) >= drop_target_limit
                    and np.sign(init_i) != np.sign(target_i)
            )
            return bool(
                abs(target_i) >= drop_target_limit
                and abs(float(rb)) <= drop_zero_limit
                and not near_zero_target
                and not sign_flip_move
            )

        while pending:
            done = []
            for mag in list(pending):
                rb = self._pv_get(f"{mag}:current", default=np.nan)
                if not np.isfinite(rb):
                    continue
                last_readback[mag] = float(rb)
                if abs(float(rb) - target[mag]) <= tol:
                    done.append(mag)
            for mag in done:
                pending.remove(mag)
            if not pending:
                return
            if time.time() >= deadline:
                drop_magnets = [
                    mag for mag in sorted(pending)
                    if mag in last_readback and _is_suspicious_zero_drop(mag, last_readback[mag])
                ]
                if drop_magnets:
                    raise CurrentDropToZeroError(
                        "Current readback dropped near 0 A after nonlinear knob apply.",
                        target={mag: float(target[mag]) for mag in drop_magnets},
                        readback={mag: float(last_readback.get(mag, np.nan)) for mag in drop_magnets},
                        magnets=drop_magnets,
                    )
                raise TimeoutError(f"Current timeout. pending={sorted(pending)}")
            time.sleep(poll)

    def apply_linear_knobs(self, knob_values: dict, tol=15, timeout=30.0, poll=0.05, settle_dt=0.5, use_trim=True,
                           base_positions=None):
        dpos = self._build_linear_deltas(knob_values)
        self._apply_positions_batch(
            dpos,
            tol=tol,
            timeout=timeout,
            poll=poll,
            settle_dt=settle_dt,
            use_trim=use_trim,
            base_positions=base_positions,
        )

    def apply_nonlinear_knobs(self, knob_values: dict, tol=0.05, timeout=15.0, poll=0.05, base_currents=None):
        dcur = self._build_nonlinear_deltas(knob_values)
        self._apply_currents_batch(dcur, tol=tol, timeout=timeout, poll=poll, base_currents=base_currents)

    def _move_corrector_knob_absolute(self, knob_name: str, target_value: float, *, timeout=30.0, poll=0.05,
                                      settle_dt=0.1):
        key = str(knob_name)
        if key not in self.corrector_knob_pvs:
            raise KeyError(f"Unknown corrector knob: {knob_name}")

        pvs = self.corrector_knob_pvs[key]
        self._pv_put(pvs["pv_set"], float(target_value))
        time.sleep(float(settle_dt))
        self._pv_put(pvs["pv_move"], 1)

        deadline = time.time() + float(timeout)
        saw_busy = False
        last_busy = np.nan
        while True:
            busy = self._pv_get(pvs["pv_busy"], default=np.nan, timeout=0.7)
            last_busy = busy
            if np.isfinite(busy):
                if int(busy) == 1:
                    saw_busy = True
                elif int(busy) == 0 and saw_busy:
                    self._corrector_knob_values[key] = float(target_value)
                    return
            if time.time() >= deadline:
                raise TimeoutError(
                    f"Corrector move timeout for {key}: target={float(target_value):+.6f}, busy={last_busy}"
                )
            time.sleep(poll)

    def apply_corrector_knobs(self, knob_values: dict, timeout=30.0, poll=0.05, settle_dt=0.1):
        for knob, value in knob_values.items():
            self._move_corrector_knob_absolute(
                knob,
                float(value),
                timeout=timeout,
                poll=poll,
                settle_dt=settle_dt,
            )

    def _move_zscan_knob_absolute(self, knob_name: str, target_value: float, *, scan_mode_label="30", timeout=30.0, poll=0.05, settle_dt=0.05, busy_seen_grace=0.2):
        key = str(knob_name)
        if key != self.zscan_knob_name:
            raise KeyError(f"Unknown Z scan knob: {knob_name}")

        pvs = self._zscan_pvs_for_mode(scan_mode_label)
        self._pv_put(pvs["pv_position"], float(target_value))
        time.sleep(float(settle_dt))
        self._pv_put(pvs["pv_state"], 2)

        t0 = time.time()
        deadline = t0 + float(timeout)
        saw_busy = False
        last_busy = np.nan
        while True:
            busy = self._pv_get(pvs["pv_busy"], default=np.nan, timeout=0.7)
            last_busy = busy
            if np.isfinite(busy):
                if int(busy) == 1:
                    saw_busy = True
                elif int(busy) == 0 and (saw_busy or (time.time() - t0) >= float(busy_seen_grace)):
                    self._zscan_knob_values[key] = float(target_value)
                    return
            if time.time() >= deadline:
                raise TimeoutError(
                    f"Z scan move timeout for {key} in mode={pvs['mode']} ({pvs['axis']}): "
                    f"target={float(target_value):+.6f}, busy={last_busy}"
                )
            time.sleep(poll)

    def _decode_ipbsm_dat(self, raw: bytes):
        offset = 0
        save_amplitude, = struct.unpack_from("<d", raw, offset);
        offset += 8
        save_eamplitude, = struct.unpack_from("<d", raw, offset);
        offset += 8
        save_beamsize, = struct.unpack_from("<d", raw, offset);
        offset += 8
        save_ebeamsize, = struct.unpack_from("<d", raw, offset);
        offset += 8
        save_average, = struct.unpack_from("<d", raw, offset);
        offset += 8
        save_phase, = struct.unpack_from("<d", raw, offset);
        offset += 8
        save_filename = raw[offset:offset + 256].split(b"\x00", 1)[0].decode(errors="ignore")
        offset += 256
        save_ict_average, = struct.unpack_from("<d", raw, offset)
        return {
            "modulation": float(save_amplitude),
            "error": abs(float(save_eamplitude)),
            "beamsize": float(save_beamsize),
            "ebeamsize": float(save_ebeamsize),
            "average": float(save_average),
            "phase": float(save_phase),
            "filename": str(save_filename),
            "ict_average": float(save_ict_average),
        }

    def get_ipbsm_full(self, timeout=300, file_wait=330.0, poll=0.1, start_return_timeout=5.0, start_to_finish_wait=1.0):
        with self._ipbsm_lock:
            try:
                prev_mtime = os.path.getmtime(self.datafile)
            except FileNotFoundError:
                prev_mtime = 0.0

            trigger_pv_name = "IPBSM:FringeScan:RemoteScanStart"
            end_pv_name = "IPBSM:FringeScan:RemoteScanFinish"

            self.pv_trigger.put(1)
            t_start_deadline = time.time() + float(start_return_timeout)
            while True:
                v_start = self._pv_get(trigger_pv_name, default=np.nan, timeout=0.7)
                if np.isfinite(v_start) and int(v_start) == 0:
                    break
                if time.time() >= t_start_deadline:
                    reset_notes = []
                    try:
                        self.pv_trigger.put(0)
                        reset_notes.append(f"{trigger_pv_name} was reset to 0 by optimizer.")
                    except Exception as exc:
                        reset_notes.append(f"Failed to reset {trigger_pv_name} to 0: {exc}")
                    try:
                        self.pv_end.put(0)
                        reset_notes.append(f"{end_pv_name} was set to 0 for cleanup.")
                    except Exception as exc:
                        reset_notes.append(f"Failed to set {end_pv_name} to 0 during cleanup: {exc}")
                    reset_msg = " ".join(reset_notes)
                    raise TimeoutError(
                        f"{trigger_pv_name} stayed at 1 for more than {float(start_return_timeout):.1f}s after trigger. "
                        "Please check that the 'Remote' checkbox on the Fringe Scan GUI is ON. "
                        f"{reset_msg}"
                    )
                time.sleep(poll)

            if float(start_to_finish_wait) > 0.0:
                time.sleep(float(start_to_finish_wait))

            t_deadline = time.time() + float(timeout)
            while True:
                v = self._pv_get(end_pv_name, default=np.nan, timeout=0.7)
                if np.isfinite(v) and int(v) == 1:
                    break
                if time.time() >= t_deadline:
                    raise TimeoutError(
                        f"IPBSM measurement timeout: {end_pv_name} never became 1 within {float(timeout):.1f}s."
                    )
                time.sleep(poll)

            self.pv_end.put(0)

            t_file_deadline = time.time() + float(file_wait)
            try:
                for _ in os.scandir(os.path.dirname(self.datafile)):
                    break
            except Exception:
                pass

            last_seen_mtime = prev_mtime
            while True:
                try:
                    mtime = os.path.getmtime(self.datafile)
                except FileNotFoundError:
                    mtime = 0.0
                last_seen_mtime = mtime
                if mtime > prev_mtime:
                    break
                if time.time() >= t_file_deadline:
                    raise TimeoutError(
                        f"IPBSM datafile not updated: prev_mtime={prev_mtime}, last_mtime={last_seen_mtime}"
                    )
                time.sleep(poll)

            with open(self.datafile, "rb") as f:
                raw = f.read()
            return self._decode_ipbsm_dat(raw)

    def _pv_get(self, pv_name, default=np.nan, timeout=0.7):
        try:
            pv = PV(pv_name)
            if not pv.wait_for_connection(timeout=timeout):
                return default
            val = pv.get(timeout=timeout)
            if val is None:
                return default
            return float(val)
        except Exception:
            return default

    def _pv_put(self, pv_name, value):
        PV(pv_name).put(value)

    def _build_qmag_pv_names(self, mag_name, slacsys_mov_notation=None):
        # slacsys_mov_notation is kept as an explicit hook for future system-specific naming.
        if slacsys_mov_notation is None:
            slacsys_mov_notation = f"{mag_name}:"

        return {
            "pv_read_lvdt1": slacsys_mov_notation + "lvdt1",
            "pv_read_lvdt2": slacsys_mov_notation + "lvdt2",
            "pv_read_lvdt3": slacsys_mov_notation + "lvdt3",
            "pv_read_lvdt_x": slacsys_mov_notation + "x:lvdt",
            "pv_read_lvdt_y": slacsys_mov_notation + "y:lvdt",
            "pv_read_lvdt_tilt": slacsys_mov_notation + "tilt:lvdt",
            "pv_read_lvdt1_raw": slacsys_mov_notation + "lvdt1:raw",
            "pv_read_lvdt2_raw": slacsys_mov_notation + "lvdt2:raw",
            "pv_read_lvdt3_raw": slacsys_mov_notation + "lvdt3:raw",
            "pv_set_x": f"{mag_name}:MAG:DES:X",
            "pv_set_y": f"{mag_name}:MAG:DES:Y",
            "pv_set_roll": f"{mag_name}:MAG:DES:R",
            "pv_read_enc_x": f"{mag_name}:ENC:X",
            "pv_read_enc_y": f"{mag_name}:ENC:Y",
            "pv_read_enc_roll": f"{mag_name}:ENC:R",
            "pv_read_enc1_mrp": f"{mag_name}:ENC1:MRP",
            "pv_read_enc2_mrp": f"{mag_name}:ENC2:MRP",
            "pv_read_enc3_mrp": f"{mag_name}:ENC3:MRP",
            "pv_read_enc1_step": f"{mag_name}:ENC1",
            "pv_read_enc2_step": f"{mag_name}:ENC2",
            "pv_read_enc3_step": f"{mag_name}:ENC3",
            "pv_read_enc1_angl": f"{mag_name}:ENC1:ANG",
            "pv_read_enc2_angl": f"{mag_name}:ENC2:ANG",
            "pv_read_enc3_angl": f"{mag_name}:ENC3:ANG",
            "pv_read_pot_x": f"{mag_name}:POT:X",
            "pv_read_pot_y": f"{mag_name}:POT:Y",
            "pv_read_pot_roll": f"{mag_name}:POT:R",
            "pv_read_pot1_mrp": f"{mag_name}:POT1:MRP",
            "pv_read_pot2_mrp": f"{mag_name}:POT2:MRP",
            "pv_read_pot3_mrp": f"{mag_name}:POT3:MRP",
            "pv_read_pot1_ang": f"{mag_name}:POT1:ANG",
            "pv_read_pot2_ang": f"{mag_name}:POT2:ANG",
            "pv_read_pot3_ang": f"{mag_name}:POT3:ANG",
            "pv_read_pot1_count": f"{mag_name}:POT1:RBV",
            "pv_read_pot2_count": f"{mag_name}:POT2:RBV",
            "pv_read_pot3_count": f"{mag_name}:POT3:RBV",
            "pv_read_pot1_volt": f"{mag_name}:POT1:INV",
            "pv_read_pot2_volt": f"{mag_name}:POT2:INV",
            "pv_read_pot3_volt": f"{mag_name}:POT3:INV",
            "pv_read_pot1_volt_mean": f"{mag_name}:POT1:INV:AVG",
            "pv_read_pot2_volt_mean": f"{mag_name}:POT2:INV:AVG",
            "pv_read_pot3_volt_mean": f"{mag_name}:POT3:INV:AVG",
            "pv_read_pot1_lim_upper": f"{mag_name}:POT1:PLM",
            "pv_read_pot1_lim_lower": f"{mag_name}:POT1:NLM",
            "pv_read_pot2_lim_upper": f"{mag_name}:POT2:PLM",
            "pv_read_pot2_lim_lower": f"{mag_name}:POT2:NLM",
            "pv_read_pot3_lim_upper": f"{mag_name}:POT3:PLM",
            "pv_read_pot3_lim_lower": f"{mag_name}:POT3:NLM",
            "pv_set_mtr1_stop": f"{mag_name}:M1:STOP",
            "pv_set_mtr2_stop": f"{mag_name}:M2:STOP",
            "pv_set_mtr3_stop": f"{mag_name}:M3:STOP",
            "pv_set_mtr1_pos_abs": f"{mag_name}:M1:MOV:ABS",
            "pv_set_mtr2_pos_abs": f"{mag_name}:M2:MOV:ABS",
            "pv_set_mtr3_pos_abs": f"{mag_name}:M3:MOV:ABS",
            "pv_set_mtr1_pos_rel": f"{mag_name}:M1:MOV:REL",
            "pv_set_mtr2_pos_rel": f"{mag_name}:M2:MOV:REL",
            "pv_set_mtr3_pos_rel": f"{mag_name}:M3:MOV:REL",
            "pv_motorsteps1": f"{mag_name}:M1:POS",
            "pv_motorsteps2": f"{mag_name}:M2:POS",
            "pv_motorsteps3": f"{mag_name}:M3:POS",
            "pv_read_mtr1_step": f"{mag_name}:M1:DES:STP",
            "pv_read_mtr2_step": f"{mag_name}:M2:DES:STP",
            "pv_read_mtr3_step": f"{mag_name}:M3:DES:STP",
            "pv_read_mtr1_angl": f"{mag_name}:M1:ANG",
            "pv_read_mtr2_angl": f"{mag_name}:M2:ANG",
            "pv_read_mtr3_angl": f"{mag_name}:M3:ANG",
            "pv_read_mtr1_busy": f"_{mag_name}:M1:BUSY",
            "pv_read_mtr2_busy": f"_{mag_name}:M2:BUSY",
            "pv_read_mtr3_busy": f"_{mag_name}:M3:BUSY",
            "pv_set_mag_move": f"{mag_name}:MOVE",
            "pv_dotrim": f"{mag_name}:TRIM",
            "pv_set_mag_x": f"{mag_name}:MAG:X",
            "pv_set_mag_y": f"{mag_name}:MAG:Y",
            "pv_set_mag_roll": f"{mag_name}:MAG:R",
            "pv_set_mag_stop": f"{mag_name}:STOP",
        }
