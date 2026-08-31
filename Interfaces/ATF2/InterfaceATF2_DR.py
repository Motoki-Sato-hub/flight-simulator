import numpy as np
import time, math,os
from Interfaces.AbstractMachineInterface import AbstractMachineInterface

from epics import PV, ca, caget

'''
TO IMPLEMENT!!!

def get_quadrupoles(self):
    pass
def set_quadrupoles(self):
    pass

'''

class InterfaceATF2_DR(AbstractMachineInterface):
    def get_name(self):
        return 'ATF2_DR'

    def __init__(self, nsamples=1, nominal_intensity=0.1, wfs_intensity=0.125):
        self.nsamples = nsamples
        self.log = print
        # Bpms and correctors in beamline order
        sequence = [
             'MB1R', 'MB2R', 'ZV1R', 'ZH1R', 'MB3R', 'MB4R', 'ZV2R', 'ZH2R',
             'MB5R', 'MB6R', 'ZV3R', 'ZH3R', 'MB7R', 'MB8R', 'ZV4R', 'ZH4R',
             'MB9R', 'MB10R', 'ZV5R', 'ZH5R', 'MB11R', 'MB12R', 'ZV6R', 'ZH6R',
             'MB13R', 'MB14R', 'ZV7R', 'ZH7R', 'MB15R', 'MB16R', 'ZV8R', 'ZH8R',
             'MB17R', 'MB18R', 'ZV9R', 'ZH9R', 'MB19R', 'MBX1', 'MBX2', 'MB21R', 'MB22R', 'ZH10R', 'ZV10R',
             'MB23R', 'ZH11R', 'MB24R','ZV11R', 'MB25R', 'ZH12R', 'MB26R','ZV12R',
             'MB27R', 'ZV13R', 'MB28R', 'ZH13R', 'MB29R', 'ZV14R', 'MB30R', 'ZH14R', 'ZV15R',
             'MB31R', 'ZV16R', 'ZH15R', 'MB32R', 'MB33R', 'ZV17R', 'ZH16R',
             'MB34R', 'ZV18R', 'MB35R', 'ZH17R', 'MB36R', 'MB37R', 'ZV19R', 'ZH18R',
             'MB38R', 'MB39R', 'ZV20R', 'ZH19R', 'MB40R', 'MB41R', 'ZV21R', 'ZH20R',
             'MB42R', 'MB43R', 'ZV22R', 'ZH21R', 'MB44R', 'MB45R', 'ZV23R', 'ZH22R',
             'MB46R', 'MB47R', 'ZV24R', 'ZH23R', 'MB48R', 'MB49R', 'ZV25R', 'ZH24R',
             'MB50R', 'MB51R', 'ZV26R', 'ZH25R', 'MB52R', 'MB53R', 'ZV27R', 'ZH26R',
             'MB54R', 'MB55R', 'ZV28R', 'ZH27R', 'MB56R', 'MB57R', 'ZV29R', 'ZH28R',
             'MB58R', 'MB59R', 'ZV30R', 'ZH29R', 'MB60R', 'MB61R', 'ZV31R', 'ZH30R',
             'MB62R', 'MB63R', 'ZV32R', 'ZH31R', 'MB64R', 'MB65R', 'ZH32R', 'ZV33R',
             'MB66R', 'ZV34R', 'MB67R', 'ZH33R', 'MB68R', 'MB69R', 'ZH34R', 'ZV35R',
             'MB70R', 'ZV36R', 'MB71R', 'ZH35R', 'ZV37R', 'MB72R', 'ZH36R', 'MB73R', 'ZV38R',
             'MB74R', 'ZH37R', 'ZV39R', 'MB76R', 'ZV40R', 'MB77R', 'ZH38R',
             'MB78R', 'ZV41R', 'MB79R', 'ZH42R', 'MB80R', 'ZV43R', 'ZH40R',
             'MB81R','MB82R', 'ZV44R', 'ZH41R', 'MB83R', 'ZV45R', 'MB84R', 'MB85R', 'MB86R',
             'ZV46R', 'ZH43R', 'MB87R', 'MB88R', 'ZV47R', 'ZH44R', 'MB89R', 'MB90R',
             'ZV48R', 'ZH45R', 'MB91R', 'MB92R', 'ZV49R', 'ZH46R', 'MB93R', 'MB94R',
             'ZV50R', 'ZH47R', 'MB95R', 'MB96R', 'ZV51R', 'ZH48R', 'MB97R', 'MB98R'
        ]

        # ATF2' BPMs Epics names
        # https://atf.kek.jp/atfbin/view/ATF/EPICS_DATABASE
        monitors = [
            'MB1R', 'MB2R', 'MB3R', 'MB4R', 'MB5R', 'MB6R', 'MB7R', 'MB8R',
            'MB9R', 'MB10R', 'MB11R', 'MB12R', 'MB13R', 'MB14R', 'MB15R', 'MB16R', 'MB17R', 'MB18R', 'MB19R','MB21R', 'MB22R',
            'MB9R', 'MB10R', 'MB11R', 'MB12R', 'MB13R', 'MB14R', 'MB15R', 'MB16R', 'MB17R', 'MB18R', 'MB19R', 'MB21R', 'MB22R',
            'MB23R', 'MB24R', 'MB25R', 'MB26R', 'MB27R', 'MB28R', 'MB29R', 'MB30R', 'MB31R', 'MB32R', 'MB33R',
            'MB34R', 'MB35R', 'MB36R', 'MB37R', 'MB38R', 'MB39R', 'MB40R', 'MB41R', 'MB42R', 'MB43R', 'MB44R', 'MB45R',
            'MB46R', 'MB47R', 'MB48R', 'MB49R', 'MB50R', 'MB51R', 'MB52R', 'MB53R', 'MB54R', 'MB55R', 'MB56R', 'MB57R',
            'MB58R', 'MB59R', 'MB60R', 'MB61R', 'MB62R', 'MB63R', 'MB64R', 'MB65R', 'MB66R', 'MB67R', 'MB68R', 'MB69R',
            'MB70R', 'MB71R', 'MB72R', 'MB73R', 'MB74R', 'MB76R', 'MB77R', 'MB78R', 'MB79R', 'MB80R', 'MB81R', 'MB82R',
            'MB83R', 'MB84R', 'MB85R', 'MB86R', 'MB87R', 'MB88R', 'MB89R', 'MB90R', 'MB91R', 'MB92R', 'MB93R', 'MB94R',
            'MB95R', 'MB96R', 'MB97R', 'MB98R'
        ]
        # monitors = ['MB1R', 'MB3R', 'MB4R', 'MB5R', 'MB7R', 'MB8R', 'MB9R', 'MB10R', 'MB11R', 'MB12R', 'MB13R', 'MB14R',
        #             'MB15R', 'MB16R', 'MB18R', 'MB21R', 'MB22R', 'MB24R', 'MB25R', 'MB26R', 'MB27R', 'MB28R', 'MB29R',
        #             'MB31R', 'MB33R', 'MB34R', 'MB35R', 'MB36R', 'MB37R', 'MB38R', 'MB40R', 'MB41R', 'MB42R', 'MB43R',
        #             'MB44R', 'MB45R', 'MB48R', 'MB49R', 'MB50R', 'MB51R', 'MB52R', 'MB53R', 'MB54R', 'MB55R', 'MB56R',
        #             'MB58R', 'MB59R', 'MB61R', 'MB62R', 'MB63R', 'MB64R', 'MB65R', 'MB66R', 'MB67R', 'MB68R', 'MB69R',
        #             'MB71R']

        self.sextupoles = [
            "SF1R.1", "SD1R.1", "SBH1R.1", "SBH1R.2", "SQF1R.1", "SQF1R.2",
            "SF1R.2", "SD1R.2", "SBH1R.3", "SBH1R.4", "SQF1R.3", "SQF1R.4",
            "SF1R.3", "SD1R.3", "SBH1R.5", "SBH1R.6", "SQF1R.5", "SQF1R.6",
            "SF1R.4", "SD1R.4", "SBH1R.7", "SBH1R.8", "SQF1R.7", "SQF1R.8",
            "SF1R.5", "SD1R.5", "SBH1R.9", "SBH1R.10", "SQF1R.9", "SQF1R.10",
            "SF1R.6", "SD1R.6", "SBH1R.11", "SBH1R.12", "SQF1R.11", "SQF1R.12",
            "SF1R.7", "SD1R.7", "SBH1R.13", "SBH1R.14", "SQF1R.13", "SQF1R.14",
            "SF1R.8", "SD1R.8", "SBH1R.15", "SBH1R.16", "SQM3R.1", "SQM3R.2",
            "SF1R.9", "SD1R.9", "SBH1R.17", "SBH1R.18", "SBH1R.19", "SBH1R.20",
            "SF1R.10", "SD1R.10", "SBH1R.21", "SBH1R.22", "SQM22R.1", "SQM22R.2",
            "SF1R.11", "SD1R.11", "SBH1R.23", "SBH1R.24", "SQF1R.15", "SQF1R.16",
            "SF1R.12", "SD1R.12", "SBH1R.25", "SBH1R.26", "SQF1R.17", "SQF1R.18",
            "SF1R.13", "SD1R.13", "SBH1R.27", "SBH1R.28", "SQF1R.19", "SQF1R.20",
            "SF1R.14", "SD1R.14", "SBH1R.29", "SBH1R.30", "SQF1R.21", "SQF1R.22",
            "SF1R.15", "SD1R.15", "SBH1R.31", "SBH1R.32", "SQF1R.23", "SQF1R.24",
            "SF1R.16", "SD1R.16", "SBH1R.33", "SBH1R.34", "SQF1R.25", "SQF1R.26",
            "SF1R.17", "SD1R.17", "SBH1R.35", "SBH1R.36", "SQF1R.27", "SQF1R.28",
            "SF1R.18", "SD1R.18", "SBH1R.37", "SBH1R.38", "SQF1R.29", "SQF1R.30",
            "SF1R.19", "SD1R.19", "SBH1R.39", "SBH1R.40", "SQF1R.31", "SQF1R.32",
            "SF1R.20", "SD1R.20", "SBH1R.41", "SBH1R.42", "SQF1R.33", "SQF1R.34",
            "SF1R.21", "SD1R.21", "SBH1R.43", "SBH1R.44", "SQF1R.35", "SQF1R.36",
            "SF1R.22", "SD1R.22", "SBH1R.45", "SBH1R.46", "SQF1R.37", "SQF1R.38",
            "SF1R.23", "SD1R.23", "SBH1R.47", "SBH1R.48", "SQF1R.39", "SQF1R.40",
            "SF1R.24", "SD1R.24", "SBH1R.49", "SBH1R.50", "SQF1R.41", "SQF1R.42",
            "SF1R.25", "SD1R.25", "SBH1R.51", "SBH1R.52", "SQM3R.3", "SQM3R.4",
            "SF1R.26", "SD1R.26", "SBH1R.53", "SBH1R.54", "SBH1R.55", "SBH1R.56",
            "SF1R.27", "SD1R.27", "SBH1R.57", "SBH1R.58", "SQM22R.3", "SQM22R.4",
            "SF1R.28", "SD1R.28", "SBH1R.59", "SBH1R.60", "SQF1R.43", "SQF1R.44",
            "SF1R.29", "SD1R.29", "SBH1R.61", "SBH1R.62", "SQF1R.45", "SQF1R.46",
            "SF1R.30", "SD1R.30", "SBH1R.63", "SBH1R.64", "SQF1R.47", "SQF1R.48",
            "SF1R.31", "SD1R.31", "SBH1R.65", "SBH1R.66", "SQF1R.49", "SQF1R.50",
            "SF1R.32", "SD1R.32", "SBH1R.67", "SBH1R.68", "SQF1R.51", "SQF1R.52",
            "SF1R.33", "SD1R.33", "SBH1R.69", "SBH1R.70", "SQF1R.53", "SQF1R.54",
            "SF1R.34", "SD1R.34", "SBH1R.71", "SBH1R.72", "SQF1R.55", "SQF1R.56",
        ]

        self.quadrupoles = [
            "QF2R.1", "QF1R.1", "QF2R.2", "QF1R.2", "QF2R.3", "QF1R.3",
            "QF2R.4", "QF1R.4", "QF2R.5", "QF1R.5", "QF2R.6", "QF1R.6",
            "QF2R.7", "QF1R.7", "QM1R.1", "QM2R.1", "QM3R.1", "QM4R.1",
            "QM5R.1", "QM6R.1", "QM7R.1", "QM8R.1", "QM9R.1", "QM10R.1",
            "QM11R.1", "QM12R.1", "QM13R.1", "QM14R.1", "QM15R.1", "QM16R.1",
            "QM17R.1", "QM18R.1", "QM19R.1", "QM20R.1", "QM21R.1", "QM22R.1",
            "QM23R.1", "QF1R.8", "QF2R.8", "QF1R.9", "QF2R.9", "QF1R.10",
            "QF2R.10", "QF1R.11", "QF2R.11", "QF1R.12", "QF2R.12", "QF1R.13",
            "QF2R.13", "QF1R.14", "QF2R.14", "QF1R.15", "QF2R.15", "QF1R.16",
            "QF2R.16", "QF1R.17", "QF2R.17", "QF1R.18", "QF2R.18", "QF1R.19",
            "QF2R.19", "QF1R.20", "QF2R.20", "QF1R.21", "QM1R.2", "QM2R.2",
            "QM3R.2", "QM4R.2", "QM5R.2", "QM6R.2", "QM7R.2", "QM8R.2",
            "QM9R.2", "QM10R.2", "QM11R.2", "QM12R.2", "QM13R.2", "QM14R.2",
            "QM15R.2", "QM16R.2", "QM17R.2", "QM18R.2", "QM19R.2", "QM20R.2",
            "QM21R.2", "QM22R.2", "QM23R.2", "QF1R.22", "QF2R.21", "QF1R.23",
            "QF2R.22", "QF1R.24", "QF2R.23", "QF1R.25", "QF2R.24", "QF1R.26",
            "QF2R.25", "QF1R.27", "QF2R.26", "QF1R.28",
        ]
        # Use list comprehension to filter out strings starting with 'Z' or 'z'
        monitors_from_sequence = [string for string in sequence if not string.lower().startswith('z')]
        # Check if the bpms in the config files are known to Epics
        bpm_ok = all(bpm in monitors for bpm in monitors_from_sequence)
        if not bpm_ok:
            bpms_unknown = [bpm for bpm in monitors_from_sequence if bpm not in monitors]
            self.log(f'Unknown bpms {bpms_unknown} removed from list')
        # Only retain BPMs in config file which are known by Epics
        sequence_filtered = [element for element in sequence if (element in monitors) or element.lower().startswith('z')]
        # Subset of BPMs and correctors from the config file
        self.sequence = sequence_filtered
        self.bpms = [string for string in self.sequence if not string.lower().startswith('z')]
        self.corrs = [string for string in self.sequence if string.lower().startswith('z')]
        self.screens = []
        # Bunch current monitors
        self.ict_names = [
            'gun:GUNcharge', 'l0:L0charge', 'linacbt:LNEcharge', 'linacbt:BTMcharge',
            'ext:EXTcharge', 'linacbt:BTEcharge', 'BIM:DR:nparticles', 'BIM:IP:nparticles'
        ]
        self.nominal_laser_intensity = nominal_intensity
        self.test_laser_intensity = wfs_intensity
        #self.laser_intensity = PV('RFGun:LasetIntensity1:Read').get()
        self.twiss_path = os.path.join(os.path.dirname(__file__), "DR_ATF2", "ATF_DR_twiss_file.tws")

        self.xsr_pvs = {
            "proj_v_sigma_ave": "XSR:PROJ_V_SIGMA_AVE",
            "proj_v_mean_ave": "XSR:PROJ_V_MEAN_AVE",
            "proj_h_sigma_ave": "XSR:PROJ_H_SIGMA_AVE",
            "proj_h_mean_ave": "XSR:PROJ_H_MEAN_AVE",
        }
        self.arc_dispersion_pv = "MONITOR:DR:ARCDISPERSION"
        self.laser_intensity1 = PV('RFGun:LaserIntensity1:Read').get()
        self.laser_intensity2 = PV('RFGun:LaserIntensity2:Read').get()
        self.machine_name = "ATF2"

    def get_beam_factors(self):
        # TO BE REPLACED WITH A PV OF REAL BEAM ENERGY
        Pref = 1.2999999e3
        gamma_rel = np.sqrt((Pref / 0.51099895) ** 2 + 1.0)
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

    def change_energy(self):
        PV('RAMP:CONTROL_ON_SW').put(1)
        self._wait_for_pv_readback('RAMP:CONTROL_ON_SW', 1)
        ### delta_freq MUST MATCH :MI2: to EPICS --> means "MINUS2"
        delta_freq = +4 # kHz
        # PV('RAMP:MI2:ONOFF_SW').put(1)
        PV('RAMP:PL4:ONOFF_SW').put(1)
        self._wait_for_pv_readback('RAMP:PL4:ONOFF_SW', 1)
        DR_freq = 714e3; # 714 MHz in kHz
        DR_momentum_compaction = 2.1e-3
        dP_P = -delta_freq / DR_freq / DR_momentum_compaction
        return dP_P

    def reset_energy(self):
        PV('RAMP:CONTROL_OFF_SW').put(0)
        self._wait_for_pv_readback('RAMP:CONTROL_ON_SW', 0)

    def change_intensity(self, intensity=0.1):
        print(f'Changing laser intensity to {intensity}...')
        laser_intensity1 = 10000 * float(intensity) / self.laser_intensity2
        PV('RFGun:LaserIntensity1:Write').put(laser_intensity1)
        self._wait_for_pv_readback('RFGun:LaserIntensity1:Read', laser_intensity1)
        return self

    def reset_intensity(self):
        print('Resetting laser intensity...')
        PV('RFGun:LaserIntensity1:Write').put(self.laser_intensity1)
        self._wait_for_pv_readback('RFGun:LaserIntensity1:Read', self.laser_intensity1)
        return self

    def get_beam_settings(self):
        settings = {"energy": {}, "intensity": {}}
        for name, pv_name in (
            ("ramp_control", "RAMP:CONTROL_ON_SW"),
            ("ramp_pl4", "RAMP:PL4:ONOFF_SW"),
        ):
            settings["energy"][name] = float(PV(pv_name).get())
            settings["intensity"]["laser_intensity1"] = float(PV('RFGun:LaserIntensity1:Read').get())
        return settings

    def restore_beam_settings(self, settings):
        settings = settings or {}
        energy = settings.get("energy", {})
        for name, pv_name in (("ramp_control", "RAMP:CONTROL_ON_SW"),
                              ("ramp_pl4", "RAMP:PL4:ONOFF_SW")):
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
        name_to_index = {string: index for index, string in enumerate(self.sequence)}
        return [name_to_index.get(name, np.nan) for name in names]

    def log_messages(self, console):
        self.log = console or print

    def get_target_dispersion(self, names=None): # for DR too
        if names is None:
            names = self.get_bpms()["names"]
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
        self.log("Reading ict's...")
        charge = []
        for ict in self.ict_names:
            pv = PV(f'{ict}')
            if 0: # Reading the icts is time consuming and unnecessary for SysID and BBA
                charge.append(pv.get())
            else:
                charge.append(1.0)
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
        self.log("Reading correctors' strengths...")
        bdes, bact = [], []
        for corrector in self.corrs:
            pv_des = PV(f'{corrector}:currentWrite')
            pv_act = PV(f'{corrector}:currentRead')
            bdes.append(pv_des.get())
            bact.append(pv_act.get())

        correctors = {
            "names": self.corrs,
            "bdes": np.array(bdes),
            "bact": np.array(bact),
        }

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
        p = PV('DR:monitors')
        x, y, tmit = [], [], []
        max_attempts = max(1, int(self.nsamples) * 3)
        sample = 0
        attempts = 0

        while sample < int(self.nsamples) and attempts < max_attempts:
            attempts += 1
            try:
                self.log(f'Sample = {sample} attempt to read bpms = {attempts}')
                data_returned = p.get()
                if data_returned is None:
                    raise RuntimeError('DR:monitors returned None')
                data_returned = np.asarray(data_returned, dtype=float)
                if data_returned.size == 0:
                    raise RuntimeError('DR:monitors returned an empty array')
                a = data_returned.reshape((-1, 10))

                if names is None:
                    names = self.bpms

                if isinstance(names, str):
                    names = [names]

                x_, y_, tmit_ = [], [], []
                for name in names:
                    # it can also be written as
                    # re.findall(r'\d+', "MB12R") --> ['12']
                    k = int(name[2:-1])-1
                    charge = a[k,3]
                    status = a[k,0]
                    valid = (status == 1) & np.isfinite(charge) & (charge > 0)
                    if valid:
                        x_.append(a[k,1])
                        y_.append(a[k,2])
                        tmit_.append(charge)
                    else:
                        x_.append(np.nan)
                        y_.append(np.nan)
                        tmit_.append(np.nan)

                x.append(np.asarray(x_))
                y.append(np.asarray(y_))
                tmit.append(np.asarray(tmit_))

                self.log(f'Interface::get_bpms() = {x}')
                sample += 1
                time.sleep(1)

            except Exception as e:
                self.log(f'An error occurred while reading DR BPM sample: {e}')
                time.sleep(1)

        x = np.asarray(x)
        y = np.asarray(y)
        tmit = np.asarray(tmit)

        bpms = {
            "names": names,
            "x": np.vstack(x) / 1e3, # mm
            "y": np.vstack(y) / 1e3, # mm
            "tmit": np.vstack(tmit),
        }

        return bpms

    # def get_sextupoles(self, names=None):
    #     if names is None:
    #         names = self.sextupoles
    #     if isinstance(names, str):
    #         names = [names]
    #     bdes, bact = [], []
    #
    #     for name in names:
    #         bdes.append(self.pv_get(f"{name}:currentWrite"))
    #         bact.append(self.pv_get(f"{name}:currentRead"))
    #     return {
    #         "names": names,
    #         "bdes": np.array(bdes, dtype=float),
    #         "bact": np.array(bact, dtype=float),
    #     }

    def set_sextupoles(self, names, values):
        if isinstance(names, str):
            names = [names]
        values = np.asarray(values, dtype=float).reshape(-1)
        if len(names) != len(values):
            raise ValueError("len(names) != len(values) in set_sextupoles")
        for name, value in zip(names, values):
            self._pv_put(f"{name}:currentWrite", float(value))
            self._wait_for_corrector_readback(name, value)

    def _wait_for_corrector_readback(self, corrector, target, tolerance=1e-4, timeout=1.0, poll_interval=0.05):
        readback_pv = PV(f'{corrector}:currentRead')
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
        self.log(
            f'Warning: {corrector}:currentRead did not reach target {float(target):.6g} '
            f'within {timeout:.2f}s. Last readback = {last_value:.6g}'
        )
        return False

    def set_correctors(self, names, corr_vals):
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        if len(names) != len(corr_vals):
            self.log('Error: len(names) != len(corr_vals) in set_correctors(names, corr_vals)')
            return
        for corrector, corr_val in zip(names, corr_vals):
            pv_des = PV(f'{corrector}:currentWrite')
            pv_des.put(corr_val)
            self._wait_for_corrector_readback(corrector, corr_val)

    def vary_correctors(self, names, corr_vals):
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        if len(names) != len(corr_vals):
            self.log('Error: len(names) != len(corr_vals) in vary_correctors(names, corr_vals)')
            return
        for corrector, corr_val in zip(names, corr_vals):
            pv_des = PV(f'{corrector}:currentWrite')
            curr_val = self.make_safe_float(pv_des.get(), default=np.nan)
            if not np.isfinite(curr_val):
                self.log(f'Warning: could not read current value for {corrector}; skipping.')
                continue
            target = curr_val + float(corr_val)
            pv_des.put(target)
            self._wait_for_corrector_readback(corrector, target)

    '''
    SATO-SAN'S METHODS:
    '''

    def pv_get(self, pv_name, default=np.nan):
        try:
            val = PV(pv_name).get()
            if val is None:
                return default
            return float(val)
        except Exception:
            return default

    def pv_put(self, pv_name, value):
        PV(pv_name).put(float(value))

    def pv_put_many(self, pv_to_value):
        for pv_name, value in pv_to_value.items():
            self.pv_put(pv_name, value)

    def read_current(self, pv_write: str, pv_read: str = None, quantize_phase: bool = True) -> float:
        pv_name = pv_read or pv_write
        return self.pv_get(pv_name)

    def get_skew_series_names(self):
        return ["SD1R", "SF1R"]

    def get_skew_names(self, series: str):
        series_name = str(series).strip().upper()
        if series_name not in self.get_skew_series_names():
            raise ValueError(f"Unknown skew series: {series}")
        return [f"{series_name}.{idx}" for idx in range(1, 35)]

    def _skew_pv_names(self, name: str):
        series, idx = str(name).strip().upper().split(".")
        base = f"{series}_{int(idx)}"
        return {
            "set": f"{base}:internalCurrentWrite",
            "read": f"{base}:current",
        }

    def get_skews(self, series: str):
        names, currents = [], []
        for name in self.get_skew_names(series):
            pvs = self._skew_pv_names(name)
            names.append(name)
            currents.append(self.pv_get(pvs["read"]))
        arr = np.asarray(currents, dtype=float)
        return {
            "series": str(series).strip().upper(),
            "names": names,
            "current": arr,
            "abs_current": np.abs(arr),
        }

    def get_top_skews_by_abs_current(self, series: str, count: int = 5):
        skews = self.get_skews(series)
        vals = np.asarray(skews["current"], dtype=float)
        valid = np.isfinite(vals)
        idx = np.where(valid)[0]
        idx = idx[np.argsort(-np.abs(vals[idx]))]
        top = idx[:max(int(count), 0)]
        return {
            "series": skews["series"],
            "names": [skews["names"][i] for i in top],
            "current": vals[top],
            "abs_current": np.abs(vals[top]),
        }

    def set_skew_current(self, name: str, value: float):
        pvs = self._skew_pv_names(name)
        self.pv_put(pvs["set"], value)

    def set_skew_currents(self, name_to_value):
        for name, value in name_to_value.items():
            self.set_skew_current(name, value)

    def get_skew_current(self, name: str):
        pvs = self._skew_pv_names(name)
        return self.pv_get(pvs["read"])

    def get_xsr(self):
        return {name: self.pv_get(pv) for name, pv in self.xsr_pvs.items()}

    def get_arc_dispersion(self):
        return self.pv_get(self.arc_dispersion_pv)
