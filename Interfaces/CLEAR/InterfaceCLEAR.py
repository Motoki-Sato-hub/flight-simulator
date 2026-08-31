from Interfaces.CLEAR.Setup_files.CLEAR_BPM_getHV import baseline_correct, find_peak, threshold_integral, plot_peak, plot_integral, change_inverted_bpm_polarity
from Interfaces.CLEAR.InterfaceCLEAR_RFTrack import InterfaceCLEAR_RFTrack
import sys, time, math, os, json
from scipy.integrate import trapezoid
from enum import Enum
from scipy.optimize import curve_fit
import numpy as np
try:
    import pyda
    import pyda_japc
except ImportError:
    pyda = None
    pyda_japc = None

try:
    from Interfaces.CLEAR import config
except ImportError:
    import config
from Interfaces.AbstractMachineInterface import AbstractMachineInterface

class BPMsMode(Enum):
    peak = "peak"
    baseline_peak = "baseline_peak"
    integral = "integral"
    integral_window = "integral_window"
    integral_threshold = "integral_threshold"

scaling_factors = {
    "BPM0530": {"H": -0.403, "V": -0.396},
    "BPM0595": {"H": -0.351, "V": -0.374},
    "BPM0690": {"H": -0.400, "V": -0.410},
    "BPM0820": {"H": -0.346, "V": -0.392},
    "BPM0890": {"H": -0.391, "V": -0.417},
}

class CLEAR_real_machine(AbstractMachineInterface):
    def get_name(self):
        return 'CLEAR'

    def __init__(self, nsamples=3, bg_shots=10.0 ):
        self.screen_backgrounds = {}
        self.chosen_ict = "CA.BCMTHZ/Acquisition#charge"
        self.steps_readback_position = 0.0
        self.bpm_mode = BPMsMode.integral_threshold
        self.nsamples = nsamples
        self.electronmass = 0.51099895 # MeV/c^2
        self.Pref = 198 # MeV/c
        self.machine_name = "CLEAR"
        self.tracking_interface = InterfaceCLEAR_RFTrack()
        self.energy_param = [
            'CA.BEAM/Acquisition#momentum',
            'CA.BEAM/Acquisition#energy',
        ]
        self.is_simulation = False
        self.context_acquisition = "SCT.USER.SETUP"
        self.context_empty = ""
        self.log = print
        self.client = pyda.SimpleClient(provider=pyda_japc.JapcProvider())
        self.rf_phase_nominal = 125 # degrees
        self.rf_phase_test = 145 # degrees

        # Bpms and correctors in beamline order
        sequence = [
            'CA.DHG0130', 'CA.DVG0130', #'CA.BPC0220',
            'CA.DHG0225', 'CA.DVG0225', #'CA.BPC0240',
            'CA.DHG0245', 'CA.DVG0245', #'CA.BPC0260',
            'CA.DHG0265', #'CA.BPC0310',
            'CA.DHG0320', 'CA.DVG0320', #'CA.SDV0340',
            'CA.QFD0350', 'CA.QDD0355', 'CA.QFD0360',
            'CA.DHG0385', 'CA.DVG0385',
            'CA.BTV0390L', 'CA.BTV0390H',
            'CA.QFD0510', 'CA.QDD0515', 'CA.QFD0520',
            'CA.BPM0530', 'CA.DHJ0540', 'CA.DVJ0540',
            'CA.DHJ0590', 'CA.DVJ0590', 'CA.BPM0595',
            'CA.BTV0620', 'CA.BPM0690', 'CA.DHJ0710', 'CA.DVJ0710',
            'CA.BTV0730', 'CA.QFD0760', 'CA.QDD0765', 'CA.QFD0770',
            'CA.DHJ0780', 'CA.DVJ0780', 'CA.BTV0810', 'CA.BPM0820',
            'CA.DHJ0840', 'CA.DVJ0840',
            'CA.QDD0870', 'CA.QFD0880', 'CA.BPM0890', 'CA.BTV0910',
        ]

        monitors = [
                     # 'CA.BPC0220', 'CA.BPC0240', 'CA.BPC0260',
                     # 'CA.BPC0310',
                    'CA.BPM0530', 'CA.BPM0595',
                     'CA.BPM0690', 'CA.BPM0820', 'CA.BPM0890',
        ]

        correctors = [
            'CA.DHG0130', 'CA.DVG0130',
            'CA.DHG0225', 'CA.DVG0225',
            'CA.DHG0245', 'CA.DVG0245',
            'CA.DHG0265', 'CA.DVG0265',
            'CA.DHG0320', 'CA.DVG0320',
            'CA.DHG0385', 'CA.DVG0385',
            'CA.DHJ0540', 'CA.DVJ0540',
            'CA.DHJ0590', 'CA.DVJ0590',
            'CA.DHJ0710', 'CA.DVJ0710',
            'CA.DHJ0780', 'CA.DVJ0780',
            'CA.DHJ0840', 'CA.DVJ0840',
            # 'CA.SDV0340',
        ]

        self.screen_status_params = {
            "CA.BTV0390L": "CA.BTV0390_CAS.BTV0420/OPSettingSystem1#positionChannel1",
            "CA.BTV0390H": "CA.BTV0390_CAS.BTV0420/OPSettingSystem1#positionChannel1",
            "CA.BTV0620":  "CAS.BTV0440_CA.BTV0620/OPSettingSystem2#positionChannel5",
            "CA.BTV0730":  "CA.BTV0730_CA.BTV0800/OPSettingSystem1#positionChannel1",
            "CA.BTV0810":  "CA.BTV0805_CA.BTV0810/OPSettingSystem2#positionChannel5",
            "CA.BTV0910":  "CA.BTV0910_CAS.BTV0930/OPSettingSystem1#positionChannel1",
        }

        self.corrector_set_params = {name: f'{name}/SettingPPM' for name in correctors}
        self.corrector_get_params = {name: f'{name}/Acquisition' for name in correctors}

        self.sextupoles = []
        self.quadrupoles = list(config.quad_names)
        monitors_from_sequence = [element for element in sequence if element in monitors]
        bpm_ok = all(bpm in monitors for bpm in monitors_from_sequence)
        if not bpm_ok:
            bpms_unknown = [bpm for bpm in monitors_from_sequence if bpm not in monitors]
            self.log(f'Unknown bpms {bpms_unknown} removed from list')

        sequence_filtered = [
            element for element in sequence
            if element in monitors or element in correctors or element in self.quadrupoles or element in config.cameras
        ]
        self.sequence = sequence_filtered
        self.bpms = [element for element in self.sequence if element in monitors]
        self.corrs = [element for element in self.sequence if element in correctors]
        self.screen_names = list(config.cameras.keys())
        self.screens = self.screen_names
        self.screen_config = config.cameras # it consists also screens that we should not use!

        # ["CA.BTV0125","CA.BTV0215","CA.BTV0235", "CA.BTV0390", "CAS.BTV0420", "CAS.BTV0440", "CAS.CAMVESPER1",
        #  "CA.BTV0545","CA.BTV0620", "CA.BTV0730","CA.BTV0800", "CA.BTV0805", "CA.BTV0810" , "CA.BTV0875",
        #  "CA.BTV0910", "CAS.BTV0930","CA.CAMAIR1", "CA.CAMAIR2", "CA.CAMAIR3","CA.CAMAIR4", "CS.BTV0120", "CS.BTV0305",
        #  "CS.BTV0420","CS.BTV0520", "CS.BTVVAC1","CS.CAMVAC1","CS.CAMAIR1","CS.CAMAIR2","CS.CAMAIR3","CS.CAMAIR4",
        #  "PHIN.BTV01","PHIN.BTV.Spectro","PHIN.VCAT"]

        self.bpm_indexes = [index for index, string in enumerate(sequence) if string in self.bpms]

        # Bunch current monitors
        self.ict_names = [
            'CA.BCMGUN/Acquisition#charge',
            'CA.BCMVESPER/Acquisition#charge',
            'CA.BCM0395/Acquisition#charge',
            'CS.BCM0620/Acquisition#charge',
            'CA.BCMTHZ/Acquisition#charge',
            'CA.BCMTHZ2/Acquisition#charge',
        ]

        self.quadrupoles = list(config.quad_names)
        self.quad_set_params = dict(zip(config.quad_names, config.current_set_params))
        self.quad_get_params = dict(zip(config.quad_names, config.current_get_params))
        self.quad_status_params = dict(zip(config.quad_names, config.current_status_params))
        self.cam_props = self.CamList() # Load camera configuration from assets/cameras.json
        self.camList = list(self.cam_props.keys())
        self.lattice = self.tracking_interface.lattice
        self.start = self.tracking_interface.start
        self.end = self.tracking_interface.end
        self.bg_shots = int(bg_shots)

    def CamList(self):
        _JSON_PATH = os.path.join(os.path.dirname(__file__), 'cameras.json')
        """Return the full device configuration dict (keyed by BTV device name)."""
        with open(_JSON_PATH) as f:
            data = json.load(f)
        return data['devices']

    def _get_charge_reference(self):
        value = self.get_icts(names = self.chosen_ict)
        return value["charge"]

    def _give_elements_to_show_beamline(self, quad_selected):
        start_quad_element_name = quad_selected
        return start_quad_element_name

    def _get_elements_positions_show_beamline(self, names=None):
        if isinstance(names, str):
            names = [names]
        all_names = [name for name in self.tracking_interface.get_sequence() if names is None or name in names]

        return {
            "names": all_names,
            "S": np.array([self._get_tracking_element(name).get_S("entrance") for name in all_names], dtype=float),
        }

    def get_beam_factors(self):
        pref = self.Pref
        # data = self.client.get("CA.BEAM/Acquisition").data
        # for field in ("momentum", "energy"):
        #     value = self.make_safe_float(data.get(field), default=np.nan)
        #     if np.isfinite(value) and value > 0:
        #         pref = value
        #         break
        pref = 195
        gamma_rel = np.sqrt((pref / self.electronmass) ** 2 + 1.0)
        beta_rel = np.sqrt(1.0 - 1.0 / gamma_rel ** 2)
        beta_gamma = gamma_rel * beta_rel
        return gamma_rel, beta_rel, beta_gamma

    def _get_twiss_s_positions(self, names):
        positions = []
        for name in names:
            element = self._get_tracking_element(name)
            positions.append(float(element.get_S("exit")) if element is not None else np.nan)
        return positions

    def _get_tracking_element(self, name):
        """Return the model element, including the BTV0390L/H machine aliases."""
        candidates = (name, str(name).rstrip("LH"))
        for candidate in dict.fromkeys(candidates):
            try:
                element = self.lattice[candidate]
            except Exception:
                continue
            if isinstance(element, list):
                if not element:
                    continue
                element = element[-1]
            return element
        return None

    @staticmethod
    def make_safe_float(value, default=np.nan):  # so even if japc address returns none, empty array or whatever, interface still works
        try:
            if value is None:
                return float(default)
            arr = np.asarray(value)
            if arr.size == 0:
                return float(default)
            return float(arr.flat[0])
        except Exception:
            return float(default)

    def change_energy(self):
        energy_readback = self.client.get('CK.LL-MKS11/Setting').data['PhaseSh_SP'] # changes value globally
        self.log(f"Value before changing energy: {energy_readback}")
        new_energy = self.rf_phase_test
        self.client.set('CK.LL-MKS11/Setting', data = {"PhaseSh_SP" : new_energy})
        self.log(f"Value after changing energy: {new_energy}")
        self._wait_for_japc_readback('CK.LL-MKS11/Setting', 'PhaseSh_SP', new_energy)
        after_energy_change = self.client.get('CK.LL-MKS11/Setting').data['PhaseSh_SP']
        print(after_energy_change)
        return after_energy_change

    def reset_energy(self):
        print(f"Resetting energy to {self.rf_phase_nominal}...")
        self.client.set('CK.LL-MKS11/Setting', data = {"PhaseSh_SP" : self.rf_phase_nominal})
        self._wait_for_japc_readback('CK.LL-MKS11/Setting', 'PhaseSh_SP', self.rf_phase_nominal)
        print(f"Energy has been reset to {self.rf_phase_nominal}...")
        after_energy_reset = self.client.get('CK.LL-MKS11/Setting').data['PhaseSh_SP']
        print(after_energy_reset)

    def change_intensity(self):
        self.steps_readback_position = self.client.get('CO.TOWB.102.UVATT2/Setting').data['position']
        self.steps_readback_position_min = self.client.get('CO.TOWB.102.UVATT2/Setting').data['position_min']
        self.steps_readback_position_max =self.client.get('CO.TOWB.102.UVATT2/Setting').data['position_max']
        print(f'Changing intensity to ...')
        nominal_settings_steps = self.steps_readback_position
        N_steps = 1000 # to be verified!
        new_laser_settings = nominal_settings_steps + N_steps
        self.log(f"The new laser settings will be set to {new_laser_settings}... Nominal value is {self.steps_readback_position}.")
        self.client.set('CO.TOWB.102.UVATT2/Setting', data={"position": new_laser_settings})
        self._wait_for_japc_readback('CO.TOWB.102.UVATT2/Setting', 'position', new_laser_settings)
        self.log(f"The new laser settings has been set to {new_laser_settings}. Nominal value was {self.steps_readback_position}.")
        after_intensity_change = self.client.get('CO.TOWB.102.UVATT2/Setting').data['position']
        self.log(f"Read after change of intensity:", after_intensity_change)
        return self

    def reset_intensity(self):
        print(f"Resetting intensity to {self.steps_readback_position}...")
        self.client.set('CO.TOWB.102.UVATT2/Setting', data = {"position" : self.steps_readback_position})
        self._wait_for_japc_readback('CO.TOWB.102.UVATT2/Setting', 'position', self.steps_readback_position)
        print(f"Intensity steps has been reset to {self.steps_readback_position}...")
        after_intensity_reset = self.client.get('CO.TOWB.102.UVATT2/Setting').data['position']
        self.log(f"Read after reset of intensity:", after_intensity_reset)

    def get_beam_settings(self):
        settings = {"energy": {}, "intensity": {}}
        settings["energy"]["mks11_phase"] = self.make_safe_float(self.client.get('CK.LL-MKS11/Setting').data['PhaseSh_SP'])
        settings["intensity"]["uvatt2_position"] = self.make_safe_float(self.client.get('CO.TOWB.102.UVATT2/Setting').data['position'])
        return settings

    def restore_beam_settings(self, settings):
        settings = settings or {}
        phase = self.make_safe_float(settings.get("energy", {}).get("mks11_phase"))
        if np.isfinite(phase):
            self.client.set('CK.LL-MKS11/Setting', data={"PhaseSh_SP": phase})
            self._wait_for_japc_readback('CK.LL-MKS11/Setting', 'PhaseSh_SP', phase)

        position = self.make_safe_float(settings.get("intensity", {}).get("uvatt2_position"))
        if np.isfinite(position):
            self.client.set('CO.TOWB.102.UVATT2/Setting', data={"position": position})
            self._wait_for_japc_readback('CO.TOWB.102.UVATT2/Setting', 'position', position)
        return True

    def get_sequence(self):
        return self.sequence

    def get_hcorrectors_names(self):
        return [name for name in self.corrs if ("DHG" in name) or ("DHJ" in name)]

    def get_vcorrectors_names(self):
        return [name for name in self.corrs if ("DVG" in name) or ("DVJ" in name) or ("SDV" in name)]

    def get_elements_indices(self, names):
        if isinstance(names, str):
            names = [names]
        name_to_index = {string: index for index, string in enumerate(self.sequence)}
        return [name_to_index.get(name, np.nan) for name in names]

    def _read_screen_status(self, screen_name):
        try:
            address = self.screen_status_params[screen_name]
            property_address, field = address.rsplit("#", 1)
            value = self.client.get(property_address, context=self.context_empty).data[field]
            return self.make_safe_float(value)
        except Exception as exc:
            self.log(f"Could not read screen status for {screen_name}: {exc}")
            return np.nan

    @staticmethod
    def _camera_frame_id(camera_data):
        timestamp = float(np.asarray(camera_data["imageTimeStamp"]).ravel()[0])
        image = np.asarray(camera_data["image2D"])
        return ("image", image.shape, image.dtype.str, hash(image.tobytes()))

    def _acquire_screen_data(self, screen_name, previous_frame_id=None, timeout=5.0):
        japc_camera = self.screen_config.get(screen_name, {}).get("japc_name", screen_name.rstrip("LH"))
        camera_config = self.screen_config.get(screen_name, {})
        selector = camera_config.get("japc_selector", self.context_empty)
        deadline = time.perf_counter() + timeout
        while time.perf_counter() < deadline:
            camera_data = self.client.get(f"{japc_camera}.DigiCam/LastImage", context=selector).data
            if previous_frame_id is None or self._camera_frame_id(camera_data) != previous_frame_id:
                return camera_data
            time.sleep(0.1)
        self.log(f"No new camera frame for {screen_name} within {timeout:.1f} s.")
        return None

    def _get_screen_pixel_calibration(self, screen_name):
        camera_config = self.screen_config.get(screen_name, {})
        japc_camera = camera_config.get("japc_name", screen_name.rstrip("LH"))
        selector = camera_config.get("japc_selector", self.context_empty)
        calibration = self.client.get(f"{japc_camera}.DigiCam/CalibrationSetting", context=selector).data
        hpixel = self.make_safe_float(calibration.get("pixelCalSet1"))
        vpixel = self.make_safe_float(calibration.get("pixelCalSet2"))
        return hpixel, vpixel

    def _orient_screen_image(self, screen_name, image, hpixel, vpixel):
        camera_config = self.screen_config.get(screen_name, {})
        japc_camera = camera_config.get("japc_name", screen_name.rstrip("LH"))
        camera_properties = self.cam_props.get(japc_camera, {})
        oriented = np.asarray(image, dtype=float)
        if camera_properties.get("flip_hor", 0): oriented = np.fliplr(oriented)
        if camera_properties.get("flip_ver", 0): oriented = np.flipud(oriented)
        rotate = int(camera_properties.get("rotate", 0)) % 4
        if rotate:
            oriented = np.rot90(oriented, rotate)
            if rotate % 2:
                hpixel, vpixel = vpixel, hpixel
        return oriented, hpixel, vpixel

    def get_icts(self, names=None):
        #BCM_THZ = ('CA.BCMTHZ/Acquisition#charge', 'SCT.USER.SETUP')
        self.log("Reading ict's...")
        if names is None:
            names = self.ict_names
        if isinstance(names, str):
            names = [names]
        charge = []
        for name in names:
            property_address, field = name.rsplit("#", 1)
            try:
                value = self.client.get(property_address, context=self.context_acquisition).data[field]
            except Exception:
                value = np.nan
            charge.append(self.make_safe_float(value))
        return {
            "names": list(names),
            "charge": np.asarray(charge, dtype=float),
        }

    def get_correctors(self, names=None):
        #{corr_name}/SettingPPM#current
        self.log("Reading correctors' strengths...")
        selected_names = self.corrs if names is None else ([names] if isinstance(names, str) else list(names))

        bdes, bact = [], []
        for corrector in selected_names:
            setting_data = self.client.get(self.corrector_set_params[corrector],context = self.context_empty).data
            acquisition_data = self.client.get(self.corrector_get_params[corrector], context = self.context_acquisition).data
            bdes.append(setting_data['current'])
            bact.append(acquisition_data['currentAverage'])

        return {
            "names": list(selected_names),
            "bdes": np.asarray(bdes, dtype=float),
            "bact": np.asarray(bact, dtype=float),
        }

    def get_bpms(self, names=None):
        window = (260, 360)
        self.log('Reading bpms...')
        selected_names = self.bpms if names is None else ([names] if isinstance(names, str) else list(names))
        x, y, tmit = [], [], []
        mode = self.bpm_mode
        for sample in range(self.nsamples):
            self.log(f'Sample = {sample}')
            x_sample, y_sample, tmit_sample = [], [], []
            for bpm in selected_names:
                hsamples = self.client.get(f"{bpm}H-SA/SamplesFromTrigger", context = self.context_acquisition).data
                vsamples = self.client.get(f"{bpm}V-SA/SamplesFromTrigger", context=self.context_acquisition).data
                ssamples = self.client.get(f"{bpm}S-SA/SamplesFromTrigger", context=self.context_acquisition).data

                H_samples = change_inverted_bpm_polarity(np.asarray(hsamples["samples"], dtype=float).ravel(), bpm)
                V_samples = change_inverted_bpm_polarity(np.asarray(vsamples["samples"], dtype=float).ravel(), bpm)
                S_samples = np.asarray(ssamples["samples"], dtype=float).ravel()

                H_b_samples = baseline_correct(H_samples)
                V_b_samples = baseline_correct(V_samples)
                S_b_samples = baseline_correct(S_samples)

                s_sum = np.sum(S_samples[320:330])

                if mode == BPMsMode.peak: # Find peak (largest magnitude, keeping the sign)
                    H, H_idx = find_peak(H_samples)
                    V, V_idx = find_peak(V_samples)
                    # plot_peak([H_samples, V_samples], [H_idx, V_idx], [H, V], ["H", "V"], BPM)

                elif mode == BPMsMode.baseline_peak: # Find peak on baseline corrected signal (largest magnitude, keeping the sign)
                    H, H_idx = find_peak(H_b_samples)
                    V, V_idx = find_peak(V_b_samples)
                    # plot_peak([H_b_samples, V_b_samples], [H_idx, V_idx], [H, V], ["H", "V"], BPM, )

                elif mode == BPMsMode.integral: # Integration of full BPM signal with baseline correction
                    H = trapezoid(H_b_samples)
                    V = trapezoid(V_b_samples)
                    # plot_integral(signals=[H_b_samples, V_b_samples], integrals=[H, V], labels=["H", "V"], BPM=BPM)

                elif mode == BPMsMode.integral_window: # Integration of fixed window BPM signal with baseline correction
                    window_start, window_end = window
                    H = trapezoid(H_b_samples[window_start:window_end])
                    V = trapezoid(V_b_samples[window_start:window_end])
                    # plot_integral(signals=[H_b_samples, V_b_samples], integrals=[H, V], labels=["H", "V"], BPM=BPM, starts=[window_start, window_start], ends=[window_end, window_end], )

                elif mode == BPMsMode.integral_threshold: # Integration between 5% of the peak threshold region with baseline correction
                    H, H_start, H_end, H_peak_idx = threshold_integral(H_b_samples)
                    V, V_start, V_end, V_peak_idx = threshold_integral(V_b_samples)
                    S, S_start, S_end, S_peak_idx = threshold_integral(S_b_samples)
                    bpm_key = bpm[3:] if bpm.startswith("CA.") else bpm

                    H = (H / S) / scaling_factors[bpm_key]["H"]
                    V = (V / S) / scaling_factors[bpm_key]["V"]

                else:
                    raise ValueError(
                        f"Unknown mode '{mode}'. " "Choose from: 'peak', 'baseline_peak', " "'integral', 'integral_window', 'integral_threshold'.")

                x_sample.append(H)
                y_sample.append(V)
                tmit_sample.append(s_sum)
            x.append(x_sample)
            y.append(y_sample)
            tmit.append(tmit_sample)
            time.sleep(1)

        return {
            "names": list(selected_names),
            "x": np.asarray(x, dtype=float),
            "y": np.asarray(y, dtype=float),
            "tmit": np.asarray(tmit, dtype=float),
        }

    def _wait_for_japc_readback(self, property_address, field, target, *, context=None, tolerance=5e-3, timeout=10.0):
        def read_value():
            data = self.client.get(property_address, context=context).data
            return self.make_safe_float(data.get(field), default=np.nan)
        return self._wait_for_readback(read_value, target, description=f"{property_address}#{field}", tolerance=tolerance, timeout=timeout)

    def _wait_for_corrector_readbacks(self, names, targets, tolerance=5e-3, timeout=10.0, poll_interval=0.05):
        targets = {name: float(target) for name, target in zip(names, targets)}
        pending = set(targets)
        last_values = {name: np.nan for name in targets}
        deadline = time.perf_counter() + timeout

        while pending and time.perf_counter() < deadline:
            for corrector in tuple(pending):
                try:
                    data = self.client.get(self.corrector_get_params[corrector], context=self.context_acquisition).data
                    value = self.make_safe_float(data.get("currentAverage"))
                except Exception:
                    value = np.nan
                last_values[corrector] = value
                if np.isfinite(value) and abs(value - targets[corrector]) <= tolerance:
                    pending.remove(corrector)
            if pending:
                time.sleep(poll_interval)

        for corrector in pending:
            self.log(
                f"Warning: {corrector} did not reach target {targets[corrector]:.6g} "
                f"within {timeout:.2f}s. Last readback = {last_values[corrector]:.6g}"
            )
        return not pending

    def _wait_for_quadrupole_readback(self, quadrupole, target, tolerance=5e-3, timeout=10.0):
        readback_param = self.quad_get_params[quadrupole]
        property_address, field = readback_param.rsplit("#", 1)
        return self._wait_for_japc_readback(property_address, field, target, context=self.context_acquisition, tolerance=tolerance, timeout=timeout)

    def set_correctors(self, names, corr_vals):
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        if len(names) != len(corr_vals):
            self.log('Error: len(names) != len(corr_vals) in set_correctors(names, corr_vals)')
            return
        for corrector, corr_val in zip(names, corr_vals):
            self.client.set(self.corrector_set_params[corrector], data={'current': corr_val})
        return self._wait_for_corrector_readbacks(names, corr_vals)

    def vary_correctors(self, names, corr_vals):
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        if len(names) != len(corr_vals):
            self.log('Error: len(names) != len(corr_vals) in vary_correctors(names, corr_vals)')
            return
        current = self.get_correctors(names)['bdes']
        target = current + np.asarray(corr_vals, dtype=float)
        self.set_correctors(names, target)

    def get_quadrupoles(self, names=None):
        if names is None:
            names = self.quadrupoles
        if isinstance(names, str):
            names = [names]

        ides = []
        iact = []

        for quadrupole in names:
            set_address = self.quad_set_params[quadrupole]
            set_property, set_field = set_address.rsplit("#", 1)

            get_address = self.quad_get_params[quadrupole]
            get_property, get_field = get_address.rsplit("#", 1)

            try:
                set_value = self.client.get(set_property, context=self.context_empty).data[set_field]
            except Exception:
                set_value = np.nan
            try:
                get_value = self.client.get(get_property, context=self.context_acquisition).data[get_field]
            except Exception:
                get_value = np.nan

            ides.append(self.make_safe_float(set_value))
            iact.append(self.make_safe_float(get_value))

        ides = np.asarray(ides, dtype=float)
        iact = np.asarray(iact, dtype=float)
        try:
            bdes = np.asarray([self.current_to_k1l(name, current) for name, current in zip(names, ides)], dtype=float)
            bact = np.asarray([self.current_to_k1l(name, current) for name, current in zip(names, iact)], dtype=float)
            self.update_tracking_model_with_japc_readback(names=names, nominal_bdes_value=bdes, nominal_bact_value=bact)

        except Exception as exc:
            self.log(f"CLEAR current-to-K1L conversion failed for {names}: {exc}")
            bdes = np.full(len(names), np.nan, dtype=float)
            bact = np.full(len(names), np.nan, dtype=float)

        return {
            "names": list(names),
            "bdes": bdes,
            "bact": bact,
            "ides": ides,
            "iact": iact,
        }

    def set_quadrupoles(self, names, k1l_values):
        if isinstance(names, str):
            names = [names]
        if not isinstance(k1l_values, (list, tuple, np.ndarray)):
            k1l_values = [k1l_values]
        if len(names) != len(k1l_values):
            raise ValueError(f"len(names)={len(names)} != len(k1l_values)={len(k1l_values)}")

        for quadrupole, k1l in zip(names, k1l_values):
            current_A = self.k1l_to_current(quadrupole, k1l)
            address = self.quad_set_params[quadrupole]
            property_address, field = address.rsplit("#", 1)
            self.client.set(property_address, data={field: current_A})
            self._wait_for_quadrupole_readback(quadrupole, current_A)

    def _get_screen_movement_info(self, screen_name):
        btv_key = screen_name.rstrip("LH")
        cam = self.cam_props.get(btv_key)
        if cam is None: raise RuntimeError(f"Camera {btv_key} not found")
        if not bool(cam.get("screenInstalled")): raise RuntimeError(f"Camera {btv_key} not installed")

        screen_props = {
            "btv_key": btv_key,
            "btvdevice": cam.get('controlDeviceName'),  # None when not configured
            "ctrl_type" : cam.get('controlDeviceType'),
            "ctrl_fields" : cam.get('controlDeviceFields', {}),
            "screen_mover_device" : cam.get('screenMoverDevice'), # None for most cameras
            "screen_mover_type" : cam.get('screenMoverType'),
            "screen_mover_fields" : cam.get('screenMoverFields') or {}
        }

        screen_props["system"] = int(screen_props["ctrl_fields"].get("system", 1))
        screen_props["has_custom_screen_mover"] = isinstance(screen_props["screen_mover_device"], str)

        if screen_props["system"] == 1:
            screen_props["set_prop"] = 'OPSettingSystem1'
            screen_props["get_prop"] = 'ExpertSettingDCSystem1'
            screen_props["get_set_field"] = 'positionChannel1'
            screen_props["description_field"] = 'dcm1DriverNames'

        elif screen_props["system"] == 2:
            screen_props["set_prop"] = 'OPSettingSystem2'
            screen_props["get_prop"] = 'ExpertSettingDCSystem2'
            screen_props["get_set_field"] = 'positionChannel5'
            screen_props["description_field"] = 'dcm3DriverNames'

        else:
            screen_props["set_prop"] = screen_props["get_prop"] = screen_props["field"] = screen_props["description_field"] = None

        return screen_props

    def insert_screen(self, screen_name):
        info = self._get_screen_movement_info(screen_name)
        current_screen_inout_status = self.client.get(f"{info['btvdevice']}/{info['set_prop']}").data[info['get_set_field']] # 0 or not == 0 means screen is out, whatever else means IN
        if current_screen_inout_status.value == 0:
            self.log(f"Inserting {screen_name}...")
            self.client.set(f"{info['btvdevice']}/{info['set_prop']}", data={f"{info['get_set_field']}": 1}) # 1, meaning INSERT the screen
            reached_target = self._wait_for_screen_target_position(screen_name, 1)
            if not reached_target: raise RuntimeError(f"Screen {screen_name} was not inserted within time.")
            self.log(f"Inserted {screen_name}!")
            current_screen_inout_status2 = self.client.get(f"{info['btvdevice']}/{info['set_prop']}").data[info['get_set_field']]
            print("Current Screen Inout Status:", current_screen_inout_status2)
        else:
            print(current_screen_inout_status.value)
            self.log(f"Screen {screen_name} already inserted")
            return
            return

    def extract_screen(self, screen_name):
        screen_name = screen_name.rstrip("LH")
        info = self._get_screen_movement_info(screen_name)
        current_screen_inout_status = self.client.get(f"{info['btvdevice']}/{info['set_prop']}").data[info['get_set_field']]  # 0 or not == 0 means screen is out, whatever else means IN
        if current_screen_inout_status.value == 0:
            self.log(f"Screen {screen_name} already extracted")
            return
        else:
            self.log(f"Extracting {screen_name}...")
            self.client.set(f"{info['btvdevice']}/{info['set_prop']}", data={f"{info['get_set_field']}": 0})  # 0, meaning EXTRACT the screen
            reached_target = self._wait_for_screen_target_position(screen_name, 0)
            if not reached_target: raise RuntimeError(f"Screen {screen_name} was not extracted within time.")
            self.log(f"Extracted {screen_name}!")

    def acquire_screen_background(self, screen_name, frames = None):
        if frames is None: frames = self.bg_shots
        previous_data = self._acquire_screen_data(screen_name)
        previous_frame_id = self._camera_frame_id(previous_data) if previous_data is not None else None
        self.extract_screen(screen_name)
        background_frames = []
        for frame in range(frames):
            self.log(f"Acquiring {frame}/{frames} background frames...")
            camera_data = self._acquire_screen_data(screen_name, previous_frame_id)
            if camera_data is None:
                continue
            previous_frame_id = self._camera_frame_id(camera_data)
            image = np.asarray(camera_data['image2D'], dtype=float)
            background_frames.append(image)
            self.log(f"Acquired {frame}/{frames} background frames...")
            # camgui.bgAction samples its background frames at 100 ms intervals.
            time.sleep(0.1)
        if not background_frames:
            raise RuntimeError(f"No background frames available for {screen_name}")
        self.log(f"Acquired {frames} background frames. Calculating the mean...")
        bg_img = np.mean(np.stack(background_frames, axis=0), axis=0)
        self.log(f"Mean calculated.")
        self.screen_backgrounds[screen_name] = bg_img
        return bg_img

    def _wait_for_screen_target_position(self, screen_name, target, timeout=10.0, poll_interval=0.05):
        info = self._get_screen_movement_info(screen_name)
        t0 = time.perf_counter()
        while time.perf_counter() - t0 < timeout:
            current_screen_inout_status = self.client.get(f"{info['btvdevice']}/{info['set_prop']}").data[info['get_set_field']]  # 0 or not == 0 means screen is out, whatever else means IN
            if current_screen_inout_status.value == 0 and target==0: return True
            if current_screen_inout_status.value > 0 and target >0: return True
            time.sleep(poll_interval)
        self.log(
            f'Warning: {screen_name} did not reach target state = {target:.6g} '
            f'within {timeout:.2f}s. Last readback = {current_screen_inout_status.value:.6g}'
        )
        return False

    def acquire_screen_image(self, screen_name):
        previous_data = self._acquire_screen_data(screen_name)
        previous_frame_id = self._camera_frame_id(previous_data) if previous_data is not None else None
        self.insert_screen(screen_name)
        camera_data = self._acquire_screen_data(screen_name, previous_frame_id)
        if camera_data is None: raise RuntimeError(f"No camera data available for {screen_name}")
        beam_img = np.asarray(camera_data['image2D'], dtype=float)
        bg_img = self.screen_backgrounds[screen_name]
        subtracted_img = beam_img - bg_img
        subtracted_img[~np.isfinite(subtracted_img)] = 0.0
        return subtracted_img, bg_img.copy(), beam_img

    def get_screens(self, names=None):
        background_images, beam_images = [], []
        self.log(f'Reading screens {names}...')
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
        hedges_all = []
        vedges_all = []
        inout_list = []

        for screen_name in selected_names:
            raw_hpixel, raw_vpixel = self._get_screen_pixel_calibration(screen_name)
            hpixel, vpixel = raw_hpixel, raw_vpixel

            for attempt in range(3):  # re-acquire the image if the Gaussian fit came back as a nan
                try:
                    if screen_name not in self.screen_backgrounds:
                        self.log(f"Acquiring background image for {screen_name}.")
                        self.acquire_screen_background(screen_name, frames = 10)
                    subtracted_img, bg_img, beam_img = self.acquire_screen_image(screen_name)
                    subtracted_img, hpixel, vpixel = self._orient_screen_image(screen_name, subtracted_img, raw_hpixel, raw_vpixel)
                    bg_img, _, _ = self._orient_screen_image(screen_name, bg_img, raw_hpixel, raw_vpixel)
                    beam_img, _, _ = self._orient_screen_image(screen_name, beam_img, raw_hpixel, raw_vpixel)
                    x_mean, y_mean, sigx, sigy, total, img, hedges, vedges = self._screen_data_from_image(subtracted_img, hpixel, vpixel)
                    if np.isfinite(sigx) and np.isfinite(sigy):
                        break

                except Exception as e:
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

            status = self._read_screen_status(screen_name) # is screen inserted or extracted?
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
            "names": list(selected_names),
            "hpixel": np.asarray(hpixel_list, dtype=float), # mm
            "vpixel": np.asarray(vpixel_list, dtype=float), # mm
            "x": np.asarray(xb_list, dtype=float),
            "y": np.asarray(yb_list, dtype=float),
            "sigx": np.asarray(sigx_list, dtype=float),
            "sigy": np.asarray(sigy_list, dtype=float),
            "sum": np.asarray(sum_list, dtype=float),
            "hedges": hedges_all, #imagePositionSet1
            "vedges": vedges_all,
            "background_images": background_images,
            "beam_images": beam_images,
            "images": images, # subtracted images
            "S": np.asarray(s_positions, dtype=float),
            "inout": np.asarray(inout_list, dtype=float),
        }

        return screens

    @staticmethod
    def _quad_sign(name):
        if "QFD" in name:
            return 1.0
        if "QDD" in name:
            return -1.0
        raise ValueError(f"Unknown interface quadrupole: {name}")

    def current_to_k1l(self, name, current_A, pref_mev_c=None):
        current_A = float(current_A)
        pref_mev_c = self.tracking_interface.Pref if pref_mev_c is None else float(pref_mev_c)
        if not np.isfinite(current_A) or not np.isfinite(pref_mev_c) or pref_mev_c <= 0:
            return np.nan

        length = float(self._get_tracking_element(name).get_length())
        k1 = self.tracking_interface.get_Quad_K_from_I(current_A, length, pref_mev_c)
        return self._quad_sign(name) * k1 * length

    def k1l_to_current(self, name, k1l, pref_mev_c=None):
        k1l = float(k1l)
        pref_mev_c = self.tracking_interface.Pref if pref_mev_c is None else float(pref_mev_c)
        if not np.isfinite(k1l) or not np.isfinite(pref_mev_c) or pref_mev_c <= 0:
            raise ValueError(f"Invalid K1L or reference momentum for {name}")
        a = float(self.tracking_interface.get_ITF(0.0))
        b = float(a - self.tracking_interface.get_ITF(1.0))
        target = self._quad_sign(name) * k1l * pref_mev_c / 299.8
        discriminant = a * a - 4.0 * b * target
        if discriminant < 0.0:
            raise ValueError(f"K1L={k1l:.6g} is outside the CLEAR calibration range for {name}")
        solutions = ((a - np.sqrt(discriminant)) / (2.0 * b),
                 (a + np.sqrt(discriminant)) / (2.0 * b))
        return float(min(solutions, key=abs))

    def update_tracking_model_with_japc_readback(self, nominal_bdes_value=None, nominal_bact_value=None, names=None):
        if isinstance(names, str):
            names = [names]
        bact = np.asarray(nominal_bact_value, dtype=float)
        self.tracking_interface.set_quadrupoles(names, bact)

    def predict_emittance_scan_response(self, *args, **kwargs):
        return self.tracking_interface.predict_emittance_scan_response(*args, **kwargs)

    def get_R_matrix_scan(self, *args, **kwargs):
        return self.tracking_interface.get_R_matrix_scan(*args, **kwargs)

    def get_phase_space_transport_to_screens(self, *args, **kwargs):
        return self.tracking_interface.get_phase_space_transport_to_screens(*args, **kwargs)

    def get_twiss_evolution(self, *args, **kwargs):
        return self.tracking_interface.get_twiss_evolution(*args, **kwargs)

    @staticmethod
    def _gaussian(x, amplitude, center, sigma, offset):
        return amplitude * np.exp(-((x - center) ** 2) / (2.0 * sigma ** 2)) + offset

    @classmethod
    def _fit_projection(cls, axis, projection):
        baseline_subtracted = projection - np.min(projection)
        normalisation = baseline_subtracted.sum()
        centre = baseline_subtracted.dot(axis) / normalisation
        rms = np.sqrt(baseline_subtracted.dot((axis - centre) ** 2) / normalisation)
        fitted, _ = curve_fit(cls._gaussian, axis, projection, p0=[np.max(projection) - np.min(projection), centre, rms, np.min(projection)]) # [amplitude, center, sigma, background]
        return fitted[1], abs(fitted[2]) # sigx, sigy

    def _screen_data_from_image(self, image, hpixel, vpixel): # better be subtracted!
        img = np.flipud(np.asarray(image, dtype=float).copy())
        img[~np.isfinite(img)] = 0.0
        ny, nx = img.shape # e.g. ny = no. of rows, nx = no. of columns
        x_pixels_positions = hpixel * np.linspace(-nx / 2, nx / 2, nx)
        y_pixels_positions = vpixel * np.linspace(-ny / 2, ny / 2, ny)
        summed_intensity = np.sum(img)
        proj_x = np.mean(img, axis=0)
        proj_y = np.mean(img, axis=1)

        x_mean_positions, sigx = self._fit_projection(x_pixels_positions, proj_x)
        y_mean_positions, sigy = self._fit_projection(y_pixels_positions, proj_y)

        hedges = np.r_[x_pixels_positions - hpixel / 2, x_pixels_positions[-1] + hpixel / 2]
        vedges = np.r_[y_pixels_positions - vpixel / 2, y_pixels_positions[-1] + vpixel / 2]

        return x_mean_positions, y_mean_positions, sigx, sigy, summed_intensity, img, hedges, vedges

    def log_messages(self, console):
        self.log = console or print

    def _read_screen_setting(self, screen_name):
        japc_camera = self.screen_config.get(screen_name, {}).get('japc_name', screen_name.rstrip('LH'))
        try:
            return self.client.get(f'{japc_camera}.DigiCam/Setting', context = self.context_empty).data
        except Exception as e:
            print(e)
            return None
