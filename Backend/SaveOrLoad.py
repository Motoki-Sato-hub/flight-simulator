import os, pickle, json, re
from datetime import datetime
try:
    from PyQt6.QtWidgets import (
        QDialog, QVBoxLayout, QDialogButtonBox,
        QRadioButton, QLabel,QFileDialog, QMessageBox
        )
    from PyQt6.QtCore import QEvent, Qt
except ImportError:
    from PyQt5.QtWidgets import (
        QDialog, QVBoxLayout, QDialogButtonBox,
        QRadioButton, QLabel,QFileDialog, QMessageBox
        )
    from PyQt5.QtCore import QEvent, Qt
import numpy as np
from Backend.State import State


def _match_exactly():
    return Qt.MatchFlag.MatchExactly if hasattr(Qt, "MatchFlag") else Qt.MatchExactly

class SaveOrLoad():

    def _refactor_names_order(self,elements_list, selected_names):
        state=getattr(self,"initial_state",None)

        if state is not None:
            if hasattr(self,"correctors_list") and elements_list is self.correctors_list:
                actuator_mode = getattr(getattr(self,"actuator_mode",None),"name", "Kicker")
                if actuator_mode == "QM" and hasattr(self, "qm_corrs"):
                    ref_list = [str(x) for x in self.qm_corrs]
                else:
                    ref_list = [str(x) for x in state.get_correctors()["names"]]
            elif hasattr(self, "bpms_list") and elements_list is self.bpms_list:
                    ref_list=[str(x) for x in state.get_bpms()["names"]]
            elif hasattr(self,"quadrupoles_list") and elements_list is self.quadrupoles_list:
                ref_list = [str(x) for x in state.get_quadrupoles()["names"]]
            elif hasattr(self,"screens_list") and elements_list is self.screens_list:
                ref_list = [str(x) for x in state.get_screens()["names"]]
            else:
                ref_list=[elements_list.item(i).text() for i in range(elements_list.count())]
        else:
            ref_list = [elements_list.item(i).text() for i in range(elements_list.count())]

        selected_set={str(x) for x in selected_names}
        return [name for name in ref_list if name in selected_set]


    def _saving_func(self, elements_list, filename, saving_name, *, use_dialog=True,base_dir=None):  # * - must be passed by keyword
        items = elements_list.selectedItems()
        if not items:
            selected_names = [elements_list.item(i).text() for i in range(elements_list.count())]
        else:
            selected_names=[it.text() for it in items]

        ordered_names=self._refactor_names_order(elements_list,selected_names)
        base = base_dir or self.dir_name
        os.makedirs(base, exist_ok=True)
        if use_dialog:
            fn, _ = QFileDialog.getSaveFileName(self, f"{saving_name}", os.path.join(base, f"{filename}"),"Text (*.txt)")
            if not fn:
                return
        else:
            fn = os.path.join(base, filename)
        with open(fn, "w") as f:
            for it in ordered_names:
                f.write(f"{it}\n")

    def _loading_func(self, elements_list, filename, loading_name, *, use_dialog=True, base_dir=None):
        base = base_dir or self.dir_name
        if use_dialog:
            fn, _ = QFileDialog.getOpenFileName(self, f"{loading_name}", os.path.join(base, f"{filename}"),"Text (*.txt)")
        else:
            fn = os.path.join(base, filename)
        selected = None
        if fn and os.path.isfile(fn):
            with open(fn, "r") as f:
                selected = [ln.strip() for ln in f]
        if not selected:
            state = getattr(self, "initial_state", None)
            if state is not None:
                if hasattr(self, "bpms_list") and elements_list is self.bpms_list:
                    selected = state.get_bpms()["names"]
                elif hasattr(self, "correctors_list") and elements_list is self.correctors_list:
                    actuator_mode = getattr(getattr(self, "actuator_mode", None), "name", "Kicker")
                    if actuator_mode == "QM" and hasattr(self, "qm_corrs"):
                        selected = self.qm_corrs
                    else:
                        selected = state.get_correctors()["names"]
                elif hasattr(self, "quadrupoles_list") and elements_list is self.quadrupoles_list:
                    selected = state.get_quadrupoles()["names"]
                elif hasattr(self, "screens_list") and elements_list is self.screens_list:
                    selected = state.get_screens()["names"]
                else:
                    selected = [elements_list.item(i).text() for i in range(elements_list.count())]
        elements_list.clearSelection()
        for name in selected:
            for it in elements_list.findItems(name, _match_exactly()):
                it.setSelected(True)

    def _load_correctors(self):
        self._loading_func(loading_name="Load Correctors", filename="correctors.txt",elements_list=self.correctors_list)

    def _load_bpms(self):
        self._loading_func(loading_name="Load BPMs", filename="bpms.txt", elements_list=self.bpms_list)

    def _load_quadrupoles(self):
        self._loading_func(loading_name="Load Quadrupoles", filename="quadrupoles.txt", elements_list=self.quadrupoles_list)

    def _load_screens(self):
        self._loading_func(loading_name="Load Screens", filename="screens.txt", elements_list=self.screens_list)

    def _pick_and_load_data_dir(self, button_ui, button_name, oper):
        default_dir = f"~/CERN-Flight_Simulator-Data/"
        default_dir = os.path.expanduser(os.path.expandvars(default_dir))
        os.makedirs(default_dir, exist_ok=True)
        folder = QFileDialog.getExistingDirectory(self, "Select data directory", default_dir)
        if not folder:
            return
        info = self._find_useful_files(folder)
        if not info["ok"]:
            QMessageBox.warning(self, "Load data", "Wrong data directory selected")
        self._data_dirs[oper] = info
        button_ui.setText(folder)

    def _pick_and_load_disp_data(self):
        self._pick_and_load_data_dir(oper="dfs", button_ui=self.dfs_response_3, button_name="DFS Data Loaded")

    def _pick_and_load_wake_data(self):
        self._pick_and_load_data_dir(oper="wfs", button_ui=self.wfs_response_3, button_name="WFS Data Loaded")

    def _pick_and_load_traj_data(self):
        self._pick_and_load_data_dir(oper="traj", button_ui=self.trajectory_response_3,button_name="Trajectory Data Loaded")

    def _pick_and_load_machine_status_file(self):
        directory = os.path.expanduser(os.path.expandvars("~/CERN-Flight_Simulator-Data/"))
        filename, _ = QFileDialog.getOpenFileName(self, "Load Machine Status", directory, "Machine status (*.pkl)")
        if not filename:
            return
        try:
            state = State(filename=filename)
        except Exception as exc:
            QMessageBox.warning(self, "Restore machine status", f"Could not load this file:\n{exc}")
            return

        interface_id = f"{type(self.interface).__module__}.{type(self.interface).__name__}"
        if state.get_interface_id() != interface_id:
            QMessageBox.warning(self, "Restore machine status", "This machine status was created on a different interface/machine. Choose appropriate one or change inetrface.")
            return
        self.log("Restoring correctores from machine status...")
        correctors = state.get_correctors()
        self.log("Correctors restored!")
        self.log("Restoring quadrupoles from machine status...")
        quadrupoles = state.get_quadrupoles()
        self.log("Quadrupoles restored!")
        self.log("Restoring sextupoles from machine status...")
        sextupoles = state.get_sextupoles()
        self.log("Sextupoles restored!")
        self.log("Restoring correctors' strengths from machine status...")
        if self.interface.restore_correctors_state(state) is False:
            self.log("Warning: not every corrector was set at its saved current.")
            QMessageBox.warning(self, "Restore machine status", "Not every corrector was set back at its saved current. Check them on the machine.")
        else:
            self.log("Correctors restored!")
        self.log("Restoring quadrupoles' state from machine status...")
        self.interface.restore_quadrupoles_state(state)
        self.log("Quadrupoles restored!")
        self.log("Restoring sextupoles' state from machine status...")
        self.interface.restore_sextupoles_state(state)
        self.log("Sextupoles restored!")
        self.log("Restoring energy and intensity from machine status...")
        self.interface.restore_beam_settings(state.get_beam_settings())
        self.log("Energy and intensity restored!")

        timestamp = state.get_timestamp()
        if getattr(timestamp, "tzinfo", None) is not None:
            saved_at = timestamp.isoformat(sep=" ", timespec="seconds")
        else:
            saved_at = timestamp.strftime("%Y-%m-%d %H:%M:%S") + " (legacy/local timezone)"
        if hasattr(self, "machine_status_file"):
            self.machine_status_file.setText(filename)
        if hasattr(self, "log"):
            self.log(f"Machine status restored from {filename}; saved at {saved_at}")
        QMessageBox.information(self, "Restore machine status", f"Machine status restored.\nSaved at: {saved_at}")

    def save_session_settings(self, w1, w2, w3, rcond, iters, gain, beta, max_horizontal_current,max_vertical_current, is_triangular,bpm_weights,Axx, Ayy,Axy,Ayx, Bx, By, is_jitter_subtraction_checked, machine_state_file=None):
        clock_now = getattr(self, "_clock_now", datetime.now)
        saved_at = clock_now()
        save_session_dir = getattr(self, "_session_dir", None)
        if save_session_dir is None:
            time_str = saved_at.strftime("%y%m%d%H%M%S")
            default_dir = os.path.expanduser(os.path.expandvars("~/CERN-Flight_Simulator-Data/"))
            save_session_dir = os.path.join(default_dir, f"BBA_{self.interface.get_name()}{time_str}_session_settings")
        os.makedirs(save_session_dir, exist_ok=True)
        self.session_database_3.setText(save_session_dir)
        if machine_state_file is not None:
            machine_state_file.save(filename=os.path.join(save_session_dir, "machine_status.pkl"))
        self._saving_func(elements_list=self.correctors_list, filename="correctors.txt", saving_name="Save Correctors",
                          use_dialog=False, base_dir=save_session_dir)
        self._saving_func(elements_list=self.bpms_list, filename="bpms.txt", saving_name="Save BPMs", use_dialog=False,
                          base_dir=save_session_dir)

        correction_settings = {
            "saved_at": saved_at.isoformat(timespec="seconds"),
            "timezone": getattr(self, "_clock_zone_name", None),
            "nsamples": int(self.interface.nsamples),
            "beam_change_values": getattr(self, "_beam_change_values", {}),
            "actuator_mode": "Kicker",
            "w1": w1,
            "w2": w2,
            "w3": w3,
            "rcond": rcond,
            "iters": iters,
            "gain": gain,
            "beta": beta,
            "max_horizontal_current" : max_horizontal_current,
            "max_vertical_current" : max_vertical_current,
            "is_triangular" : is_triangular,
            "is_jitter_subtraction_checked": is_jitter_subtraction_checked,
            "bpm_weights" : bpm_weights,
            "data_dirs": {k: (self._expand_data_path(v["dir"]) if v else None)
                          for k, v in self._data_dirs.items()},
        }

        with open(os.path.join(save_session_dir, "correction_settings.json"), "w") as f:
            json.dump(correction_settings, f, indent=2)

        def __save_graph_data(path, series):
            with open(path, "w") as f:
                f.write("Iteration\tvalue\n")
                for i, v in enumerate(series, start=1):
                    f.write(f"{i}\t{v}\n")

        __save_graph_data(os.path.join(save_session_dir, "trajectory_x_after_correction.txt"), self._hist_orbit_x)
        __save_graph_data(os.path.join(save_session_dir, "dispersion_x_after_correction.txt"), self._hist_disp_x)
        __save_graph_data(os.path.join(save_session_dir, "wakefield_x_after_correction.txt"), self._hist_wake_x)
        __save_graph_data(os.path.join(save_session_dir, "trajectory_y_after_correction.txt"), self._hist_orbit_y)
        __save_graph_data(os.path.join(save_session_dir, "dispersion_y_after_correction.txt"), self._hist_disp_y)
        __save_graph_data(os.path.join(save_session_dir, "wakefield_y_after_correction.txt"), self._hist_wake_y)

        corrs, bpms = self._get_selection()
        R0xx, R0yy, R0xy, R0yx, B0x, B0y, R1xx, R1yy, R1xy, R1yx, B1x, B1y, R2xx, R2yy, R2xy, R2yx, B2x, B2y, hcorrs0, vcorrs0, hcorrs1, vcorrs1, hcorrs2, vcorrs2, bpms0, bpms1, bpms2 = self._get_data_from_loaded_directories(
            selected_corrs=corrs, selected_bpms=bpms)
        correction_matrices = {
            "Axx": Axx, "Ayy": Ayy, "Axy":Axy, "Ayx":Ayx, "B0x": B0x, "B0y": B0y,
            "R0xx": R0xx, "R1xx": R1xx, "R2xx": R2xx,
            "R0yy": R0yy, "R1yy": R1yy, "R2yy": R2yy,
            "R0xy": R0xy, "R0yx": R0yx, "R1xy": R1xy, "R1yx": R1yx,
            "R2xy": R2xy, "R2yx": R2yx,
            "Bx": Bx, "By": By,
        }
        with open(os.path.join(save_session_dir, "correction_matrices.pkl"), "wb") as f:
            pickle.dump(correction_matrices, f)

    def _find_matching_emittance_settings(self, states_dir):
        states_dir = os.path.normpath(states_dir)
        parent_dir = os.path.dirname(states_dir)
        settings_path = os.path.join(parent_dir, "emittance_settings.json")
        if not os.path.isfile(settings_path):
            return None
        try:
            with open(settings_path, "r") as f:
                settings = json.load(f)
        except Exception:
            return None
        saved_states_dir = settings.get("states_dir") or settings.get("data_session")
        if not saved_states_dir:
            return None
        saved_states_dir = os.path.normpath(os.path.expanduser(os.path.expandvars(str(saved_states_dir))))
        if saved_states_dir != states_dir and os.path.basename(saved_states_dir) != os.path.basename(states_dir):
            return None
        return settings

    def load_screens_data(self, folder=None):
        if folder is None:
            default_dir = f"~/CERN-Flight_Simulator-Data/"
            default_dir = os.path.expanduser(os.path.expandvars(default_dir))
            os.makedirs(default_dir, exist_ok=True)
            folder = QFileDialog.getExistingDirectory(self, "Select database", default_dir)
            if not folder:
                return
        folder = os.path.abspath(os.path.expanduser(os.path.expandvars(str(folder))))
        if not os.path.isdir(folder):
            QMessageBox.warning(self, "Load session", "Selected data directory does not exist.")
            return
        quad_selected = None
        screens = []
        folder_base_name = os.path.basename(os.path.normpath(folder))
        state_file_pattern = re.compile(r"screen_(\d+)_step_(\d+)_shot_(\d+)\.pkl")

        if folder_base_name.startswith(("states_", "screens_data_")) or any(state_file_pattern.fullmatch(name) for name in os.listdir(folder)):
            state_folder_path = folder
        else:
            states_folders = [name for name in sorted(os.listdir(folder)) if name.startswith(("states_", "screens_data_")) and os.path.isdir(os.path.join(folder, name))]
            if not states_folders:
                QMessageBox.information(self, "Load session", "No states or screen data folder found.")
                return
            state_folder_path = os.path.join(folder, states_folders[0])

        self.load_screens_data_database.setText(state_folder_path)
        print(f"Loading scan data from {state_folder_path}")
        state_folder_name = os.path.basename(os.path.normpath(state_folder_path))
        saved_scan_settings = self._find_matching_emittance_settings(state_folder_path)
        is_quad_scan = bool(saved_scan_settings.get("is_quad_scan", not state_folder_name.startswith("screens_data"))) if saved_scan_settings else not state_folder_name.startswith("screens_data")
        if is_quad_scan and state_folder_name.startswith("states_"):
            quad_selected = state_folder_name.removeprefix("states_").rsplit("_", 2)[0]
        if is_quad_scan and saved_scan_settings:
            quad_selected = saved_scan_settings.get("quad_name") or quad_selected

        self.quadrupoles_list.clearSelection()
        if quad_selected:
            for it in self.quadrupoles_list.findItems(quad_selected, _match_exactly()):
                it.setSelected(True)
        state_files = []
        for filename in sorted(os.listdir(state_folder_path)):
            if state_file_pattern.fullmatch(filename):
                state_files.append(os.path.join(state_folder_path, filename))
        self.loaded_state_files = []
        self.loaded_states_from_scan = []
        for state_file in state_files:
            try:
                self.loaded_states_from_scan.append(State(filename=state_file))
                self.loaded_state_files.append(state_file)
            except Exception as e:
                print(f"Couldn't load {state_file}, because {e}")
        print(f"Loaded {len(self.loaded_states_from_scan)} states")

        screens_by_index = {}
        for path, state in zip(self.loaded_state_files, self.loaded_states_from_scan):
            filename = os.path.basename(path)
            parts = filename.replace(".pkl", "").split("_")
            screen_i = int(parts[1]) # screen_0000_step_0003_shot_0001.pkl -> 0000
            if screen_i in screens_by_index: continue
            screen_data = state.get_screens()
            state_screen_names = list(screen_data.get("names", []))
            if state_screen_names:
                screens_by_index[screen_i] = str(state_screen_names[0])

        screens = [screens_by_index[index] for index in sorted(screens_by_index)]
        add_missing_screens = getattr(self, "_add_missing_screens_to_list", None)
        if callable(add_missing_screens):
            add_missing_screens(screens)
        self.screens_list.clearSelection()
        list_names = [str(self.screens_list.item(i).data(Qt.ItemDataRole.UserRole) or self.screens_list.item(i).text())
                      for i in range(self.screens_list.count())]
        for screen in screens:
            listed = self.interface.match_screen_name(screen, list_names)
            if listed is None:
                continue
            for it in self.screens_list.findItems(listed, _match_exactly()):
                it.setSelected(True)
        if quad_selected:
            self.quadrupoles_list.blockSignals(True)
            self.quadrupoles_list.clearSelection()

            for i in range(self.quadrupoles_list.count()):
                item = self.quadrupoles_list.item(i)
                item_name = str(item.data(Qt.ItemDataRole.UserRole) or item.text())
                if item_name == quad_selected:
                    item.setSelected(True)
                    break

            self.quadrupoles_list.blockSignals(False)
            self._last_selected_quadrupoles = [quad_selected]

        read_filenames = []
        for path in self.loaded_state_files:
            filename = os.path.basename(path)
            parts = filename.replace(".pkl", "").split("_")
            screen_i = int(parts[1]) # screen_0000_step_0003_shot_0001.pkl -> 0000
            step_i = int(parts[3]) # screen_0000_step_0003_shot_0001.pkl -> 0003
            shot_i = int(parts[5]) # screen_0000_step_0003_shot_0001.pkl -> 0001
            read_filenames.append((screen_i, step_i, shot_i))

        if not read_filenames:
            QMessageBox.warning(self, "Load session", "Couldn't find names like screen_0000_step_0003_shot_0001.pkl.")
            return
        nshots = max(shot_i for screen_i, step_i, shot_i in read_filenames)+1
        nscreens = max(screen_i for screen_i, step_i, shot_i in read_filenames)+1


        screen_current_ranges = {}
        if is_quad_scan:
            current_A_min = float(self.minimum_current.value())
            current_A_max = float(self.maximum_current.value())
            steps_per_screen = {}
            for screen_i, step_i, shot_i in read_filenames:
                steps_per_screen.setdefault(screen_i, set()).add(step_i)
            scan_steps = max(len(indices) for indices in steps_per_screen.values())
            if saved_scan_settings:
                if saved_scan_settings.get("scan_steps") is not None:
                    scan_steps = int(saved_scan_settings["scan_steps"])
                if saved_scan_settings.get("current_A_min") is not None:
                    current_A_min = float(saved_scan_settings["current_A_min"])
                if saved_scan_settings.get("current_A_max") is not None:
                    current_A_max = float(saved_scan_settings["current_A_max"])
                screen_current_ranges = dict(saved_scan_settings.get("screen_current_ranges") or {})
                self.minimum_current.setValue(current_A_min)
                self.maximum_current.setValue(current_A_max)
        else:
            current_A_min, current_A_max, scan_steps = 0.0, 0.0, 0.0
            quad_selected = None

        self._screen_current_ranges = {str(screen): (float(values[0]), float(values[1])) for screen, values in screen_current_ranges.items()}
        self._update_per_screen_ranges_button()

        print(f"Nshots: {nshots}, Scan steps per screen: {scan_steps}")

        is_fit_quad_strength_checked = bool(self.fit_quadrupole_strength_checkbox.isChecked())

        self.emittance_settings = {
            "current_A_min": current_A_min,
            "current_A_max": current_A_max,
            "screen_current_ranges": screen_current_ranges,
            "scan_steps": scan_steps,
            "nshots": nshots,
            "nscreens": nscreens,
            "is_quad_scan": is_quad_scan,
            "is_fit_quad_strength_checked": is_fit_quad_strength_checked,
            "bounds": self._get_bounds_from_gui(),
            "screens": screens if screens is not None else [],
            "quad_name": quad_selected if quad_selected else None,
        }

        self.steps_settings.setValue(int(scan_steps))
        self.meas_per_step.setValue(int(nshots))

        return self.loaded_states_from_scan

    def load_scan_and_optimization_settings(self):
        default_dir = "~/CERN-Flight_Simulator-Data/"
        default_dir = os.path.expanduser(os.path.expandvars(default_dir))
        os.makedirs(default_dir, exist_ok=True)
        folder = QFileDialog.getExistingDirectory(self, "Select session directory.", default_dir)
        if not folder:
            return
        self.session_directory.setText(folder)
        emittance_settings_path = os.path.join(folder, "emittance_settings.json")
        if not os.path.isfile(emittance_settings_path):
            QMessageBox.warning(self, "Load session settings", "Selected directory does not contain emittance_settings.json.")
            return
        try:
            with open(emittance_settings_path, "r") as f:
                saved_settings = json.load(f)
            quad_txt = os.path.join(folder, "quadrupoles.txt")
            if os.path.isfile(quad_txt):
                self._loading_func(elements_list=self.quadrupoles_list, filename="quadrupoles.txt", loading_name="Load Quadrupoles", use_dialog=False, base_dir=folder)

            screens_txt = os.path.join(folder, "screens.txt")
            if os.path.isfile(screens_txt):
                self._loading_func(elements_list=self.screens_list, filename="screens.txt", loading_name="Load Screens", use_dialog=False, base_dir=folder)

            states_dir = (saved_settings.get("states_dir") or saved_settings.get("data_session"))

            if states_dir:
                states_dir = os.path.expanduser(os.path.expandvars(str(states_dir)))
            if not states_dir or not os.path.isdir(states_dir):
                states_dir = None
                for name in os.listdir(folder):
                    candidate = os.path.join(folder, name)
                    if (os.path.isdir(candidate) and (name.startswith("states_") or name.startswith("screens_data_"))):
                        states_dir = candidate
                        break
            if not states_dir or not os.path.isdir(states_dir):
                QMessageBox.warning(self, "Load session settings", "Could not find the saved screen State files for this session.")
                return

            loaded_states = self.load_screens_data(states_dir)
            if not loaded_states:
                QMessageBox.warning(self,"Load session settings", "No screen State files could be loaded for this session.")
                return
            self.emittance_settings = saved_settings
            self.load_screens_data_database.setText(states_dir)

            if "is_quad_scan" not in self.emittance_settings:
                data_folder_name = os.path.basename(os.path.normpath(states_dir))
                self.emittance_settings["is_quad_scan"] = not data_folder_name.startswith("screens_data")
            if self.emittance_settings.get("current_A_min") is not None:
                self.minimum_current.setValue(float(self.emittance_settings["current_A_min"]))
            if self.emittance_settings.get("current_A_max") is not None:
                self.maximum_current.setValue(float(self.emittance_settings["current_A_max"]))
            if "screen_current_ranges" in self.emittance_settings:
                self._screen_current_ranges = {str(screen): (float(values[0]), float(values[1])) for screen, values in dict(self.emittance_settings["screen_current_ranges"]).items()}
                self._update_per_screen_ranges_button()
            if "scan_steps" in self.emittance_settings:
                self.steps_settings.setValue(int(self.emittance_settings["scan_steps"]))
            if "nshots" in self.emittance_settings:
                self.meas_per_step.setValue(int(self.emittance_settings["nshots"]))
            if "is_fit_quad_strength_checked" in self.emittance_settings:
                self.fit_quadrupole_strength_checkbox.setChecked(bool(self.emittance_settings["is_fit_quad_strength_checked"]))
            saved_bounds = self.emittance_settings.get("bounds")
            if saved_bounds:
                self._set_bounds_from_saved_settings(saved_bounds)
            self.session = self._get_session_data_from_database()

            if self.session is None:
                QMessageBox.warning(self,"Load session settings", "Could not reconstruct the scan session from the saved State files.")
                return
            if not saved_bounds or "quad_k1l_0" not in saved_bounds:
                self._set_default_quad_strength_bounds_from_session(self.session)

            self._refresh_plot_comboboxes_from_session(self.session)
            self._draw_live_scan(self.session)
            self._clear_fit_panel()

            self.minimum_current.setEnabled(False)
            self.maximum_current.setEnabled(False)
            self.per_screen_ranges_button.setEnabled(False)
            self.steps_settings.setEnabled(False)
            self.meas_per_step.setEnabled(False)
            self.quadrupoles_list.setEnabled(False)

            QMessageBox.information(self, "Data directory selected", "Loaded session. Run optimization to calculate a new fit.")

        except Exception as e:
            QMessageBox.warning(self, "Load session settings", f"Could not load session:\n{e}")
            return

    def load_session_settings(self):
        default_dir = f"~/CERN-Flight_Simulator-Data/"
        default_dir = os.path.expanduser(os.path.expandvars(default_dir))
        os.makedirs(default_dir, exist_ok=True)
        folder = QFileDialog.getExistingDirectory(self, "Select database", default_dir)
        if not folder:
            return
        if hasattr(self, "session_database_3"):
            self.session_database_3.setText(folder)
        correction_settings_path = os.path.join(folder, "correction_settings.json")
        if not os.path.isfile(correction_settings_path):
            QMessageBox.warning(self, "Load session", "Selected folder doesn't contain proper correction settings.")
            return
        try:
            with open(correction_settings_path, "r") as f:
                settings = json.load(f)
                actuator_mode = settings.get("actuator_mode", "Kicker")
        except Exception as e:
            QMessageBox.warning(self,"Load session",f"Couldn't read correction_settings.json: {e}")
            return

        if hasattr(self, "actuator_mode_combo"):
            mode_text = "Quadrupole movers" if actuator_mode == "QM" else "Correctors"
            idx = self.actuator_mode_combo.findText(mode_text)
            if idx >= 0:
                self.actuator_mode_combo.setCurrentIndex(idx)

        if hasattr(self, "_refresh_corrector_list"):
            self._refresh_corrector_list()
        if hasattr(self, "_update_qm_widgets_visibility"):
            self._update_qm_widgets_visibility()
        if hasattr(self, "_refresh_metric_plots_for_mode"):
            self._refresh_metric_plots_for_mode()
        if hasattr(self, "_refresh_specific_bpm_candidates"):
            self._refresh_specific_bpm_candidates()

        if actuator_mode == "QM" and os.path.isfile(os.path.join(folder, "quadrupole_movers.txt")):
            actuator_filename = "quadrupole_movers.txt"
        else:
            actuator_filename = "correctors.txt"

        self._loading_func(elements_list=self.correctors_list, filename=actuator_filename,
                           loading_name="Load Actuators", use_dialog=False, base_dir=folder)
        self._loading_func(elements_list=self.bpms_list, filename="bpms.txt", loading_name="Load BPMs",
                           use_dialog=False, base_dir=folder)

        paths = settings.get("data_dirs") or {}
        for key, path in paths.items():
            if not path:
                continue
            path = os.path.expanduser(os.path.expandvars(path))
            info = self._find_useful_files(path)
            if info["ok"]:
                self._data_dirs[key] = info
            else:
                QMessageBox.warning(self, "Load session", "Data directory not found")
                return
        saved_at = settings.get("saved_at")
        message = "Loaded session"
        if saved_at:
            message += f"\nSaved at: {saved_at}"
        QMessageBox.information(self.session_database_3, "Data directory selected", message)

        if "w1" in settings: self.lineEdit.setText(str(settings["w1"]))
        if "w2" in settings: self.lineEdit_2.setText(str(settings["w2"]))
        if "w3" in settings: self.lineEdit_3.setText(str(settings["w3"]))
        if "rcond" in settings: self.lineEdit_4.setText(str(settings["rcond"]))
        if "iters" in settings:  self.lineEdit_5.setText(str(settings["iters"]))
        if "gain" in settings: self.lineEdit_6.setText(str(settings["gain"]))
        if "beta" in settings: self.lineEdit_beta.setText(str(settings["beta"]))
        if "nsamples" in settings and hasattr(self, "nsamples_input"):
            self.nsamples_input.setText(str(settings["nsamples"]))
        if settings.get("beam_change_values") and hasattr(self, "_load_beam_change_values"):
            self._load_beam_change_values(settings["beam_change_values"])
        if "is_triangular" in settings: self.triangular_checkbox.setChecked(settings["is_triangular"])
        if "is_jitter_subtraction_checked" in settings: self.subtract_jitter_checkbox.setChecked(settings["is_jitter_subtraction_checked"])
        if "bpm_weights" in settings:
            for bpm_name, weights in settings["bpm_weights"].items():
                try:
                    w1_bpm, w2_bpm, w3_bpm = weights
                    self.bpm_weights[str(bpm_name)] = (float(w1_bpm), float(w2_bpm), float(w3_bpm))
                except Exception:
                    continue
            for i in range(self.bpms_list.count()):
                item = self.bpms_list.item(i)
                self._update_bpm_weights(item)

        data_dirs = settings.get("data_dirs") or {}
        if hasattr(self, "trajectory_response_3"):
            self.trajectory_response_3.setText(data_dirs.get("traj") or "")

        if actuator_mode == "Kicker":
            if "max_horizontal_current" in settings:
                self.max_horizontal_current_spinbox.setValue(settings["max_horizontal_current"])
            if "max_vertical_current" in settings:
                self.max_vertical_current_spinbox.setValue(settings["max_vertical_current"])
            if hasattr(self, "dfs_response_3"):
                self.dfs_response_3.setText(data_dirs.get("dfs") or "")
            if hasattr(self, "wfs_response_3"):
                self.wfs_response_3.setText(data_dirs.get("wfs") or "")
        else:
            if hasattr(self, "_refresh_specific_bpm_candidates"):
                self._refresh_specific_bpm_candidates()
            if "specific_bpm" in settings and hasattr(self, "specific_bpm_combo"):
                idx = self.specific_bpm_combo.findText(str(settings["specific_bpm"]))
                if idx >= 0:
                    self.specific_bpm_combo.setCurrentIndex(idx)
            if "max_horizontal_range" in settings:
                self.max_horizontal_current_spinbox.setValue(settings["max_horizontal_range"])
            if "max_vertical_range" in settings:
                self.max_vertical_current_spinbox.setValue(settings["max_vertical_range"])

    def save_emittance_measurement_session(self, session=None, is_fit_quad_strength_checked=None, bounds=None, target_dir=None):
        save_session_dir = target_dir or getattr(self, "dir_name", None) or self.session_directory.text()
        os.makedirs(save_session_dir, exist_ok=True)
        preserve_status = getattr(self, "_preserve_quadrupole_status_files", None)
        if callable(preserve_status):
            preserve_status(save_session_dir)
        self.session_directory.setText(save_session_dir)
        self._saving_func(elements_list=self.quadrupoles_list, filename="quadrupoles.txt", saving_name="Save quadrupoles", use_dialog=False, base_dir=save_session_dir)
        self._saving_func(elements_list=self.screens_list, filename="screens.txt", saving_name="Save screens", use_dialog=False, base_dir=save_session_dir)
        settings_path = os.path.join(save_session_dir, "emittance_settings.json")

        if os.path.isfile(settings_path):
            try:
                with open(settings_path, "r") as f:
                    self.emittance_settings = json.load(f)
            except Exception:
                self.emittance_settings = {}
        else:
            self.emittance_settings = {}
        self.emittance_settings.pop("ls_steps", None)
        for status_key in ("quadrupoles_status", "model_quadrupoles_status"):
            status = getattr(self, status_key, None)
            if status is not None:
                self.emittance_settings[status_key] = status

        if session is not None:
            state_files = []
            for step in session.get("scan_steps", []):
                state_files.extend(step.get("state_files", []))
            state_files = list(dict.fromkeys(str(path) for path in state_files))
            self.emittance_settings.update({
                "current_A_min": session.get("current_A_min"),
                "current_A_max": session.get("current_A_max"),
                "screen_current_ranges": session.get("screen_current_ranges", {}),
                "scan_steps": session.get("steps"),
                "nshots": session.get("nshots"),
                "data_session": self.load_screens_data_database.text(),
                "states_dir": (session.get("states_dir") or self.load_screens_data_database.text()),
                "state_files": state_files,
                "quad_name": session.get("quad_name"),
                "screens": session.get("screens", []),
                "nscreens": session.get("nscreens", len(session.get("screens", []))),
                "is_quad_scan": session.get("is_quad_scan"),
                "reference_screen": session.get("reference_screen"),
                "K1L_0": session.get("K1L_0"),
                "K1L_values": session.get("K1L_values"),
                "current_values": session.get("current_values"),
                "quad_value_unit": session.get("quad_value_unit"),
                "sigma_unit": session.get("sigma_unit", "mm"),
            })

            if is_fit_quad_strength_checked is not None:
                self.emittance_settings["is_fit_quad_strength_checked"] = bool(is_fit_quad_strength_checked)

            if bounds is not None:
                self.emittance_settings["bounds"] = dict(bounds)

            if "optimization_result" in session:
                optimization_result_path = os.path.join(save_session_dir, "optimization_result.pkl")
                with open(optimization_result_path, "wb") as f:
                    pickle.dump(session["optimization_result"], f)
                self.emittance_settings["optimization_result_file"] = "optimization_result.pkl"

        with open(settings_path, "w") as f:
            json.dump(self.emittance_settings, f, indent=2)
        return save_session_dir

    def save_session_settings_qm_correction(self, w1, w2, w3, specific_bpm, rcond, iters, gain, beta, max_horizontal_range, max_vertical_range, is_triangular, bpm_weights, response, is_jitter_subtraction_checked):
        save_session_dir = getattr(self, "_session_dir", None)
        if save_session_dir is None:
            time_str = datetime.now().strftime("%y%m%d%H%M%S")
            default_dir = os.path.expanduser(os.path.expandvars("~/CERN-Flight_Simulator-Data/"))
            save_session_dir = os.path.join(default_dir, f"BBA_{self.interface.get_name()}{time_str}_QM_session_settings")
        os.makedirs(save_session_dir, exist_ok=True)
        self.session_database_3.setText(save_session_dir)
        self._saving_func(elements_list=self.correctors_list, filename="quadrupole_movers.txt", saving_name="Save Quadrupole Movers", use_dialog=False, base_dir=save_session_dir)
        self._saving_func(elements_list=self.bpms_list, filename="bpms.txt", saving_name="Save BPMs", use_dialog=False, base_dir=save_session_dir)
        correction_settings = {
            "actuator_mode": "QM",
            "w1": w1,
            "w2": w2,
            "w3": w3,
            "specific_bpm": specific_bpm,
            "rcond": rcond,
            "iters": iters,
            "gain": gain,
            "beta": beta,
            "max_horizontal_range": max_horizontal_range,
            "max_vertical_range": max_vertical_range,
            "is_triangular": is_triangular,
            "is_jitter_subtraction_checked": is_jitter_subtraction_checked,
            "bpm_weights": bpm_weights,
            "data_dirs": {
                k: (self._expand_data_path(v["dir"]) if v else None)
                for k, v in self._data_dirs.items()
            },
        }

        with open(os.path.join(save_session_dir, "correction_settings.json"), "w") as f:
            json.dump(correction_settings, f, indent=2)

        def __save_graph_data(path, series):
            with open(path, "w") as f:
                f.write("Iteration\tvalue\n")
                for i, v in enumerate(series, start=1):
                    f.write(f"{i}\t{v}\n")

        __save_graph_data(os.path.join(save_session_dir, "trajectory_x_after_correction.txt"), self._hist_orbit_x)
        __save_graph_data(os.path.join(save_session_dir, "trajectory_y_after_correction.txt"), self._hist_orbit_y)
        __save_graph_data(os.path.join(save_session_dir, "trajectory_combined_after_correction.txt"), self._hist_orbit)

        qm_matrices = {}
        if response is not None:
            for key in ("qcorrs", "bpms", "R_xx", "R_xy", "R_yx", "R_yy", "T_xx", "T_yy"):
                if key in response:
                    qm_matrices[key] = response[key]

        with open(os.path.join(save_session_dir, "qm_correction_matrices.pkl"), "wb") as f:
            pickle.dump(qm_matrices, f)

        return save_session_dir

    def _save_bumps_machine_status(self, machine_state=None):
        if self._session_dir is None:
            saved_at = self._clock_now()
            time_str = saved_at.strftime("%y%m%d%H%M%S")
            default_dir = os.path.expanduser(os.path.expandvars("~/CERN-Flight_Simulator-Data/"))
            self._session_dir = os.path.join(default_dir, f"Bumps_{self.interface.get_name()}{time_str}_session_settings")
            os.makedirs(self._session_dir, exist_ok=True)
        machine_status_file = os.path.join(self._session_dir, "machine_status.pkl")
        if not os.path.isfile(machine_status_file):
            machine_state = machine_state if machine_state is not None else self.interface.get_state()
            machine_state.timestamp = self._clock_now()
            machine_state.save(filename=machine_status_file)
        self.loadsave_session_edit.setText(self._session_dir)
        return machine_state

    def _bumps_session_settings(self):
        try:
            rcond = float(self.pinv_edit.text())
        except ValueError:
            rcond = self.pinv_value
        try:
            beta = float(self.beta_edit.text())
        except ValueError:
            beta = self.beta_value
        return {
            "saved_at": self._clock_now().isoformat(timespec="seconds"),
            "timezone": self._clock_zone_name,
            "interface_id": f"{type(self.interface).__module__}.{type(self.interface).__name__}",
            "actuator_mode": "Kicker",
            "response_directory": self._expand_path(self.response_dir_edit.text()) if self.response_dir_edit.text().strip() else None,
            "nsamples": int(self.interface.nsamples),
            "rcond": rcond,
            "beta": beta,
            "max_horizontal_current": self.max_horizontal_current_spinbox.value(),
            "max_vertical_current": self.max_vertical_current_spinbox.value(),
            "is_triangular": bool(self.triangular_checkbox.isChecked()),
            "targets": {
                str(name): {plane: None if value is None else float(value) for plane, value in values.items()}
                for name, values in self.targets.items()
            },
        }

    def _save_bumps_reference_orbit(self):
        if self.reference is None:
            return
        reference = {
            "names": [str(name) for name in self.reference["names"]],
            "x": np.asarray(self.reference["x"], dtype=float).tolist(),
            "y": np.asarray(self.reference["y"], dtype=float).tolist(),
        }
        with open(os.path.join(self._session_dir, "reference_orbit.json"), "w") as file:
            json.dump(reference, file, indent=2)

    def _save_bumps_result_orbit(self):
        if self._result_data is None:
            return
        result = {key: value for key, value in self._result_data.items() if value is not None}
        result["bpms"] = np.asarray(result["bpms"], dtype=str)
        np.savez(os.path.join(self._session_dir, "result_orbit.npz"), **result)

    def _save_bumps_session(self, machine_state=None):
        self._save_bumps_machine_status(machine_state)
        self._saving_func(self.correctors_list, "correctors.txt", "Save Correctors", use_dialog=False,
                          base_dir=self._session_dir)
        self._saving_func(self.bpms_list, "bpms.txt", "Save BPMs", use_dialog=False, base_dir=self._session_dir)
        with open(os.path.join(self._session_dir, "correction_settings.json"), "w") as file:
            json.dump(self._bumps_session_settings(), file, indent=2)
        self._save_bumps_reference_orbit()
        self._save_bumps_result_orbit()

    def _append_bumps_kicks(self, before, after):
        kicks_path = os.path.join(self._session_dir, "kicks.txt")
        write_header = not os.path.exists(kicks_path)
        with open(kicks_path, "a") as file:
            if write_header:
                file.write("time\titeration\tcorrector\tbdes_before\tapplied_kick\tbdes_after\n")
            timestamp = self._clock_now().isoformat(timespec="seconds")
            for corrector, previous, current in zip(self.corrector_names, before, after):
                applied_kick = float(current) - float(previous)
                file.write(f"{timestamp}\t{self._session_kick_count}\t{corrector}\t{float(previous):.12g}\t{applied_kick:.12g}\t{float(current):.12g}\n")

    @staticmethod
    def _bumps_state_matches_interface(state, interface):
        interface_id = f"{type(interface).__module__}.{type(interface).__name__}"
        return state.get_interface_id() == interface_id

    def _pick_and_load_bumps_machine_status_file(self):
        directory = os.path.expanduser(os.path.expandvars("~/CERN-Flight_Simulator-Data/"))
        filename, _ = QFileDialog.getOpenFileName(self, "Load Machine Status", directory, "Machine status (*.pkl)")
        if filename:
            self._load_bumps_machine_status_file(filename)

    def _load_bumps_machine_status_file(self, filename):
        try:
            state = State(filename=filename)
        except Exception as exc:
            QMessageBox.warning(self, "Restore machine status", f"Could not load this file:\n{exc}")
            return False
        if not self._bumps_state_matches_interface(state, self.interface):
            QMessageBox.warning(self, "Restore machine status", "This machine status was created on a different interface/machine. Choose appropriate one or change interface.")
            return False
        if self.interface.restore_correctors_state(state) is False:
            QMessageBox.warning(self, "Restore machine status", "Not every corrector was set back at its saved current. Check them on the machine.")
        self.interface.restore_quadrupoles_state(state)
        self.interface.restore_sextupoles_state(state)
        self.interface.restore_beam_settings(state.get_beam_settings())
        self.corr_status_edit.setText(filename)
        timestamp = state.get_timestamp()
        saved_at = timestamp.isoformat(sep=" ", timespec="seconds") if getattr(timestamp, "tzinfo", None) is not None else timestamp.strftime("%Y-%m-%d %H:%M:%S") + " (legacy/local timezone)"
        QMessageBox.information(self, "Restore machine status", f"Machine status restored.\nSaved at: {saved_at}")
        return True

    def _pick_and_load_bumps_session(self):
        directory = os.path.expanduser(os.path.expandvars("~/CERN-Flight_Simulator-Data/"))
        folder = QFileDialog.getExistingDirectory(self, "Select Bumps session directory", directory)
        if folder:
            self._load_bumps_session(folder)

    def _load_bumps_session(self, folder):
        folder = self._expand_path(folder)
        settings_path = os.path.join(folder, "correction_settings.json")
        if not os.path.isfile(settings_path):
            QMessageBox.warning(self, "Load session", "Selected folder doesn't contain proper correction settings.")
            return
        try:
            with open(settings_path, "r") as file:
                settings = json.load(file)
        except Exception as exc:
            QMessageBox.warning(self, "Load session", f"Couldn't read correction_settings.json: {exc}")
            return
        saved_interface_id = settings.get("interface_id")
        current_interface_id = f"{type(self.interface).__module__}.{type(self.interface).__name__}"
        if saved_interface_id is not None and saved_interface_id != current_interface_id:
            QMessageBox.warning(self, "Load session", "This session was created on a different interface/machine. Choose appropriate one or change interface.")
            return
        response_directory = settings.get("response_directory")
        if response_directory:
            self.response_matrix_data = None
            self._load_response_directory(response_directory)
            if self.response_matrix_data is None:
                return
        self.reference = None
        self.restore_state = None
        self._loading_func(self.correctors_list, "correctors.txt", "Load Correctors", use_dialog=False, base_dir=folder)
        self._loading_func(self.bpms_list, "bpms.txt", "Load BPMs", use_dialog=False, base_dir=folder)
        if "rcond" in settings:
            self.pinv_edit.setText(str(settings["rcond"]))
        if "beta" in settings:
            self.beta_edit.setText(str(settings["beta"]))
        if "nsamples" in settings:
            self.bpm_samples_edit.setText(str(settings["nsamples"]))
        if "max_horizontal_current" in settings:
            self.max_horizontal_current_spinbox.setValue(float(settings["max_horizontal_current"]))
        if "max_vertical_current" in settings:
            self.max_vertical_current_spinbox.setValue(float(settings["max_vertical_current"]))
        if "is_triangular" in settings:
            self.triangular_checkbox.setChecked(bool(settings["is_triangular"]))
        self.targets = {
            str(name): {plane: values.get(plane) for plane in ("x", "y")}
            for name, values in settings.get("targets", {}).items()
        }
        reference_path = os.path.join(folder, "reference_orbit.json")
        if os.path.isfile(reference_path):
            try:
                with open(reference_path, "r") as file:
                    reference = json.load(file)
                self.reference = {
                    "names": [str(name) for name in reference["names"]],
                    "x": np.asarray(reference["x"], dtype=float),
                    "y": np.asarray(reference["y"], dtype=float),
                }
            except Exception:
                self.reference = None
        reference_state_path = os.path.join(folder, "reference_machine_status.pkl")
        if os.path.isfile(reference_state_path):
            try:
                reference_state = State(filename=reference_state_path)
                if self._bumps_state_matches_interface(reference_state, self.interface):
                    self.restore_state = reference_state
            except Exception:
                self.restore_state = None
        self.restore_button.setEnabled(self.restore_state is not None)
        result_path = os.path.join(folder, "result_orbit.npz")
        if os.path.isfile(result_path):
            try:
                with np.load(result_path, allow_pickle=False) as result:
                    measured_x = result["measured_x"] if "measured_x" in result.files else None
                    measured_y = result["measured_y"] if "measured_y" in result.files else None
                    self._draw_result_plot(
                        result["bpms"].astype(str).tolist(), result["desired_x"], result["desired_y"],
                        result["predicted_x"], result["predicted_y"], measured_x, measured_y,
                    )
            except Exception:
                self._result_data = None
                self._draw_result_placeholder()
        self._session_dir = folder
        self._session_kick_count = 0
        kicks_path = os.path.join(folder, "kicks.txt")
        if os.path.isfile(kicks_path):
            with open(kicks_path, "r") as file:
                for line in file:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) > 1:
                        try:
                            self._session_kick_count = max(self._session_kick_count, int(fields[1]))
                        except ValueError:
                            pass
        self.loadsave_session_edit.setText(folder)
        self.apply_button.setEnabled(False)
        self._refresh_desired_plot()
        saved_at = settings.get("saved_at")
        QMessageBox.information(self, "Load session", f"Loaded session" + (f"\nSaved at: {saved_at}" if saved_at else ""))
