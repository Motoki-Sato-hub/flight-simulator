import os, pickle, json
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
            for it in elements_list.findItems(name, Qt.MatchFlag.MatchExactly):
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
        self.interface.restore_correctors_state(state)
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
        if saved_states_dir != states_dir:
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
        else:
            folder = folder
        self.load_screens_data_database.setText(folder)
        quad_selected = None
        screens = []
        folder_base_name = os.path.basename(os.path.normpath(folder))
        if folder_base_name.startswith("states_"):
            quad_selected = folder_base_name.removeprefix("states_").split("_")[0]
        else:
            quad_selected = None

        self.quadrupoles_list.clearSelection()
        if quad_selected:
            for it in self.quadrupoles_list.findItems(quad_selected, Qt.MatchFlag.MatchExactly):
                it.setSelected(True)
        else:
            for name in os.listdir(folder):
                if name.startswith("states_"):
                    quad_selected = name.removeprefix("states_")
                    for it in self.quadrupoles_list.findItems(quad_selected, Qt.MatchFlag.MatchExactly):
                        it.setSelected(True)
                    break
                else:
                    self.quadrupoles_list.setEnabled(False)
                    quad_selected = None

        state_files = []
        state_folder_path = folder
        for filename in sorted(os.listdir(state_folder_path)):
            if filename.endswith(".pkl"):
                state_files.append(os.path.join(state_folder_path, filename))
        self.loaded_state_files = state_files
        self.loaded_states_from_scan = []
        for state_file in state_files:
            try:
                self.loaded_states_from_scan.append(State(filename=state_file))
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
        self.screens_list.clearSelection()
        for screen in screens:
            for it in self.screens_list.findItems(screen, Qt.MatchFlag.MatchExactly):
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

        folder_name = os.path.basename(os.path.normpath(folder))
        is_quad_scan = not folder_name.startswith("screens_data")

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


        if is_quad_scan:
            delta_min = float(self.delta_min_scan.value())
            delta_max = float(self.delta_max_scan.value())
            scan_steps = max(step_i for screen_i, step_i, shot_i in read_filenames)+1
            saved_scan_settings = self._find_matching_emittance_settings(folder)
            if saved_scan_settings:
                if "delta_min" in saved_scan_settings:
                    delta_min = float(saved_scan_settings["delta_min"])
                if "delta_max" in saved_scan_settings:
                    delta_max = float(saved_scan_settings["delta_max"])
                self.delta_min_scan.setValue(delta_min)
                self.delta_max_scan.setValue(delta_max)
        else:
            delta_min, delta_max, scan_steps = 0.0, 0.0, 0.0
            quad_selected = None

        print(f"Nshots: {nshots}, Scan steps: {scan_steps}")

        is_fit_quad_strength_checked = bool(self.fit_quadrupole_strength_checkbox.isChecked())

        self.emittance_settings = {
            "delta_min": delta_min,
            "delta_max": delta_max,
            "scan_steps": scan_steps,
            "nshots": nshots,
            "nscreens": nscreens,
            "is_quad_scan": is_quad_scan,
            "is_fit_quad_strength_checked": is_fit_quad_strength_checked,
            "bounds": self._get_bounds_from_gui(),
            "screens": screens if screens is not None else [],
            "quad_name": quad_selected if quad_selected else None,
        }

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
            if "delta_min" in self.emittance_settings:
                self.delta_min_scan.setValue(float(self.emittance_settings["delta_min"]))
            if "delta_max" in self.emittance_settings:
                self.delta_max_scan.setValue(float(self.emittance_settings["delta_max"]))
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

            self.delta_min_scan.setEnabled(False)
            self.delta_max_scan.setEnabled(False)
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

        if session is not None:
            state_files = []
            for step in session.get("scan_steps", []):
                state_files.extend(step.get("state_files", []))
            state_files = list(dict.fromkeys(str(path) for path in state_files))
            self.emittance_settings.update({
                "delta_min": session.get("delta_min"),
                "delta_max": session.get("delta_max"),
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
        __save_graph_data(os.path.join(save_session_dir, "orbit_rms_x_after_correction.txt"), self._hist_abs_rms_x)
        __save_graph_data(os.path.join(save_session_dir, "orbit_rms_y_after_correction.txt"), self._hist_abs_rms_y)
        __save_graph_data(os.path.join(save_session_dir, "orbit_rms_xy_after_correction.txt"), self._hist_abs_rms_xy)

        qm_matrices = {}
        if response is not None:
            for key in ("qcorrs", "bpms", "R_xx", "R_xy", "R_yx", "R_yy", "T_xx", "T_yy"):
                if key in response:
                    qm_matrices[key] = response[key]

        with open(os.path.join(save_session_dir, "qm_correction_matrices.pkl"), "wb") as f:
            pickle.dump(qm_matrices, f)

        return save_session_dir
