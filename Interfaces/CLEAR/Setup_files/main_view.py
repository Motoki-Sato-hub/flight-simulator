import pyda
import pyda_japc
import json
import time
import os
import threading
import numpy as np
import scipy.io as sio
import scipy.ndimage as ndi
import pylogbook
from collections import deque

import PyQt6.QtCore as QtCore
from PyQt6.QtGui import QIcon, QIntValidator
from PyQt6.QtWidgets import (
    QWidget,
    QLabel,
    QLineEdit,
    QMessageBox,
    QTabWidget,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QSizePolicy,
    QStatusBar,
)

from view import initBGBox, initBMBox, initCentBox, initDevBox, initBTVBox
from view import initFiltBox, initImageBox, initProfBox, initPropBox, initZoomBox
from view.canvas import ImageView
from assets import CamList

class View(QWidget):
    """Init Self"""

    camera_frame_ready = QtCore.pyqtSignal()

    def __init__(self):
        super().__init__()

        #self.GUI = "TSG41.AWAKE-GUI-SUPPORT"
        #self.GUI2 = "TSG41.AWAKE-GUI-SUPPORT2"
        #self.GUI = ""
        #self.GUI2 = ""
        self.initJAPC()
        self._frame_lock = threading.Lock()
        self._frame_queue = deque(maxlen=1)
        self._subscription_generation = 0
        self._frame_update_pending = False
        self._dropped_frames = 0
        self._last_drop_report = 0
        self._beam_marks_path = os.path.join(
            os.path.expanduser("~"),
            ".config",
            "camgui",
            "beam_marks.json",
        )
        self.beam_marks = self._load_beam_marks()
        self._multi_save_total = 0
        self._multi_save_remaining = 0
        self._multi_save_camera = None
        self._multi_save_dir = None
        self._multi_save_paths = []
        self.can_connect = False
        self.dev_name = ""
        self.cam_device = ""
        self.dev_prop = ""
        self.dev_sett = ""
        self.dev_settBI = ""
        self.dev_settBI2 = ""
        self.dev_roi = ""
        self.camera_frame_ready.connect(
            self._process_latest_cam_data,
            QtCore.Qt.ConnectionType.QueuedConnection,
        )
        
        self.initUI()
        self.initLogbook()

    def initJAPC(self):
        """pyda initialization"""
        provider = pyda_japc.JapcProvider()
        self.client_get = pyda.SimpleClient(provider=provider)
        self.client_callback = pyda.CallbackClient(provider=provider)
        self.client_asynch = pyda.AsyncIOClient(provider=pyda_japc.JapcProvider())

    def initLogbook(self):
        """Create logbook instance"""

        # Adapting to pylogbook 3.5.0
        # For this to work, run `eval $(rbac-authenticate -l)` before launching
        # the GUI, or use cam_gui_launcher.sh which does this for location login.
        # For this to work, need to run from a control room machine

        self.elogs = {}
        self.elog = None
        try:
            self.elog = self._logbook_client("CLEAR")

        except Exception as e:
            if os.environ.get("RBAC_TOKEN_SERIALIZED"):
                print(f"initLogbook: {e}")
            self.statusBar.showMessage("Failed to setup logbook")
            time.sleep(1)

    def _logbook_name(self):
        return "CTF2" if self.dev_name.startswith("PHIN") else "CLEAR"

    def _logbook_client(self, name):
        if name not in self.elogs:
            self.elogs[name] = pylogbook.ActivitiesClient(name)
        return self.elogs[name]

    def _ensure_logbook(self, name=None):
        name = name or self._logbook_name()
        if name in self.elogs:
            return self.elogs[name]

        if not os.environ.get("RBAC_TOKEN_SERIALIZED"):
            self.statusBar.showMessage("No RBAC token; logbook unavailable")
            return None

        try:
            client = self._logbook_client(name)
            if name == "CLEAR":
                self.elog = client
            return client
        except Exception as e:
            print(f"_ensure_logbook: {e}")
            self.statusBar.showMessage(f"Failed to setup {name} logbook")
            return None

    def _tab_with_widgets(self, widgets):
        layout = QVBoxLayout()
        for widget in widgets:
            layout.addWidget(widget)
        layout.addStretch(1)

        tab = QWidget()
        tab.setLayout(layout)
        return tab

    def _create_save_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()

        row = QHBoxLayout()
        row.addWidget(QLabel("Images"))

        self.save_n_input = QLineEdit("10", self)
        self.save_n_input.setValidator(QIntValidator(1, 10000, self))
        self.save_n_input.setMaximumWidth(80)
        row.addWidget(self.save_n_input)
        row.addStretch(1)
        layout.addLayout(row)

        self.save_n_button = QPushButton("Save N Images", self)
        self.save_n_button.setStyleSheet("background-color:#999999")
        self.save_n_button.clicked[bool].connect(self.saveNImagesAction)
        layout.addWidget(self.save_n_button)

        self.save_n_status = QLabel("No multi-image save active", self)
        self.save_n_status.setWordWrap(True)
        layout.addWidget(self.save_n_status)

        layout.addStretch(1)
        tab.setLayout(layout)
        return tab

    def _reset_start_button(self):
        self.btn.setChecked(False)
        self.btn.setText("Start")

    def _stop_and_reset_start_button(self):
        self.stop_subs()
        self._reset_start_button()

    def _run_with_subscription_paused(self, action):
        was_running = self.sub_state
        if was_running:
            self.stop_subs()

        try:
            action()
        finally:
            if was_running:
                self.start_subs()

    def _ask_to_move_current_screen_out(self):
        if not self.ask_screen_out.isChecked():
            return

        camera = self.cam_props.get(self.dev_name, {})
        if not camera.get("screenInstalled"):
            return

        answer = QMessageBox.question(
            self,
            "Remove screen?",
            f"Move the {self.dev_name} screen to out before changing camera?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return

        labels = [
            self.position_dropdown.itemText(index).strip().lower()
            for index in range(self.position_dropdown.count())
        ]
        out_index = next(
            (index for index, label in enumerate(labels) if label == "out"),
            next(
                (index for index, label in enumerate(labels) if "out" in label),
                -1,
            ),
        )
        if out_index < 0:
            QMessageBox.warning(
                self,
                "No out position",
                f"No 'out' screen position is configured for {self.dev_name}.",
            )
            return

        self.position_dropdown.setCurrentIndex(out_index)
        self.positionScreenAction(self.position_dropdown.currentText())

    def initUI(self):
        """Initialize GUI"""

        # Load camera configuration from assets/cameras.json
        self.cam_props = CamList.CamList()
        self.camList = list(self.cam_props.keys())

        # Tab1 PropBoxes
        dev_box = initDevBox.initDevBox(self)
        btv_box = initBTVBox.initBTVBox(self)
        prop_box = initPropBox.initPropBox(self)
        image_box = initImageBox.initImageBox(self)


        # Tab2 PropBoxes
        bg_box = initBGBox.initBGBox(self)
        filt_box = initFiltBox.initFiltBox(self)
        zoom_box = initZoomBox.initZoomBox(self)

        # Tab3 PropBoxes
        bm_box = initBMBox.initBMBox(self)
        cent_box = initCentBox.initCentBox(self)

        # Tab4 PropBoxes
        prof_box = initProfBox.initProfBox(self)

        # Add tabs
        self.tabs = QTabWidget()
        self.tabs.addTab(
            self._tab_with_widgets([dev_box, btv_box, prop_box, image_box]),
            "Select",
        )
        self.tabs.addTab(
            self._tab_with_widgets([bg_box, filt_box, zoom_box]),
            "Process",
        )
        self.tabs.addTab(
            self._tab_with_widgets([bm_box, cent_box]),
            "Target",
        )
        self.tabs.addTab(
            self._tab_with_widgets([prof_box]),
            "Profiles",
        )
        self.tabs.addTab(self._create_save_tab(), "Save")

        # Create logbook button
        self.logbook = QPushButton("Save and Log", self)
        self.logbook.setStyleSheet("background-color:#999999")
        self.logbook.clicked[bool].connect(self.print2elog)

        # Create a plotting window
        self.view = ImageView(self, width=5, height=4, dpi=100)
        self.view.title = "Image"

        # Create Status Bar
        self.statusBar = QStatusBar(self)
        self.statusBar.setSizeGripEnabled(False)

        # Create Layout
        # Left pane: tabs + always-visible logbook button + status bar
        left_vbox = QVBoxLayout()
        left_vbox.addWidget(self.tabs)
        left_vbox.addWidget(self.logbook)
        left_vbox.addWidget(self.statusBar)

        left_widget = QWidget()
        left_widget.setLayout(left_vbox)
        left_widget.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Preferred)

        # Right pane: canvas expands to fill remaining space
        self.view.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        hbox = QHBoxLayout()
        hbox.setSpacing(6)
        hbox.addWidget(left_widget)      # stretch=0 by default → fixed width
        hbox.addWidget(self.view, 1)     # stretch=1 → takes all extra space
        self.setLayout(hbox)
        self.setGeometry(1600, 300, 1400, 700)
        self.setStyleSheet("QWidget { font-family: Sans Serif; font-size: 12pt; }")

        # Make a window
        self.setWindowTitle("Image Viewer")
        self.setWindowIcon(
            QIcon(os.path.join(os.path.dirname(__file__), "clearlogo.png"))
        )

        # Setup a timer to stop subscription after inactivity
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.stop_subs_timeout)
        self.timer.setSingleShot(True)

        # Start the show
        self.show()
        self.sub = None        # no subscription exists yet
        self.sub_state = False
        if self._select_current_camera():
            self.statusBar.showMessage("Ready")

    def cam_select(self, i):
        """Select Cam — stop current subscription, reconfigure, restart only if was running."""
        self._ask_to_move_current_screen_out()
        was_running = self.sub_state
        self.stop_subs()

        camera_ready = self._select_current_camera()
        if camera_ready:
            self.unzoomAction()
        self._reset_frame_mailbox()
        # Only auto-restart if a subscription was already active before the change.
        # If the user had manually stopped, changing camera should not silently restart.
        if was_running and camera_ready:
            self.start_subs()


    def trig_select(self, i):
        """Select Trig"""
        self._stop_and_reset_start_button()
        camera_ready = self.set_dev_prop()
        self.btn.setEnabled(camera_ready)
        self.fExtract.setEnabled(camera_ready)

    def rateAction(self):
        """Rate Action"""
        self._stop_and_reset_start_button()
        camera_ready = self.set_dev_prop()
        self.btn.setEnabled(camera_ready)
        self.fExtract.setEnabled(camera_ready)

    def _select_current_camera(self):
        """Refresh camera and movement state without letting failed devices crash."""
        camera_ready = self.set_dev_prop()
        try:
            self.readBTVScreenFilterOptions()
        except Exception as e:
            print(f"_select_current_camera movement options: {e}")
            self._set_movement_unavailable("Unavailable")

        self.btn.setEnabled(camera_ready)
        self.fExtract.setEnabled(camera_ready)
        if not camera_ready:
            self._reset_start_button()
        return camera_ready

    def set_dev_prop(self):
        """Set Device/Property"""
        try:
            self.dev_name = self.camList[self.cams.currentIndex()]
            trig_ind = self.trigs.currentIndex()
            try:
                self.fRate = int(self.rate_input.text())
            except ValueError:
                self.fRate = 1
                self.rate_input.setText("1")

            if trig_ind == 0:
                self.rate_input.setEnabled(True)
                prop_name = "LastImage"
                #if self.fRate > 4:
                #    self.japc.setSelector("")
                #else:
                #    str_rate = "{:0.0f}".format(1000 / self.fRate)
                #    self.japc.setSelector("PULL." + str_rate + "")
            elif trig_ind == 1 or trig_ind == 2:
                self.rate_input.setEnabled(False)
                prop_name = "ExtractionImage"
            elif trig_ind == 3:
                self.rate_input.setEnabled(False)
                prop_name = "OneHertzImage"
            else:
                raise RuntimeError(f"unknown trigger index {trig_ind}")

            # The JSON key is the BTV device name (e.g. 'CA.BTV0810').
            # Image acquisition uses the associated DigiCam FESA device.
            self.cam_device = self.cam_props[self.dev_name]["digiCamDevice"]
            self.dev_prop    = self.cam_device + "/" + prop_name
            self.dev_sett    = self.cam_device + "/PublishedSettings"
            self.dev_settBI  = self.cam_device + "/CameraSetting"
            self.dev_settBI2 = self.cam_device + "/CalibrationSetting"
            self.dev_roi     = self.cam_device + "/Setting"
        except Exception as e:
            print(f"set_dev_prop: {e}")
            self.can_connect = False
            self.statusBar.showMessage(f"Invalid camera configuration: {e}")
            return False

        return self.get_dev_info()

    def get_dev_info(self):
        """Get Information about device"""

        try:
            data_obj = self.async_get(self.dev_prop)
            settBI_obj = self.async_get(self.dev_settBI)
            settBI2_obj = self.async_get(self.dev_settBI2)
 
            flip_hor = self.cam_props[self.dev_name].get("flip_hor", 0)
            flip_ver = self.cam_props[self.dev_name].get("flip_ver", 0)
            rotate = self.cam_props[self.dev_name].get("rotate", 0)

            print("Looking at ", self.dev_prop)
            print("Setting rotate, flip_hor, flip_ver to: ", rotate, flip_hor, flip_ver)

            settings = settBI_obj.data
            calibration = settBI2_obj.data
            image_data = data_obj.data
            im = np.asarray(image_data["image2D"])
            if im.ndim != 2 or im.size == 0:
                raise RuntimeError(f"invalid image2D shape {im.shape}")

            self.exp_input.setText(str(settings["cameraExposureTimeUs"]))

            self.del_input.setText(str(settings["cameraTriggerDelay"]))

            self.gain_input.setText(str(settings["cameraGain"]))

            self.px_sz1 = float(calibration["pixelCalSet1"])
            self.px_sz2 = float(calibration["pixelCalSet2"])

            self.pix_input1.setText(str(self.px_sz1))
            self.pix_input2.setText(str(self.px_sz2))

            height = int(im.shape[0])
            width = int(im.shape[1])
            self.camera_image_shape = im.shape

            self.ts = 0

            t_string = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(self.ts))
            self.view.title = t_string

            ax1, ax2 = self._image_axes(width, height)
            
            self.view.ax1 = ax1
            self.view.ax2 = ax2
            self.view.roi = [ax1[0], ax1[-1], ax2[0], ax2[-1]]
            self.view.proj_roi = [ax1[0], ax1[-1], ax2[0], ax2[-1]]
            self.view.proj_roi_manual = False
            self.view.im = im
            self.view.camera_name = self.dev_name
            self.view.bg = np.zeros(im.shape)
            self._load_current_beam_mark()

            self.view.fliplr = flip_hor
            self.view.flipud = flip_ver
            self.view.rotate = rotate
            #self.view.GUI2 = self.GUI2
            self._sync_view_options(sub_bg=False)

            if self.view.rotate == 1 or self.view.rotate == 3:
                print("rotate axes and roi")
                self.view.ax1 = ax2
                self.view.ax2 = ax1
                self.view.roi = [ax2[0], ax2[-1], ax1[0], ax1[-1]]
                self.view.frame_ana.frame_analyzed = False  # Remember to rotate

            print("done rotating axes and roi")

            self.view.update_figure()

            self.can_connect = True
            print("can_connect: ", self.can_connect)
            self.btn.setEnabled(True)
            self.fExtract.setEnabled(True)
            return True

        except Exception as e:
            print(f"get_dev_info: {e}")
            self.can_connect = False
            self.btn.setEnabled(False)
            self.fExtract.setEnabled(False)
            self._reset_start_button()
            self.statusBar.showMessage(
                "Unable to connect to device " + self.cams.currentText()
            )
            return False

    def _image_axes(self, width, height):
        """Build physical image axes from current pixel calibration."""
        ax1 = self.px_sz1 * np.linspace(-width / 2, width / 2, width)
        ax2 = self.px_sz2 * np.linspace(-height / 2, height / 2, height)
        return ax1, ax2

    def _apply_pixel_axes_to_view(self):
        """Rebuild displayed physical axes after pixel calibration changes."""
        if not hasattr(self, "camera_image_shape"):
            return

        height, width = self.camera_image_shape
        ax1, ax2 = self._image_axes(width, height)

        if self.view.rotate == 1 or self.view.rotate == 3:
            ax1, ax2 = ax2, ax1

        self.view.ax1 = ax1
        self.view.ax2 = ax2
        full_roi = self._full_image_roi()
        self.view.roi = full_roi
        self.view.proj_roi = full_roi
        self.view.proj_roi_manual = False
        self.view.frame_ana.frame_analyzed = True
        self.view.update_figure()

    def _axes_ready(self):
        return (
            np.asarray(getattr(self.view, "ax1", [])).size > 0
            and np.asarray(getattr(self.view, "ax2", [])).size > 0
        )

    def _rebuild_axes_from_current_image(self):
        image = np.asarray(getattr(self.view, "im", []))
        if image.ndim != 2 or image.size == 0:
            return False

        height, width = image.shape
        self.camera_image_shape = image.shape
        ax1, ax2 = self._image_axes(width, height)
        if self.view.rotate == 1 or self.view.rotate == 3:
            ax1, ax2 = ax2, ax1
        self.view.ax1 = ax1
        self.view.ax2 = ax2
        return self._axes_ready()

    def _full_image_roi(self):
        if not self._axes_ready() and not self._rebuild_axes_from_current_image():
            raise RuntimeError("image axes are not available")
        return [self.view.ax1[0], self.view.ax1[-1], self.view.ax2[0], self.view.ax2[-1]]

    def _set_movement_unavailable(self, text="None"):
        if hasattr(self, "position_readback"):
            self.position_readback.setText(text)
        if hasattr(self, "filter_readback"):
            self.filter_readback.setText(text)
        if hasattr(self, "position_dropdown"):
            self.position_dropdown.blockSignals(True)
            self.position_dropdown.clear()
            self.position_dropdown.addItems(["No options"])
            self.position_dropdown.blockSignals(False)
        if hasattr(self, "filter_dropdown"):
            self.filter_dropdown.blockSignals(True)
            self.filter_dropdown.clear()
            self.filter_dropdown.addItems(["No options"])
            self.filter_dropdown.blockSignals(False)

    def _sync_view_options(self, sub_bg=None):
        """Copy current UI options into the ImageView."""
        self.view.autoscale = self.autoState
        self.view.cax = [self.cmin, self.cmax]
        self.view.cmap = self.cmap_combo.currentText()
        self.view.appMedFilt = self.filt_use.isChecked()
        self.view.appFitGauss = self.gauss_use.isChecked()
        self.view.sub_bg = self.bg_use.isChecked() if sub_bg is None else sub_bg

        self.view.show_bm = self.bm_show.isChecked()
        self.view.show_bmline = self.bmline_show.isChecked()
        self.view.show_bm_avg = self.bm_show_avg.isChecked()
        self.view.text_bm = self.bm_text.isChecked()
        self.view.textdelta_bm = self.bm_textdelta.isChecked()
        self.view.deltaalarm_bm = self.bm_deltaalarm.isChecked()
        self.view.show_ct = self.ct_show.isChecked()
        self.view.show_pj = self.pj_show.isChecked()
        self.view.show_pjbox = self.pjbox_show.isChecked()
        self.view.check_flatness = self.flatness_check.isChecked()
        try:
            self.view.flatness_half_width_px = max(1, int(self.flatness_w_input.text()))
        except ValueError:
            self.view.flatness_half_width_px = 5

        self.view.pj_yaxis_manual = self.pj_yaxis_manual.isChecked()
        if self.view.pj_yaxis_manual:
            try:
                self.view.pj_ymin = float(self.pj_ymin_input.text())
                self.view.pj_ymax = float(self.pj_ymax_input.text())
            except ValueError:
                self.view.pj_yaxis_manual = False

        if self.axis_equal.isChecked():
            self.view.aspect = "equal"
        elif self.axis_auto.isChecked():
            self.view.aspect = "auto"

        self.view.LogScale = self.scale_combo.currentText() == "Log"

    def cRangeAction(self):
        """Change caxis range"""
        try:
            self.cmin = int(self.cmin_input.text())
            self.cmax = int(self.cmax_input.text())
        except ValueError:
            self.statusBar.showMessage("Invalid color range")

    def autoscale(self):
        """Change autoscaling"""
        self.autoState = self.autoOn.isChecked()
        if self.autoState:
            self.cmin_input.setEnabled(False)
            self.cmax_input.setEnabled(False)
        else:
            self.cmin_input.setEnabled(True)
            self.cmax_input.setEnabled(True)
            self.cRangeAction()

    def async_set(self, prop, val):
        print(f'SET on {prop} with {val} called')
        """Async Set"""
        attempts = [
            (prop, val),
        ]
        if '#' in prop:
            prop_devname, prop_fieldname = prop.split("#", 1)
            attempts.append((prop_devname, {prop_fieldname: val}))

        errors = []
        for param, data in attempts:
            try:
                return self.client_get.set(param, data=data)
            except TypeError as e:
                errors.append(e)
                try:
                    return self.client_get.set(param, data)
                except Exception as pos_e:
                    errors.append(pos_e)
            except Exception as e:
                errors.append(e)

        raise errors[-1]

    def async_get(self, prop):
        """Async Get"""
        val = self.client_get.get(prop)
        return val

    def _property_data_dict(self, data):
        """Return pyda property data as a mutable plain dict."""
        if hasattr(data, "items"):
            return dict(data.items())
        return dict(data)

    def _set_pixel_calibration(self, px_sz1, px_sz2):
        current = self._property_data_dict(self.async_get(self.dev_settBI2).data)
        current["pixelCalSet1"] = px_sz1
        current["pixelCalSet2"] = px_sz2

        attempts = [
            lambda: self.async_set(self.dev_settBI2, current),
            lambda: self.async_set(
                self.dev_settBI2,
                {"pixelCalSet1": px_sz1, "pixelCalSet2": px_sz2},
            ),
            lambda: (
                self.async_set(self.dev_settBI2 + "#pixelCalSet1", px_sz1),
                self.async_set(self.dev_settBI2 + "#pixelCalSet2", px_sz2),
            ),
        ]

        errors = []
        for attempt in attempts:
            try:
                attempt()
                return
            except Exception as e:
                errors.append(e)

        raise errors[-1]

    def bgAction(self):
        """Get background"""
        sub_state = self.sub_state
        if sub_state:
            self.stop_subs()

        try:
            n_bg = int(self.avg_input.text())
            if n_bg <= 0:
                raise ValueError("number of background images must be positive")
            self.view.bg.fill(0)
            for i in range(n_bg):
                data_obj = self.async_get(self.dev_prop)
                bg = data_obj.data["image2D"]
                self.view.bg += bg
                time.sleep(0.1)

            self.view.bg = self.view.bg / n_bg

        except Exception as e:
            print(f"bgAction: {e}")
            self.statusBar.showMessage(
                "Unable to acquire background from " + self.cams.currentText()
            )

        self.bg_btn.setChecked(False)

        if sub_state:
            self.start_subs()

    def use_bg(self, state=None):
        """Use background"""
        if not self.bg_use.isChecked():
            return
        self.bg_show.blockSignals(True)
        self.bg_show.setChecked(False)
        self.bg_show.blockSignals(False)
        if np.sum(self.view.bg) == 0:
            self.bg_use.blockSignals(True)
            self.bg_use.setChecked(False)
            self.bg_use.blockSignals(False)
            self.statusBar.showMessage("Acquire BG first")

    def show_bg(self, state=None):
        """Show background"""
        if not self.bg_show.isChecked():
            return

        try:
            self._stop_and_reset_start_button()
            bg = np.asarray(self.view.bg)
            if bg.ndim != 2 or bg.size == 0 or np.sum(bg) == 0:
                self.bg_show.blockSignals(True)
                self.bg_show.setChecked(False)
                self.bg_show.blockSignals(False)
                self.statusBar.showMessage("Acquire BG first")
                return

            self.bg_use.blockSignals(True)
            self.bg_use.setChecked(False)
            self.bg_use.blockSignals(False)

            if self.camera_image_shape != bg.shape:
                self.camera_image_shape = bg.shape
                height, width = bg.shape
                ax1, ax2 = self._image_axes(width, height)
                if self.view.rotate == 1 or self.view.rotate == 3:
                    ax1, ax2 = ax2, ax1
                self.view.ax1 = ax1
                self.view.ax2 = ax2
                full_roi = self._full_image_roi()
                self.view.roi = full_roi
                self.view.proj_roi = full_roi
                self.view.proj_roi_manual = False

            self.view.im = bg.copy()
            self.view.sub_bg = False
            self.view.frame_ana.frame_analyzed = False
            self.view.update_figure()
        except Exception as e:
            print(f"show_bg: {e}")
            self.bg_show.blockSignals(True)
            self.bg_show.setChecked(False)
            self.bg_show.blockSignals(False)
            self.statusBar.showMessage(f"Show background failed: {e}")

    def zoomAction(self):
        """Set view zoom: user draws a box → image display is cropped to that region"""
        try:
            self._full_image_roi()
            self.view.drawRect(mode='zoom')
        except Exception as e:
            print(f"zoomAction: {e}")
            self.statusBar.showMessage(f"Cannot start zoom selection: {e}")

    def setProjRegion(self):
        """Set projection region: user draws a box → boxed projections and Gaussian fit use that region"""
        if self.flatness_check.isChecked():
            self.proj_btn.setChecked(False)
            self.statusBar.showMessage("Disable Check flatness before setting a projection region")
            return

        try:
            self._full_image_roi()
            self.view.drawRect(mode='proj')
        except Exception as e:
            print(f"setProjRegion: {e}")
            self.statusBar.showMessage(f"Cannot start projection selection: {e}")


    def bmAction(self):
        """Beam mark action"""
        try:
            self._full_image_roi()
            self.setCursor(QtCore.Qt.CursorShape.CrossCursor)
            self.view.placeBeamMark()
        except Exception as e:
            print(f"bmAction: {e}")
            self.setCursor(QtCore.Qt.CursorShape.ArrowCursor)
            self.statusBar.showMessage(f"Cannot place marker: {e}")

    def resetAverageAction(self):
        """Beam mark average reset action"""
        try:
            avg_array_size = int(self.bm_avg_input.text())
            if avg_array_size <= 0:
                raise ValueError("average size must be positive")
            self.view.init_avg_bm(size=avg_array_size)
        except Exception as e:
            print(f"resetAverageAction: {e}")
            self.statusBar.showMessage(f"Average reset failed: {e}")

    def bmCentroid(self):
        """Beam mark centroid"""
        try:
            self.view.centerBeamMark()
        except Exception as e:
            print(f"bmCentroid: {e}")
            self.statusBar.showMessage(f"Centroid marker failed: {e}")

    def markerOptionsAction(self, checked=False):
        """Apply marker visibility or geometry changes immediately."""
        try:
            self._sync_view_options()
            if getattr(self.view, "im", np.array([])).size:
                self.view.update_figure()
        except Exception as e:
            print(f"markerOptionsAction: {e}")
            self.statusBar.showMessage(f"Marker update failed: {e}")

    def _load_beam_marks(self):
        try:
            with open(self._beam_marks_path, "r", encoding="utf-8") as stream:
                data = json.load(stream)
            return data if isinstance(data, dict) else {}
        except FileNotFoundError:
            return {}
        except Exception as e:
            print(f"Unable to load beam marks: {e}")
            return {}

    def _load_current_beam_mark(self):
        mark = self.beam_marks.get(self.dev_name, {})
        self.view.x_mark = float(mark.get("x", 0.0))
        self.view.y_mark = float(mark.get("y", 0.0))
        self.view.angle_mark = float(mark.get("angle", 0.0))
        self.bm_x_input.setText(f"{self.view.x_mark:.3f}")
        self.bm_y_input.setText(f"{self.view.y_mark:.3f}")
        self.bm_a_input.setText(f"{self.view.angle_mark:.3f}")
        if mark:
            self.bm_show.blockSignals(True)
            self.bm_show.setChecked(True)
            self.bm_show.blockSignals(False)

    def _save_current_beam_mark(self):
        self.beam_marks[self.dev_name] = {
            "x": float(self.view.x_mark),
            "y": float(self.view.y_mark),
            "angle": float(self.view.angle_mark),
        }
        try:
            directory = os.path.dirname(self._beam_marks_path)
            os.makedirs(directory, exist_ok=True)
            temporary_path = self._beam_marks_path + ".tmp"
            with open(temporary_path, "w", encoding="utf-8") as stream:
                json.dump(self.beam_marks, stream, indent=2, sort_keys=True)
            os.replace(temporary_path, self._beam_marks_path)
        except Exception as e:
            print(f"Unable to save beam marks: {e}")
            self.statusBar.showMessage(f"Unable to save beam mark: {e}")

    def update_bm_x(self):
        """Beam mark action X"""
        try:
            self.view.x_mark = np.double(self.bm_x_input.text())
            self._save_current_beam_mark()
            self.markerOptionsAction()
        except Exception as e:
            print(f"update_bm_x: {e}")
            self.statusBar.showMessage(f"Marker x update failed: {e}")

    def update_bm_y(self):
        """Beam mark action Y"""
        try:
            self.view.y_mark = np.double(self.bm_y_input.text())
            self._save_current_beam_mark()
            self.markerOptionsAction()
        except Exception as e:
            print(f"update_bm_y: {e}")
            self.statusBar.showMessage(f"Marker y update failed: {e}")

    def update_bm_angle(self):
        """Beam mark action Angle"""
        try:
            self.view.angle_mark = np.double(self.bm_a_input.text())
            self._save_current_beam_mark()
            self.markerOptionsAction()
        except Exception as e:
            print(f"update_bm_angle: {e}")
            self.statusBar.showMessage(f"Marker angle update failed: {e}")

    def update_bm(self, x_mark, y_mark):
        """Beam mark action"""
        self.bm_show.setChecked(True)
        self.setCursor(QtCore.Qt.CursorShape.ArrowCursor)
        self.view.x_mark = x_mark
        self.view.y_mark = y_mark
        self.bm_x_input.setText(str(round(x_mark, 3)))
        self.bm_y_input.setText(str(round(y_mark, 3)))
        self._save_current_beam_mark()
        self.markerOptionsAction()

    def unzoomAction(self):
        """Reset view zoom to the full image (projection roi unchanged)"""
        self.zoom_btn.setChecked(False)

        def reset_zoom():
            self.view.roi = self._full_image_roi()
            self.view.frame_ana.frame_analyzed = True  # Don't rotate it again
            self.view.update_figure()

        try:
            self._run_with_subscription_paused(reset_zoom)
        except Exception as e:
            print(f"unzoomAction: {e}")
            self.statusBar.showMessage(f"Cannot reset zoom: {e}")

    def resetProjRegion(self):
        """Reset projection roi to the full image (view zoom unchanged)"""
        if self.flatness_check.isChecked():
            self.statusBar.showMessage("Disable Check flatness before resetting projection region")
            return

        self.proj_btn.setChecked(False)

        def reset_projection_region():
            self.view.proj_roi = self._full_image_roi()
            self.view.proj_roi_manual = False
            self.view.frame_ana.frame_analyzed = True  # Don't rotate it again
            self.view.update_figure()

        try:
            self._run_with_subscription_paused(reset_projection_region)
        except Exception as e:
            print(f"resetProjRegion: {e}")
            self.statusBar.showMessage(f"Cannot reset projection: {e}")

    def profilesYAxisAction(self):
        """Apply profile y-axis mode changes immediately."""
        try:
            self._sync_view_options()
            if getattr(self.view, "im", np.array([])).size:
                self.view.update_figure()
        except Exception as e:
            print(f"profilesYAxisAction: {e}")
            self.statusBar.showMessage(f"Profile y-axis update failed: {e}")

    def _set_projection_region_buttons_enabled(self, enabled):
        if hasattr(self, "proj_btn"):
            self.proj_btn.setEnabled(enabled)
            if not enabled:
                self.proj_btn.setChecked(False)
        if hasattr(self, "reset_proj_btn"):
            self.reset_proj_btn.setEnabled(enabled)

    def flatnessAction(self, checked=False):
        """Apply flatness-profile mode around the current beam-marker position."""
        try:
            if self.flatness_check.isChecked():
                self.pjbox_show.setChecked(True)

            self._set_projection_region_buttons_enabled(
                not self.flatness_check.isChecked()
            )
            self._sync_view_options()
            if getattr(self.view, "im", np.array([])).size:
                self.view.update_figure()
        except Exception as e:
            print(f"flatnessAction: {e}")
            self.statusBar.showMessage(f"Flatness check update failed: {e}")


    def setROIAction(self):
        """Let the user draw a rectangle and apply it as a hardware camera ROI."""
        try:
            self._full_image_roi()
            self._stop_and_reset_start_button()
            self.statusBar.showMessage("Draw a box to select the camera ROI")
            self.view.drawRect(mode='camera_roi')
        except Exception as e:
            print(f"setROIAction: {e}")
            self.statusBar.showMessage(f"Cannot start camera ROI selection: {e}")

    def _axis_value_to_index(self, axis, value):
        return int(np.argmin(np.abs(np.asarray(axis) - value)))

    def _camera_roi_to_indices(self, roi):
        x_min, x_max, y_min, y_max = roi
        col_min = self._axis_value_to_index(self.view.ax1, x_min)
        col_max = self._axis_value_to_index(self.view.ax1, x_max)
        row_min = self._axis_value_to_index(self.view.ax2, y_min)
        row_max = self._axis_value_to_index(self.view.ax2, y_max)
        col_min, col_max = sorted((col_min, col_max))
        row_min, row_max = sorted((row_min, row_max))

        height, width = self.camera_image_shape
        col_min = max(0, min(col_min, width - 1))
        col_max = max(0, min(col_max, width - 1))
        row_min = max(0, min(row_min, height - 1))
        row_max = max(0, min(row_max, height - 1))

        return [
            int(col_min),
            int(row_min),
            int(col_max - col_min + 1),
            int(row_max - row_min + 1),
        ]

    def _camera_setting_dict(self):
        return self._property_data_dict(self.async_get(self.dev_roi).data)

    def _set_camera_hardware_roi(self, image_roi, enabled):
        image_roi = np.asarray(image_roi, dtype=np.int32)
        setting = self._camera_setting_dict()
        setting["imageROI"] = image_roi
        setting["imageROIEnable"] = bool(enabled)

        attempts = [
            lambda: self.async_set(self.dev_roi, setting),
            lambda: self.async_set(
                self.dev_roi,
                {"imageROI": image_roi, "imageROIEnable": bool(enabled)},
            ),
            lambda: (
                self.async_set(self.dev_roi + "#imageROI", image_roi),
                self.async_set(self.dev_roi + "#imageROIEnable", bool(enabled)),
            ),
        ]

        errors = []
        for attempt in attempts:
            try:
                attempt()
                return
            except Exception as e:
                errors.append(e)

        raise errors[-1]

    def _full_camera_hardware_roi(self):
        setting = self._camera_setting_dict()
        if "imageWindow" in setting:
            window = np.asarray(setting["imageWindow"], dtype=np.int32)
            if window.size == 4:
                return [int(v) for v in window]

        height, width = self.camera_image_shape
        return [0, 0, int(width), int(height)]

    def resetROIAction(self):
        """Disable the hardware camera ROI and restore the full image window."""
        try:
            self._stop_and_reset_start_button()
            full_roi = self._full_camera_hardware_roi()
            self._set_camera_hardware_roi(full_roi, False)
            self.view.camera_roi = None
            self.statusBar.showMessage(f"Camera ROI reset to {full_roi}")
        except Exception as e:
            print(f"resetROIAction: {e}")
            self.statusBar.showMessage(f"Camera ROI reset failed: {e}")

    def cameraROICallback(self, roi):
        image_roi = self._camera_roi_to_indices(roi)
        try:
            self._set_camera_hardware_roi(image_roi, True)
            self.view.camera_roi = roi
            self.statusBar.showMessage(f"Camera ROI set to {image_roi}")
        except Exception as e:
            print(f"cameraROICallback: {e}")
            self.statusBar.showMessage(f"Camera ROI set failed: {e}")
    

    def start_subs(self):
        """Start subscription for the current device.

        Always stops any existing subscription first so that self.sub is never
        silently orphaned (which would leave the old callback firing indefinitely).
        """
        if not self.can_connect:
            self._reset_start_button()
            self.statusBar.showMessage(
                "Unable to connect to device " + self.cams.currentText()
            )
            return False

        # Defensive stop: if a previous subscription exists (e.g. start_subs
        # was called twice, or stop_subs previously failed silently), kill it
        # before creating a new one so we never lose the reference.
        if self.sub is not None:
            self._do_stop_sub(self.sub)
            self.sub = None

        self.sub_state = False
        self._subscription_generation += 1
        subscription_generation = self._subscription_generation
        self._reset_frame_mailbox()
        try:
            self.sub = self.client_callback.subscribe(
                self.dev_prop,
                context='',
                callback=lambda jval, generation=subscription_generation: (
                    self.proc_cam_data(jval, generation)
                ),
            )
            self.sub.start()
        except Exception as e:
            print(f"start_subs: {e}")
            if self.sub is not None:
                self._do_stop_sub(self.sub)
                self.sub = None
            self._subscription_generation += 1
            self._reset_frame_mailbox()
            self._reset_start_button()
            self.statusBar.showMessage(
                "Unable to start subscription for " + self.cams.currentText()
            )
            return False

        # TODO (pyda docs): confirm that sub.start() / sub.stop() is sufficient
        # to fully start/stop server-side delivery, or whether an additional
        # sub.close() / sub.unsubscribe() call is needed to release resources
        # and guarantee that no further callbacks can fire after stop().

        self.sub_state = True
        self.zoom_btn.setChecked(False)
        self.statusBar.showMessage("Subscriptions Started")
        # TODO not sure we need this timer feature....
        self.timer.start(24 * 60 * 60 * 1000)
        self.count_shots = 0
        return True



    def proc_cam_data(self, jval, generation):
        """Receive one incoming camera frame.

        Guard against stale callbacks: after stop_subs() the subscription object
        may still deliver one or more in-flight frames before the server acknowledges
        the stop. Dropping callbacks from old subscription generations prevents
        displaying data from the previous camera after a camera switch.

        Keep only the newest waiting frame. Rendering every queued frame would
        make a high-rate camera fall progressively behind live data.
        """
        if not self.sub_state or generation != self._subscription_generation:
            return

        received_at = time.time()
        should_emit = False
        with self._frame_lock:
            if generation != self._subscription_generation:
                return
            if len(self._frame_queue) == self._frame_queue.maxlen:
                self._dropped_frames += 1
            self._frame_queue.append((jval, generation, received_at))
            if self._frame_update_pending:
                pass
            else:
                self._frame_update_pending = True
                should_emit = True

        if should_emit:
            self.camera_frame_ready.emit()

    def _process_latest_cam_data(self):
        """Render one queued frame on the Qt GUI thread."""
        with self._frame_lock:
            if self._frame_queue:
                jval, generation, received_at = self._frame_queue.popleft()
            else:
                jval = None
                generation = self._subscription_generation
                received_at = None

        if (
            not self.sub_state
            or generation != self._subscription_generation
            or jval is None
        ):
            self._reset_frame_mailbox()
            return

        try:
            self._render_cam_data(jval, received_at)
        except Exception as e:
            print(f"_process_latest_cam_data: {e}")
            self.statusBar.showMessage(f"Camera frame failed: {e}")
            self.stop_subs()
            self._reset_start_button()
        finally:
            emit_again = False
            dropped_report = None
            with self._frame_lock:
                if (
                    self.sub_state
                    and self._frame_queue
                    and self._frame_queue[0][1] == self._subscription_generation
                ):
                    emit_again = True
                else:
                    self._frame_update_pending = False

                if self._dropped_frames - self._last_drop_report >= 100:
                    self._last_drop_report = self._dropped_frames
                    dropped_report = self._dropped_frames

            if dropped_report is not None:
                print(f"Dropped {dropped_report} stale camera frames while rendering.")

            if emit_again:
                self.camera_frame_ready.emit()

    def _reset_frame_mailbox(self):
        """Clear queued camera frames and pending GUI updates."""
        with self._frame_lock:
            self._frame_queue.clear()
            self._frame_update_pending = False
            self._dropped_frames = 0
            self._last_drop_report = 0

    def _render_cam_data(self, jval, received_at):
        """Process and display one camera frame on the Qt GUI thread."""
        data = getattr(jval, "data", None)
        if data is None:
            raise RuntimeError("camera frame has no image2D field")

        try:
            im = np.asarray(data["image2D"])
        except Exception as e:
            raise RuntimeError("camera frame has no image2D field") from e
        if im.ndim != 2 or im.size == 0:
            raise RuntimeError(f"invalid image2D shape {im.shape}")

        self.view.image_saturated = self._has_saturated_cluster(im)
        shape_changed = (
            hasattr(self, "camera_image_shape")
            and self.camera_image_shape != im.shape
        )
        self.camera_image_shape = im.shape

        # Reject Extraction images in other cycles
        # TODO: AWAKE features to be cross-checked/separated form this application, 
        # or made general enough ....
        if self.trigs.currentIndex() == 2:
            try:
                cycle_name = str(data["cycleNameGlob"])
                if "AWAKE1" not in cycle_name:
                    print("Skip non-AWAKE Extraction", cycle_name)
                    self.statusBar.showMessage("Skip non-AWAKE Extraction")
                    return
            except Exception as e:
                print("Bad cycle name")
                print(e)

        try:
            image_timestamp = data["imageTimeStamp"]
        except Exception:
            image_timestamp = time.time() * 1e9
        self.ts = float(image_timestamp) / 1e9
        self.view.image_timestamp = self.ts
        self.view.frame_received_timestamp = received_at

        self.view.im = im
        self._sync_view_options()
        if shape_changed:
            self._apply_pixel_axes_to_view()

        self.view.frame_ana.frame_analyzed = False

        full_update_interval = 5
        need_full_update = (
            self.count_shots % full_update_interval == 0
            or self.view.appMedFilt
            or self.view.sub_bg
            or self.view.appFitGauss
            or self.view.show_pj
            or self.view.show_pjbox
            or self.view.show_bm_avg
            or self.view.text_bm
            or self.view.show_ct
        )
        if need_full_update:
            self.view.update_figure()
        else:
            self.view.update_image_only()
        self.count_shots += 1
        self._save_requested_frame()

    def _has_saturated_cluster(self, image):
        mask = np.asarray(image) > 4050
        if np.count_nonzero(mask) <= 100:
            return False

        labels, _ = ndi.label(mask, structure=np.ones((3, 3), dtype=bool))
        cluster_sizes = np.bincount(labels.ravel())
        return cluster_sizes.size > 1 and np.max(cluster_sizes[1:]) > 100

    def _do_stop_sub(self, sub):
        """Low-level: stop a single subscription object and log errors.

        Separated so it can be called both from stop_subs() and from the
        defensive check at the top of start_subs().

        TODO (pyda docs): check whether sub.stop() alone is enough to
        - cancel server-side delivery immediately
        - prevent any further callbacks from firing (thread-safety)
        - release server resources
        or whether sub.close() / sub.unsubscribe() / similar must also be called.
        """
        try:
            sub.stop()
            # TODO: sub.close() here if pyda requires it
        except Exception as e:
            print(f'stop_subs: error stopping subscription: {e}')

    def stop_subs(self):
        """Stop the current subscription and mark state as inactive."""
        self._subscription_generation += 1
        if self.sub is not None:
            self._do_stop_sub(self.sub)
            self.sub = None   # drop reference so it can't be accidentally reused
        self.sub_state = False
        self._reset_frame_mailbox()
        self.timer.stop()
        self.statusBar.showMessage("Subscriptions Stopped")

    def stop_subs_timeout(self):
        """Called by the inactivity timer — same as stop_subs but updates the button."""
        self.stop_subs()
        self._reset_start_button()
        self.statusBar.showMessage("Subscriptions timed out")

    """ Button Action """

    def doAction(self, pressed):
        if pressed:
            self.btn.setText("Stop")
            self.start_subs()
        else:
            self.btn.setText("Start")
            self.stop_subs()

    """Read all settings shown in the Select tab."""

    def extrAction(self):
        errors = []

        try:
            self.readBTVScreenFilterOptions()
            self.readBTVLampSettings()
        except Exception as e:
            print(f"extrAction movement settings: {e}")
            errors.append("device movement")

        try:
            self.readCameraProperties()
        except Exception as e:
            print(f"extrAction camera properties: {e}")
            errors.append("camera properties")

        if errors:
            self.statusBar.showMessage("Read failed for " + ", ".join(errors))
        else:
            self.statusBar.showMessage("All settings read")

    def readCameraProperties(self):
        """Read Camera Properties fields from FESA and update the GUI."""
        settBI_obj = self.async_get(self.dev_settBI)
        settBI2_obj = self.async_get(self.dev_settBI2)

        settings = settBI_obj.data
        calibration = settBI2_obj.data

        self.exp_input.setText(str(settings["cameraExposureTimeUs"]))
        self.del_input.setText(str(settings["cameraTriggerDelay"]))
        self.gain_input.setText(str(settings["cameraGain"]))

        self.px_sz1 = float(calibration["pixelCalSet1"])
        self.px_sz2 = float(calibration["pixelCalSet2"])
        self.pix_input1.setText(str(self.px_sz1))
        self.pix_input2.setText(str(self.px_sz2))
        self._apply_pixel_axes_to_view()

    def readBTVLampSettings(self):
        """Read voltage/output fields shown in Device Movement."""
        self.getBTVDevInfo()
        fields = self._btv_lamp_fields("ExpertSettingDCSystem")
        if fields is None:
            self.v1_input.setText("0.0")
            self.v2_input.setText("0.0")
            self.onoff_combo.blockSignals(True)
            self.onoff_combo.setCurrentText("OFF")
            self.onoff_combo.blockSignals(False)
            return

        prop, v1_field, v2_field, enable_fields = fields
        data = self.async_get(self.btvdevice + "/" + prop).data
        self.v1_input.setText(str(self._pyda_value(data[v1_field])))
        self.v2_input.setText(str(self._pyda_value(data[v2_field])))

        outputs = [
            self._pyda_value(data[field])
            for field in enable_fields
            if field in data
        ]
        if outputs:
            text = "ON" if all(self._btv_output_enabled(value) for value in outputs) else "OFF"
            self.onoff_combo.blockSignals(True)
            self.onoff_combo.setCurrentText(text)
            self.onoff_combo.blockSignals(False)

    def _btv_lamp_fields(self, prefix):
        if self.BTVDevPos == 1:
            return (
                f"{prefix}1",
                "powerChannel3",
                "powerChannel4",
                ("enableChannel3", "enableChannel4"),
            )
        if self.BTVDevPos == 2:
            return (
                f"{prefix}2",
                "powerChannel7",
                "powerChannel8",
                ("enableChannel7", "enableChannel8"),
            )
        return None

    def _pyda_value(self, value):
        if hasattr(value, "available_bit_names") and hasattr(value, "is_bit_high"):
            return {
                name: value.is_bit_high(name)
                for name in value.available_bit_names()
            }
        if hasattr(value, "items"):
            return dict(value.items())
        return getattr(value, "value", value)

    def _btv_output_enabled(self, value):
        if isinstance(value, dict):
            if "ENABLE_MOVEMENT" in value:
                return bool(value["ENABLE_MOVEMENT"])
            return any(
                bool(state)
                for name, state in value.items()
                if "ENABLE" in name and "DISABLE" not in name
            )
        try:
            return float(value) != 0
        except Exception:
            return bool(value)

    def _btv_output_value(self, current, enabled):
        if isinstance(current, dict):
            # Treat the combo as a simple lamp output on/off, not regulation mode.
            # CA.BTV0730 readback shows lamp output ON as ENABLE_MOVEMENT=True.
            if enabled:
                return {
                    key: key == "ENABLE_MOVEMENT"
                    for key in current
                }
            return {key: False for key in current}
        return 4 if enabled else 0

    def _set_btv_lamp_power(self, field, value):
        fields = self._btv_lamp_fields("OPSettingSystem")
        if fields is None:
            raise RuntimeError("no BTVCTRL controller configured")

        prop = fields[0]
        endpoint = self.btvdevice + "/" + prop
        data = self._property_data_dict(self.async_get(endpoint).data)
        min_value = self._pyda_value(data.get(f"{field}_min", 0))
        max_value = self._pyda_value(data.get(f"{field}_max", 100))
        if value < float(min_value) or value > float(max_value):
            raise ValueError(f"{field} must be between {min_value} and {max_value}")

        errors = []
        for attempt in (
            lambda: self.client_get.set(endpoint, data={field: value}),
            lambda: self.async_set(endpoint + "#" + field, value),
            lambda: self.async_set(
                endpoint,
                {key: self._pyda_value(val) for key, val in data.items()} | {field: value},
            ),
        ):
            try:
                attempt()
                break
            except Exception as e:
                errors.append(e)
        else:
            raise errors[-1]

        check = self.async_get(endpoint).data
        return self._pyda_value(check[field])

    def _set_btv_lamp_output(self, enabled):
        fields = self._btv_lamp_fields("OPSettingSystem")
        if fields is None:
            raise RuntimeError("no BTVCTRL controller configured")

        prop, _, _, enable_fields = fields
        endpoint = self.btvdevice + "/" + prop
        data = self._property_data_dict(self.async_get(endpoint).data)
        values = {
            field: self._btv_output_value(self._pyda_value(data[field]), enabled)
            for field in enable_fields
            if field in data
        }
        if not values:
            raise RuntimeError("no lamp output enable fields found")

        errors = []
        for attempt in (
            lambda: self.client_get.set(endpoint, data=values),
            lambda: [
                self.async_set(endpoint + "#" + field, value)
                for field, value in values.items()
            ],
            lambda: self.async_set(
                endpoint,
                {key: self._pyda_value(val) for key, val in data.items()} | values,
            ),
        ):
            try:
                attempt()
                break
            except Exception as e:
                errors.append(e)
        else:
            raw_value = 4 if enabled else 0
            for attempt in (
                lambda: self.client_get.set(
                    endpoint,
                    data={field: raw_value for field in values},
                ),
                lambda: [
                    self.async_set(endpoint + "#" + field, raw_value)
                    for field in values
                ],
            ):
                try:
                    attempt()
                    break
                except Exception as e:
                    errors.append(e)
            else:
                raise errors[-1]

        check = self.async_get(endpoint).data
        return [
            self._btv_output_enabled(self._pyda_value(check[field]))
            for field in values
        ]

    """ Exposure Time Action """

    def expAction(self):
        try:
            val = np.double(self.exp_input.text())
            self.async_set(self.dev_settBI + "#cameraExposureTimeUs", val)
            check = str(self.async_get(self.dev_settBI).data["cameraExposureTimeUs"])
            self.statusBar.showMessage(str("Exposure Changed to " + check))
        except Exception as e:
            print(f"expAction: {e}")
            self.statusBar.showMessage(
                "Error. Invalid exposure. Check logs in /tmp/"
            )

    """ Trigger Delay Action """

    def delAction(self):
        try:
            val = np.double(self.del_input.text())
            self.async_set(self.dev_settBI + "#cameraTriggerDelay", val)
            check = str(self.async_get(self.dev_settBI).data["cameraTriggerDelay"])
            self.statusBar.showMessage(str("Delay Changed to " + check))
        except Exception as e:
            print(f"delAction: {e}")
            self.statusBar.showMessage(f"Delay change failed: {e}")

    """ Pixel Calibration Action """

    def pixAction(self):
        try:
            self.px_sz1 = float(self.pix_input1.text())
            self.px_sz2 = float(self.pix_input2.text())

            self._set_pixel_calibration(self.px_sz1, self.px_sz2)

            check = self.async_get(self.dev_settBI2).data
            self.px_sz1 = float(check["pixelCalSet1"])
            self.px_sz2 = float(check["pixelCalSet2"])
            self.pix_input1.setText(str(self.px_sz1))
            self.pix_input2.setText(str(self.px_sz2))
            self._apply_pixel_axes_to_view()
            self.statusBar.showMessage(
                f"Pixel Size Changed to {self.px_sz1}, {self.px_sz2}"
            )
 

        except Exception as e:
            print(f"pixAction: {e}")
            self.statusBar.showMessage(f"Pixel size change failed: {e}")

    """ Gain Action """

    def gainAction(self):
        try:
            val = np.double(self.gain_input.text())
            self.async_set(self.dev_settBI + "#cameraGain", val)
            check = str(self.async_get(self.dev_settBI).data["cameraGain"])
            self.statusBar.showMessage(str("Gain Changed to " + check))
        except Exception as e:
            print(f"gainAction: {e}")
            self.statusBar.showMessage(f"Gain change failed: {e}")

    """ Camera Reboot Action """

    def _reboot_timeout_expected(self, error):
        text = str(error).lower()
        return "timed out" in text or "timeout" in text

    def _set_start_button_running(self):
        self.btn.blockSignals(True)
        self.btn.setChecked(True)
        self.btn.setText("Stop")
        self.btn.blockSignals(False)

    def _restart_after_camera_reboot(self, camera_device, attempt=1, max_attempts=10):
        if camera_device != self.cam_device:
            return

        self.statusBar.showMessage(
            f"Waiting for {self.cam_device} after reboot ({attempt}/{max_attempts})"
        )
        if self.start_subs():
            self._set_start_button_running()
            self.statusBar.showMessage(f"{self.cam_device} image restored after reboot")
            return

        if attempt < max_attempts:
            QtCore.QTimer.singleShot(
                3000,
                lambda: self._restart_after_camera_reboot(
                    camera_device,
                    attempt + 1,
                    max_attempts,
                ),
            )
        else:
            self._reset_start_button()
            self.statusBar.showMessage(
                f"{self.cam_device} did not reconnect after reboot"
            )

    def rebootCameraAction(self, pressed=False):
        if not self.cam_device:
            self.statusBar.showMessage("No camera selected")
            return

        endpoint = self.cam_device + "/RebootCamera"
        duration_s = 3
        restart_after_reboot = self.sub_state
        camera_device = self.cam_device
        try:
            if self.sub_state:
                self.stop_subs()
                self._reset_start_button()
            self.statusBar.showMessage(f"Rebooting {self.cam_device} for {duration_s} s")

            response = self.async_set(
                endpoint,
                {"rebootDurationSec": duration_s, "powerSwitch": True},
            )
            exception = getattr(response, "exception", None)
            if exception is not None:
                raise RuntimeError(str(exception))

            self.statusBar.showMessage(
                f"Reboot requested for {self.cam_device} ({duration_s} s)"
            )
        except Exception as e:
            if not self._reboot_timeout_expected(e):
                print(f"rebootCameraAction: {e}")
                self.statusBar.showMessage(f"Camera reboot failed: {e}")
                return

            print(f"rebootCameraAction: expected timeout during reboot: {e}")
            self.statusBar.showMessage(
                f"Reboot request sent to {self.cam_device}; waiting for reconnect"
            )

        if restart_after_reboot:
            QtCore.QTimer.singleShot(
                (duration_s + 2) * 1000,
                lambda: self._restart_after_camera_reboot(camera_device),
            )

    """ V1 Voltage Action — set lamp 1 voltage """

    def v1Action(self):
        self.getBTVDevInfo()
        cam = self.cam_props.get(self.dev_name, {})
        if cam.get('fixValueLamps'):
            self.statusBar.showMessage("V1: lamp voltage is fixed for this device")
            return
        try:
            val = float(self.v1_input.text())
            fields = self._btv_lamp_fields("OPSettingSystem")
            if fields is None:
                self.statusBar.showMessage("V1: no controller configured")
                return
            check = self._set_btv_lamp_power(fields[1], val)
            self.v1_input.setText(str(check))
            self.statusBar.showMessage(f"V1 set to {float(check):.2f} %")
        except Exception as e:
            print(f'v1Action: {e}')
            self.statusBar.showMessage(f"V1 set failed: {e}")

    """ V2 Voltage Action — set lamp 2 voltage """

    def v2Action(self):
        self.getBTVDevInfo()
        cam = self.cam_props.get(self.dev_name, {})
        if cam.get('fixValueLamps'):
            self.statusBar.showMessage("V2: lamp voltage is fixed for this device")
            return
        try:
            val = float(self.v2_input.text())
            fields = self._btv_lamp_fields("OPSettingSystem")
            if fields is None:
                self.statusBar.showMessage("V2: no controller configured")
                return
            check = self._set_btv_lamp_power(fields[2], val)
            self.v2_input.setText(str(check))
            self.statusBar.showMessage(f"V2 set to {float(check):.2f} %")
        except Exception as e:
            print(f'v2Action: {e}')
            self.statusBar.showMessage(f"V2 set failed: {e}")

    """ ON/OFF Action — enable or disable lamp outputs """

    def onoffAction(self, text):
        self.getBTVDevInfo()
        try:
            enabled = text == "ON"
            states = self._set_btv_lamp_output(enabled)
            readback = "ON" if states and all(states) else "OFF"
            self.onoff_combo.blockSignals(True)
            self.onoff_combo.setCurrentText(readback)
            self.onoff_combo.blockSignals(False)
            self.statusBar.showMessage(f"Lamp output set to {readback}")
        except Exception as e:
            print(f'onoffAction: {e}')
            self.statusBar.showMessage(f"ON/OFF set failed: {e}")

    """ Button Action """

    def reboot(self, pressed):
        self.statusBar.showMessage("New FESA. Cannot. (actually never used)")

    """ Print to elogbook """

    def _image_save_dir(self):
        return "/clear/data/Python/operation/camgui/data/" + time.strftime(
            "%Y/%m/%d/", time.localtime(self.ts)
        )

    def _image_file_stem(self):
        t_string = time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime(self.ts))
        return self.dev_name.replace(".", "_") + "_" + t_string

    def _logbook_snapshot_lines(self, image_dir):
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(self.ts))
        image_shape = getattr(self.view.frame_ana, "frame", self.view.im).shape
        trigger = self.trigs.currentText() if hasattr(self, "trigs") else "unknown"
        exposure = self.exp_input.text() if hasattr(self, "exp_input") else "unknown"
        delay = self.del_input.text() if hasattr(self, "del_input") else "unknown"
        gain = self.gain_input.text() if hasattr(self, "gain_input") else "unknown"
        px_x = getattr(self, "px_sz1", "unknown")
        px_y = getattr(self, "px_sz2", "unknown")

        return [
            f"Timestamp: {timestamp}",
            f"Camera: {self.dev_name}",
            f"Logbook: {self._logbook_name()}",
            f"DigiCam device: {self.cam_device}",
            f"Image property: {self.dev_prop}",
            f"Trigger: {trigger}; exposure: {exposure} us; delay: {delay}; gain: {gain}",
            f"Pixel size: x={px_x} mm, y={px_y} mm; image shape: {image_shape[1]} x {image_shape[0]} px",
            f"Image directory: {image_dir}",
        ]

    def _save_logbook_snapshot(self, f_name, image_dir):
        fig = self.view.axes.get_figure()
        subplotpars = fig.subplotpars
        old_subplotpars = {
            "left": subplotpars.left,
            "right": subplotpars.right,
            "bottom": subplotpars.bottom,
            "top": subplotpars.top,
            "wspace": subplotpars.wspace,
            "hspace": subplotpars.hspace,
        }
        old_frameon = fig.get_frameon()
        old_facecolor = fig.get_facecolor()
        old_edgecolor = fig.get_edgecolor()
        old_fig_alpha = fig.patch.get_alpha()
        old_fig_visible = fig.patch.get_visible()
        old_axes_state = [
            (ax, ax.get_facecolor(), ax.patch.get_alpha())
            for ax in fig.axes
        ]

        lines = self._logbook_snapshot_lines(image_dir)

        header_text = fig.text(
            0.02,
            0.205,
            "Camera log snapshot",
            ha="left",
            va="top",
            fontsize=10,
            fontweight="bold",
            color="black",
        )
        left_text = fig.text(
            0.02,
            0.165,
            "\n".join(lines[:5]),
            ha="left",
            va="top",
            fontsize=8,
            color="black",
        )
        right_text = fig.text(
            0.50,
            0.165,
            "\n".join(lines[5:7]),
            ha="left",
            va="top",
            fontsize=8,
            color="black",
        )
        path_text = fig.text(
            0.02,
            0.04,
            lines[7].replace(": ", ":\n", 1),
            ha="left",
            va="bottom",
            fontsize=7,
            color="black",
        )
        metadata_artists = [header_text, left_text, right_text, path_text]

        try:
            fig.set_frameon(True)
            fig.patch.set_visible(True)
            fig.patch.set_facecolor("white")
            fig.patch.set_edgecolor("white")
            fig.patch.set_alpha(1.0)
            for ax in fig.axes:
                ax.set_facecolor("white")
                ax.patch.set_alpha(1.0)
            fig.subplots_adjust(bottom=0.28)
            if hasattr(self.view, "_sync_colorbar_position"):
                self.view._sync_colorbar_position()
            fig.savefig(
                f_name,
                dpi=150,
                facecolor="white",
                edgecolor="white",
                transparent=False,
            )
            from PIL import Image
            with Image.open(f_name) as image:
                if image.mode in ("RGBA", "LA"):
                    flattened = Image.new("RGB", image.size, "white")
                    flattened.paste(image.convert("RGB"), mask=image.getchannel("A"))
                    flattened.save(f_name)
                elif image.mode != "RGB" or "transparency" in image.info:
                    image.convert("RGB").save(f_name)
        finally:
            for artist in metadata_artists:
                artist.remove()
            fig.subplots_adjust(**old_subplotpars)
            fig.set_frameon(old_frameon)
            fig.patch.set_visible(old_fig_visible)
            fig.patch.set_facecolor(old_facecolor)
            fig.patch.set_edgecolor(old_edgecolor)
            fig.patch.set_alpha(old_fig_alpha)
            for ax, facecolor, alpha in old_axes_state:
                ax.set_facecolor(facecolor)
                ax.patch.set_alpha(alpha)
            fig.canvas.draw_idle()

    def print2elog(self, pressed):
        sub_state = self.sub_state
        self.stop_subs()
        self.statusBar.showMessage("Sending to logbook")
        image_dir = self._image_save_dir()
        f_name = "/tmp/" + self._image_file_stem() + ".png"
        logbook_name = self._logbook_name()
        sent_to_logbook = False
        try:
            try:
                self._save_logbook_snapshot(f_name, image_dir)
                try:
                    elog = self._ensure_logbook(logbook_name)
                    if elog is not None:
                        entry = elog.add_event(self.dev_name)
                        entry.attach_file(f_name)
                        sent_to_logbook = True
                except Exception as e:
                    print(f"print2elog: {e}")
                    time.sleep(0.5)
                    self.statusBar.showMessage("Failed to send to logbook")
                    time.sleep(1)
            except Exception as e:
                print(f"print2elog: {e}")
                time.sleep(0.5)
                self.statusBar.showMessage("Failed to send to logbook")
                time.sleep(1)
            finally:
                if os.path.exists(f_name):
                    os.remove(f_name)

            try:
                self.save_image()
            except Exception as e:
                print(f"print2elog save_image: {e}")
                self.statusBar.showMessage(f"Logbook sent; image save failed: {e}")
                return

            if sent_to_logbook:
                self.statusBar.showMessage(f"Saved and sent to {logbook_name}!")
            else:
                self.statusBar.showMessage("Saved image; logbook failed")
        finally:
            if sub_state:
                self.start_subs()

    """ Save to .mat file """

    def _unique_image_path(self, image_dir, suffix=""):
        stem = self._image_file_stem()
        candidate = os.path.join(image_dir, f"{stem}{suffix}.mat")
        index = 1
        while os.path.exists(candidate):
            candidate = os.path.join(image_dir, f"{stem}{suffix}_{index:03d}.mat")
            index += 1
        return candidate

    def _current_image_data(self):
        if hasattr(self.view, "_analysis_input_frame"):
            image = self.view._analysis_input_frame()[0]
        else:
            image = self.view.frame_ana.frame
        return {
            "image": image,
            "x_ax": self.view.ax1,
            "y_ax": self.view.ax2,
        }

    def _write_image_mat(self, m_name):
        sio.savemat(m_name, {"data_struct": self._current_image_data()})

    def save_image(self, image_dir=None, suffix=""):
        # TODO: this must be some configuration external to the code
        path = image_dir or self._image_save_dir()
        if not os.path.exists(path):
            os.makedirs(path)

        m_name = self._unique_image_path(path, suffix)
        self._write_image_mat(m_name)
        return m_name

    def saveNImagesAction(self, pressed=False):
        try:
            total = int(self.save_n_input.text())
        except ValueError:
            self.statusBar.showMessage("Enter a valid number of images")
            return

        if total <= 0:
            self.statusBar.showMessage("Enter at least one image")
            return

        self._multi_save_total = total
        self._multi_save_remaining = total
        self._multi_save_camera = self.dev_name
        self._multi_save_dir = None
        self._multi_save_paths = []
        self.save_n_status.setText(
            f"Saving {total} images from {self._multi_save_camera}..."
        )
        self.statusBar.showMessage(
            f"Saving next {total} images from {self._multi_save_camera}"
        )

        if not self.sub_state:
            self.btn.setChecked(True)
            self.btn.setText("Stop")
            if not self.start_subs():
                self._cancel_multi_save("camera not connected")

    def _save_requested_frame(self):
        if self._multi_save_remaining <= 0:
            return
        if self.dev_name != self._multi_save_camera:
            self._cancel_multi_save("camera changed")
            return

        try:
            if self._multi_save_dir is None:
                self._multi_save_dir = self._image_save_dir()
            saved_path = self.save_image(self._multi_save_dir, suffix="_seq")
        except Exception as e:
            print(f"_save_requested_frame: {e}")
            self._cancel_multi_save(f"save failed: {e}")
            return

        self._multi_save_paths.append(saved_path)
        self._multi_save_remaining -= 1
        saved = self._multi_save_total - self._multi_save_remaining
        self.save_n_status.setText(
            f"Saved {saved}/{self._multi_save_total} images from {self._multi_save_camera}"
        )

        if self._multi_save_remaining == 0:
            self._finish_multi_save()

    def _cancel_multi_save(self, reason):
        total = self._multi_save_total
        saved = len(self._multi_save_paths)
        self._multi_save_total = 0
        self._multi_save_remaining = 0
        self._multi_save_camera = None
        self._multi_save_dir = None
        self._multi_save_paths = []
        self.save_n_status.setText(f"Save stopped: {reason}")
        self.statusBar.showMessage(f"Saved {saved}/{total} images; {reason}")

    def _finish_multi_save(self):
        total = self._multi_save_total
        camera = self._multi_save_camera
        image_dir = self._multi_save_dir
        saved_paths = list(self._multi_save_paths)

        self._multi_save_total = 0
        self._multi_save_remaining = 0
        self._multi_save_camera = None
        self._multi_save_dir = None
        self._multi_save_paths = []

        self.save_n_status.setText(
            f"Saved {total} images from {camera} in {image_dir}"
        )
        self.statusBar.showMessage(f"Saved {total} images; sending logbook note")
        self._log_multi_save_note(camera, total, image_dir, saved_paths)

    def _log_multi_save_note(self, camera, total, image_dir, saved_paths):
        logbook_name = self._logbook_name()
        note = f"{total} images of {camera} were saved in folder {image_dir}"
        try:
            elog = self._ensure_logbook(logbook_name)
            if elog is None:
                self.statusBar.showMessage(
                    f"Saved {total} images; logbook unavailable"
                )
                return
            entry = elog.add_event(note)
            if saved_paths:
                note_path = os.path.join(
                    "/tmp",
                    f"{camera.replace('.', '_')}_save_note.txt",
                )
                with open(note_path, "w", encoding="utf-8") as note_file:
                    note_file.write(note + "\n")
                    note_file.write("Files:\n")
                    for path in saved_paths:
                        note_file.write(path + "\n")
                try:
                    entry.attach_file(note_path)
                finally:
                    if os.path.exists(note_path):
                        os.remove(note_path)
            self.statusBar.showMessage(
                f"Saved {total} images and logged to {logbook_name}"
            )
        except Exception as e:
            print(f"_log_multi_save_note: {e}")
            self.statusBar.showMessage(f"Saved {total} images; logbook note failed")

    """ GTFO """

    def closeEvent(self, event):
        self.stop_subs()   # ensure subscription is stopped before the window closes
        QWidget.closeEvent(self, event)

    """ Move the BTV screen """

    def positionScreenAction(self, text):
        self.getBTVDevInfo()

        if not self.screenInstalled:
            self.position_readback.setText('None')
            return

        label = self.position_dropdown.currentText()
        if label == 'No options':
            self.position_readback.setText('None')
            return

        # ---- Custom stepper mover (e.g. CA.BTV0875 / BStepMotorVME) ----
        if self.has_custom_screen_mover:
            # Custom mover: read raw position value and reverse-map to label
            if self.screenMoverType in ['BStepMotorVME', 'NewFocusPicomotor']:
                try: 
                    setpoints = self.screenMoverFields.get('setpoints', {})
                    if label not in setpoints:
                        self.position_readback.setText('No setpoint')
                        self.statusBar.showMessage(f"No setpoint configured for {label}")
                        return
                    if self.screenMoverType == 'BStepMotorVME':
                        self.async_set(self.screenMoverDevice + '/Move', {'mode': int(2), 'value': int(setpoints[label]), 'units': int(2)})
                    elif self.screenMoverType == 'NewFocusPicomotor':
                        self.async_set(self.screenMoverDevice + '/Setting#position', int(setpoints[label]))
                    else:
                        print(f'positionScreenAction: unknown mover type {self.screenMoverType!r}')
                        self.position_readback.setText('Unknown type')

                except Exception as e:
                    print(f'positionScreenAction (mover): {e}')
                    self.position_readback.setText('Error')
            return

        # ---- Standard BTVCTRL screen command ----
        if self.BTVDevPos == 1:
            prop = 'OPSettingSystem1#positionChannel1'
        elif self.BTVDevPos == 2:
            prop = 'OPSettingSystem2#positionChannel5'
        else:
            self.position_readback.setText('No ctrl')
            return

        if not self.btvdevice:
            self.position_readback.setText('No ctrl')
            return

        index = self.position_dropdown.currentIndex()
        try:
            self.async_set(self.btvdevice + '/' + prop, int(index))
            self.position_readback.setText(label)
        except Exception as e:
            print(f'positionScreenAction (BTVCTRL): {e}')
            self.position_readback.setText('None')


    """ Move the BTV filter """

    def filterScreenAction(self, text):
        # Filter is always driven by BTVCTRL, even on cameras with a custom screen mover.
        self.getBTVDevInfo()
        if not self.btvdevice:
            self.filter_readback.setText('No ctrl')
            return
        if self.BTVDevPos == 1:
            prop = 'OPSettingSystem1#positionChannel2'
        elif self.BTVDevPos == 2:
            prop = 'OPSettingSystem2#positionChannel6'
        else:
            self.filter_readback.setText('No ctrl')
            return

        index = self.filter_dropdown.currentIndex()
        try:
            self.async_set(self.btvdevice + '/' + prop, int(index))
            self.filter_readback.setText(self.filter_dropdown.itemText(index))
        except Exception as e:
            print(f'filterScreenAction: {e}')
            self.filter_readback.setText('None')


    """ Get all the BTV parameters (screen position + filter position) """

    def getBTVdevpos(self):
        self.getBTVDevInfo()

        # ---- Screen readback ----
        if not self.screenInstalled:
            self.position_readback.setText('None')
        elif self.has_custom_screen_mover:
            # Custom mover: read raw position value and reverse-map to label
            if self.screenMoverType in ['BStepMotorVME', 'NewFocusPicomotor']:
                try:
                    raw_pos  = int(self.async_get(self.screenMoverDevice + '/Acquisition').data['position'])
                    setpoints = self.screenMoverFields.get('setpoints', {})
                    # Build reverse map: integer value → label
                    rev = {v: k for k, v in setpoints.items()}
                    label = rev.get(raw_pos, str(raw_pos))
                    self.position_readback.setText(label)
                    # Sync the dropdown to the current position if the label exists
                    idx = self.position_dropdown.findText(label)
                    if idx >= 0:
                        self.position_dropdown.blockSignals(True)
                        self.position_dropdown.setCurrentIndex(idx)
                        self.position_dropdown.blockSignals(False)
                except Exception as e:
                    print(f'getBTVdevpos (mover screen): {e}')
                    self.position_readback.setText('None')

        elif self.screenInstalled and self.btvdevice:
            # Standard BTVCTRL screen readback
            try:
                if self.BTVDevPos == 1:
                    status_screen = self.async_get(self.btvdevice + '/ExpertSettingDCSystem1').data['positionChannel1']
                    self.position_readback.setText(self._combo_text(self.position_dropdown, int(status_screen.value)))

                elif self.BTVDevPos == 2:
                    status_screen = self.async_get(self.btvdevice + '/ExpertSettingDCSystem2').data['positionChannel5']
                    self.position_readback.setText(self._combo_text(self.position_dropdown, int(status_screen.value)))

                else:
                    self.position_readback.setText('None')

            except Exception as e:
                print(f'getBTVdevpos (BTVCTRL screen): {e}')
                self.position_readback.setText('None')
        else:
            self.position_readback.setText('None')
 

        # ---- Filter readback (always BTVCTRL) ----
        if self.filterInstalled and self.btvdevice:
            try:
                if self.BTVDevPos == 1:
                    status_filter = self.async_get(self.btvdevice + '/ExpertSettingDCSystem1').data['positionChannel2']
                    self.filter_readback.setText(self._combo_text(self.filter_dropdown, int(status_filter.value)))

                elif self.BTVDevPos == 2:
                    status_filter = self.async_get(self.btvdevice + '/ExpertSettingDCSystem2').data['positionChannel6']
                    self.filter_readback.setText(self._combo_text(self.filter_dropdown, int(status_filter.value)))
                else:
                    self.filter_readback.setText('None')
            except Exception as e:
                print(f'getBTVdevpos (filter): {e}')
                self.filter_readback.setText('None')

        else:
            self.filter_readback.setText('None')


    """ See what screens/filters are available on a given BTV """

    def _option_names(self, values):
        if values is None:
            return []
        return [str(value).strip() for value in list(values) if str(value).strip()]

    def _combo_text(self, combo, index):
        if 0 <= index < combo.count():
            return combo.itemText(index)
        return 'None'

    def readBTVScreenFilterOptions(self):
        self.getBTVDevInfo()

        # ---- Screen options ----
        try:
            if not self.screenInstalled:
                self.screen_names = []
            elif self.has_custom_screen_mover:
                if self.screenMoverFields:
                    setpoints = self.screenMoverFields.get('setpoints', {})
                    self.screen_names = self._option_names(setpoints.keys())
                else:
                    self.screen_names = []
            elif self.screenInstalled and self.btvdevice and self.BTVDevPos == 1:
                self.screen_names = self._option_names(
                    self.async_get(self.btvdevice + '/Description').data['dcm1DriverNames']
                )
            elif self.screenInstalled and self.btvdevice and self.BTVDevPos == 2:
                self.screen_names = self._option_names(
                    self.async_get(self.btvdevice + '/Description').data['dcm3DriverNames']
                )
            else:
                self.screen_names = []
        except Exception as e:
            print(f"readBTVScreenFilterOptions (screen): {e}")
            self.screen_names = []

        # ---- Filter options (always BTVCTRL) ----
        try:
            if self.filterInstalled and self.btvdevice:
                if self.BTVDevPos == 1:
                    self.filter_names = self._option_names(
                        self.async_get(self.btvdevice + '/Description').data['dcm2DriverNames']
                    )
                elif self.BTVDevPos == 2:
                    self.filter_names = self._option_names(
                        self.async_get(self.btvdevice + '/Description').data['dcm4DriverNames']
                    )
                else:
                    self.filter_names = []
            else:
                self.filter_names = []
        except Exception as e:
            print(f"readBTVScreenFilterOptions (filter): {e}")
            self.filter_names = []

        if not self.screen_names:
            self.screen_names = ['No options']
        if not self.filter_names:
            self.filter_names = ['No options']

        self.position_dropdown.blockSignals(True)
        self.position_dropdown.clear()
        self.position_dropdown.addItems(self.screen_names)
        self.position_dropdown.blockSignals(False)

        self.filter_dropdown.blockSignals(True)
        self.filter_dropdown.clear()
        self.filter_dropdown.addItems(self.filter_names)
        self.filter_dropdown.blockSignals(False)
        self.getBTVdevpos()

    def getBTVDevInfo(self):
        """Look up the BTV controller and screen-mover info for the current camera."""
        cam = self.cam_props.get(self.dev_name, {})

        # ---- Filter/lamp controller (BTVCTRL) ----
        self.btvdevice  = cam.get('controlDeviceName')   # None when not configured
        ctrl_type       = cam.get('controlDeviceType')
        ctrl_fields     = cam.get('controlDeviceFields', {})
        if ctrl_type == 'BTVCTRL':
            self.BTVDevPos = ctrl_fields.get('system', 0)
        else:
            # Future controller types will be handled here
            self.BTVDevPos = 0

        # ---- Custom screen mover (e.g. BStepMotorVME) ----
        self.screenMoverDevice = cam.get('screenMoverDevice')   # None for most cameras
        self.screenMoverType   = cam.get('screenMoverType')
        self.screenMoverFields = cam.get('screenMoverFields') or {}
        self.screenInstalled = bool(cam.get('screenInstalled'))
        self.has_custom_screen_mover = (
            self.screenInstalled
            and isinstance(self.screenMoverDevice, str)
        )
        self.filterInstalled = bool(cam.get('filterInstalled'))

        def _move_screen(self, screen_name, requested_position):
            control_device, prop, positions = (
                self._get_screen_movement_info(screen_name)
            )

            normalized_positions = [
                position.lower()
                for position in positions
            ]

            requested_position = requested_position.lower()

            if requested_position not in normalized_positions:
                raise ValueError(
                    f"Position {requested_position!r} not available for "
                    f"{screen_name}. Available positions: {positions}"
                )

            index = normalized_positions.index(requested_position)

            self.client.set(
                f"{control_device}/{prop}",
                int(index),
                context=self.context_empty,
            )

            self.log(
                f"Moving {screen_name} to {positions[index]!r}, "
                f"index {index}"
            )

            return positions[index]









