try:
    from PyQt6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QDialogButtonBox, QRadioButton, QLabel,QMessageBox,QPushButton, QButtonGroup
    from PyQt6.QtCore import QEvent, Qt
    from PyQt6.QtGui import QPixmap
except ImportError:
    from PyQt5.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QDialogButtonBox, QRadioButton, QLabel,QMessageBox, QPushButton, QButtonGroup
    from PyQt5.QtCore import QEvent, Qt
    from PyQt5.QtGui import QPixmap

import importlib
import importlib.util
from pathlib import Path
from Interfaces.interface_setup import INTERFACE_SETUP


ASSETS_DIR = Path(__file__).resolve().parent.parent / "UI files" / "Assets"
ACCELERATOR_BRANDING = {
    "ATF2": ("KEK", "atf2_logo.png", "kek_logo.png"),
    "FACET2": ("SLAC", "facet2_logo.png", "slac_logo.png"),
    "CLEAR": ("CERN", "clear_logo.png", "CERN_logo.png"),
}

def _no_focus_policy():
    return Qt.FocusPolicy.NoFocus if hasattr(Qt, "FocusPolicy") else Qt.NoFocus

def _strong_focus_policy():
    return Qt.FocusPolicy.StrongFocus if hasattr(Qt, "FocusPolicy") else Qt.StrongFocus


def _keep_aspect_ratio():
    return Qt.AspectRatioMode.KeepAspectRatio if hasattr(Qt, "AspectRatioMode") else Qt.KeepAspectRatio


def _smooth_transformation():
    return Qt.TransformationMode.SmoothTransformation if hasattr(Qt, "TransformationMode") else Qt.SmoothTransformation


def _branding_header(accelerator):
    branding = ACCELERATOR_BRANDING.get(accelerator)
    if branding is None:
        return None

    institute, accelerator_logo, institute_logo = branding
    header = QHBoxLayout()
    header.setContentsMargins(0, 0, 0, 6)
    header.setSpacing(12)
    header.addStretch()

    for logo_name in (accelerator_logo, institute_logo):
        pixmap = QPixmap(str(ASSETS_DIR / logo_name))
        if pixmap.isNull():
            continue
        logo = QLabel()
        logo.setPixmap(pixmap.scaled(140, 58, _keep_aspect_ratio(), _smooth_transformation()))
        logo.setAlignment(Qt.AlignmentFlag.AlignCenter if hasattr(Qt, "AlignmentFlag") else Qt.AlignCenter)
        header.addWidget(logo)

    header.addStretch()
    return header

def get_available_interface_setup():
    available_setup = {}
    for machine_name, entries in INTERFACE_SETUP.items():
        available_entries = []
        for entry in entries:
            module_name = entry.get("module", "")
            class_name = entry.get("class_name", "")
            if not module_name:
                continue
            try:
                module_spec = importlib.util.find_spec(module_name)
            except (ModuleNotFoundError, ValueError, ImportError):
                continue

            if module_spec is None:
                continue
            try:
                module = importlib.import_module(module_name)
            except Exception:
                continue
            if not hasattr(module, class_name):
                continue
            available_entries.append(entry)
        if available_entries:
            available_setup[machine_name] = available_entries
    return available_setup

class SelectAcc(QDialog):
    def __init__(self,machines,parent=None):
        super().__init__(parent)
        self.setFocusPolicy(_strong_focus_policy())
        self.setWindowTitle("Select a machine")
        self.selected_machine = None

        layout = QVBoxLayout(self)

        prompt_label = QLabel("Choose one of the following accelerators:")
        prompt_label.setAlignment(
            Qt.AlignmentFlag.AlignCenter
            if hasattr(Qt, "AlignmentFlag")
            else Qt.AlignCenter
        )
        layout.addWidget(prompt_label)

        self.radio_buttons = []
        self.button_group = QButtonGroup(self)
        self.button_group.setExclusive(True)

        radio_column = QVBoxLayout()

        for index, acc in enumerate(machines):
            rb = QRadioButton(acc)
            rb.setFocusPolicy(_no_focus_policy())
            self.radio_buttons.append(rb)
            self.button_group.addButton(rb, index)

            radio_column.addWidget(rb)

        centered_radio_block = QHBoxLayout()
        centered_radio_block.addStretch()
        centered_radio_block.addLayout(radio_column)
        centered_radio_block.addStretch()

        layout.addLayout(centered_radio_block)

        if self.radio_buttons:
            self.radio_buttons[0].setChecked(True) #default ATF2

        buttons = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        self.button_box = QDialogButtonBox(buttons)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        self.button_box.button(QDialogButtonBox.StandardButton.Ok).setDefault(True)
        self.button_box.setFocusPolicy(_no_focus_policy())
        for button in self.button_box.buttons():
            button.setFocusPolicy(_no_focus_policy())
        layout.addWidget(self.button_box)
        self.setFocus()

        # self.installEventFilter(self)

    def accept(self):
        for rb in self.radio_buttons:
            if rb.isChecked():
                self.selected_machine = rb.text()
                break
        super().accept()

    def _move_selection(self, step):
        if not self.radio_buttons:
            return
        current_index = next((i for i, rb in enumerate(self.radio_buttons) if rb.isChecked()), 0)
        new_index = (current_index + int(step)) % len(self.radio_buttons)
        self.radio_buttons[new_index].setChecked(True)

    def keyPressEvent(self, event):
        key = event.key()
        enter_keys = {getattr(Qt.Key, "Key_Return", None), getattr(Qt.Key, "Key_Enter", None)}
        down_keys = {getattr(Qt.Key, "Key_Down", None), getattr(Qt.Key, "Key_Right", None)}
        up_keys = {getattr(Qt.Key, "Key_Up", None), getattr(Qt.Key, "Key_Left", None)}

        if key in enter_keys:
            ok_button = self.button_box.button(QDialogButtonBox.StandardButton.Ok)
            if ok_button.isEnabled():
                ok_button.click()
                return
        if key in down_keys:
            self._move_selection(1)
            return
        if key in up_keys:
            self._move_selection(-1)
            return

        super().keyPressEvent(event)

class InterfaceSelectionDialog(QDialog):
    def __init__(self, selected_acc,parent=None):
        super().__init__(parent)
        self.setFocusPolicy(_strong_focus_policy())
        self.setWindowTitle("Select an Interface")
        self.selected_interface = None
        self.selected_interface_name = None
        self.selected_acc = selected_acc
        self.go_back=False
        self.available_interface_setup = get_available_interface_setup()
        self.are_more_machines = len(self.available_interface_setup.keys()) > 1
        self.entries = self.available_interface_setup.get(selected_acc, [])

        layout = QVBoxLayout(self)

        branding_header = _branding_header(selected_acc)
        if branding_header is not None:
            layout.addLayout(branding_header)
        prompt_label = QLabel("Choose one of the following Interfaces:")
        prompt_label.setAlignment(Qt.AlignmentFlag.AlignCenter
            if hasattr(Qt, "AlignmentFlag")
            else Qt.AlignCenter
        )
        layout.addWidget(prompt_label)
        self.radio_buttons = []
        self.button_group = QButtonGroup(self)
        self.button_group.setExclusive(True)

        radio_column = QVBoxLayout()

        for index, entry in enumerate(self.entries):
            rb = QRadioButton(entry["display_name"])
            rb.setFocusPolicy(_no_focus_policy())
            self.radio_buttons.append(rb)
            self.button_group.addButton(rb, index)

            radio_column.addWidget(rb)

        centered_radio_block = QHBoxLayout()
        centered_radio_block.addStretch()
        centered_radio_block.addLayout(radio_column)
        centered_radio_block.addStretch()

        layout.addLayout(centered_radio_block)

        if self.radio_buttons:
            self.radio_buttons[0].setChecked(True)

        buttons = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        self.button_box = QDialogButtonBox(buttons)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        if self.are_more_machines:
            self.back_button=QPushButton("Back")
            self.button_box.addButton(self.back_button,QDialogButtonBox.ButtonRole.ActionRole)
            self.back_button.clicked.connect(self._go_back)
        self.button_box.button(QDialogButtonBox.StandardButton.Ok).setDefault(True)
        self.button_box.setFocusPolicy(_no_focus_policy())
        for button in self.button_box.buttons():
            button.setFocusPolicy(_no_focus_policy())
        layout.addWidget(self.button_box)
        self.setFocus()

        # self.installEventFilter(self)

    def _go_back(self):
        self.go_back=True
        self.reject()

    def accept(self):
        selected_entry = None
        for rb, entry in zip(self.radio_buttons, self.entries):
            if rb.isChecked():
                selected_entry = entry
                break

        if selected_entry is None:
            super().accept()
            return

        module_name = selected_entry["module"]
        class_name = selected_entry["class_name"]
        settings = dict(selected_entry.get("settings", {}))

        # try:
        #     module = importlib.import_module(module_name)
        #     cls = getattr(module, class_name)
        #     self.selected_interface = cls(**settings)
        #
        #     super().accept()
        #
        # except Exception as e:
        #     QMessageBox.critical(self,"Interface unavailable",f"This interface is unavailable. {e}")


        module = importlib.import_module(module_name)
        cls = getattr(module, class_name)
        self.selected_interface = cls(**settings)

        super().accept()


    def _move_selection(self, step):
        if not self.radio_buttons:
            return
        current_index = next((i for i, rb in enumerate(self.radio_buttons) if rb.isChecked()), 0)
        new_index = (current_index + int(step)) % len(self.radio_buttons)
        self.radio_buttons[new_index].setChecked(True)

    def keyPressEvent(self, event):
        key = event.key()
        enter_keys = {getattr(Qt.Key, "Key_Return", None), getattr(Qt.Key, "Key_Enter", None)}
        down_keys = {getattr(Qt.Key, "Key_Down", None), getattr(Qt.Key, "Key_Right", None)}
        up_keys = {getattr(Qt.Key, "Key_Up", None), getattr(Qt.Key, "Key_Left", None)}
        back_keys = {getattr(Qt.Key, "Key_Escape", None), getattr(Qt.Key, "Key_Backspace", None)}

        if key in enter_keys:
            ok_button = self.button_box.button(QDialogButtonBox.StandardButton.Ok)
            if ok_button.isEnabled():
                ok_button.click()
                return

        if key in back_keys and self.are_more_machines:
            self._go_back()
            return

        if key in down_keys:
            self._move_selection(1)
            return

        if key in up_keys:
            self._move_selection(-1)
            return

        super().keyPressEvent(event)


def choose_acc_and_interface(parent=None):
    available_interface_setup = get_available_interface_setup()
    machines = sorted(available_interface_setup.keys())
    if not machines:
        QMessageBox.critical(parent, "Interface unavailable", "No available interfaces were found.")
        return None

    if len(machines)==1:
        accelerator=machines[0]
        interface_dialog=InterfaceSelectionDialog(accelerator,parent=parent)
        if interface_dialog.exec():
            return interface_dialog.selected_interface
        return None
    while True:
        acc_dialog=SelectAcc(machines,parent=parent)
        if acc_dialog.exec():
            accelerator=acc_dialog.selected_machine
        else:
            return None
        interface_dialog=InterfaceSelectionDialog(accelerator,parent=None)
        if interface_dialog.exec():
            return interface_dialog.selected_interface
        if getattr(interface_dialog,"go_back",False):
            continue
        return None
