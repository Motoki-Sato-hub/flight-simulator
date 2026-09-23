import os
try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (QVBoxLayout, QDialog, QLabel, QPlainTextEdit)
    from PyQt6 import uic
except ImportError:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (QVBoxLayout, QDialog, QLabel, QPlainTextEdit)
    from PyQt5 import uic

BOUNDS_POPUP_UI = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "UI files", "BoundsPopup.ui"))

class BoundsForParameter(QDialog):
    def __init__(self, parameter_name, parent=None, unit="", hint=""):
        super().__init__(parent)
        uic.loadUi(BOUNDS_POPUP_UI, self)
        self.setWindowTitle(f"Fit bounds - {parameter_name}")
        self.parameter_label.setText(f"{parameter_name} [{unit}]" if unit else str(parameter_name))
        self.hint_label.setText(hint)
        for spinbox in (self.lower_bound, self.higher_bound):
            spinbox.setSuffix(f" {unit}" if unit else "")
        self.button_apply.clicked.connect(self.accept)
        self.button_cancel.clicked.connect(self.reject)

    def set_values(self,bound_low,bound_high):
        self.lower_bound.setValue(bound_low)
        self.higher_bound.setValue(bound_high)

    def get_values(self):
        low = float(self.lower_bound.value())
        high = float(self.higher_bound.value())
        return min(low, high), max(low, high)
