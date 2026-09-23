try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel,
                                 QDoubleSpinBox, QPushButton, QScrollArea, QWidget, QFrame)
except ImportError:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel,
                                 QDoubleSpinBox, QPushButton, QScrollArea, QWidget, QFrame)

class ScanCurrentRanges(QDialog):
    def __init__(self, screens, default_range, ranges=None, steps=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Quadrupole scan current range per screen")
        self.screens = [str(screen) for screen in screens]
        self.default_range = (float(default_range[0]), float(default_range[1]))
        ranges = dict(ranges or {})

        layout = QVBoxLayout(self)
        hint = "Screen is scanned between minimum and maximum range."
        if steps:
            hint = f"{hint} Every screen uses {int(steps)} steps."
        hint_label = QLabel(hint)
        hint_label.setWordWrap(True)
        layout.addWidget(hint_label)

        content = QWidget(self)
        grid = QGridLayout(content)
        grid.addWidget(QLabel("Screen"), 0, 0)
        grid.addWidget(QLabel("Minimum [A]"), 0, 1)
        grid.addWidget(QLabel("Maximum [A]"), 0, 2)

        self.spinboxes = {}
        for row, screen in enumerate(self.screens, start=1):
            low, high = ranges.get(screen, self.default_range)
            minimum_spinbox = self._make_spinbox(float(low))
            maximum_spinbox = self._make_spinbox(float(high))
            grid.addWidget(QLabel(screen), row, 0)
            grid.addWidget(minimum_spinbox, row, 1)
            grid.addWidget(maximum_spinbox, row, 2)
            self.spinboxes[screen] = (minimum_spinbox, maximum_spinbox)

        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame if hasattr(QFrame, "Shape") else QFrame.NoFrame)
        scroll.setWidget(content)
        layout.addWidget(scroll)

        buttons = QHBoxLayout()
        reset_button = QPushButton("Same range for all screens")
        reset_button.clicked.connect(self._apply_default_to_all)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        apply_button = QPushButton("Apply")
        apply_button.setDefault(True)
        apply_button.clicked.connect(self.accept)
        buttons.addWidget(reset_button)
        buttons.addStretch(1)
        buttons.addWidget(cancel_button)
        buttons.addWidget(apply_button)
        layout.addLayout(buttons)

    def _make_spinbox(self, value):
        spinbox = QDoubleSpinBox(self)
        spinbox.setDecimals(4)
        spinbox.setRange(-10000.0, 10000.0)
        spinbox.setSuffix(" A")
        spinbox.setValue(value)
        return spinbox

    def _apply_default_to_all(self):
        for minimum_spinbox, maximum_spinbox in self.spinboxes.values():
            minimum_spinbox.setValue(self.default_range[0])
            maximum_spinbox.setValue(self.default_range[1])

    def get_values(self):
        ranges = {}
        for screen, (minimum_spinbox, maximum_spinbox) in self.spinboxes.items():
            low = float(minimum_spinbox.value())
            high = float(maximum_spinbox.value())
            ranges[screen] = (min(low, high), max(low, high))
        return ranges
