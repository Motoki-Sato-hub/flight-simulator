import numpy as np
try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QListWidget,
                                 QListWidgetItem, QPushButton)
except ImportError:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QListWidget,
                                 QListWidgetItem, QPushButton)

class ScanPointSelection(QDialog):
    POINT_ROLE = int(Qt.ItemDataRole.UserRole) + 1

    def __init__(self, points, excluded, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Scan points used in the fit")
        self.resize(560, 460)
        excluded = set(excluded or ())

        layout = QVBoxLayout(self)
        self.point_list = QListWidget(self)
        for point in points:
            item = QListWidgetItem(self._label(point))
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            key = (int(point["step"]), int(point["screen_index"]))
            item.setCheckState(Qt.CheckState.Unchecked if key in excluded else Qt.CheckState.Checked)
            item.setData(self.POINT_ROLE, key)
            self.point_list.addItem(item)
        layout.addWidget(self.point_list)

        buttons = QHBoxLayout()
        use_all_button = QPushButton("Use all points")
        use_all_button.clicked.connect(self._use_all)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        apply_button = QPushButton("Apply")
        apply_button.setDefault(True)
        apply_button.clicked.connect(self.accept)
        buttons.addWidget(use_all_button)
        buttons.addStretch(1)
        buttons.addWidget(cancel_button)
        buttons.addWidget(apply_button)
        layout.addLayout(buttons)

    @staticmethod
    def _label(point):
        def fmt(value):
            value = float(value)
            return f"{value:.3f}" if np.isfinite(value) else "-"
        quad_text = f"{float(point['quad_value']):.4g} {point['quad_unit']}"
        return (f"step {int(point['step']):>3}   {quad_text:>14}   {point['screen']:<14}"
                f"   sigx = {fmt(point['sigx']):>9}   sigy = {fmt(point['sigy']):>9}")

    def _use_all(self):
        for row in range(self.point_list.count()):
            self.point_list.item(row).setCheckState(Qt.CheckState.Checked)

    def get_excluded(self):
        excluded = set()
        for row in range(self.point_list.count()):
            item = self.point_list.item(row)
            if item.checkState() != Qt.CheckState.Checked:
                excluded.add(tuple(item.data(self.POINT_ROLE)))
        return excluded
