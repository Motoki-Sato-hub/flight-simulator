try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QVBoxLayout, QHBoxLayout, QDialog, QLabel, QWidget, QComboBox, QSizePolicy
except ImportError:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QDialog, QLabel, QWidget, QComboBox, QSizePolicy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from Backend.draw_beamline_script import drawBeamline

class ShowBeamline(QDialog):
    def __init__(self, interface, parent=None, quad_selected=None, screens = None):
        super().__init__(parent)
        self.interface = interface
        quadrupoles = list(getattr(self.interface, "quadrupoles", []))
        try:
            self.quad_selected = quad_selected[0]
        except:
            self.quad_selected = quadrupoles[0]
        self.screens = [str(screen) for screen in screens]
        self.lattice = getattr(self.interface, "lattice", None)
        self.setWindowTitle("Beamline View")
        self.setWindowFlag(Qt.WindowType.WindowContextHelpButtonHint, False)
        self.setMinimumSize(800, 500)
        self.resize(1000, 650)
        self.setSizeGripEnabled(True)
        self.figure = Figure(figsize=(10, 6), constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        plot_widget = QWidget(self)
        plot_widget.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        plot_layout = QVBoxLayout(plot_widget)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        plot_layout.setSpacing(0)
        plot_layout.addWidget(self.canvas, 1)
        header = QWidget(self)
        header_layout = QHBoxLayout(header)
        header_layout.setContentsMargins(0, 0, 0, 0)
        header_layout.addWidget(QLabel("Start element:", self))
        self.start_element_combobox = QComboBox(self)
        header_layout.addWidget(self.start_element_combobox)
        header_layout.addWidget(QLabel("Last element:", self))
        self.last_element_combobox = QComboBox(self)
        header_layout.addWidget(self.last_element_combobox)
        header_layout.addStretch(1)
        self.last_screen = self.screens[-1]


        start_quad_element = self.interface._give_elements_to_show_beamline(self.quad_selected)
        try:
            self.start_element_name = start_quad_element.get_name()
        except:
            self.start_element_name = start_quad_element
        start_positions = self.interface._get_elements_positions_show_beamline(names=self.start_element_name)["S"]
        self.first_element_position = float(start_positions[0])
        self.last_element_position = self.interface._get_elements_positions_show_beamline(names = self.last_screen)['S'][0]

        self.lattice_names = self.interface._get_elements_positions_show_beamline()["names"]
        self.start_element_combobox.addItems([str(e) for e in self.lattice_names])
        self.last_element_combobox.addItems([str(e) for e in self.lattice_names])

        self.start_element_combobox.setCurrentText(self.start_element_name)
        self.last_element_combobox.setCurrentText(self.last_screen)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)
        layout.addWidget(header, 0)
        layout.addWidget(plot_widget, 1)
        self.start_element_combobox.currentTextChanged.connect(self._display_beamline_view)
        self.last_element_combobox.currentTextChanged.connect(self._display_beamline_view)


    def _display_beamline_view(self):

        start_element = self.start_element_combobox.currentText()
        last_element = self.last_element_combobox.currentText()
        self.figure.clear()
        ax = self.figure.subplots(1, 1)
        drawer = drawBeamline(ax)
        start_s = float(self.interface._get_elements_positions_show_beamline(names=start_element)["S"][0])
        end_s = float(self.interface._get_elements_positions_show_beamline(names=last_element)["S"][0])
        ax.set_xlim(start_s, end_s)
        self.lattice.accept(drawer)
        ax_bottom = ax.secondary_xaxis("bottom")
        ax_bottom.set_xlabel("S [m]")
        ax_bottom.set_xticks(np.linspace(start_s, end_s, 6))
        ax_bottom.tick_params(axis="x", labelsize=9)
        self.canvas.draw()

