try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QVBoxLayout, QHBoxLayout, QDialog, QLabel, QWidget, QComboBox, QSizePolicy
except ImportError:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QDialog, QLabel, QWidget, QComboBox, QSizePolicy
import numpy as np
import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

class DisplayScreenImages(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Screen Images")
        self.setWindowFlag(Qt.WindowType.WindowContextHelpButtonHint, False)
        self.setMinimumSize(800, 500)
        self.resize(1000, 650)
        self.setSizeGripEnabled(True)
        self.session = None
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
        header_layout.addWidget(QLabel("Screen:", self))
        self.screen_combobox = QComboBox(self)
        self.screen_combobox.currentIndexChanged.connect(self._update_screen_image)
        header_layout.addWidget(self.screen_combobox)
        header_layout.addWidget(QLabel("Step:", self))
        self.step_combobox = QComboBox(self)
        self.step_combobox.currentIndexChanged.connect(self._update_screen_image)
        header_layout.addWidget(self.step_combobox)
        header_layout.addStretch(1)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)
        layout.addWidget(header, 0)
        layout.addWidget(plot_widget, 1)

    def _plot_screen_image(self, session):
        self.session = session
        screens = list(session.get("screens", []))
        steps = session.get("steps", [])
        self.screen_combobox.blockSignals(True)
        self.screen_combobox.clear()
        self.screen_combobox.addItems([str(s) for s in screens])
        self.screen_combobox.blockSignals(False)

        self.step_combobox.blockSignals(True)
        self.step_combobox.clear()
        if steps == 0:
            self.step_combobox.addItems(["0"])
        else:
            self.step_combobox.addItems([str(s) for s in range(1, steps+1)])
        self.step_combobox.blockSignals(False)
        self._update_screen_image()

    def _update_screen_image(self, index=None):
        if self.session is None:
            return
        screen_index = self.screen_combobox.currentIndex()
        step_index = self.step_combobox.currentIndex()
        images = self.session.get("images", [])
        shot_images = images[step_index][screen_index]
        shot_images = [np.asarray(image, dtype=float) for image in shot_images if image is not None]
        if not shot_images:
            fig = self.canvas.figure
            fig.clear()
            ax = fig.add_subplot(111)
            screen_name = self.screen_combobox.currentText()
            step_name = self.step_combobox.currentText()
            ax.text(0.5, 0.5, f"No image saved for {screen_name} at step {step_name}.", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            self.canvas.draw_idle()
            return
        image = np.nanmean(np.stack(shot_images, axis=0), axis=0)
        hedges = self._get_edges("hedges", step_index, screen_index)
        vedges = self._get_edges("vedges", step_index, screen_index)

        if image.shape == (hedges.size - 1, vedges.size - 1):
            image = image.T
        intensity_x = np.nansum(image, axis=0)
        intensity_y = np.nansum(image, axis=1)

        fig = self.canvas.figure
        fig.clear()
        gridspec=fig.add_gridspec(2,2,width_ratios=(4,1), height_ratios=(1,4), hspace=0.05, wspace=0.05)
        ax = fig.add_subplot(gridspec[1,0])
        ax_x = fig.add_subplot(gridspec[0,0],sharex=ax)
        ax_y = fig.add_subplot(gridspec[1,1],sharey=ax)
        x_coordinates = 0.5 * (hedges[:-1] + hedges[1:])
        y_coordinates = 0.5 * (vedges[:-1] + vedges[1:])
        ax.imshow(image, origin="lower", extent=[hedges[0], hedges[-1], vedges[0], vedges[-1]], aspect="auto", cmap="jet")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax_x.plot(x_coordinates, intensity_x)
        ax_y.plot(intensity_y, y_coordinates)
        ax_x.tick_params(labelbottom=False)
        ax_y.tick_params(labelleft=False)
        self.canvas.draw_idle()

    def _get_edges(self, key, step_index, screen_index):
        values = self.session.get(key, [])
        if step_index >= len(values) or screen_index >= len(values[step_index]):
            return None
        for shot_edges in values[step_index][screen_index]:
            if shot_edges is not None:
                edges = np.asarray(shot_edges, dtype=float)
                if edges.ndim == 1 and np.all(np.isfinite(edges)):
                    return edges
        return None







