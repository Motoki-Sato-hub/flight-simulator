from Backend.State import State
from Backend.Response import Response
import pprint
from Backend.ResponseMatrix_DFS_WFS import ResponseMatrix_DFS_WFS
try:
    from PyQt6 import uic
    from PyQt6.QtWidgets import (
        QApplication, QMainWindow, QVBoxLayout,
        QLineEdit, QListWidget, QPushButton,
        QCheckBox, QFileDialog, QSizePolicy,QMessageBox,
        QWidget,
        )
    from PyQt6.QtCore import Qt,QTimer
    from PyQt6.QtGui import QPixmap
    pyqt_version = 6

except ImportError:
    from PyQt5 import uic
    from PyQt5.QtWidgets import (
        QApplication, QMainWindow, QVBoxLayout,
        QLineEdit, QListWidget, QPushButton,
        QCheckBox, QFileDialog, QSizePolicy,QMessageBox,
        QWidget,
        )
    from PyQt5.QtCore import Qt,QTimer
    from PyQt5.QtGui import QPixmap
    pyqt_version = 5
import numpy as np
import glob,sys,os,argparse,matplotlib
from Backend.SaveOrLoad import SaveOrLoad
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

class PlotPopup(QMainWindow):
    def __init__(self, title="Response Matrix", parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(1100, 850)
        central = QWidget(self)
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        layout.setContentsMargins(6, 6, 6, 6)
        self.plot = MatplotlibWidget(central)
        self.plot.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.plot)

class MatplotlibWidget(FigureCanvas):
    def __init__(self, parent=None, title='', orbit=None):
        fig = Figure(tight_layout=True)
        super().__init__(fig)
        self.setParent(parent)
        self.axes = fig.add_subplot(111)

class MainWindow(QMainWindow, SaveOrLoad, ResponseMatrix_DFS_WFS):
    """Response-matrix tool for steering-corrector SysID data."""

    def __init__(self,data_dir_1=None,data_dir_2=None,comp_difference=False,auto_click_compute=False):
        super().__init__()
        uic.loadUi("UI files/ComputeResponseMatrix_GUI.ui", self)
        self._load_logo()
        self.cwd = os.getcwd()
        self.R=None
        self.response_matrix_popup = None
        self.data_dir_1=data_dir_1
        self.data_dir_2=data_dir_2
        self.comp_difference=comp_difference
        self.auto_click_compute=auto_click_compute
        self.setWindowTitle("Compute Response Matrix Tool")
        self.setGeometry(100, 100, 600, 700)
        self.data_directory_1.setText(self.cwd)
        self.choose_directory_1.clicked.connect(lambda: self._pick_directory_with_data(self.data_directory_1))
        self.choose_directory_2.clicked.connect(lambda: self._pick_directory_with_data(self.data_directory_2))
        self.compute_button.clicked.connect(self.__compute_button_clicked)
        self.save_as_button.clicked.connect(self.__save_as_button_clicked)
        self.diff_checkbox.toggled.connect(self._compute_difference_clicked)
        self._compute_difference_clicked(self.diff_checkbox.isChecked())
        self.label_actuator_type.setVisible(False)
        self.actuator_type_combo.setVisible(False)
        self.plot_singular_values_button.clicked.connect(self._plot_singular_values)
        self.load_correctors_button.clicked.connect(self.__load_correctors_button_clicked)
        self.load_bpms_button.clicked.connect(self.__load_bpms_button_clicked)
        # self.hcorrector_prefixes=("zh", "zx")
        # self.vcorrector_prefixes=("zv")

        layout = self.plot_widget.layout()
        if layout is None:
            layout = QVBoxLayout(self.plot_widget)
        self.plot = MatplotlibWidget(self.plot_widget)
        self.plot.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.plot)
        self.plot.mpl_connect("button_press_event", self._handle_plot_double_click)

        if self.data_dir_1:
            self.data_dir_1=self._expand_path(self.data_dir_1)
            self.data_directory_1.setText(self.data_dir_1)
            self._load_lists_from_directory(self.data_dir_1)
        if self.comp_difference==True:
            self.diff_checkbox.setChecked(True)
            self._compute_difference_clicked(checked=True)
        if self.data_dir_2:
            self.data_dir_2=self._expand_path(self.data_dir_2)
            self.data_directory_2.setText(self.data_dir_2)
            self._load_lists_from_directory(self.data_dir_2)
        if self.auto_click_compute:
            QTimer.singleShot(0, self.__compute_button_clicked)

    def _handle_plot_double_click(self, event):
        if event is None:
            return
        if getattr(event, "dblclick", False) and getattr(event, "button", None) == 1:
            self._show_response_matrix_popup()

    def _show_response_matrix_popup(self):
        if self.R is None:
            QMessageBox.information(self, "No response matrix", "Compute a response matrix first.")
            return
        if self.response_matrix_popup is None:
            self.response_matrix_popup = PlotPopup("Response Matrix", parent=self)
        self._draw_response_matrix(self.response_matrix_popup.plot, self.R)
        self.response_matrix_popup.show()
        self.response_matrix_popup.raise_()
        self.response_matrix_popup.activateWindow()

    def _draw_response_matrix(self, plot, R):
        plot.figure.clf()
        fig1 = plot.figure

        ax1 = fig1.add_subplot(2, 2, 1, projection='3d')
        ax3 = fig1.add_subplot(2, 2, 3, projection='3d')
        ax2 = fig1.add_subplot(2, 2, 2, projection='3d')
        ax4 = fig1.add_subplot(2, 2, 4, projection='3d')

        x = np.array(range(len(R.hcorrs)))
        y = np.array(range(len(R.bpms)))
        X, Y = np.meshgrid(x, y)

        ax1.plot_surface(X, Y, R.Rxx, cmap='viridis')
        ax1.set_title('$R_{xx}$')
        ax1.set_xlabel(f'{self._actuator_label()} [#]')
        ax1.set_ylabel('BPM [#]')

        ax3.plot_surface(X, Y, R.Ryx, cmap='viridis')
        ax3.set_title('$R_{yx}$')
        ax3.set_xlabel(f'{self._actuator_label()} [#]')
        ax3.set_ylabel('BPM [#]')

        x = np.array(range(len(R.vcorrs)))
        X, Y = np.meshgrid(x, y)

        ax2.plot_surface(X, Y, R.Rxy, cmap='viridis')
        ax2.set_title('$R_{xy}$')
        ax2.set_xlabel(f'{self._actuator_label()} [#]')
        ax2.set_ylabel('BPM [#]')

        ax4.plot_surface(X, Y, R.Ryy, cmap='viridis')
        ax4.set_title('$R_{yy}$')
        ax4.set_xlabel(f'{self._actuator_label()} [#]')
        ax4.set_ylabel('BPM [#]')

        fig1.tight_layout()
        plot.draw()
        plot.flush_events()
        plot.repaint()

    def _compute_response_of_one_data_directory(self,directory):
        directory=self._expand_path(directory)
        datafiles=sorted(glob.glob(os.path.join(directory,"DATA*.pkl")))
        if not datafiles:
            QMessageBox.warning(self, "Error", "No data files found")
            return
        S=State(filename=datafiles[0])
        self.sequence=S.get_sequence()
        correctors = [self.correctors_list.item(i).text() for i in range(self.correctors_list.count()) if self.correctors_list.item(i).isSelected()]
        bpms = [self.bpms_list.item(i).text() for i in range(self.bpms_list.count()) if self.bpms_list.item(i).isSelected()]

        if not correctors:
            for i in range(self.correctors_list.count()):
                self.correctors_list.item(i).setSelected(True)
            correctors=self.correctors

        if not bpms:
            for i in range(self.bpms_list.count()):
                self.bpms_list.item(i).setSelected(True)
            bpms=self.bpms

        Rxx, Ryy, Rxy, Ryx, Bx, By, hcorrs, vcorrs, bpms=self._compute_response_matrix_from_directory(directory=directory, correctors=correctors, bpms=bpms, triangular=bool(self.triangular_checkbox.isChecked()), actuator_mode = self._current_actuator_mode())
        R = Response()
        R.bpms = bpms
        R.hcorrs = hcorrs
        R.vcorrs = vcorrs
        R.Rxx = Rxx
        R.Rxy = Rxy
        R.Ryx = Ryx
        R.Ryy = Ryy
        print(R.Ryy)
        R.Bx = Bx
        R.By = By
        return R

    def _expand_path(self,path):
        expanded_path=(path or "").strip()
        expanded_path=os.path.expandvars(os.path.expanduser(expanded_path))
        return os.path.abspath(os.path.normpath(expanded_path))

    def _current_actuator_mode(self):
        return "correctors"

    def _actuator_label(self):
        return "Corrector"

    def _refresh_device_lists_for_current_mode(self):
        folder = self._expand_path(self.data_directory_1.text())
        if folder:
            self._load_lists_from_directory(folder)

    def _compute_difference_clicked(self,checked):
        checked=bool(checked)
        self.label_dir_2.setEnabled(checked)
        self.choose_directory_2.setEnabled(checked)
        self.data_directory_2.setEnabled(checked)

    def _plot_singular_values(self):
        if self.R is None:
            QMessageBox.warning(self, "Error", "Load directory data and compute response matrix first.")
            return

        def get_SV(R):
            R = R.copy()
            R[np.isnan(R)] = 0
            U, S, Vh = np.linalg.svd(R)
            return S
        plt.semilogy(get_SV(self.R.Rxx), label='Rxx')
        plt.semilogy(get_SV(self.R.Rxy), label='Rxy', linestyle='dashed')
        plt.semilogy(get_SV(self.R.Ryx), label='Ryx', linestyle='dashed')
        plt.semilogy(get_SV(self.R.Ryy), label='Ryy')
        plt.xlabel('Singular Value')
        plt.ylabel('Value')
        plt.legend()
        plt.show()

    def _pick_directory_with_data(self,line_edit):
        # base=self._expand_path(line_edit.text() or self.cwd) or self.cwd
        default_dir = f"~/CERN-Flight_Simulator-Data/"
        default_dir = os.path.expanduser(os.path.expandvars(default_dir))
        folder=QFileDialog.getExistingDirectory(self,"Select data directory",default_dir)
        if not folder:
            return
        line_edit.setText(folder)
        if line_edit is self.data_directory_1:
            self._load_lists_from_directory(folder)

    def _load_lists_from_directory(self,folder):
        folder=self._expand_path(folder)
        datafiles=sorted(glob.glob(os.path.join(folder,"DATA*.pkl")))
        if not datafiles:
            return
        S=State(filename=datafiles[0])
        self.sequence=S.get_sequence()
        self.bpms=list(S.get_bpms()["names"])

        self.correctors = list(S.get_correctors()["names"])
        self.correctorsGroup.setTitle("Correctors")

        self.correctors_list.clear()
        self.correctors_list.addItems([str(c) for c in self.correctors])

        self.bpms_list.clear()
        self.bpms_list.addItems([str(b) for b in self.bpms])

    def _substract_matrices(self,R1,R2):
        if R1.Rxx.shape != R2.Rxx.shape or R1.Ryy.shape != R2.Ryy.shape:
            raise RuntimeError("Response matrices must have same shape")
        R = Response()
        R.bpms = R1.bpms
        R.hcorrs = R1.hcorrs
        R.vcorrs = R1.vcorrs
        R.Rxx = R1.Rxx-R2.Rxx
        R.Rxy = R1.Rxy-R2.Rxy
        R.Ryx = R1.Ryx-R2.Ryx
        R.Ryy = R1.Ryy-R2.Ryy
        R.Bx = R1.Bx-R2.Bx
        R.By = R1.By-R2.By
        return R

    def __save_correctors_button_clicked(self):
        self._save_correctors()

    def __load_correctors_button_clicked(self):
        self._load_correctors()

    def __clear_correctors_button_clicked(self):
        self.correctors_list.clearSelection()

    def __save_bpms_button_clicked(self):
        self._save_bpms()

    def __load_bpms_button_clicked(self):
        self._load_bpms()

    def __clear_bpms_button_clicked(self):
        self.bpms_list.clearSelection()

    def _plot_response_matrix(self,R):
        self._draw_response_matrix(self.plot, R)
        if self.response_matrix_popup is not None and self.response_matrix_popup.isVisible():
            self._draw_response_matrix(self.response_matrix_popup.plot, R)

    def __compute_button_clicked(self):
        try:
            directory_1=self._expand_path(self.data_directory_1.text())
            if not directory_1:
                QMessageBox.warning(self, "Error", "No data directory specified")
                return
            R1=self._compute_response_of_one_data_directory(directory_1)
            if R1 is None:
                QMessageBox.warning(self, "Error", "No valid DATA pairs found in the first directory")
                return
            if self.diff_checkbox.isChecked():
                directory_2 = self._expand_path(self.data_directory_2.text())
                if not directory_2:
                    QMessageBox.warning(self, "Error", "No second data directory specified")
                    return
                R2 = self._compute_response_of_one_data_directory(directory_2)
                if R2 is None:
                    QMessageBox.warning(self, "Error", "No valid DATA pairs found in the second directory")
                    return
                self.R=self._substract_matrices(R1=R1,R2=R2)
            else:
                self.R=R1
            self._plot_response_matrix(R=self.R)

        except Exception as e:
            from traceback import print_exception
            QMessageBox.warning(self, "Error", str(e))
            print_exception(e)

    def __save_as_button_clicked(self):
        default_dir = f"~/CERN-Flight_Simulator-Data/"
        default_dir = os.path.expanduser(os.path.expandvars(default_dir))
        dir_name = default_dir + '/response.pkl'
        os.chdir (self.cwd)
        filename, _ = QFileDialog.getSaveFileName(None, "Save Response Matrix", dir_name, "Pickle Files (*.pkl)")
        if filename:
            self.R.save(filename)

    def _load_logo(self):
        if not hasattr(self, "logo_label"):
            return

        self.logo_label.setText("")
        self.logo_label.setScaledContents(False)

        transform_mode = (
            Qt.TransformationMode.SmoothTransformation
            if pyqt_version == 6
            else Qt.SmoothTransformation
        )

        base_dir = os.path.dirname(os.path.abspath(__file__))
        logo_path = os.path.join(base_dir, "UI files", "Assets", "CERN_logo.png")

        if not os.path.isfile(logo_path):
            return

        pixmap = QPixmap(logo_path)
        if pixmap.isNull():
            return

        scaled = pixmap.scaledToHeight(80, transform_mode)
        self.logo_label.setPixmap(scaled)
        self.logo_label.setToolTip(logo_path)

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description='Compute Response Matrix GUI')
    parser.add_argument("--dir1",default=None,help="First data directory")
    parser.add_argument("--dir2",default=None,help="Second data directory")
    parser.add_argument("--diff",action="store_true",help="Difference between responses")
    parser.add_argument("--compute",action="store_true",help="Auto-click Compute button")

    args=parser.parse_args()

    dir1=args.dir1
    dir2=args.dir2

    app = QApplication(sys.argv)
    window = MainWindow(data_dir_1=args.dir1, data_dir_2=args.dir2, comp_difference=args.diff, auto_click_compute=args.compute)
    window.show()
    sys.exit(app.exec())
