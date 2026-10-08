from Backend.State import State
from Backend.Response import Response
from Backend.ResponseMatrix_DFS_WFS_rnd import ResponseMatrix_DFS_WFS
import numpy as np
import glob,sys,os,argparse,matplotlib
try:
    from PyQt6 import uic
    from PyQt6.QtWidgets import QApplication, QMainWindow, QFileDialog, QListWidget, QVBoxLayout, QSizePolicy, \
        QMessageBox, QComboBox, QLabel
    from PyQt6.QtCore import Qt, QThread, QObject, pyqtSignal, pyqtSlot, QTimer
    from PyQt6.QtTest import QTest
    pyqt_version = 6
except ImportError:
    from PyQt5 import uic
    from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog, QListWidget, QVBoxLayout, QSizePolicy, \
        QMessageBox, QComboBox, QLabel
    from PyQt5.QtCore import Qt, QThread, QObject, pyqtSignal, pyqtSlot, QTimer
    from PyQt5.QtTest import QTest
    pyqt_version = 5
from Backend.SaveOrLoad import SaveOrLoad
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class MatplotlibWidget(FigureCanvas):
    def __init__(self, parent=None, title='', orbit=None):
        fig = Figure(tight_layout=True)
        super().__init__(fig)
        self.setParent(parent)
        self.axes = fig.add_subplot(111)

class MainWindow(QMainWindow, SaveOrLoad):
    def __init__(self,data_dir_1=None,data_dir_2=None,comp_difference=False,auto_click_compute=False):
        super().__init__()
        self._rnd_response = ResponseMatrix_DFS_WFS()
        uic.loadUi("UI files/ComputeResponseMatrix_GUI.ui", self)
        self.cwd = os.getcwd()
        self.R=None
        self.correctors=[]
        self.bpms=[]
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
        self.load_correctors_button.clicked.connect(self.__load_correctors_button_clicked)
        self.load_bpms_button.clicked.connect(self.__load_bpms_button_clicked)

        layout = self.plot_widget.layout()
        if layout is None:
            layout = QVBoxLayout(self.plot_widget)
        self.plot = MatplotlibWidget(self.plot_widget)
        self.plot.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.plot)
        self._add_sample_kind_control()

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

    def _add_sample_kind_control(self):
        self.sample_kind_label = QLabel("State:", self)
        self.sample_kind_combo = QComboBox(self)
        self.sample_kind_combo.addItem("nominal", "nominal")
        self.sample_kind_combo.addItem("energy", "energy")
        self.sample_kind_combo.currentIndexChanged.connect(self._reload_lists_for_selected_kind)
        self.dirA_layout.addWidget(self.sample_kind_label)
        self.dirA_layout.addWidget(self.sample_kind_combo)

    def _reload_lists_for_selected_kind(self):
        folder = self._expand_path(self.data_directory_1.text())
        if folder:
            self._load_lists_from_directory(folder)

    def _expand_path(self,path):
        expanded_path=(path or "").strip()
        expanded_path=os.path.expandvars(os.path.expanduser(expanded_path))
        return os.path.abspath(os.path.normpath(expanded_path))

    def _compute_difference_clicked(self,checked):
        checked=bool(checked)
        self.label_dir_2.setEnabled(checked)
        self.choose_directory_2.setEnabled(checked)
        self.data_directory_2.setEnabled(checked)

    def _pick_directory_with_data(self,line_edit):
        base=self._expand_path(line_edit.text() or self.cwd) or self.cwd
        folder=QFileDialog.getExistingDirectory(self,"Select data directory",base)
        if not folder:
            return
        line_edit.setText(folder)
        if line_edit is self.data_directory_1:
            self._load_lists_from_directory(folder)

    def _load_lists_from_directory(self,folder):
        folder=self._expand_path(folder)
        info = self._rnd_response._find_useful_files(
            folder, sample_kind=self.sample_kind_combo.currentData())
        if not info["ok"]:
            return
        S=State(filename=info["pairs"][0])
        self.sequence=S.sequence
        self.correctors=list(S.get_correctors()["names"])
        self.bpms=list(S.get_bpms()["names"])

        self.correctors_list.clear()
        self.correctors_list.addItems([str(c) for c in self.correctors])

        self.bpms_list.clear()
        self.bpms_list.addItems([str(b) for b in self.bpms])

    def _compute_response_of_one_data_directory(self, directory):
        directory = self._expand_path(directory)
        correctors = [item.text() for item in self.correctors_list.selectedItems()] or list(self.correctors)
        bpms = [item.text() for item in self.bpms_list.selectedItems()] or list(self.bpms)

        Rxx, Ryy, Rxy, Ryx, Bx, By, hcorrs, vcorrs, fitted_bpms = \
            self._rnd_response._compute_response_matrix_from_directory(
                directory=directory,
                correctors=correctors,
                bpms=bpms,
                triangular=bool(self.triangular_checkbox.isChecked()),
                sample_kind=self.sample_kind_combo.currentData(),
            )

        R = Response()
        R.bpms = fitted_bpms
        R.hcorrs = hcorrs
        R.vcorrs = vcorrs
        R.Rxx = Rxx
        R.Rxy = Rxy
        R.Ryx = Ryx
        R.Ryy = Ryy
        R.Bx = Bx
        R.By = By
        return R

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
        dir_name = self.data_directory_1.text() + '/correctors.txt'
        filename, _ = QFileDialog.getOpenFileName(None, "Open File", dir_name, "Text Files (*.txt)")
        if filename:
            with open(filename, 'r') as f:
                selected_correctors = [line.strip() for line in f]

        self.correctors_list.clearSelection()
        for item in selected_correctors:
            items = self.correctors_list.findItems(item, Qt.MatchFlag.MatchExactly)
            for item in items:
                item.setSelected(True)

    def __clear_correctors_button_clicked(self):
        self.correctors_list.clearSelection()

    def __save_bpms_button_clicked(self):
        self._save_bpms()

    def __load_bpms_button_clicked(self):
        dir_name = self.data_directory_1.text() + '/bpms.txt'
        filename, _ = QFileDialog.getOpenFileName(None, "Open File", dir_name, "Text Files (*.txt)")
        if filename:
            with open(filename, 'r') as f:
                selected_bpms = [line.strip() for line in f]

        self.bpms_list.clearSelection()
        for item in selected_bpms:
            items = self.bpms_list.findItems(item, Qt.MatchFlag.MatchExactly)
            for item in items:
                item.setSelected(True)

    def __clear_bpms_button_clicked(self):
        self.bpms_list.clearSelection()

    def _plot_response_matrix(self,R):
        # Clear existing figure
        self.plot.figure.clf()

        # === 3D surface plots ===
        # If you want to show both sets (2D + 3D), you need two canvases.
        # But if you want to reuse the same canvas, just do this after the first draw:

        # Clear and re-plot the 3D figure
        self.plot.figure.clf()
        fig1 = self.plot.figure
        # fig1.set_size_inches(12, 10)

        ax1 = fig1.add_subplot(2, 2, 1, projection='3d')
        ax3 = fig1.add_subplot(2, 2, 3, projection='3d')
        ax2 = fig1.add_subplot(2, 2, 2, projection='3d')
        ax4 = fig1.add_subplot(2, 2, 4, projection='3d')

        x = np.array(range(len(R.hcorrs)))
        y = np.array(range(len(R.bpms)))
        X, Y = np.meshgrid(x, y)

        ax1.plot_surface(X, Y, R.Rxx, cmap='viridis')
        ax1.set_title('$R_{xx}$')
        ax1.set_xlabel('Corrector [#]')
        ax1.set_ylabel('BPM [#]')
        ax3.plot_surface(X, Y, R.Ryx, cmap='viridis')
        ax3.set_title('$R_{yx}$')
        ax3.set_xlabel('Corrector [#]')
        ax3.set_ylabel('BPM [#]')

        x = np.array(range(len(R.vcorrs)))
        X, Y = np.meshgrid(x, y)

        ax2.plot_surface(X, Y, R.Rxy, cmap='viridis')
        ax2.set_title('$R_{xy}$')
        ax2.set_xlabel('Corrector [#]')
        ax2.set_ylabel('BPM [#]')

        ax4.plot_surface(X, Y, R.Ryy, cmap='viridis')
        ax4.set_title('$R_{yy}$')
        ax4.set_xlabel('Corrector [#]')
        ax4.set_ylabel('BPM [#]')

        fig1.tight_layout()

        self.plot.draw()
        self.plot.flush_events()
        self.plot.repaint()

    def __compute_button_clicked(self):
        try:
            directory_1=self._expand_path(self.data_directory_1.text())
            if not directory_1:
                QMessageBox.warning(self, "Error", "No data directory specified")
                return
            R1=self._compute_response_of_one_data_directory(directory_1)

            if self.diff_checkbox.isChecked():
                directory_2 = self._expand_path(self.data_directory_2.text())
                if not directory_2:
                    QMessageBox.warning(self, "Error", "No second data directory specified")
                    return
                R2 = self._compute_response_of_one_data_directory(directory_2)
                self.R=self._substract_matrices(R1=R1,R2=R2)
            else:
                self.R=R1
            self._plot_response_matrix(R=self.R)

        except Exception as e:
            from traceback import print_exception
            QMessageBox.warning(self, "Error", str(e))
            print_exception(e)

    def __save_as_button_clicked(self):
        dir_name = self.cwd + '/response2.pkl'
        os.chdir (self.cwd)
        filename, _ = QFileDialog.getSaveFileName(None, "Save Response Matrix", dir_name, "Pickle Files (*.pkl)")
        if filename:
            self.R.save(filename)

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
    window = MainWindow(data_dir_1=args.dir1,
        data_dir_2=args.dir2,
        comp_difference=args.diff,
        auto_click_compute=args.compute,)
    window.show()
    sys.exit(app.exec())
