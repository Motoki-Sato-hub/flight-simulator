try:
    from PyQt6 import uic
    from PyQt6.QtCore import Qt,QProcess,QProcessEnvironment
    from PyQt6.QtWidgets import (QApplication, QRadioButton,QSizePolicy, QMainWindow, QFileDialog, QListWidget, QListWidgetItem,QMessageBox,QProgressDialog, QVBoxLayout, QPushButton, QDialog, QLabel,QStyledItemDelegate)
    from PyQt6.QtGui import QPainter
    pyqt_version = 6
except ImportError:
    from PyQt5 import uic
    from PyQt5.QtCore import Qt,QProcess,QProcessEnvironment
    from PyQt5.QtWidgets import (QApplication, QRadioButton,QSizePolicy, QMainWindow, QFileDialog, QListWidget, QListWidgetItem,QMessageBox,QProgressDialog, QVBoxLayout, QPushButton, QDialog, QLabel,QStyledItemDelegate)
    from PyQt5.QtGui import QPainter
    pyqt_version = 5
from datetime import datetime
import sys, time,os
from pathlib import Path

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.cwd = os.getcwd()
        self._cancel = False
        ui_path = os.path.join(os.path.dirname(__file__), "UI files/MainPanel.ui")
        uic.loadUi(ui_path, self)
        self.application_choice=None
        self.compute_matrix_button.clicked.connect(self.handle_compute_response_matrix_click)
        self.compute_matrix_qm_button.clicked.connect(self.handle_compute_response_matrix_qm_click)
        self.bba_interface_button.clicked.connect(self.handle_bba_click)
        self.bba_qm_interface_button.clicked.connect(self.handle_bba_qm_click)
        self.emittance_interface_button.clicked.connect(self.handle_emittance_click)
        self.knobs_interface_button.clicked.connect(self.handle_knobs_click)
        self.sysid_interface_button.clicked.connect(self.handle_sysid_click)
        self.sysid_qm_interface_button.clicked.connect(self.handle_sysid_qm_click)
        # self.dispersion_interface_button.clicked.connect(self.handle_dispersion_click)
        # self.linac_interface_button.clicked.connect(self.handle_linac_click)
        self._procs=[]
        self.setWindowTitle("Choose the application")

    def handling(self, app_name,cwd=None, args=None):
        path = os.path.join(os.path.dirname(__file__), app_name)

        proc = QProcess(self)
        proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        proc.setWorkingDirectory(cwd or self.cwd)

        env = QProcessEnvironment.systemEnvironment()
        proc.setProcessEnvironment(env)

        argv = [path] + (args or [])
        proc.start(sys.executable, argv)

        proc.readyReadStandardOutput.connect(lambda p=proc: print(bytes(p.readAllStandardOutput()).decode(errors="ignore"), end=""))
        self._procs.append(proc)

    def handle_sysid_click(self):
        self.handling("SysID_GUI.py")

    def handle_sysid_qm_click(self):
        self.handling("SysID_QM_GUI.py")

    def handle_bba_click(self):
        self.handling("BBA_GUI.py")

    def handle_bba_qm_click(self):
        self.handling("BBA_QM_GUI.py")

    def handle_emittance_click(self):
        self.handling("Emittance_Measurement_GUI.py")

    def handle_knobs_click(self):
        self.handling("Knobs/GUIs/IPBSM_Opt_GUI.py")

    def handle_compute_response_matrix_click(self):
        self.handling("ComputeResponseMatrix_GUI.py")

    def handle_compute_response_matrix_qm_click(self):
        self.handling("ComputeResponseMatrix_QM_GUI.py")

    # def handle_linac_click(self):
    #     self.handling("Knobs/Linac_Opt.py")
    #
    # def handle_dispersion_click(self):
    #     self.handling("Tantative_combined_GUI.py", args=["--only-tab", "Dispersion"])

if __name__ == "__main__":
    app = QApplication([])
    app.setQuitOnLastWindowClosed(False)

    time_str = datetime.now().strftime("%Y%m%d_%H%M%S")

    w = MainWindow()
    w.show()
    sys.exit(app.exec())
