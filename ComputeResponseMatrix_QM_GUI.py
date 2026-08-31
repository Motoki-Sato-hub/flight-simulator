"""Response-matrix tool dedicated to quadrupole-mover SysID data."""

import argparse
import os
import re
import sys

from Backend.State import State
from ComputeResponseMatrix_GUI import MainWindow as CorrectorResponseMatrixWindow


class MainWindow(CorrectorResponseMatrixWindow):
    """Compute a response matrix from ``DATA_<quad>_<axis>_*.pkl`` files."""

    def __init__(self, data_dir_1=None, data_dir_2=None, comp_difference=False, auto_click_compute=False):
        super().__init__(data_dir_1, data_dir_2, comp_difference, auto_click_compute)
        self.setWindowTitle("Compute QM Response Matrix Tool")

    def _current_actuator_mode(self):
        return "quadrupole_movers"

    def _actuator_label(self):
        return "Quadrupole mover"

    def _load_lists_from_directory(self, folder):
        folder = self._expand_path(folder)
        datafiles = sorted(self._datafiles(folder))
        if not datafiles:
            return
        state = State(filename=datafiles[0])
        self.sequence = state.get_sequence()
        self.bpms = list(state.get_bpms()["names"])
        self.correctors = self._quadrupole_movers_from_datafiles(datafiles)
        self.correctorsGroup.setTitle("Quadrupole movers")

        self.correctors_list.clear()
        self.correctors_list.addItems([str(name) for name in self.correctors])
        self.bpms_list.clear()
        self.bpms_list.addItems([str(name) for name in self.bpms])

    @staticmethod
    def _datafiles(folder):
        import glob

        return glob.glob(os.path.join(folder, "DATA*.pkl"))

    def _quadrupole_movers_from_datafiles(self, datafiles):
        pattern = re.compile(r"DATA_(.+)_(p|m)(\d+)\.pkl$")
        names = []
        seen = set()
        sequence_index = {str(name): i for i, name in enumerate(self.sequence)}
        for datafile in datafiles:
            match = pattern.search(os.path.basename(datafile))
            if not match:
                continue
            tag = match.group(1)
            if not tag.endswith(("_x", "_y")):
                continue
            name = tag[:-2]
            if name not in seen:
                seen.add(name)
                names.append(name)
        return sorted(names, key=lambda name: sequence_index.get(name, 10**9))


def main():
    parser = argparse.ArgumentParser(description="Compute quadrupole-mover response matrix")
    parser.add_argument("--dir1", default=None, help="First data directory")
    parser.add_argument("--dir2", default=None, help="Second data directory")
    parser.add_argument("--diff", action="store_true", help="Difference between responses")
    parser.add_argument("--compute", action="store_true", help="Auto-click Compute button")
    args = parser.parse_args()

    try:
        from PyQt6.QtWidgets import QApplication
    except ImportError:
        from PyQt5.QtWidgets import QApplication

    app = QApplication(sys.argv)
    window = MainWindow(args.dir1, args.dir2, args.diff, args.compute)
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
