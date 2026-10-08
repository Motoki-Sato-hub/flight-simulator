import RF_Track as rft
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class drawBeamline(rft.UserVisitor):
    def __init__(self, ax, latticeframe=False, volumeframe=False):
        super().__init__()
        self.ax = ax
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_xlabel('')
        self.ax.set_ylabel('')
        self.ax.axis('off')
        self.volumeframe = volumeframe
        self.latticeframe = latticeframe
        self.color_map = {
            "dipole": "#6699FF",  # soft blue
            "quadrupole": "#99CC99",  # light green
            "sextupole": "#FFCC99",  # pale orange
            "multipole": "#FFCC99",  # pale orange, same as sextupole
            "solenoid": "#99CCCC",  # light teal
            "rf_cavity": "#CC99FF",  # lavender
            "diagnostic": "#E66567",  # soft red
        }

    def visit(self, e):
        print(type(e), e.get_name())
        s0 = e.get_S('entrance')
        s1 = e.get_S('exit')
        cm = self.color_map

        if type(e) == rft.Drift:
            if s1 > s0:
                self.ax.plot([s0, s1], [0, 0], color='black', linewidth=1)
            else:
                self.ax.plot([s1, s0], [0, 0], color='white', linewidth=1)

        elif type(e) == rft.Quadrupole:
            a, b = (0, 0.75) if e.get_gradient() > 0 else (-0.75, 0)
            rect = patches.Rectangle((s0, a), s1 - s0, b - a, facecolor=cm['quadrupole'])
            self.ax.add_patch(rect)

        elif type(e) == rft.Solenoid:
            a, b = (-0.75, 0.75)
            rect = patches.Rectangle((s0, a), s1 - s0, b - a, facecolor=cm['solenoid'])
            self.ax.add_patch(rect)

        elif type(e) == rft.Multipole:
            a, b = (-0.25, 0.25)
            rect = patches.Rectangle((s0, a), s1 - s0, b - a, facecolor=cm['multipole'])
            self.ax.add_patch(rect)
            #pass

        elif type(e) == rft.Sextupole:
            a, b = (-0.75, 0.75)
            rect = patches.Rectangle((s0, a), s1 - s0, b - a, facecolor=cm['sextupole'])
            self.ax.add_patch(rect)

        elif type(e) in [rft.RF_FieldMap, rft.RF_FieldMap_1d, rft.RF_FieldMap_1d_CINT,
                         rft.RF_FieldMap_2d, rft.RF_FieldMap_2d_CINT]:
            a, b = (-0.75, 0.75)
            rect = patches.Rectangle((s0, a), s1 - s0, b - a, facecolor=cm['rf_cavity'])
            self.ax.add_patch(rect)

        elif type(e) == rft.SBend:
            e1, e2 = e.get_E1(), e.get_E2()
            points = [
                (s0 - np.tan(e1), -1),  # lower-left
                (s1 + np.tan(e2), -1),  # lower-right
                (s1 - np.tan(e2), +1),  # upper-right
                (s0 + np.tan(e1), +1)  # upper-left
            ]
            polygon = patches.Polygon(points, closed=True, facecolor=cm['dipole'])
            self.ax.add_patch(polygon)

        elif type(e) in [rft.Bpm, rft.Screen]:
            x = 0.5 * (s0 + s1)
            self.ax.scatter([x], [0], marker="v", s=250, color=cm["diagnostic"], zorder=5)

        elif type(e) == rft.Lattice and self.latticeframe:
            a, b = -1.1, 1.1
            rect = patches.Rectangle((s0, a), s1 - s0, b - a, linewidth=1, linestyle='dashed', edgecolor='black',
                                     facecolor='None', clip_on=False)
            self.ax.add_patch(rect)

        elif type(e) == rft.Volume and self.volumeframe:
            a, b = -1.1, 1.1
            rect = patches.Rectangle((s0, a), s1 - s0, b - a, linewidth=1, linestyle='dotted', edgecolor='black',
                                     facecolor='None', clip_on=False)
            self.ax.add_patch(rect)

        if not isinstance(e, (rft.Drift)):
            x = 0.5 * (s0 + s1)
            if type(e) == rft.Quadrupole:
                if e.get_gradient() > 0:
                    y = 0.85
                    va = "bottom"
                else:
                    y = -0.85
                    va = "top"
            elif type(e) in [rft.Bpm, rft.Screen]:
                y = 0.1
                va = "bottom"
            elif type(e)  == rft.Multipole:
                y = 0.35
                va = "bottom"
            else:
                y = 0.1
                va = "bottom"
            self.ax.text(x, y, e.get_name(), rotation=90, ha="center", va=va, fontsize=9, clip_on=True)