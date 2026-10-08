"""Compatibility name for the SAD-derived ATF BT RF-Track interface.

Historically this module tried to append an undefined BT lattice to the old
Linac RF-Track class.  The BT has a distinct 1.3 GeV/c transport model now;
the physical Linac-to-BT injection match remains a separate calibration task.
"""

from Interfaces.ATF2.InterfaceATF2_BT_RFTrack import InterfaceATF2_BT_RFTrack


class InterfaceATF2_LinacBT_RFTrack(InterfaceATF2_BT_RFTrack):
    def get_name(self):
        return "ATF2_LinacBT_RFT"
