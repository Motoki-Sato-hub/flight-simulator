from Interfaces.ATF2.InterfaceATF2_Linac import InterfaceATF2_Linac, LINAC_SEQUENCE, _build_beamline

BT_SEQUENCE = LINAC_SEQUENCE + [
    # Transport-line BPM/corrector order from Interfaces/ATF2/Linac_ATF2/atfbt199912_daihon.sad.
    # CELLST defines ML10T / ML11T in SAD, but LINAC:monitors currently exposes the end monitors
    # as MB10T / MB11T, so keep the existing EPICS-facing names here.
    "ZX10T", "ZX11T", "ML1T", "ZH10T", "ZV11T", "MB1T", "ZX12T",
    "ML2T", "ZY20T", "ZY21T", "ML101T", "ML102T", "ZY22T", "ZY23T", "ML103T",
    "ZX30T", "ML3T", "ZX31T", "ZV30T", "ZH30T", "ML104T", "ZX32T", "ML4T",
    "ML105T", "ZV40T", "ZH40T", "ML5T", "ML6T",
    "ZX50T", "ML106T", "ZX51T", "ML7T", "ZX52T", "ZV50T", "ML8T", "ZH50T", "ZV51T", "ML9T",
    "MB10T", "MB11T",
]


class InterfaceATF2_LinacBT(InterfaceATF2_Linac):
    """Linac + BT machine view used by LinacOpt."""

    def get_name(self):
        return 'ATF2_LinacBT' # Beam transport

    def __init__(self, nsamples=1):
        super().__init__(nsamples=nsamples)
        self.sequence, self.bpms, self.corrs, self.bpm_indexes = _build_beamline(sequence=BT_SEQUENCE, monitors=self.monitors)

    def get_hcorrectors_names(self):
        return [string for string in self.corrs if string.lower().startswith(('zh', 'zx'))]

    def get_vcorrectors_names(self):
        return [string for string in self.corrs if string.lower().startswith(('zv', 'zy'))]
