from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

USE_REAL_MACHINE = False

if USE_REAL_MACHINE:
    from Interfaces.CLEAR.InterfaceCLEAR import CLEAR_real_machine

    machine = CLEAR_real_machine(nsamples=3)
else:
    from Interfaces.CLEAR.InterfaceCLEAR_RFTrack import InterfaceCLEAR_RFTrack

    machine = InterfaceCLEAR_RFTrack(nsamples=1)


def show(label, data):
    print(f"\n{label}")
    for key, value in data.items():
        if key not in {"images", "beam_images", "background_images", "hedges", "vedges"}:
            print(f"  {key}: {value}")


print(f"Connected to: {machine.get_name()}")
print("Beamline:", machine.get_sequence()[:5], "...")

corrector = machine.get_hcorrectors_names()[0]
bpm = machine.bpms[0]
quad = machine.quadrupoles[0]
screen = machine.screens[0]

show("Corrector settings/readbacks", machine.get_correctors([corrector]))
show("BPM position and intensity", machine.get_bpms([bpm]))
show("Charge monitor", machine.get_icts())
show("Quadrupole strength", machine.get_quadrupoles([quad]))
if not USE_REAL_MACHINE:
    show("Screen centroid, size and image", machine.get_screens([screen]))

machine.set_correctors([corrector], [0.0])     # set absolute corrector value
machine.vary_correctors([corrector], [0.01])   # add a corrector increment
machine.set_quadrupoles([quad], [0.0])         # set K1L quadrupole strength
machine.change_energy(); machine.reset_energy()       # temporary energy change
machine.change_intensity(); machine.reset_intensity() # temporary charge change
show("Screen data", machine.get_screens([screen]))   # inserts/moves BTV on real CLEAR

state = machine.get_state()

machine.restore_correctors_state(state)
machine.restore_quadrupoles_state(state)
machine.restore_beam_settings(state.get_beam_settings())
