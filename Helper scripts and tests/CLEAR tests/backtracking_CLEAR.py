import RF_Track as rft
from Interfaces.CLEAR.InterfaceCLEAR_RFTrack import InterfaceCLEAR_RFTrack

interface = InterfaceCLEAR_RFTrack()

Twiss = rft.Bunch6d_twiss()
Twiss.beta_x = 17.9        # m
Twiss.beta_y = 8.07 # 2.73        # m
Twiss.alpha_x = 0.236
Twiss.alpha_y = -1.72 # 0.339
Twiss.emitt_x = 12.7     # mm.mrad normalised emittance
Twiss.emitt_y = 4.34 # 0.727      # mm.mrad
Twiss.sigma_t = 0
Twiss.sigma_pt = 0
Twiss.mean_xp = 0.0
Twiss.mean_yp = 0

B_qfd350 = rft.Bunch6d_QR(rft.electronmass, interface.population, interface.Q, interface.Pref, Twiss, interface.nparticles)

end_index = interface.sequence.index("CA.QS0350")
L = rft.Lattice()

for element in list(interface.lattice["*"])[:end_index + 1]:
    L.append(element)

B_stline = L.btrack(B_qfd350)

info = B_stline.get_info()
print(f"beta_x  = {info.beta_x:.4f} m")
print(f"alpha_x = {info.alpha_x:.4f}")
print(f"beta_y  = {info.beta_y:.4f} m")
print(f"alpha_y = {info.alpha_y:.4f}")
print(f"emit_x  = {info.emitt_x:.4f} m")
print(f"emit_y  = {info.emitt_y:.4f} m")

B_check = L.track(B_stline)