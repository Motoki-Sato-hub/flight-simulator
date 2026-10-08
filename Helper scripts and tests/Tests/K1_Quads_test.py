import RF_Track as rft
import numpy as np
from Interfaces.ATF2.InterfaceATF2_Ext_RFTrack import InterfaceATF2_Ext_RFTrack

I = InterfaceATF2_Ext_RFTrack()
quad_name = "QF17X"
quadrupole_data = I.get_quadrupoles([quad_name])
K1L_0 = float(quadrupole_data["bdes"][0])
deltas = np.linspace(-0.10, 0.10, 5)
K1L_values = K1L_0 * (1 + deltas)
screens = ["OTR0X", "OTR1X", "OTR2X", "OTR3X"]

for K1L in K1L_values:
    I.set_quadrupoles([quad_name],[K1L])
    actual_K1L = I.get_quadrupoles(names=[quad_name])["bdes"][0]
    screen_data = I.get_screens(names=screens)

    print("requested K1L:", K1L)
    print("actual K1L:", actual_K1L)
    print("sigx:" , screen_data['sigx'])
    print("sigy:" , screen_data['sigy'])
print("Result from RF-Track:")

print(f"K1L_0: {K1L_0}")


from Interfaces.ATF2.InterfaceATF2_Ext import InterfaceATF2_Ext

I2 = InterfaceATF2_Ext()
#K1L_0_2 = I2.current_to_k1('QD18X', 31.7475)
K1L_0_2 = I2.current_to_k1('QF17X', 31.7475)

print("Reading from clibmagnet:")
print(f"K1L_0: {K1L_0_2}")





