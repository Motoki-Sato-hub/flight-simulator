print("================================================================================")
print("Testing BPM CA.BPM0530, reading its values directly from interface, using get_bpms method... ")

result = I.get_bpms(names=["CA.BPM0530"])
print(f"names from method get_bpms: {result['names']}")
print(f"bdes  from method get_bpms: {result['bdes']}")
print(f"bact  from method get_bpms: {result['bact']}")


# print("Methods to test Orbit Measurement.")
# state = State(correctors=correctors)
# print("Reading method 'get_orbit'...")
# bpms = I.get_bpms(names=["CA.BPM0530H-AS"])
# orbit = state.get_orbit()
# print("Successfully run method 'get_orbit'!")
# print("Reading method 'get_state'...")
# I.get_state()
# print("Successfully run method 'get_state'!")
# print("Reading method 'get_bpms'...")
# I.get_bpms()
# print("Successfully run method 'get_bpms'!")
# print("Orbit Measurement is able to run.")

print("Trying to read corrector current...")
print(japc.getParam("CA.DHG0130/SettingPPM#current")) # this works!, A
print(japc.getParam("CA.DHG0130/Status"))
print("Reading method 'get_correctors'...")
I.get_correctors()
print("Successfully run method 'get_correctors'!")

print("Methods to test reading screens.")
print("Reading method 'get_screens' with inserted screen...")
#I.get_screens()
print("Last image reading...")
#pixelCalSet1 or 2
#print(japc.getParam("CA.BTV0125.DigiCam/LastImage#image2D"))
print(japc.getParam("CA.BTV0390_CAS.BTV0420/OPSettingSystem1#positionChannel1")) # anything else than 0 meanssscreen is in

# 0 out, 1 in
print("Successfully run method 'get_screens' with inserted screen!")

print("Reading method 'get_screens' without inserted screen (no beam)...")
#I.get_screens()
print("Successfully run method 'get_screens' without inserted screen (no beam)!")

print(japc.getParam("CA.QFD0350/SettingPPM#current")) # this works!, A


