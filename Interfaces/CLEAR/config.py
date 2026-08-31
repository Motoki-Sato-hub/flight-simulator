# If true, no setpoints will be changed on the machine
no_set = False

# Number of images for background measurement
n_background_images = 5

# Set of cameras available in the camera widget
cameras = {
    'CA.BTV0390L': {
        'japc_name': 'CA.BTV0390',
        'japc_selector': '',

        # Image (width, height) before taking into account image rotation
        'image_size': (1936, 1216),

        's_x_res': 0.045, 
        'delta_s_x_res': 0.0376,

        's_y_res': 0.045,
        'delta_s_y_res': 0.0347,
    },

    'CA.BTV0390H': {
        'japc_name': 'CA.BTV0390',
        'japc_selector': '',

        # Image (width, height) before taking into account image rotation
        'image_size': (1936, 1216),

        's_x_res': 0.0376, 
        'delta_s_x_res': 0.0376,

        's_y_res': 0.0347,
        'delta_s_y_res': 0.0347,
    },

    'CA.BTV0620': {
        'japc_name': 'CA.BTV0620',
        'japc_selector': '',
        
        # Image (width, height) before taking into account image rotation
        'image_size': (1936, 1216),

        's_x_res': 0.08, 
        'delta_s_x_res': 0.08,

        's_y_res': 0.08,
        'delta_s_y_res': 0.08,
    },

    'CA.BTV0730': {
        'japc_name': 'CA.BTV0730',
        'japc_selector': '',

        # Image (width, height) before taking into account image rotation
        'image_size': (1936, 1216),

        's_x_res': 0.220,  #0.220? Anto/Giac to be checked 
        'delta_s_x_res': 0.050,

        's_y_res': 0.080, #0.080? Anto/Giac to be checked 
        'delta_s_y_res': 0.03,
    },

    #'CA.BTV0805': {
    #    'japc_name': 'CA.BTV0805',
    #    'japc_selector': '',
    #
    #    # Image (width, height) before taking into account image rotation
    #    'image_size': (258, 231),
    #
    #    # Pixels, (from, to)
    #    # Whole screen, including screws
    #    'AOI_0': (0, 0),
    #    'AOI_1': (258, 231),
    #
    #    's_x_res': 0.0134, 
    #    'delta_s_x_res': 0.050,
    #
    #   's_y_res': 0.0128,
    #    'delta_s_y_res': 0.03,
    #},

    'CA.BTV0810': {
        'japc_name': 'CA.BTV0810',
        'japc_selector': '',

        # Image (width, height) before taking into account image rotation
        'image_size': (1936, 1216),

        # Pixels, (from, to)
        # Whole screen, including screws
        'AOI_0': (129, 1100),
        'AOI_1': (509, 1785),

        's_x_res': 0.220, 
        'delta_s_x_res': 0.050,

        's_y_res': 0.080,
        'delta_s_y_res': 0.03,
    },

    'CA.BTV0910': {
        'japc_name': 'CA.BTV0910',
        'japc_selector': '',

        # Image (width, height) before taking into account image rotation
        'image_size': (1936, 1216),

        's_x_res': 0.220, 
        'delta_s_x_res': 0.050,

        's_y_res': 0.060,
        'delta_s_y_res': 0.03,
        }
}

# Set of cameras available in the calibration widget
calibrated_cameras = [
    'CA.BTV0125',
    'CA.BTV0215',
    'CA.BTV0390L',
    'CA.BTV0390H',
    'CAS.BTV0420',
    'CA.BTV0620',
    'CA.BTV0730',
    'CA.BTV0805',
    'CA.BTV0810',
    'CA.BTV0910'
    # CA.BTV0235, 545, 800, 805, 875,

]

# Quadrupoole magnet calibration
quad_length = 0.226

# Integrated transfer function
def get_ITF(I):
    return 1.29404711e-2  - 2.59458259e-07*I # T/A

# A value of 0.05 means that the current should not differ from the setpoint by more than 5% of the smallest scan step size
current_tolerance = 0.05

# Number of extra measurements at the start and end of the scan range
margin_measurements = 2

current_set_params = [
    'CA.QFD0350/SettingPPM#current',
    'CA.QDD0355/SettingPPM#current',
    'CA.QFD0360/SettingPPM#current',

    'CA.QFD0510/SettingPPM#current',
    'CA.QDD0515/SettingPPM#current',
    'CA.QFD0520/SettingPPM#current',

    'CA.QFD0760/SettingPPM#current',
    'CA.QDD0765/SettingPPM#current',
    'CA.QFD0770/SettingPPM#current',

    'CA.QDD0870/SettingPPM#current',
    'CA.QFD0880/SettingPPM#current'
]


current_get_params = [
    'CA.QFD0350/Acquisition#currentAverage',
    'CA.QDD0355/Acquisition#currentAverage',
    'CA.QFD0360/Acquisition#currentAverage',

    'CA.QFD0510/Acquisition#currentAverage',
    'CA.QDD0515/Acquisition#currentAverage',
    'CA.QFD0520/Acquisition#currentAverage',

    'CA.QFD0760/Acquisition#currentAverage',
    'CA.QDD0765/Acquisition#currentAverage',
    'CA.QFD0770/Acquisition#currentAverage',

    'CA.QDD0870/Acquisition#currentAverage',
    'CA.QFD0880/Acquisition#currentAverage'
]

# Selector ''
current_status_params = [
    'CA.QFD0350/Status#mode',
    'CA.QDD0355/Status#mode',
    'CA.QFD0360/Status#mode',

    'CA.QFD0510/Status#mode',
    'CA.QDD0515/Status#mode',
    'CA.QFD0520/Status#mode',

    'CA.QFD0760/Status#mode',
    'CA.QDD0765/Status#mode',
    'CA.QFD0770/Status#mode',

    'CA.QDD0870/Status#mode',
    'CA.QFD0880/Status#mode'
]

quad_names = [
    'CA.QFD0350', 
    'CA.QDD0355', 
    'CA.QFD0360',

    'CA.QFD0510', 
    'CA.QDD0515', 
    'CA.QFD0520',

    'CA.QFD0760', 
    'CA.QDD0765', 
    'CA.QFD0770',

    'CA.QDD0870', 
    'CA.QFD0880'
]