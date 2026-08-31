INTERFACE_SETUP = {
    "ATF2": [
        {
            "display_name": "ATF2 Damping Ring",
            "module": "Interfaces.ATF2.InterfaceATF2_DR",
            "class_name": "InterfaceATF2_DR",
            "clock_timezone": "Asia/Tokyo",
            "settings": {"nsamples": 3},
            "units":
            {
                "corrector_strength":"T*mm",
                "bpm_position": "mm",
                "sysid_corrector_kick": 0.01
            },
            "bounds": {
                "emit_x_norm": [1.0, 10.0],
                "beta_x0": [1.0, 12.0],
                "alpha_x0": [-10.0, 3.0],
                "emit_y_norm": [0.005, 0.15],
                "beta_y0": [0.5, 8.0],
                "alpha_y0": [-4.0, 6.0],
            }
        },
        {
            "display_name": "ATF2 Damping Ring RFTrack",
            "module": "Interfaces.ATF2.InterfaceATF2_DR_RFTrack",
            "class_name": "InterfaceATF2_DR_RFTrack",
            "clock_timezone": "Europe/Zurich",
            "settings": {"jitter":0.0, "bpm_resolution":0.0, "nsamples":1},
            "units":
                {
                    "corrector_strength": "T*mm",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 0.01
                },
            "bounds": {
                "emit_x_norm": [1.0, 10.0],
                "beta_x0": [1.0, 12.0],
                "alpha_x0": [-10.0, 3.0],
                "emit_y_norm": [0.005, 0.15],
                "beta_y0": [0.5, 8.0],
                "alpha_y0": [-4.0, 6.0],
            }
        },
        {
            "display_name": "ATF2 Extraction Line",
            "module": "Interfaces.ATF2.InterfaceATF2_Ext",
            "class_name": "InterfaceATF2_Ext",
            "clock_timezone": "Asia/Tokyo",
            "settings": {"nsamples": 3},
            "units":
                {
                    "corrector_strength": "T*mm",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 0.01,
                    "em_sigma_unit": "mm"
                },
            "bounds":
                {
                    "emit_x_norm": [0.5, 8.0],
                    "beta_x0": [0.2, 5.0],
                    "alpha_x0": [-4.0, 2.0],
                    "emit_y_norm": [0.005, 0.12],
                    "beta_y0": [2.0, 20.0],
                    "alpha_y0": [-8.0, 2.0],
                }
        },
        {
            "display_name": "ATF2 Extraction Line RFTrack",
            "module": "Interfaces.ATF2.InterfaceATF2_Ext_RFTrack",
            "class_name": "InterfaceATF2_Ext_RFTrack",
            "clock_timezone": "Europe/Zurich",
            "settings": {"jitter":0.0, "bpm_resolution":0.00, "nsamples":1},
            "units":
                {
                    "corrector_strength": "T*mm",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 0.01,
                    "em_sigma_unit": "mm"

                },

        "bounds": {
            "emit_x_norm": [5, 10.0],
            "beta_x0": [0.001, 60.0],
            "alpha_x0": [-15.0, 15.0],
            "emit_y_norm": [0.01, 5.0],
            "beta_y0": [0.001, 60.0],
            "alpha_y0": [-15, 15.0],
        }

        },

        {
            "display_name": "ATF2 Linac",
            "module": "Interfaces.ATF2.InterfaceATF2_Linac",
            "class_name": "InterfaceATF2_Linac",
            "clock_timezone": "Asia/Tokyo",
            "settings": {"nsamples":3},
            "units":
                {
                    "corrector_strength": "T*mm",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 0.01
                },
            "bounds": {
                "emit_x_norm": [1e-4, 0.05],
                "beta_x0": [0.05, 5.0],
                "alpha_x0": [-6.0, 2.0],
                "emit_y_norm": [1e-4, 0.05],
                "beta_y0": [0.2, 8.0],
                "alpha_y0": [-15.0, 2.0],
            }
        },
        {
            "display_name": "ATF2 Linac Beam Transport",
            "module": "Interfaces.ATF2.InterfaceATF2_LinacBT",
            "class_name": "InterfaceATF2_LinacBT",
            "clock_timezone": "Asia/Tokyo",
            "settings": {"nsamples": 3},
            "units":
                {
                    "corrector_strength": "T*mm",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 0.01
                },
            "bounds": {
                "emit_x_norm": [1e-4, 0.05],
                "beta_x0": [0.05, 5.0],
                "alpha_x0": [-6.0, 2.0],
                "emit_y_norm": [1e-4, 0.05],
                "beta_y0": [0.2, 8.0],
                "alpha_y0": [-15.0, 2.0],
            }
        },
        {
            "display_name": "ATF2 Linac RFTrack",
            "module": "Interfaces.ATF2.InterfaceATF2_Linac_RFTrack",
            "class_name": "InterfaceATF2_Linac_RFTrack",
            "clock_timezone": "Europe/Zurich",
            "settings": {"jitter": 0.00, "bpm_resolution": 0.00},
            "units":
                {
                    "corrector_strength": "T*mm",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 0.01
                },
            "bounds": {
                "emit_x_norm": [1e-4, 0.05],
                "beta_x0": [0.05, 5.0],
                "alpha_x0": [-6.0, 2.0],
                "emit_y_norm": [1e-4, 0.05],
                "beta_y0": [0.2, 8.0],
                "alpha_y0": [-15.0, 2.0],
            }
        },
    ],

    "FACET2": [
        {
            "display_name": "FACET2 Linac RFTrack",
            "module": "Interfaces.FACET2.InterfaceFACET2_Linac_RFTrack",
            "class_name": "InterfaceFACET2_Linac_RFTrack",
            "clock_timezone": "Europe/Zurich",
            "settings": {"jitter":0.0, "bpm_resolution":0.0, "nsamples":1},
            "units":
                {
                    "corrector_strength": "gauss*m",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 0.0001
                },
            "bounds": {
                "emit_x_norm": [0.5, 8.0],
                "beta_x0": [0.5, 12.0],
                "alpha_x0": [-5.0, 3.0],
                "emit_y_norm": [0.5, 8.0],
                "beta_y0": [0.5, 12.0],
                "alpha_y0": [-5.0, 5.0],
            }
        },
        {
            "display_name": "FACET2 Linac",
            "module": "Interfaces.FACET2.InterfaceFACET2_Linac",
            "class_name": "InterfaceFACET2_Linac",
            "clock_timezone": "America/Los_Angeles",
            "settings": {"nsamples":10},
            "units":
                {
                    "corrector_strength": "gauss*m",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 0.0001
                },
            "bounds": {
                "emit_x_norm": [0.5, 8.0],
                "beta_x0": [0.5, 12.0],
                "alpha_x0": [-5.0, 3.0],
                "emit_y_norm": [0.5, 8.0],
                "beta_y0": [0.5, 12.0],
                "alpha_y0": [-5.0, 5.0],
            }
        },],

    "CLEAR": [
        {
            "display_name": "CLEAR",
            "module": "Interfaces.CLEAR.InterfaceCLEAR",
            "class_name": "CLEAR_real_machine",
            "clock_timezone": "Europe/Zurich",
            "settings": {"nsamples":3},
            "units":
                {
                    "corrector_strength": "gauss*m",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 1
                },
            "bounds":
                {
                    "emit_x_norm": [0.5, 10.0],
                    "beta_x0": [0.2, 20.0],
                    "alpha_x0": [-4.0, 2.0],
                    "emit_y_norm": [0.5, 10.0],
                    "beta_y0": [10.0, 30.0],
                    "alpha_y0": [-8.0, 2.0],
                }
        },
        {
            "display_name": "CLEAR RFTrack",
            "module": "Interfaces.CLEAR.InterfaceCLEAR_RFTrack",
            "class_name": "InterfaceCLEAR_RFTrack",
            "clock_timezone": "Europe/Zurich",
            "settings": {"jitter": 0.0, "bpm_resolution": 0.00, "nsamples": 1},
            "units":
                {
                    "corrector_strength": "gauss*m",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 0.01
                },
            "bounds": {
                "emit_x_norm": [0.0, 20.0],
                "beta_x0": [0.001, 200.0],
                "alpha_x0": [-50.0, 50.0],
                "emit_y_norm": [0.0, 20.0],
                "beta_y0": [0.001, 200.0],
                "alpha_y0": [-50.0, 50.0],
            }
        },
    ],
}
