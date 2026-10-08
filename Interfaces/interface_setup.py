INTERFACE_SETUP = {
    "ATF2": [
        {
            "display_name": "ATF2 Damping Ring",
            "module": "Interfaces.ATF2.InterfaceATF2_DR",
            "class_name": "InterfaceATF2_DR",
            "clock_timezone": "Asia/Tokyo",
            "beam_change": {
                "energy": {
                    "label": "Change energy",
                    "tooltip": "The frequency offset of RAMP/PL4 is used to calculate dP/P.",
                    "test": {"label": "Δf [kHz]", "attribute": "energy_frequency_offset_khz", "default": 4.0},
                },
                "intensity": {
                    "label": "Change intensity",
                    "tooltip": "Sets the test laser-intensity factor. Reset restores the laser setting captured when the interface was opened.",
                    "test": {"label": "Test factor", "attribute": "test_laser_intensity", "default": 0.125},
                },
            },
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
                "energy_pref": [1000.0, 2000.0],
            }
        },
        {
            "display_name": "ATF2 Damping Ring RFTrack",
            "module": "Interfaces.ATF2.InterfaceATF2_DR_RFTrack",
            "class_name": "InterfaceATF2_DR_RFTrack",
            "clock_timezone": "Europe/Zurich",
            "beam_change": {
                "energy": {"label": "Change energy", "tooltip": "Multiplies RFTrack reference momentum Pref for the DFS measurement.", "test": {"label": "Pref factor", "attribute": "dfs_test_energy", "default": 0.98}},
                "intensity": {"label": "Change intensity", "tooltip": "Multiplies the RFTrack bunch charge for the WFS measurement.", "test": {"label": "Charge factor", "attribute": "wfs_test_charge", "default": 0.90}},
            },
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
                "energy_pref": [1000.0, 2000.0],
            }
        },
        {
            "display_name": "ATF2 Extraction Line",
            "module": "Interfaces.ATF2.InterfaceATF2_Ext",
            "class_name": "InterfaceATF2_Ext",
            "clock_timezone": "Asia/Tokyo",
            "beam_change": {
                "energy": {
                    "label": "Change energy",
                    "tooltip": "The frequency offset of RAMP/PL4 is used to calculate dP/P.",
                    "test": {"label": "Δf [kHz]", "attribute": "energy_frequency_offset_khz", "default": 4.0},
                },
                "intensity": {
                    "label": "Change intensity",
                    "tooltip": "Sets the nominal and test laser-intensity factors used for WFS.",
                    "nominal": {"label": "Nominal factor", "attribute": "nominal_laser_intensity", "default": 0.1},
                    "test": {"label": "Test factor", "attribute": "test_laser_intensity", "default": 0.125},
                },
            },
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
                    "energy_pref": [1000.0, 2000.0],
                }
        },
        {
            "display_name": "ATF2 Extraction Line RFTrack",
            "module": "Interfaces.ATF2.InterfaceATF2_Ext_RFTrack",
            "class_name": "InterfaceATF2_Ext_RFTrack",
            "clock_timezone": "Europe/Zurich",
            "beam_change": {
                "energy": {"label": "Change energy", "tooltip": "Multiplies RFTrack reference momentum Pref for the DFS measurement.", "test": {"label": "Pref factor", "attribute": "dfs_test_energy", "default": 0.98}},
                "intensity": {"label": "Change intensity", "tooltip": "Multiplies the RFTrack bunch charge for the WFS measurement.", "test": {"label": "Charge factor", "attribute": "wfs_test_charge", "default": 0.90}},
            },
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
            "energy_pref": [1000.0, 2000.0],
        }

        },

        {
            "display_name": "ATF2 Linac",
            "module": "Interfaces.ATF2.InterfaceATF2_Linac",
            "class_name": "InterfaceATF2_Linac",
            "clock_timezone": "Asia/Tokyo",
            "beam_change": {
                "energy": {
                    "label": "Change energy",
                    "tooltip": "Writes CM1L RF phase. Nominal is restored after the DFS measurement.",
                    "nominal": {"label": "Nominal [deg]", "attribute": "phase_kl1", "default": 0.0},
                    "test": {"label": "Test [deg]", "attribute": "cm1l_test_phase", "default": 5.0},
                },
                "intensity": {"label": "Change intensity (RF gun laser)", "tooltip": "Sets the test laser-intensity factor. Reset restores the captured laser setting.", "test": {"label": "Test factor", "attribute": "test_laser_intensity", "default": 0.15}},
            },
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
                "energy_pref": [1000.0, 2000.0],
            }
        },
        {
            "display_name": "ATF2 Linac Beam Transport",
            "module": "Interfaces.ATF2.InterfaceATF2_LinacBT",
            "class_name": "InterfaceATF2_LinacBT",
            "clock_timezone": "Asia/Tokyo",
            "beam_change": {},
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
                "energy_pref": [1000.0, 2000.0],
            }
        },
        {
            "display_name": "ATF2 Linac RFTrack",
            "module": "Interfaces.ATF2.InterfaceATF2_Linac_RFTrack",
            "class_name": "InterfaceATF2_Linac_RFTrack",
            "clock_timezone": "Europe/Zurich",
            "beam_change": {
                "energy": {"label": "Change energy", "tooltip": "Multiplies RFTrack reference momentum Pref for the DFS measurement.", "test": {"label": "Pref factor", "attribute": "dfs_test_energy", "default": 0.98}},
                "intensity": {"label": "Change intensity", "tooltip": "Multiplies the RFTrack bunch charge for the WFS measurement.", "test": {"label": "Charge factor", "attribute": "wfs_test_charge", "default": 0.90}},
            },
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
                "energy_pref": [1000.0, 2000.0],
            }
        },
        {
            "display_name": "ATF2 Beam Transport RFTrack (SAD, 1.3 GeV)",
            "module": "Interfaces.ATF2.InterfaceATF2_BT_RFTrack",
            "class_name": "InterfaceATF2_BT_RFTrack",
            "clock_timezone": "Asia/Tokyo",
            "settings": {"jitter": 0.00, "bpm_resolution": 0.00, "nsamples": 1, "momentum_mev_c": 1300.0},
            "units":
                {
                    "corrector_strength": "T*mm",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 0.01
                },
            "bounds": {
                "emit_x_norm": [1e-4, 10.0],
                "beta_x0": [0.05, 60.0],
                "alpha_x0": [-15.0, 15.0],
                "emit_y_norm": [1e-4, 10.0],
                "beta_y0": [0.05, 60.0],
                "alpha_y0": [-15.0, 15.0],
            }
        },
    ],

    "FACET2": [
        {
            "display_name": "FACET2 Linac RFTrack",
            "module": "Interfaces.FACET2.InterfaceFACET2_Linac_RFTrack",
            "class_name": "InterfaceFACET2_Linac_RFTrack",
            "clock_timezone": "Europe/Zurich",
            "beam_change": {
                "energy": {"label": "Change energy", "tooltip": "Multiplies RFTrack reference momentum Pref for the DFS measurement.", "test": {"label": "Pref factor", "attribute": "dfs_test_energy", "default": 0.98}},
                "intensity": {"label": "Change intensity", "tooltip": "Multiplies the RFTrack bunch charge for the WFS measurement.", "test": {"label": "Charge factor", "attribute": "wfs_test_charge", "default": 0.90}},
            },
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
            "beam_change": {
                "energy": {
                    "label": "Change energy",
                    "tooltip": "Writes FACET2 BC11, BC14 and BC20 feedback-vernier setpoints. Reset returns all three to zero.",
                    "nominal": {"label": "BC11 [MeV]", "attribute": "bba_bc11_energy_offset_mev", "default": -3.0},
                    "test": {"label": "BC14/20 [MeV]", "attribute": "bba_downstream_energy_offset_mev", "default": -40.0},
                },
                "intensity": {"label": "Change intensity (UV waveplate)", "tooltip": "Offsets the FACET2 UV waveplate; charge is then measured and used as the new charge setpoint.", "test": {"label": "UVWP Δ [deg]", "attribute": "bba_uvwp_offset_deg", "default": -2.5}},
            },
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
            "beam_change": {
                "energy": {
                    "label": "Change energy",
                    "tooltip": "Writes MKS11 PhaseSh_SP. Nominal is restored after the DFS measurement.",
                    "nominal": {"label": "Nominal [deg]", "attribute": "rf_phase_nominal", "default": 125.0},
                    "test": {"label": "Test [deg]", "attribute": "rf_phase_test", "default": 145.0},
                },
                "intensity": {
                    "label": "Scale quadrupoles",
                    "tooltip": "Enter the scaling factor of quadrupoles. It will scale the current quadrupole current to a scale * currents.",
                    "test": {"label": "Scaling factor", "attribute": "quadrupole_scaling_factor", "default": 0.5},
                },
            },
            "settings": {"nsamples":3},
            "units":
                {
                    "corrector_strength": "A",
                    "quadrupole_strength": "A",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 1,
                    "bba_max_h_strength": 9.99,
                    "bba_max_v_strength": 9.99,

                },
            "bounds": {
                "emit_x_norm": [0.0, 50.0],
                "beta_x0": [0.001, 200.0],
                "alpha_x0": [-100.0, 100.0],
                "emit_y_norm": [0.0, 50.0],
                "beta_y0": [0.001, 200.0],
                "alpha_y0": [-100.0, 100.0],
                "energy_pref": [150.0, 250.0],
            }
        },
        {
            "display_name": "CLEAR RFTrack",
            "module": "Interfaces.CLEAR.InterfaceCLEAR_RFTrack",
            "class_name": "InterfaceCLEAR_RFTrack",
            "clock_timezone": "Europe/Zurich",
            "beam_change": {
                "energy": {"label": "Change energy", "tooltip": "Multiplies RFTrack reference momentum Pref for the DFS measurement.", "test": {"label": "Pref factor", "attribute": "dfs_test_energy", "default": 0.90}},
                "intensity": {"label": "Change intensity", "tooltip": "Multiplies the RFTrack bunch charge for the WFS measurement.", "test": {"label": "Charge factor", "attribute": "wfs_test_charge", "default": 0.90}},
            },
            "settings": {"jitter": 0.0, "bpm_resolution": 0.0, "nsamples": 1},
            "units":
                {
                    "corrector_strength": "gauss*m",
                    "bpm_position": "mm",
                    "sysid_corrector_kick": 0.01,
                    "bba_max_h_strength": 9.99,
                    "bba_max_v_strength": 9.99,
                },
            "bounds": {
                "emit_x_norm": [0.0, 50.0],
                "beta_x0": [0.001, 200.0],
                "alpha_x0": [-100.0, 100.0],
                "emit_y_norm": [0.0, 50.0],
                "beta_y0": [0.001, 200.0],
                "alpha_y0": [-100.0, 100.0],
                "energy_pref": [150.0, 250.0],
            }
        },
    ],
}
