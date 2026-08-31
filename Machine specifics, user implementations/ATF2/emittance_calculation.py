
import numpy as np
import math
import os
import glob

class Emitt_Meas_Simulation:
    def __init__(self, filename='ATF2_EXT_FF_v5.2.twiss'):
        self.Pref = 1.2999999e3
        self.filename = filename

    def get_data_from_twiss_file(self):
        with open(self.filename, "r") as file:
            lines = [line.strip() for line in file if line.strip()]
        star_symbol = next(i for i, line in enumerate(lines) if line.startswith("*"))
        dollar_sign = next(i for i, line in enumerate(lines) if line.startswith("$") and i > star_symbol)
        columns = lines[star_symbol].lstrip("*").split()
        BETA_X_column = columns.index("BETX")  # m
        BETA_Y_column = columns.index("BETY")  # m
        ALPHA_X_column = columns.index("ALFX")
        ALPHA_Y_column = columns.index("ALFY")
        MU_X_column = columns.index("MUX")
        MU_Y_column = columns.index("MUY")
        NAME_column = columns.index("NAME")
        names = []
        twiss = {}

        for line in lines[dollar_sign + 1:]:
            data = line.split()
            name = data[NAME_column].strip('"')
            names.append(name)

            twiss[name] = {
                "betx": float(data[BETA_X_column]),
                "bety": float(data[BETA_Y_column]),
                "alpx": float(data[ALPHA_X_column]),
                "alpy": float(data[ALPHA_Y_column]),
                "mux": float(data[MU_X_column]),
                "muy": float(data[MU_Y_column]),
            }

        entrance_index = names.index("ATF2$START")
        entrance_name = names[entrance_index]
        entrance = twiss[entrance_name]
        otr_names = ['OTR0X', 'OTR1X', 'OTR2X', 'OTR3X']
        otrs = {name: twiss[name] for name in otr_names}
        return entrance_name, entrance, otrs

    def compute_transport_matrix(self):
        """Computes transport matrices for 2D (Mx, My) and Coupling (Mxy)."""
        entrance_name, entrance, otrs = self.get_data_from_twiss_file()

        beta_x_0 = entrance["betx"]
        beta_y_0 = entrance["bety"]

        Mx_matrix_rows = []
        My_matrix_rows = []
        Mxy_matrix_rows = []

        for otr in otrs:
            delta_mux = (otrs[otr]["mux"] - entrance["mux"]) * 2 * math.pi
            delta_muy = (otrs[otr]["muy"] - entrance["muy"]) * 2 * math.pi

            R11_x = np.sqrt((otrs[otr]["betx"]) / beta_x_0) * (math.cos(delta_mux) + entrance["alpx"] * math.sin(delta_mux))
            R12_x = np.sqrt(beta_x_0 * otrs[otr]["betx"]) * math.sin(delta_mux)

            R33_y = np.sqrt((otrs[otr]["bety"]) / beta_y_0) * (math.cos(delta_muy) + entrance["alpy"] * math.sin(delta_muy))
            R34_y = np.sqrt(beta_y_0 * otrs[otr]["bety"]) * math.sin(delta_muy)

            Mx_matrix_rows.append([R11_x**2, 2 * R11_x * R12_x, R12_x**2])
            My_matrix_rows.append([R33_y**2, 2 * R33_y * R34_y, R34_y**2])
            Mxy_matrix_rows.append([R11_x * R33_y, R11_x * R34_y, R12_x * R33_y, R12_x * R34_y])

        return np.array(Mx_matrix_rows), np.array(My_matrix_rows), np.array(Mxy_matrix_rows)
    
    def solve_least_squares(self, Mx, My, sigma_x_i, sigma_y_i):

        sigma_x_meas, *_ = np.linalg.lstsq(Mx, sigma_x_i, rcond=None)
        sigma_y_meas, *_ = np.linalg.lstsq(My, sigma_y_i, rcond=None)

        sigma1, sigma2, sigma5 = sigma_x_meas
        sigma8, sigma9, sigma10 = sigma_y_meas

        emittance_x = np.sqrt(np.abs(sigma1 * sigma5 - sigma2**2))
        emittance_y = np.sqrt(np.abs(sigma8 * sigma10 - sigma9**2))

        return emittance_x, emittance_y

    def solve_4d_least_squares(self, Mx, My, Mxy, sigma_x_i, sigma_y_i, sigma_xy_i):        
        s_x_meas, *_ = np.linalg.lstsq(Mx, sigma_x_i, rcond=None)
        s1, s2, s5 = s_x_meas
        
        s_y_meas, *_ = np.linalg.lstsq(My, sigma_y_i, rcond=None)
        s8, s9, s10 = s_y_meas
        
        s_xy_meas, *_ = np.linalg.lstsq(Mxy, sigma_xy_i, rcond=None)
        s3, s4, s6, s7 = s_xy_meas

        S = np.array([
            [s1, s2, s3, s4],
            [s2, s5, s6, s7],
            [s3, s6, s8, s9],
            [s4, s7, s9, s10]
        ])

        J = np.array([
            [0, -1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, -1],
            [0, 0, 1, 0]
        ])

        JS = J @ S
        JS2 = np.linalg.matrix_power(JS, 2)
        TJS2 = np.trace(JS2)
        DS = np.linalg.det(S)
        
        discriminant = TJS2**2 - 16 * DS
        if discriminant < 0 and abs(discriminant) < 1e-20:
            discriminant = 0.0

        sqrt_discriminant = np.sqrt(discriminant)
        term_inner_1 = -TJS2 + sqrt_discriminant
        term_inner_2 = -TJS2 - sqrt_discriminant
        
        e1 = 0.5 * np.sqrt(abs(term_inner_1))
        e2 = 0.5 * np.sqrt(abs(term_inner_2))
        
        return e1, e2, S


def load_sigmas_from_npz(filepath):
    """Extracts sigmas from the saved OTR .npz file and formats them in m^2 for tracking matrices."""
    data = np.load(filepath)
    
    sigma_x_m2 = []
    sigma_y_m2 = []
    sigma_xy_m2 = []
    
    for i in range(4):
        sig_h_um = data[f'OTR{i}_SigmaH_um']
        sig_v_um = data[f'OTR{i}_SigmaV_um']
        
        sigma_x_m2.append((sig_h_um * 1e-6)**2)
        sigma_y_m2.append((sig_v_um * 1e-6)**2)
        
        sig_13_m2 = data[f'OTR{i}_Sigma13_m2']
        sigma_xy_m2.append(sig_13_m2)
        
    return np.array(sigma_x_m2), np.array(sigma_y_m2), np.array(sigma_xy_m2)


if __name__ == "__main__":
    sim = Emitt_Meas_Simulation(filename='ATF2_EXT_FF_v5.2.twiss')
    save_folder = 'Data_mOTR'
    search_path = os.path.join(save_folder, '*.npz')
    npz_files = glob.glob(search_path)
    
    if not npz_files:
        print(f"Error: No .npz files found in {save_folder}/. Please run the OTR image analysis script first.")
    else:
        latest_file = max(npz_files, key=os.path.getctime)
        print(f"Loading empirical beam sizes from: {latest_file}")
        
        try:
            sigma_x_input, sigma_y_input, sigma_xy_input = load_sigmas_from_npz(latest_file)
            Mx, My, Mxy = sim.compute_transport_matrix()
            
            emitt_x, emitt_y = sim.solve_least_squares(Mx, My, sigma_x_input, sigma_y_input)
            print(f"2D Emittance X: {emitt_x:.4e} m-rad, Y: {emitt_y:.4e} m-rad")
            
            e1, e2, full_S = sim.solve_4d_least_squares(Mx, My, Mxy, sigma_x_input, sigma_y_input, sigma_xy_input)
            print(f"4D Eigen-emittances -> e1: {e1:.4e}, e2: {e2:.4e}")
            
            if emitt_x * emitt_y - e1 * e2 < 1e-20:
                print("Warning: non-physical solutions. Check the input.")

        except FileNotFoundError:
            print("Error: Could not find the .twiss file.")
        except KeyError as e:
            print(f"Error: Missing expected key in .npz file: {e}")
        except ValueError as e:
            print(f"Data processing error: {e}")