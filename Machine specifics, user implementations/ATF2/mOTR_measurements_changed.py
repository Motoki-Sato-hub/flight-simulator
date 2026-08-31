import sys
import numpy as np
import time, math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
# from epics import PV, ca 
import os
import matplotlib.patches as patches


# --- Gaussian Fit ---
def gaussian(x, amplitude, mean, stddev, offset):
    """1D Gaussian function for curve fitting."""
    return amplitude * np.exp(-((x - mean) / (2 * stddev))**2) + offset


def plot_otr_analysis_inset(ax_main, img_data, title_str, h_factor_um_px=1.0, v_factor_um_px=1.0):
      
    im = ax_main.imshow(img_data, cmap='gray', origin='lower') 
    
    h, w = img_data.shape
    proj_x = np.sum(img_data, axis=0) 
    proj_y = np.sum(img_data, axis=1) 
    x_coords = np.arange(w)
    y_coords = np.arange(h)
    
    p0_x = [np.max(proj_x), np.argmax(proj_x), w/10, np.min(proj_x)]
    p0_y = [np.max(proj_y), np.argmax(proj_y), h/10, np.min(proj_y)]
    
    sigma_h_um = None
    sigma_v_um = None
    
    cx, cy = w // 2, h // 2 
    radius_x_px, radius_y_px = 150, 150 

    try:
        popt_x, pcov_x = curve_fit(gaussian, x_coords, proj_x, p0=p0_x)
        amp_x, mean_x, stddev_x, offset_x = popt_x
        fit_x = gaussian(x_coords, *popt_x)
        beam_size_x_px = abs(stddev_x) # why 2-sigma value?
        sigma_h_um = beam_size_x_px * h_factor_um_px 
        sigma_h2_m2 = (sigma_h_um * 1e-6) ** 2
        cx = int(mean_x)
        radius_x_px = int(abs(stddev_x) * 3.0)

        popt_y, pcov_y = curve_fit(gaussian, y_coords, proj_y, p0=p0_y)
        amp_y, mean_y, stddev_y, offset_y = popt_y
        fit_y = gaussian(y_coords, *popt_y)
        beam_size_y_px = abs(stddev_y) # 2-sigma
        sigma_v_um = beam_size_y_px * v_factor_um_px 
        sigma_v2_m2 = (sigma_v_um * 1e-6) ** 2
        
        cy = int(mean_y) 
        radius_y_px = int(abs(stddev_y) * 3.0)
        # PLOTS
        ax_hist_x = ax_main.inset_axes([0.0, 1.05, 1.0, 0.2], transform=ax_main.transAxes)
        ax_hist_x.plot(x_coords, proj_x, label='H Projection')
        ax_hist_x.plot(x_coords, fit_x, 'r--', label=f'Size: {sigma_h2_m2:.2e} $m^2$')
        ax_hist_x.legend(fontsize='x-small', loc='upper right')
        ax_hist_x.axis('off')
        
        ax_hist_y = ax_main.inset_axes([1.05, 0.0, 0.2, 1.0], transform=ax_main.transAxes)
        ax_hist_y.plot(proj_y, y_coords, label='V Projection') 
        ax_hist_y.plot(fit_y, y_coords, 'r--', label=f'Size: {sigma_v2_m2:.2e} $m^2$')
        ax_hist_y.legend(fontsize='x-small', loc='lower right')
        ax_hist_y.axis('off')
        
    except RuntimeError:
        print(f"Could not fit Gaussian for {title_str}")
        ax_hist_x = ax_main.inset_axes([0.0, 1.05, 1.0, 0.2], transform=ax_main.transAxes)
        ax_hist_y = ax_main.inset_axes([1.05, 0.0, 0.2, 1.0], transform=ax_main.transAxes)
        ax_hist_x.axis('off')
        ax_hist_y.axis('off')

    # Dynamic ROI Cropping & Sigma_13 Calculation
    clean_img = np.array(img_data, dtype=np.float64)
    
    y_min = max(0, cy - radius_y_px)
    y_max = min(h, cy + radius_y_px)
    x_min = max(0, cx - radius_x_px)
    x_max = min(w, cx + radius_x_px)
    
    roi_img = clean_img[y_min:y_max, x_min:x_max]
    
    rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, 
                             linewidth=1, edgecolor='r', facecolor='none', alpha=0.5)
    ax_main.add_patch(rect)
    
    local_thresh = np.max(roi_img) * 0.05
    roi_img[roi_img < local_thresh] = 0
    
    tot_intensity = np.sum(roi_img)
    sigma_13_m2 = None
    
    if tot_intensity > 0:
        roi_y_idx, roi_x_idx = np.indices(roi_img.shape)
        
        x_bar = np.sum(roi_x_idx * roi_img) / tot_intensity
        y_bar = np.sum(roi_y_idx * roi_img) / tot_intensity
        
        sig13_px2 = np.sum((roi_x_idx - x_bar) * (roi_y_idx - y_bar) * roi_img) / tot_intensity
        
        sigma_13_um2 = sig13_px2 * h_factor_um_px * v_factor_um_px
        sigma_13_m2 = sigma_13_um2 * 1e-12
    if sigma_13_m2 is not None:
        ax_main.set_title(f"{title_str} | $\\sigma_{{13}}$: {sigma_13_m2:.2e} m$^2$")
    else:
        ax_main.set_title(title_str)
        
    # ax_main.axis('off') 
    # plt.colorbar(im, ax=ax_main, fraction=0.046, pad=0.04, label='Pixel Intensity (Counts)')
    
    return proj_x, proj_y, sigma_h_um, sigma_v_um, sigma_13_m2



def get_pixel_calibrations(otr_id_str):
    h_pv_name = f'mOTR{otr_id_str}:H:x1:Calibration:Factor'
    v_pv_name = f'mOTR{otr_id_str}:V:x1:Calibration:Factor' 
    h_factor = PV(h_pv_name).get()
    v_factor = PV(v_pv_name).get()
    if h_factor is None or v_factor is None:
        h_factor = 1.0
        v_factor = 1.0
    return h_factor, v_factor

def acquire_otr_image(otr_id_str, min_total_intensity=130000, max_retries=3):
    """
    Acquires an OTR image with background subtraction. 
    Takes the background once, then loops the beam acquisition if the 
    total integrated intensity is below min_total_intensity.
    """
    print(f"\nAcquiring Background and Beam for OTR{otr_id_str}")    
    pv_in_name = f'mOTR{otr_id_str}:Target:WRITE:IN'
    pv_out_name = f'mOTR{otr_id_str}:Target:WRITE:OUT'
    pv_img_data_name = f'mOTR{otr_id_str}:IMAGE:ArrayData'
    pv_acquire_name = f'mOTR{otr_id_str}:CAMERA:Acquire'
    
    otr_in_pv = PV(pv_in_name)
    otr_out_pv = PV(pv_out_name)
    image_data_pv = PV(pv_img_data_name)
    image_acquire_pv = PV(pv_acquire_name)
    
    print(" -> Retracting screen for background...")
    otr_out_pv.put(1)
    time.sleep(5) 
    
    print(" -> Capturing BACKGROUND image...")
    image_acquire_pv.put(1) 
    time.sleep(3) 
    bg_data = image_data_pv.get() 
    image_acquire_pv.put(0) 
    bg_img = bg_data.reshape(960, 1280).astype(np.float64)
    
    print(" -> Inserting screen for beam...")
    otr_in_pv.put(1)
    time.sleep(5) 
    
    subtracted_img = None
    for attempt in range(1, max_retries + 1):
        image_acquire_pv.put(1) 
        time.sleep(3) 
        beam_data = image_data_pv.get() 
        image_acquire_pv.put(0) 
        beam_img = beam_data.reshape(960, 1280).astype(np.float64)

        subtracted_img = beam_img - bg_img
        
        subtracted_img[subtracted_img < 0] = 0
        
        tot_intensity = np.sum(subtracted_img)
        
        if tot_intensity >= min_total_intensity:
            break 
        else:
            print(f" -> WARNING: Total intensity too low ({tot_intensity:,.0f} < threshold of {min_total_intensity:,.0f}).")
            if attempt < max_retries:
                print(" -> Retaking beam image (keeping screen inserted)...\n")
            else:
                print(" -> Max retries reached. Proceeding with the empty/low-intensity frame.")

    print(" -> Retracting screen...")
    otr_out_pv.put(1) 
    
    return subtracted_img

if __name__ == '__main__':
    img_OTR0_reshaped = acquire_otr_image(0, min_total_intensity=130000, max_retries=3)
    img_OTR1_reshaped = acquire_otr_image(1, min_total_intensity=130000, max_retries=3)
    img_OTR2_reshaped = acquire_otr_image(2, min_total_intensity=130000, max_retries=3)
    img_OTR3_reshaped = acquire_otr_image(3, min_total_intensity=130000, max_retries=3)    
    h_fact_0, v_fact_0 = get_pixel_calibrations(0)
    h_fact_1, v_fact_1 = get_pixel_calibrations(1)
    h_fact_2, v_fact_2 = get_pixel_calibrations(2)
    h_fact_3, v_fact_3 = get_pixel_calibrations(3)    

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(8, 6), constrained_layout=True)

    projX0, projY0, sigH0, sigV0, sig13_0 = plot_otr_analysis_inset(ax1, img_OTR0_reshaped, 'OTR0', h_fact_0, v_fact_0)
    projX1, projY1, sigH1, sigV1, sig13_1 = plot_otr_analysis_inset(ax2, img_OTR1_reshaped, 'OTR1', h_fact_1, v_fact_1)
    projX2, projY2, sigH2, sigV2, sig13_2 = plot_otr_analysis_inset(ax3, img_OTR2_reshaped, 'OTR2', h_fact_2, v_fact_2)
    projX3, projY3, sigH3, sigV3, sig13_3 = plot_otr_analysis_inset(ax4, img_OTR3_reshaped, 'OTR3', h_fact_3, v_fact_3) 

    plt.show() 

    SAVE_FOLDER = 'Data_mOTR'
    FILE_NAME = f'OTR_Beam_Data_{time.strftime("%Y%m%d-%H%M%S")}.npz'
    SAVE_PATH = os.path.join(SAVE_FOLDER, FILE_NAME)

    os.makedirs(SAVE_FOLDER, exist_ok=True)
    
    np.savez_compressed(SAVE_PATH, 
        timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
        OTR0_Image_Raw=img_OTR0_reshaped, OTR0_ProjX=projX0, OTR0_ProjY=projY0,
        OTR0_SigmaH_um=sigH0, OTR0_SigmaV_um=sigV0, OTR0_Sigma13_m2=sig13_0,
        
        OTR1_Image_Raw=img_OTR1_reshaped, OTR1_ProjX=projX1, OTR1_ProjY=projY1,
        OTR1_SigmaH_um=sigH1, OTR1_SigmaV_um=sigV1, OTR1_Sigma13_m2=sig13_1,
        
        OTR2_Image_Raw=img_OTR2_reshaped, OTR2_ProjX=projX2, OTR2_ProjY=projY2,
        OTR2_SigmaH_um=sigH2, OTR2_SigmaV_um=sigV2, OTR2_Sigma13_m2=sig13_2,
        
        OTR3_Image_Raw=img_OTR3_reshaped, OTR3_ProjX=projX3, OTR3_ProjY=projY3,
        OTR3_SigmaH_um=sigH3, OTR3_SigmaV_um=sigV3, OTR3_Sigma13_m2=sig13_3
    )
    print(f"\nSuccessfully saved data to: {SAVE_PATH}")

