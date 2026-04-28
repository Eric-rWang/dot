import numpy as np
import scipy.io as sio
import os
import re
import pandas as pd
from scipy.optimize import least_squares

class BloodSpectra:
    """Handles blood absorption coefficients based on Jacques (2013)."""
    def __init__(self):
        # Molar extinction coefficients (cm^-1 M^-1) - Prahl/Jacques standard
        # wavelengths: 540, 630, 660, 850, 940
        self.wvs = np.array([540, 630, 660, 850, 940])
        self.eps_hbO2 = np.array([39564, 520, 320, 1058, 1214])
        self.eps_hbR = np.array([49512, 3340, 3226, 692, 693])
        self.mua_water = np.array([0.0004, 0.002, 0.004, 0.04, 0.25])
        
    def get_mua_blood(self, saturation, hb_conc_g_dl=15.0):
        """Calculates blood mua (cm^-1)."""
        hb_molar = (hb_conc_g_dl / 10.0) / 64500.0 * 1000.0
        mua_hb = np.log(10) * hb_molar * (saturation * self.eps_hbO2 + (1 - saturation) * self.eps_hbR)
        return mua_hb + (self.mua_water * 0.95)

def get_geometry_masks(nx=100, ny=100, nz=100, Lx=2.5, Ly=2.5, Lz=2.5):
    """
    Mathematical replica of the MATLAB geometryDefinition.
   
    """
    # Create the coordinate grids (X, Y, Z matrices)
    x = np.linspace(-Lx/2, Lx/2, nx)
    y = np.linspace(-Ly/2, Ly/2, ny)
    z = np.linspace(0, Lz, nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Vessel Geometry (in cm)
    v_radius = 0.05          # 1 mm diameter / 2
    v_depth_top = 0.30       # 3 mm below the surface
    v_z_center = v_depth_top + v_radius # 0.35 cm
    v_x_center = 0.0         # Centered at X = 0

    # Logical mask for a cylinder running along the Y-axis
    # (X - Xc)^2 + (Z - Zc)^2 <= R^2
    v_mask = (X - v_x_center)**2 + (Z - v_z_center)**2 <= v_radius**2
    m_mask = ~v_mask # Everywhere else is muscle
    
    return m_mask, v_mask

def parse_x_location(folder_name):
    """Converts naming convention 'neg9p6mm' to -0.96 cm."""
    match = re.search(r'X_(.*?)mm', folder_name)
    if not match: return None
    val_str = match.group(1)
    numeric_str = val_str.replace('neg', '-').replace('p', '.')
    return float(numeric_str) / 10.0 

def extract_and_export(root_path, output_dir="reports/results", wavelengths=[660, 850, 940]):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 1. Initialize Geometry Masks once
    m_mask, v_mask = get_geometry_masks()
    voxel_vol = (2.5 / 100)**3
    
    all_data = []
    
    for wv in wavelengths:
        wv_dir = os.path.join(root_path, str(wv))
        if not os.path.exists(wv_dir):
            print(f"Skipping {wv}nm - directory not found.")
            continue
            
        sim_folders = [f for f in os.listdir(wv_dir) if f.startswith(f'MC_Sim_{wv}nm')]
        
        for folder in sim_folders:
            x_cm = parse_x_location(folder)
            folder_path = os.path.join(wv_dir, folder)
            
            try:
                # 2. Load only the relevant numerical matrices (No Class Objects)
                i_det = sio.loadmat(os.path.join(folder_path, "detectedIntensity_0.mat"))['collectedIntensity'][0,0]
                phi = sio.loadmat(os.path.join(folder_path, "detectedFluence_0.mat"))['detectedFluence']
            except FileNotFoundError:
                continue
            
            # 3. Calculate Mean Pathlengths (Equation 23)
            L_v = (np.sum(phi[v_mask]) * voxel_vol) / i_det
            L_m = (np.sum(phi[m_mask]) * voxel_vol) / i_det
            
            all_data.append({
                'wavelength_nm': wv,
                'x_location_cm': x_cm,
                'i_detected': i_det,
                'L_vessel_cm': L_v,
                'L_muscle_cm': L_m
            })

    # Export to structured CSV
    df = pd.DataFrame(all_data)
    df = df.sort_values(by=['wavelength_nm', 'x_location_cm'])
    output_path = os.path.join(output_dir, "transverse_sweep_pathlengths.csv")
    df.to_csv(output_path, index=False)
    
    print(f"--- Successfully Exported results to {output_path} ---")
    return df

# --- EXECUTION ---
if __name__ == "__main__":
    # Updated path to your specific simulation batch
    data_root = "/Volumes/Erics Goods/Optical Tests/DOT/Simulation Data/Muscle with single vessel S75 r0p0575mm"
    results_df = extract_and_export(data_root)