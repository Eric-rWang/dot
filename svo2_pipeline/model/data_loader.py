"""
data_loader.py
==============
Load MCmatlab simulation outputs and compute the derived quantities needed
for the inverse problem.
"""

import os
import numpy as np
from scipy.io import loadmat
import sys

sys.path.insert(0, './Simulation Data')

def folder_name(wavelength_nm, x_mm):
    x_str = f"{x_mm:.1f}"
    x_str = x_str.replace('-', 'neg').replace('.', 'p')
    return f"MC_Sim_{int(wavelength_nm)}nm_X_{x_str}mm"

def build_vessel_mask(nx, ny, nz, Lx, Ly, Lz,
                      vessel_radius_cm, vessel_depth_top_cm,
                      vessel_x_center_cm=0.0):
    dx = Lx / nx
    dz = Lz / nz
    x = np.linspace(-Lx/2 + dx/2, Lx/2 - dx/2, nx)
    z = np.linspace(dz/2, Lz - dz/2, nz)

    vessel_z_center = vessel_depth_top_cm + vessel_radius_cm

    X, Z = np.meshgrid(x, z, indexing='ij')
    cross_section = (X - vessel_x_center_cm)**2 + (Z - vessel_z_center)**2 <= vessel_radius_cm**2

    vessel_mask = np.repeat(cross_section[:, np.newaxis, :], ny, axis=1)
    muscle_mask = ~vessel_mask

    return vessel_mask, muscle_mask

def load_intensity(sim_folder, run_index=0):
    filepath = os.path.join(sim_folder, f"detectedIntensity_{run_index}.mat")
    data = loadmat(filepath)
    val = data['collectedIntensity']
    return float(np.squeeze(val))

def load_fluence(sim_folder, run_index=0):
    filepath = os.path.join(sim_folder, f"detectedFluence_{run_index}.mat")
    data = loadmat(filepath)
    return np.asarray(data['detectedFluence'], dtype=np.float64)

def compute_pathlengths(fluence, intensity, vessel_mask, muscle_mask, voxel_volume):
    if intensity <= 0:
        return 0.0, 0.0

    Phi_muscle = np.sum(fluence[muscle_mask]) * voxel_volume
    Phi_vein = np.sum(fluence[vessel_mask]) * voxel_volume

    L_muscle = Phi_muscle / intensity
    L_vein = Phi_vein / intensity

    return L_muscle, L_vein

def load_dataset(base_dir, wavelengths, x_positions_mm,
                 nx=100, ny=100, nz=100,
                 Lx=2.5, Ly=2.5, Lz=2.5,
                 vessel_radius_cm=0.05,
                 vessel_depth_top_cm=0.30,
                 vessel_x_center_cm=0.0,
                 run_index=0,
                 compute_fluence=True):
    
    n_wl = len(wavelengths)
    n_pos = len(x_positions_mm)
    voxel_vol = (Lx / nx) * (Ly / ny) * (Lz / nz)

    if compute_fluence:
        vessel_mask, muscle_mask = build_vessel_mask(
            nx, ny, nz, Lx, Ly, Lz,
            vessel_radius_cm, vessel_depth_top_cm, vessel_x_center_cm
        )
        n_vessel_vox = np.sum(vessel_mask)
        print(f"  Geometry: vessel radius = {vessel_radius_cm} cm, "
              f"{n_vessel_vox} vessel voxels "
              f"(volume = {n_vessel_vox * voxel_vol:.6f} cm³)")

    I = np.zeros((n_wl, n_pos))
    L_muscle_arr = np.zeros((n_wl, n_pos)) if compute_fluence else None
    L_vein_arr = np.zeros((n_wl, n_pos)) if compute_fluence else None

    for k, wl in enumerate(wavelengths):
        for i, x_mm in enumerate(x_positions_mm):
            fname = folder_name(int(wl), x_mm)
            sim_dir = os.path.join(base_dir, fname)

            if not os.path.isdir(sim_dir):
                continue

            I[k, i] = load_intensity(sim_dir, run_index)

            if compute_fluence:
                fluence = load_fluence(sim_dir, run_index)
                Lm, Lv = compute_pathlengths(
                    fluence, I[k, i], vessel_mask, muscle_mask, voxel_vol
                )
                L_muscle_arr[k, i] = Lm
                L_vein_arr[k, i] = Lv

        if compute_fluence:
            idx_max = np.argmax(L_vein_arr[k, :])
            print(f"  {int(wl)} nm: max <L_v> = {L_vein_arr[k, idx_max]:.6f} cm "
                  f"at x = {x_positions_mm[idx_max]:.1f} mm, "
                  f"mean <L_m> = {np.mean(L_muscle_arr[k, :]):.4f} cm")
        else:
            print(f"  {int(wl)} nm: intensities loaded "
                  f"(range: {I[k,:].min():.4e} to {I[k,:].max():.4e})")

    return {
        'I': I,
        'L_muscle': L_muscle_arr,
        'L_vein': L_vein_arr,
    }

def load_bulk_muscle(base_dir, wavelengths,
                     nx=100, ny=100, nz=100,
                     Lx=2.5, Ly=2.5, Lz=2.5,
                     run_index=0):
    
    n_wl = len(wavelengths)
    voxel_vol = (Lx / nx) * (Ly / ny) * (Lz / nz)

    I_bulk = np.zeros(n_wl)
    L_bulk = np.zeros(n_wl)

    for k, wl in enumerate(wavelengths):
        candidates = [folder_name(int(wl), 0.0)]

        sim_dir = None
        for cand in candidates:
            path = os.path.join(base_dir, cand)
            if os.path.isdir(path):
                sim_dir = path
                break

        if sim_dir is None:
            all_dirs = [d for d in os.listdir(base_dir)
                        if d.startswith(f"MC_Sim_{int(wl)}nm") and
                        os.path.isdir(os.path.join(base_dir, d))]
            if all_dirs:
                sim_dir = os.path.join(base_dir, all_dirs[0])
            else:
                continue

        I_bulk[k] = load_intensity(sim_dir, run_index)
        fluence = load_fluence(sim_dir, run_index)
        Phi_total = np.sum(fluence) * voxel_vol
        L_bulk[k] = Phi_total / I_bulk[k]

    return {
        'I_bulk': I_bulk,
        'L_bulk_muscle': L_bulk,
    }

def compute_modulation_from_data(base_data, pert_data):
    """Computes measured modulation M_resp = ln(I_perturbed / I_baseline)"""
    with np.errstate(divide='ignore', invalid='ignore'):
        M_resp = np.log(pert_data['I'] / base_data['I'])
    # Clean up NaNs from missing data
    M_resp[np.isnan(M_resp)] = 0.0
    M_resp[np.isinf(M_resp)] = 0.0
    return M_resp