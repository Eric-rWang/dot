"""
data_loader.py
==============

Load MCmatlab simulation outputs and compute the derived quantities needed
for the inverse problem.

What this module does
---------------------
1. Reads the .mat files saved by the MCmatlab simulation script.
2. Reconstructs the tissue geometry mask (which voxels are muscle,
   which are vessel) from the known vessel parameters.
3. Computes the compartment-integrated fluence and mean pathlength
   for each tissue compartment at each (wavelength, PD position) channel.

What comes out
--------------
Three key arrays, each of shape (n_wavelengths, n_positions):

    I_collected[k, i]  — the detected intensity at wavelength k, position i
    L_muscle[k, i]     — mean pathlength in muscle for detected photons [cm]
    L_vein[k, i]       — mean pathlength in the vessel for detected photons [cm]

These are the inputs to the forward model and the inverse solver.

File structure expected
-----------------------
The MCmatlab script saves data in folders named:
    MC_Sim_{wavelength}nm_X_{x_str}mm/
        detectedIntensity_0.mat   — contains variable 'collectedIntensity'
        detectedFluence_0.mat     — contains variable 'detectedFluence'
"""

import os
import numpy as np
from scipy.io import loadmat


# ─────────────────────────────────────────────────────────────────────────────
# Folder Naming Convention
# ─────────────────────────────────────────────────────────────────────────────

def folder_name(wavelength_nm, x_mm):
    """
    Generate the folder name that MCmatlab created for a given simulation.

    The MCmatlab script replaces '-' with 'neg' and '.' with 'p' in the
    position string to create filesystem-safe folder names.

    Parameters
    ----------
    wavelength_nm : int
        Wavelength in nm (e.g., 660).
    x_mm : float
        PD position along the array in mm (e.g., -2.4).

    Returns
    -------
    str
        Folder name, e.g. 'MC_Sim_660nm_X_neg2p4mm'.

    Examples
    --------
    >>> folder_name(660, -2.4)
    'MC_Sim_660nm_X_neg2p4mm'
    >>> folder_name(850, 4.8)
    'MC_Sim_850nm_X_4p8mm'
    >>> folder_name(940, 0.0)
    'MC_Sim_940nm_X_0p0mm'
    """
    x_str = f"{x_mm:.1f}"
    x_str = x_str.replace('-', 'neg').replace('.', 'p')
    return f"MC_Sim_{int(wavelength_nm)}nm_X_{x_str}mm"


# ─────────────────────────────────────────────────────────────────────────────
# Geometry Mask Construction
# ─────────────────────────────────────────────────────────────────────────────

def build_vessel_mask(nx, ny, nz, Lx, Ly, Lz,
                      vessel_radius_cm, vessel_depth_top_cm,
                      vessel_x_center_cm=0.0):
    """
    Construct a boolean mask identifying which voxels are inside the vessel.

    This reproduces the geometryDefinition function from the MCmatlab script.
    The vessel is an infinite cylinder running along the y-axis.

    Parameters
    ----------
    nx, ny, nz : int
        Number of voxels in each dimension.
    Lx, Ly, Lz : float
        Size of the simulation cuboid in cm.
    vessel_radius_cm : float
        Vessel radius in cm (e.g., 0.05 for 1 mm diameter).
    vessel_depth_top_cm : float
        Depth from the tissue surface to the TOP of the vessel in cm.
        The vessel center is at depth = vessel_depth_top + vessel_radius.
    vessel_x_center_cm : float
        Lateral position of vessel center in cm.

    Returns
    -------
    vessel_mask : ndarray of bool, shape (nx, ny, nz)
        True for voxels inside the vessel.
    muscle_mask : ndarray of bool, shape (nx, ny, nz)
        True for voxels outside the vessel (muscle).

    Notes
    -----
    MCmatlab uses a coordinate system where:
        x: [-Lx/2, +Lx/2], centered on the cuboid
        y: [-Ly/2, +Ly/2], centered on the cuboid
        z: [0, Lz], with z=0 at the tissue surface
    """
    # Voxel center coordinates
    dx = Lx / nx
    dz = Lz / nz
    x = np.linspace(-Lx/2 + dx/2, Lx/2 - dx/2, nx)
    z = np.linspace(dz/2, Lz - dz/2, nz)

    # Vessel center depth
    vessel_z_center = vessel_depth_top_cm + vessel_radius_cm

    # 2D cross-section mask (the cylinder is uniform along y)
    # Shape: (nx, nz)
    X, Z = np.meshgrid(x, z, indexing='ij')
    cross_section = (X - vessel_x_center_cm)**2 + (Z - vessel_z_center)**2 \
                    <= vessel_radius_cm**2

    # Expand to 3D by repeating along y
    # Shape: (nx, ny, nz)
    vessel_mask = np.repeat(cross_section[:, np.newaxis, :], ny, axis=1)
    muscle_mask = ~vessel_mask

    return vessel_mask, muscle_mask


# ─────────────────────────────────────────────────────────────────────────────
# Loading Individual Simulation Outputs
# ─────────────────────────────────────────────────────────────────────────────

def load_intensity(sim_folder, run_index=0):
    """
    Load the collected intensity from one MC simulation run.

    Parameters
    ----------
    sim_folder : str
        Path to the simulation output folder.
    run_index : int
        Which run to load (default 0).

    Returns
    -------
    float
        Collected intensity, normalized per launched photon.
        This is the total power reaching the detector divided by the
        number of photons launched.
    """
    filepath = os.path.join(sim_folder, f"detectedIntensity_{run_index}.mat")
    data = loadmat(filepath)
    val = data['collectedIntensity']
    # The value may be wrapped in a 1x1 array
    return float(np.squeeze(val))


def load_fluence(sim_folder, run_index=0):
    """
    Load the 3D detected fluence rate from one MC simulation run.

    Parameters
    ----------
    sim_folder : str
        Path to the simulation output folder.
    run_index : int
        Which run to load (default 0).

    Returns
    -------
    ndarray, shape (nx, ny, nz)
        Fluence rate at each voxel [cm^-2 per launched photon],
        restricted to photons that reached the detector
        (because onlyCollected = true in the simulation).

    Notes
    -----
    The fluence at voxel r, multiplied by the voxel volume, represents
    the total path length spent by detected photons in that voxel.
    This is the fundamental quantity from which pathlengths are derived.
    """
    filepath = os.path.join(sim_folder, f"detectedFluence_{run_index}.mat")
    data = loadmat(filepath)
    return np.asarray(data['detectedFluence'], dtype=np.float64)


# ─────────────────────────────────────────────────────────────────────────────
# Computing Compartment Mean Pathlengths
# ─────────────────────────────────────────────────────────────────────────────

def compute_pathlengths(fluence, intensity, vessel_mask, muscle_mask, voxel_volume):
    """
    Compute the mean photon pathlength in each tissue compartment.

    This is the core calculation that connects the MC simulation to the
    modified Beer-Lambert law. The mean pathlength in compartment j is:

        <L_j> = (sum of fluence in compartment j * voxel_volume) / I_detected

    Physically, this is the average distance that a detected photon
    traveled through compartment j.

    Parameters
    ----------
    fluence : ndarray, shape (nx, ny, nz)
        Detected fluence rate from MC (collected-only mode).
    intensity : float
        Collected intensity from the same MC run.
    vessel_mask : ndarray of bool, shape (nx, ny, nz)
        True for vessel voxels.
    muscle_mask : ndarray of bool, shape (nx, ny, nz)
        True for muscle voxels.
    voxel_volume : float
        Volume of a single voxel in cm^3.

    Returns
    -------
    L_muscle : float
        Mean pathlength in muscle [cm].
    L_vein : float
        Mean pathlength in the vessel [cm].

    Notes
    -----
    - L_muscle is typically 0.3-1.0 cm (detected photons spend most of
      their path in the bulk tissue).
    - L_vein is typically 0.001-0.1 cm (the vessel is small, so only a
      small fraction of the path goes through it).
    - The ratio L_vein / L_muscle tells you what fraction of the
      measurement's sensitivity is allocated to the blood compartment.
    - These pathlengths are evaluated at the absorption coefficients used
      in the MC simulation (the reference state). They are NOT the same
      as DPF * rho.
    """
    if intensity <= 0:
        return 0.0, 0.0

    # Compartment-integrated fluence [cm]
    Phi_muscle = np.sum(fluence[muscle_mask]) * voxel_volume
    Phi_vein = np.sum(fluence[vessel_mask]) * voxel_volume

    # Mean pathlengths [cm]
    L_muscle = Phi_muscle / intensity
    L_vein = Phi_vein / intensity

    return L_muscle, L_vein


# ─────────────────────────────────────────────────────────────────────────────
# Loading an Entire Dataset
# ─────────────────────────────────────────────────────────────────────────────

def load_dataset(base_dir, wavelengths, x_positions_mm,
                 nx=100, ny=100, nz=100,
                 Lx=2.5, Ly=2.5, Lz=2.5,
                 vessel_radius_cm=0.05,
                 vessel_depth_top_cm=0.30,
                 vessel_x_center_cm=0.0,
                 run_index=0,
                 compute_fluence=True):
    """
    Load all simulation data for one dataset and compute pathlengths.

    This is the main entry point for loading data. It iterates over all
    wavelengths and PD positions, loads the intensity and fluence from
    each simulation folder, and computes the compartment pathlengths.

    Parameters
    ----------
    base_dir : str
        Directory containing all simulation output folders for this dataset.
    wavelengths : list of int
        Wavelengths in nm, e.g. [660, 850, 940].
    x_positions_mm : list of float
        PD positions along the array in mm.
    nx, ny, nz : int
        Voxel grid dimensions.
    Lx, Ly, Lz : float
        Cuboid dimensions in cm.
    vessel_radius_cm : float
        Vessel radius used in this simulation.
    vessel_depth_top_cm : float
        Depth to top of vessel in cm.
    vessel_x_center_cm : float
        Lateral position of vessel center in cm.
    run_index : int
        Simulation run index.
    compute_fluence : bool
        If True, load fluence maps and compute pathlengths.
        If False, only load intensities (faster, for perturbation data).

    Returns
    -------
    dict with keys:
        'I'        : ndarray (n_wl, n_pos) — collected intensities
        'L_muscle' : ndarray (n_wl, n_pos) — muscle pathlengths [cm]
                     (None if compute_fluence=False)
        'L_vein'   : ndarray (n_wl, n_pos) — vessel pathlengths [cm]
                     (None if compute_fluence=False)
    """
    n_wl = len(wavelengths)
    n_pos = len(x_positions_mm)
    voxel_vol = (Lx / nx) * (Ly / ny) * (Lz / nz)

    # Build geometry masks (same for all runs in this dataset)
    if compute_fluence:
        vessel_mask, muscle_mask = build_vessel_mask(
            nx, ny, nz, Lx, Ly, Lz,
            vessel_radius_cm, vessel_depth_top_cm, vessel_x_center_cm
        )
        n_vessel_vox = np.sum(vessel_mask)
        print(f"  Geometry: vessel radius = {vessel_radius_cm} cm, "
              f"{n_vessel_vox} vessel voxels "
              f"(volume = {n_vessel_vox * voxel_vol:.6f} cm³)")

    # Allocate output arrays
    I = np.zeros((n_wl, n_pos))
    L_muscle_arr = np.zeros((n_wl, n_pos)) if compute_fluence else None
    L_vein_arr = np.zeros((n_wl, n_pos)) if compute_fluence else None

    # Load each simulation
    for k, wl in enumerate(wavelengths):
        for i, x_mm in enumerate(x_positions_mm):
            fname = folder_name(int(wl), x_mm)
            sim_dir = os.path.join(base_dir, fname)

            if not os.path.isdir(sim_dir):
                print(f"  WARNING: folder not found: {sim_dir}")
                continue

            # Load intensity
            I[k, i] = load_intensity(sim_dir, run_index)

            # Load fluence and compute pathlengths
            if compute_fluence:
                fluence = load_fluence(sim_dir, run_index)
                Lm, Lv = compute_pathlengths(
                    fluence, I[k, i], vessel_mask, muscle_mask, voxel_vol
                )
                L_muscle_arr[k, i] = Lm
                L_vein_arr[k, i] = Lv

        # Print summary for this wavelength
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
    """
    Load bulk muscle simulation data (homogeneous, no vessel).

    For bulk muscle, only one position per wavelength is needed since
    the medium is homogeneous and all positions give the same result.

    Parameters
    ----------
    base_dir : str
        Directory containing bulk muscle simulation folders.
    wavelengths : list of int
        Wavelengths in nm.
    nx, ny, nz : int
        Grid dimensions.
    Lx, Ly, Lz : float
        Cuboid dimensions in cm.
    run_index : int
        Simulation run index.

    Returns
    -------
    dict with keys:
        'I_bulk'        : ndarray (n_wl,) — detected intensities
        'L_bulk_muscle' : ndarray (n_wl,) — total mean pathlength [cm]
    """
    n_wl = len(wavelengths)
    voxel_vol = (Lx / nx) * (Ly / ny) * (Lz / nz)

    I_bulk = np.zeros(n_wl)
    L_bulk = np.zeros(n_wl)

    for k, wl in enumerate(wavelengths):
        # Try to find a folder for this wavelength
        # The bulk muscle sim may have been run at x=0 or another position
        candidates = [folder_name(int(wl), 0.0)]

        sim_dir = None
        for cand in candidates:
            path = os.path.join(base_dir, cand)
            if os.path.isdir(path):
                sim_dir = path
                break

        if sim_dir is None:
            # Search for any folder matching this wavelength
            all_dirs = [d for d in os.listdir(base_dir)
                        if d.startswith(f"MC_Sim_{int(wl)}nm") and
                        os.path.isdir(os.path.join(base_dir, d))]
            if all_dirs:
                sim_dir = os.path.join(base_dir, all_dirs[0])
                print(f"  Using {all_dirs[0]} for {wl} nm bulk muscle")
            else:
                raise FileNotFoundError(
                    f"No bulk muscle folder found for {wl} nm in {base_dir}")

        # Load intensity
        I_bulk[k] = load_intensity(sim_dir, run_index)

        # Load fluence — for bulk muscle, all voxels are muscle
        fluence = load_fluence(sim_dir, run_index)
        Phi_total = np.sum(fluence) * voxel_vol
        L_bulk[k] = Phi_total / I_bulk[k]

        print(f"  {int(wl)} nm: I_bulk = {I_bulk[k]:.6e}, "
              f"<L> = {L_bulk[k]:.4f} cm")

    return {
        'I_bulk': I_bulk,
        'L_bulk_muscle': L_bulk,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Self-Test
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("Testing geometry mask construction...")
    vessel_mask, muscle_mask = build_vessel_mask(
        100, 100, 100, 2.5, 2.5, 2.5,
        vessel_radius_cm=0.05,
        vessel_depth_top_cm=0.30,
    )
    voxvol = (2.5/100)**3
    print(f"  Muscle voxels: {np.sum(muscle_mask)}")
    print(f"  Vessel voxels: {np.sum(vessel_mask)}")
    print(f"  Vessel volume: {np.sum(vessel_mask) * voxvol:.6f} cm³")
    print(f"  Expected (pi*r²*L): {np.pi * 0.05**2 * 2.5:.6f} cm³")

    print("\nTesting folder name generation...")
    tests = [(660, -2.4), (850, 0.0), (940, 9.6)]
    for wl, x in tests:
        print(f"  ({wl}, {x}) -> {folder_name(wl, x)}")
