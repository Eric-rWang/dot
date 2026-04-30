"""
run_pipeline.py
===============
Main script for the SvO2 recovery pipeline.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(0, 'svo2_pipeline')

import optical_properties as op
from model.data_loader import (load_dataset, load_bulk_muscle, 
                         compute_modulation_from_data)
from forward_model import solve_svo2_corrected

# =========================================================================
# CONFIGURATION
# =========================================================================
 
BULK_DIR = "./Simulation Data/Bulk Muscle"
BASELINE_DIR = "./Simulation Data/Muscle with single vessel S75 r0p05mm"
PERTURBED_DIR = "./Simulation Data/Perturbed Vessel S75 r0p0575mm"
OUTPUT_DIR = "./reports/results/pipeline_results"
 
WAVELENGTHS = [660, 850, 940]
X_POSITIONS_MM = [-4.8, -2.4, 0.0, 2.4, 4.8]
 
NX, NY, NZ = 100, 100, 100
LX, LY, LZ = 2.5, 2.5, 2.5
 
VESSEL_RADIUS_BASELINE = 0.05
VESSEL_RADIUS_PERTURBED = 0.0575
VESSEL_DEPTH_TOP = 0.30
VESSEL_X_CENTER = 0.0
 
SVO2_TRUE = 0.75
DELTA_R_TRUE = (VESSEL_RADIUS_PERTURBED - VESSEL_RADIUS_BASELINE) / VESSEL_RADIUS_BASELINE

# =========================================================================
# PIPELINE STEPS
# =========================================================================
 
def step1_load_data():
    print("\n" + "=" * 70)
    print("  STEP 1: LOAD SIMULATION DATA")
    print("=" * 70)
 
    print("\n  Loading Dataset A (bulk muscle)...")
    bulk = load_bulk_muscle(
        BULK_DIR, WAVELENGTHS, NX, NY, NZ, LX, LY, LZ)
 
    print("\n  Loading Dataset B (vessel, r = 0.05 cm)...")
    baseline = load_dataset(
        BASELINE_DIR, WAVELENGTHS, X_POSITIONS_MM,
        NX, NY, NZ, LX, LY, LZ,
        VESSEL_RADIUS_BASELINE, VESSEL_DEPTH_TOP, VESSEL_X_CENTER,
        compute_fluence=True)
 
    print("\n  Loading Dataset C (vessel, r = 0.0575 cm)...")
    # For the perturbed dataset, we only need intensities, not fluence. 
    # This avoids the noisy dL/dr calculations entirely.
    perturbed = load_dataset(
        PERTURBED_DIR, WAVELENGTHS, X_POSITIONS_MM,
        NX, NY, NZ, LX, LY, LZ,
        VESSEL_RADIUS_PERTURBED, VESSEL_DEPTH_TOP, VESSEL_X_CENTER,
        compute_fluence=False) 
 
    return bulk, baseline, perturbed

def step4_modulation(baseline, perturbed):
    print("\n" + "=" * 70)
    print("  STEP 4: RESPIRATORY-AC MODULATION RATIOS")
    print("=" * 70)
 
    wl = np.array(WAVELENGTHS)
    x = np.array(X_POSITIONS_MM)
 
    M_resp = compute_modulation_from_data(baseline, perturbed)
 
    print(f"\n  Modulation Log-Ratio (NaN = missing data):")
    print(f"  {'x [mm]':>10s}", end="")
    for xi in x:
        print(f"{xi:>8.1f}", end="")
    print()
    for k in range(len(wl)):
        print(f"  {wl[k]:>8.0f} nm", end="")
        for i in range(len(x)):
            print(f"{M_resp[k,i]:>8.5f}", end="")
        print()
 
    return M_resp

def step5_solve_pipeline(baseline, perturbed):
    print("\n" + "=" * 70)
    print("  STEP 5: ANALYTICAL INVERSE SOLVER")
    print("=" * 70)
 
    wl = np.array(WAVELENGTHS, dtype=float)
 
    # Get Measured Log Modulation
    M_resp_measured = compute_modulation_from_data(baseline, perturbed)
 
    # Solve using the clean baseline vessel pathlengths
    result_joint = solve_svo2_corrected(M_resp_measured, baseline['L_vein'], wl)
    
    return result_joint

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
 
    print("\n" + "=" * 70)
    print("  SvO2 RECOVERY PIPELINE (Analytical)")
    print("  True SvO2 = 0.75, vessel r = 0.05 cm")
    print("=" * 70)
 
    bulk, baseline, perturbed = step1_load_data()
    M_resp = step4_modulation(baseline, perturbed)
    result_joint = step5_solve_pipeline(baseline, perturbed)
 
    print("\n" + "=" * 70)
    print("  PIPELINE COMPLETE")
    print("=" * 70)
    print(f"  AC-Analytical SvO2: {result_joint['SvO2']:.4f} "
          f"(error: {(result_joint['SvO2']-SVO2_TRUE)*100:+.2f} pp)")
    print(f"  True SvO2:          {SVO2_TRUE:.4f}")
    print(f"  Recovered dr/r:     {result_joint['dr_over_r']:.4f}")
    print(f"  True dr/r:          {DELTA_R_TRUE:.4f}")
 
if __name__ == '__main__':
    main()