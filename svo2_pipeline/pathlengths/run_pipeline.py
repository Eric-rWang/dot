"""
run_pipeline.py
===============

Main script for the SvO2 recovery pipeline.

This script executes the complete process described in the technical
documentation, step by step:

    Step 1: Load simulation data (Datasets A, B, C)
    Step 2: Compute compartment mean pathlengths from Dataset B
    Step 3: Compute vessel contrast (diagnostic check)
    Step 4: Compute respiratory-AC modulation ratios from Datasets B & C
    Step 5: Run the inverse solver (DC + respiratory-AC)
    Step 6: Diagnostic checks (initial guess sensitivity, noise robustness)

USAGE
-----
1. Edit the CONFIGURATION section to set your data folder paths.
2. Run:  python run_pipeline.py

PREREQUISITES
-------------
    pip install numpy scipy matplotlib

DATA REQUIREMENTS
-----------------
Three simulation datasets, each in its own directory:

    BULK_DIR/
        MC_Sim_660nm_X_0p0mm/
            detectedIntensity_0.mat
            detectedFluence_0.mat
        MC_Sim_850nm_X_0p0mm/
            ...
        MC_Sim_940nm_X_0p0mm/
            ...

    BASELINE_DIR/
        MC_Sim_660nm_X_neg9p6mm/
            detectedIntensity_0.mat
            detectedFluence_0.mat
        MC_Sim_660nm_X_neg7p2mm/
            ...
        ... (all 10 positions × 3 wavelengths = 30 folders)

    PERTURBED_DIR/
        (same folder structure as BASELINE_DIR, but simulated with
         vessel_radius = 0.0575 cm instead of 0.05 cm)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for saving plots
import matplotlib.pyplot as plt
import os

# Local modules
import optical_properties as op
from svo2_pipeline.model.data_loader import load_dataset, load_bulk_muscle
from svo2_pipeline.pathlengths.forward_model import solve_svo2, print_results, add_noise


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — EDIT THESE TO MATCH DATA                                 ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Directory containing bulk muscle simulation folders (Dataset A)
BULK_DIR = "./Simulation Data/Bulk Muscle"

# Directory containing baseline vessel simulation folders (Dataset B)
# Vessel radius = 0.05 cm, SvO2 = 0.75
BASELINE_DIR = "./Simulation Data/Muscle with single vessel S75 r0p05mm"

# Directory containing perturbed vessel simulation folders (Dataset C)
# Vessel radius = 0.0575 cm, SvO2 = 0.75
PERTURBED_DIR = "./Simulation Data/Perturbed Vessel S75 r0p0575mm"

# Output directory for plots and results
OUTPUT_DIR = "./reports/results/pipeline_results"

# Wavelengths used in the simulation [nm]
WAVELENGTHS = [660, 850, 940]

# PD positions along the x-axis [mm] — the 10 simulated locations
X_POSITIONS_MM = [-9.6, -7.2, -4.8, -2.4, 0.0, 2.4, 4.8, 7.2, 9.6, 12.0]

# Simulation grid parameters (must match MCmatlab settings)
NX, NY, NZ = 100, 100, 100
LX, LY, LZ = 2.5, 2.5, 2.5  # cm

# Vessel geometry (Dataset B — baseline)
VESSEL_RADIUS_BASELINE = 0.05     # cm
VESSEL_RADIUS_PERTURBED = 0.0575  # cm (15% expansion)
VESSEL_DEPTH_TOP = 0.30           # cm (surface to top of vessel)
VESSEL_X_CENTER = 0.0             # cm

# True values used in the simulation (for validation)
SVO2_TRUE = 0.75
HB_TRUE = 2.3    # mM — check mediaPropertiesFunc
DELTA_R_TRUE = (VESSEL_RADIUS_PERTURBED - VESSEL_RADIUS_BASELINE) / \
                VESSEL_RADIUS_BASELINE  # = 0.15


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 1: Load Simulation Data                                             ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

def step1_load_data():
    """
    Load all three simulation datasets.

    Returns
    -------
    bulk : dict with 'I_bulk', 'L_bulk_muscle'
    baseline : dict with 'I', 'L_muscle', 'L_vein'
    perturbed : dict with 'I'
    """
    print("\n" + "=" * 70)
    print("  STEP 1: LOAD SIMULATION DATA")
    print("=" * 70)

    # --- Dataset A: Bulk muscle ---
    print("\n  Loading Dataset A (bulk muscle, no vessel)...")
    bulk = load_bulk_muscle(
        BULK_DIR, WAVELENGTHS,
        NX, NY, NZ, LX, LY, LZ,
    )

    # --- Dataset B: Baseline vessel (r = 0.05 cm) ---
    print("\n  Loading Dataset B (vessel, r = 0.05 cm, SvO2 = 0.75)...")
    baseline = load_dataset(
        BASELINE_DIR, WAVELENGTHS, X_POSITIONS_MM,
        NX, NY, NZ, LX, LY, LZ,
        vessel_radius_cm=VESSEL_RADIUS_BASELINE,
        vessel_depth_top_cm=VESSEL_DEPTH_TOP,
        vessel_x_center_cm=VESSEL_X_CENTER,
        compute_fluence=True,  # Need fluence for pathlengths
    )

    # --- Dataset C: Perturbed vessel (r = 0.0575 cm) ---
    print("\n  Loading Dataset C (vessel, r = 0.0575 cm, SvO2 = 0.75)...")
    perturbed = load_dataset(
        PERTURBED_DIR, WAVELENGTHS, X_POSITIONS_MM,
        NX, NY, NZ, LX, LY, LZ,
        vessel_radius_cm=VESSEL_RADIUS_PERTURBED,
        vessel_depth_top_cm=VESSEL_DEPTH_TOP,
        vessel_x_center_cm=VESSEL_X_CENTER,
        compute_fluence=False,  # Only need intensities
    )

    return bulk, baseline, perturbed


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 2: Display Pathlengths (already computed in Step 1)                 ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

def step2_pathlengths(baseline):
    """
    Display the compartment pathlengths computed during data loading.

    The pathlengths L_muscle and L_vein were computed in load_dataset()
    by integrating the MC fluence over each compartment. Here we just
    display them.
    """
    print("\n" + "=" * 70)
    print("  STEP 2: COMPARTMENT MEAN PATHLENGTHS")
    print("=" * 70)

    wl = np.array(WAVELENGTHS)
    x = np.array(X_POSITIONS_MM)

    print("\n  Vessel pathlength <L_v> [mm] at each (wavelength, position):")
    print(f"  {'x [mm]':>10s}", end="")
    for xi in x:
        print(f"{xi:>8.1f}", end="")
    print()
    for k in range(len(wl)):
        print(f"  {wl[k]:>8.0f} nm", end="")
        for i in range(len(x)):
            print(f"{baseline['L_vein'][k, i]*10:>8.4f}", end="")
        print()

    print(f"\n  Muscle pathlength <L_m> [cm] (mean across positions):")
    for k in range(len(wl)):
        mean_Lm = np.mean(baseline['L_muscle'][k, :])
        print(f"    {wl[k]:.0f} nm: {mean_Lm:.4f} cm")

    print(f"\n  Pathlength ratio <L_v> / <L_m> at vessel center (x=0):")
    idx_center = np.argmin(np.abs(x))
    for k in range(len(wl)):
        ratio = baseline['L_vein'][k, idx_center] / baseline['L_muscle'][k, idx_center]
        print(f"    {wl[k]:.0f} nm: {ratio:.4f}  "
              f"({ratio*100:.2f}% of path is in vessel)")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 3: Vessel Contrast Analysis                                         ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

def step3_contrast(bulk, baseline):
    """
    Compute and display vessel contrast.

    Contrast = (I_vessel - I_bulk) / I_bulk

    This tells us how much the vessel changes the detected signal.
    """
    print("\n" + "=" * 70)
    print("  STEP 3: VESSEL CONTRAST ANALYSIS")
    print("=" * 70)

    wl = np.array(WAVELENGTHS)
    x = np.array(X_POSITIONS_MM)

    contrast = np.zeros((len(wl), len(x)))
    for k in range(len(wl)):
        contrast[k, :] = (baseline['I'][k, :] - bulk['I_bulk'][k]) / bulk['I_bulk'][k]

    # Print table
    print(f"\n  Vessel contrast [%]:")
    print(f"  {'x [mm]':>10s}", end="")
    for xi in x:
        print(f"{xi:>8.1f}", end="")
    print()
    for k in range(len(wl)):
        print(f"  {wl[k]:>8.0f} nm", end="")
        for i in range(len(x)):
            print(f"{contrast[k, i]*100:>8.2f}", end="")
        print()

    # Assessment
    print(f"\n  Peak contrast (most negative = strongest vessel signal):")
    for k in range(len(wl)):
        peak = np.min(contrast[k, :]) * 100
        idx = np.argmin(contrast[k, :])
        print(f"    {wl[k]:.0f} nm: {peak:.2f}% at x = {x[idx]:.1f} mm")
        if abs(peak) < 0.5:
            print(f"           ⚠ WARNING: contrast < 0.5% — may be below noise floor")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    colors = ['#d62728', '#2ca02c', '#1f77b4']

    ax = axes[0]
    for k in range(len(wl)):
        ax.plot(x, contrast[k, :] * 100, '-o', color=colors[k],
                linewidth=2, markersize=6, label=f'{wl[k]:.0f} nm')
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    ax.axvline(VESSEL_X_CENTER * 10, color='gray', linestyle=':',
               linewidth=1, alpha=0.5)
    ax.set_xlabel('PD Position x [mm]')
    ax.set_ylabel('Vessel Contrast [%]')
    ax.set_title('Vessel-Induced Contrast')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    for k in range(len(wl)):
        ax.plot(x, baseline['L_vein'][k, :] * 10, '-o', color=colors[k],
                linewidth=2, markersize=6, label=f'{wl[k]:.0f} nm')
    ax.axvline(VESSEL_X_CENTER * 10, color='gray', linestyle=':',
               linewidth=1, alpha=0.5)
    ax.set_xlabel('PD Position x [mm]')
    ax.set_ylabel('⟨L_v⟩ [mm]')
    ax.set_title('Vessel Pathlength')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'step3_contrast.png'), dpi=150)
    print(f"\n  Plot saved: {OUTPUT_DIR}/step3_contrast.png")

    return contrast


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 4: Compute Respiratory-AC Modulation Ratios                         ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

def step4_modulation_ratios(baseline, perturbed):
    """
    Compute the respiratory-AC modulation ratio from Datasets B and C.

    M_resp = (I_C - I_B) / I_B

    This is the fractional change in intensity when the vessel expands
    from its baseline radius to the respiratory-perturbed radius.
    """
    print("\n" + "=" * 70)
    print("  STEP 4: RESPIRATORY-AC MODULATION RATIOS")
    print("=" * 70)

    wl = np.array(WAVELENGTHS)
    x = np.array(X_POSITIONS_MM)

    M_resp = (perturbed['I'] - baseline['I']) / baseline['I']

    print(f"\n  Modulation ratio M_resp [%]:")
    print(f"  {'x [mm]':>10s}", end="")
    for xi in x:
        print(f"{xi:>8.1f}", end="")
    print()
    for k in range(len(wl)):
        print(f"  {wl[k]:>8.0f} nm", end="")
        for i in range(len(x)):
            print(f"{M_resp[k, i]*100:>8.3f}", end="")
        print()

    # The ratio of modulation ratios at two wavelengths should encode SvO2
    idx_center = np.argmin(np.abs(x))
    print(f"\n  Modulation ratio at vessel center (x = {x[idx_center]:.1f} mm):")
    for k in range(len(wl)):
        print(f"    {wl[k]:.0f} nm: {M_resp[k, idx_center]*100:.4f}%")

    if M_resp[2, idx_center] != 0:
        ratio_660_940 = M_resp[0, idx_center] / M_resp[2, idx_center]
        print(f"\n  Ratio M_resp(660) / M_resp(940) = {ratio_660_940:.3f}")
        print(f"  (This ratio is a direct function of SvO2)")

    return M_resp


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 5: Run the Inverse Solver                                           ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

def step5_solve(bulk, baseline, perturbed):
    """
    Run the inverse solver to recover SvO2.

    Performs both DC-only and DC + respiratory-AC fits for comparison.
    """
    print("\n" + "=" * 70)
    print("  STEP 5: INVERSE SOLVER — SvO2 RECOVERY")
    print("=" * 70)

    wl = np.array(WAVELENGTHS, dtype=float)

    # ── 5a: DC-only fit ───────────────────────────────────────────────────

    print("\n  --- 5a: DC-only fit ---")
    result_dc = solve_svo2(
        wl, baseline['L_muscle'], baseline['L_vein'],
        I_baseline=baseline['I'],
        I_perturbed=None,  # No AC data
        I_bulk=bulk['I_bulk'],
        L_bulk_muscle=bulk['L_bulk_muscle'],
        SvO2_init=0.60,
        Hb_init=2.0,
        verbose=False,
    )
    print_results(result_dc, wl,
                  true_SvO2=SVO2_TRUE, true_Hb=HB_TRUE)

    # ── 5b: DC + respiratory-AC fit ───────────────────────────────────────

    print("\n  --- 5b: DC + respiratory-AC fit ---")
    result_joint = solve_svo2(
        wl, baseline['L_muscle'], baseline['L_vein'],
        I_baseline=baseline['I'],
        I_perturbed=perturbed['I'],
        I_bulk=bulk['I_bulk'],
        L_bulk_muscle=bulk['L_bulk_muscle'],
        SvO2_init=0.60,
        Hb_init=2.0,
        verbose=False,
    )
    print_results(result_joint, wl,
                  true_SvO2=SVO2_TRUE, true_Hb=HB_TRUE,
                  true_delta_r=DELTA_R_TRUE)

    return result_dc, result_joint


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 6: Diagnostic Checks                                                ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

def step6_diagnostics(bulk, baseline, perturbed):
    """
    Run diagnostic tests on the solver.

    6a: Sensitivity to initial guess (do all starting points converge?)
    6b: Noise robustness (how does accuracy degrade with noise?)
    """
    print("\n" + "=" * 70)
    print("  STEP 6: DIAGNOSTIC CHECKS")
    print("=" * 70)

    wl = np.array(WAVELENGTHS, dtype=float)

    # ── 6a: Initial guess sensitivity ─────────────────────────────────────

    print("\n  --- 6a: Initial guess sensitivity (DC + resp-AC) ---")
    print(f"  {'SvO2_init':>10s} {'SvO2_recovered':>15s} {'Error':>10s} "
          f"{'Converged':>10s}")
    print("  " + "-" * 50)

    for s0 in [0.50, 0.60, 0.70, 0.80, 0.90]:
        r = solve_svo2(
            wl, baseline['L_muscle'], baseline['L_vein'],
            baseline['I'], perturbed['I'],
            bulk['I_bulk'], bulk['L_bulk_muscle'],
            SvO2_init=s0, verbose=False,
        )
        err = (r['SvO2'] - SVO2_TRUE) * 100
        print(f"  {s0:>10.2f} {r['SvO2']:>15.4f} {err:>+10.2f}pp "
              f"{str(r['success']):>10s}")

    # ── 6b: Noise robustness ─────────────────────────────────────────────

    print("\n  --- 6b: Noise robustness (DC + resp-AC) ---")
    print(f"  {'Noise σ':>8s} {'Mean SvO2':>10s} {'Std SvO2':>10s} "
          f"{'Mean err':>10s} {'Mean [Hb]':>10s} {'Mean Δr/r':>10s}")
    print("  " + "-" * 62)

    noise_levels = [0.005, 0.01, 0.02, 0.05]
    n_trials = 20
    all_noise_results = {}

    for sigma in noise_levels:
        svo2_list = []
        hb_list = []
        dr_list = []

        for trial in range(n_trials):
            I_noisy_b = add_noise(baseline['I'], sigma, seed=100 + trial)
            I_noisy_c = add_noise(perturbed['I'], sigma, seed=200 + trial)

            r = solve_svo2(
                wl, baseline['L_muscle'], baseline['L_vein'],
                I_noisy_b, I_noisy_c,
                bulk['I_bulk'], bulk['L_bulk_muscle'],
                SvO2_init=0.70, verbose=False,
            )
            svo2_list.append(r['SvO2'])
            hb_list.append(r['Hb_mM'])
            dr_list.append(r['delta_r_over_r'])

        svo2_arr = np.array(svo2_list)
        hb_arr = np.array(hb_list)
        dr_arr = np.array(dr_list)
        err_arr = svo2_arr - SVO2_TRUE

        print(f"  {sigma*100:>7.1f}% {np.mean(svo2_arr):>10.4f} "
              f"{np.std(svo2_arr):>10.4f} {np.mean(err_arr):>+10.4f} "
              f"{np.mean(hb_arr):>10.3f} {np.mean(dr_arr):>10.4f}")

        all_noise_results[sigma] = svo2_arr

    # Plot noise results
    fig, ax = plt.subplots(figsize=(8, 5))
    data = [all_noise_results[s] for s in noise_levels]
    bp = ax.boxplot(data, positions=range(len(noise_levels)), widths=0.5)
    ax.axhline(SVO2_TRUE, color='r', linestyle='--', linewidth=1.5,
               label=f'True SvO₂ = {SVO2_TRUE}')
    ax.set_xticks(range(len(noise_levels)))
    ax.set_xticklabels([f'{s*100:.1f}%' for s in noise_levels])
    ax.set_xlabel('Noise Level (σ)')
    ax.set_ylabel('Recovered SvO₂')
    ax.set_title('SvO₂ Recovery vs. Noise Level (DC + Resp-AC)')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'step6_noise_robustness.png'), dpi=150)
    print(f"\n  Plot saved: {OUTPUT_DIR}/step6_noise_robustness.png")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  MAIN                                                                     ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

def main():
    """Run the full SvO2 recovery pipeline."""

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("\n" + "╔" + "═" * 68 + "╗")
    print("║  SvO2 RECOVERY PIPELINE — Two-Compartment MBLL                    ║")
    print("║  Single vessel (r = 0.05 cm) in bulk muscle                       ║")
    print("║  True SvO2 = 0.75                                                 ║")
    print("╚" + "═" * 68 + "╝")

    # Step 1: Load data
    bulk, baseline, perturbed = step1_load_data()

    # Step 2: Display pathlengths
    step2_pathlengths(baseline)

    # Step 3: Vessel contrast
    contrast = step3_contrast(bulk, baseline)

    # Step 4: Modulation ratios
    M_resp = step4_modulation_ratios(baseline, perturbed)

    # Step 5: Solve for SvO2
    result_dc, result_joint = step5_solve(bulk, baseline, perturbed)

    # Step 6: Diagnostics
    step6_diagnostics(bulk, baseline, perturbed)

    # ── Summary ───────────────────────────────────────────────────────────

    print("\n" + "╔" + "═" * 68 + "╗")
    print("║  PIPELINE COMPLETE                                                ║")
    print("╚" + "═" * 68 + "╝")
    print(f"\n  DC-only SvO2:          {result_dc['SvO2']:.4f} "
          f"(error: {(result_dc['SvO2'] - SVO2_TRUE)*100:+.2f} pp)")
    print(f"  DC + Resp-AC SvO2:     {result_joint['SvO2']:.4f} "
          f"(error: {(result_joint['SvO2'] - SVO2_TRUE)*100:+.2f} pp)")
    print(f"  True SvO2:             {SVO2_TRUE:.4f}")
    print(f"\n  Outputs saved to: {OUTPUT_DIR}/")

    # Save numerical results
    np.savez(
        os.path.join(OUTPUT_DIR, 'pipeline_results.npz'),
        wavelengths=np.array(WAVELENGTHS),
        x_positions_mm=np.array(X_POSITIONS_MM),
        I_baseline=baseline['I'],
        I_perturbed=perturbed['I'],
        I_bulk=bulk['I_bulk'],
        L_muscle=baseline['L_muscle'],
        L_vein=baseline['L_vein'],
        contrast=contrast,
        M_resp=M_resp,
        SvO2_dc=result_dc['SvO2'],
        SvO2_joint=result_joint['SvO2'],
        Hb_joint=result_joint['Hb_mM'],
        delta_r_joint=result_joint['delta_r_over_r'],
        mua_musc_joint=result_joint['mua_musc'],
    )
    print(f"  Numerical results: {OUTPUT_DIR}/pipeline_results.npz")


if __name__ == '__main__':
    main()
