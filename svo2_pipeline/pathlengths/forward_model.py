"""
forward_model.py
================

Forward model (two-compartment MBLL) and inverse solver for SvO2 recovery.

This module contains two things:

1. The FORWARD MODEL: given tissue parameters (SvO2, [Hb], muscle mua, etc.)
   and the MC-derived pathlengths, predict what the detector should measure.

2. The INVERSE SOLVER: given what the detector actually measured, find the
   tissue parameters that best explain the data.

The forward model has two components:
    - DC: predicts the time-averaged log-intensity
    - Respiratory-AC: predicts the modulation ratio due to venous pulsation

Both are combined into a single optimization problem.
"""

import numpy as np
from scipy.optimize import least_squares
import optical_properties as op


# ─────────────────────────────────────────────────────────────────────────────
# Forward Model: DC Component
# ─────────────────────────────────────────────────────────────────────────────

def predict_dc(SvO2, Hb_mM, mua_musc, lnI0, wavelengths, L_muscle, L_vein):
    """
    Predict the DC (time-averaged) log-intensity at each measurement channel.

    This implements the two-compartment modified Beer-Lambert law:

        ln I(λ, x) = ln I₀(λ) - μₐ_musc(λ) * <L_m>(λ,x) - μₐ_blood(λ) * <L_v>(λ,x)

    Each term represents one source of attenuation:
        - ln I₀: the light that would be detected with zero absorption
        - μₐ_musc * <L_m>: attenuation due to muscle
        - μₐ_blood * <L_v>: attenuation due to the blood vessel

    Parameters
    ----------
    SvO2 : float
        Venous oxygen saturation (0-1). THIS IS WHAT WE WANT TO RECOVER.
    Hb_mM : float
        Hemoglobin concentration in mM.
    mua_musc : ndarray, shape (n_wl,)
        Muscle absorption coefficient at each wavelength [cm^-1].
    lnI0 : ndarray, shape (n_wl,)
        Log of the effective source intensity at each wavelength.
        This absorbs LED power, PD responsivity, and geometric coupling.
    wavelengths : ndarray, shape (n_wl,)
        Wavelengths in nm.
    L_muscle : ndarray, shape (n_wl, n_pos)
        Mean photon pathlength in muscle [cm], from MC simulation.
    L_vein : ndarray, shape (n_wl, n_pos)
        Mean photon pathlength in the vessel [cm], from MC simulation.

    Returns
    -------
    lnI_pred : ndarray, shape (n_wl, n_pos)
        Predicted log-intensity at each channel.
    """
    n_wl, n_pos = L_muscle.shape

    # Blood absorption at each wavelength — depends on SvO2 and [Hb]
    mua_blood = op.mua_blood_spectrum(wavelengths, SvO2, Hb_mM)  # (n_wl,)

    # Predict log-intensity at each channel
    lnI_pred = np.zeros((n_wl, n_pos))
    for k in range(n_wl):
        #   ln I = ln I0 - mua_muscle * L_muscle - mua_blood * L_vein
        lnI_pred[k, :] = (lnI0[k]
                          - mua_musc[k] * L_muscle[k, :]
                          - mua_blood[k] * L_vein[k, :])

    return lnI_pred


# ─────────────────────────────────────────────────────────────────────────────
# Forward Model: Respiratory-AC Component
# ─────────────────────────────────────────────────────────────────────────────

def predict_resp_ac(SvO2, Hb_mM, delta_r_over_r, wavelengths, L_vein):
    """
    Predict the respiratory-AC modulation ratio at each measurement channel.

    The respiratory modulation ratio is:

        M_resp(λ, x) = μₐ_blood(λ) * <L_v>(λ, x) * (Δr_v / r_v)

    This comes from the first-order perturbation of the Beer-Lambert law
    when the vessel radius changes by Δr_v during the respiratory cycle.

    Key property: the modulation ratio does NOT depend on I₀ (the
    calibration constant cancels when dividing AC by DC). This makes
    it self-calibrating.

    Parameters
    ----------
    SvO2 : float
        Venous oxygen saturation (0-1).
    Hb_mM : float
        Hemoglobin concentration in mM.
    delta_r_over_r : float
        Fractional vessel radius change during respiration.
        Simulated value: (0.0575 - 0.05) / 0.05 = 0.15.
    wavelengths : ndarray, shape (n_wl,)
        Wavelengths in nm.
    L_vein : ndarray, shape (n_wl, n_pos)
        Mean photon pathlength in the vessel [cm], from MC simulation.

    Returns
    -------
    M_pred : ndarray, shape (n_wl, n_pos)
        Predicted modulation ratio at each channel.

    Notes
    -----
    The modulation ratio encodes SvO2 through the spectral dependence
    of μₐ_blood. Taking the ratio of M_resp at two wavelengths cancels
    Δr/r and [Hb], leaving only the blood absorption ratio — which is
    a function of SvO2 alone.
    """
    n_wl, n_pos = L_vein.shape

    # Blood absorption spectrum
    mua_blood = op.mua_blood_spectrum(wavelengths, SvO2, Hb_mM)  # (n_wl,)

    # Modulation ratio at each channel
    M_pred = np.zeros((n_wl, n_pos))
    for k in range(n_wl):
        M_pred[k, :] = mua_blood[k] * L_vein[k, :] * delta_r_over_r

    return M_pred


# ─────────────────────────────────────────────────────────────────────────────
# Residual Function (what the optimizer minimizes)
# ─────────────────────────────────────────────────────────────────────────────

def residual_function(p_vec, wavelengths, L_muscle, L_vein,
                      d_dc, d_resp, w_resp):
    """
    Compute the residual vector for the joint DC + respiratory-AC fit.

    The optimizer minimizes sum(residuals^2). This function computes
    the difference between the forward model predictions and the
    measured data for all channels.

    The residual vector is:
        r = [ (ln I_pred - ln I_meas)  for all DC channels,
              w * (M_pred - M_meas)    for all AC channels ]

    where w is a weight that balances the DC and AC contributions.

    Parameters
    ----------
    p_vec : ndarray, shape (n_params,)
        Parameter vector. Layout:
            p_vec[0]          = SvO2
            p_vec[1]          = Hb_mM
            p_vec[2:2+n_wl]   = mua_musc at each wavelength
            p_vec[2+n_wl:2+2*n_wl] = lnI0 at each wavelength
            p_vec[2+2*n_wl]   = delta_r_v_resp / r_v

    wavelengths : ndarray, shape (n_wl,)
        Wavelengths in nm.

    L_muscle : ndarray, shape (n_wl, n_pos)
        Muscle pathlengths from MC.

    L_vein : ndarray, shape (n_wl, n_pos)
        Vessel pathlengths from MC.

    d_dc : ndarray, shape (n_wl * n_pos,)
        Measured DC data: ln(I_B) flattened.

    d_resp : ndarray, shape (n_wl * n_pos,)
        Measured respiratory-AC data: M_resp flattened.

    w_resp : float
        Weight applied to respiratory-AC residuals to balance
        their contribution against DC residuals.

    Returns
    -------
    residuals : ndarray, shape (2 * n_wl * n_pos,)
        Concatenated DC and weighted AC residuals.
    """
    n_wl = len(wavelengths)

    # Unpack parameter vector
    SvO2 = p_vec[0]
    Hb_mM = p_vec[1]
    mua_musc = p_vec[2: 2 + n_wl]
    lnI0 = p_vec[2 + n_wl: 2 + 2 * n_wl]
    delta_r = p_vec[2 + 2 * n_wl]

    # --- DC residuals ---
    lnI_pred = predict_dc(SvO2, Hb_mM, mua_musc, lnI0,
                          wavelengths, L_muscle, L_vein)
    res_dc = lnI_pred.ravel() - d_dc

    # --- Respiratory-AC residuals ---
    M_pred = predict_resp_ac(SvO2, Hb_mM, delta_r,
                             wavelengths, L_vein)
    res_resp = w_resp * (M_pred.ravel() - d_resp)

    # Stack into one residual vector
    return np.concatenate([res_dc, res_resp])


def residual_dc_only(p_vec, wavelengths, L_muscle, L_vein, d_dc):
    """
    Residual function for DC-only fit (no AC data).

    Same as residual_function but without the respiratory-AC component.
    Use this if the perturbation simulation data is not yet available.

    Parameters
    ----------
    p_vec : ndarray, shape (2 + 2*n_wl,)
        [SvO2, Hb_mM, mua_musc(n_wl), lnI0(n_wl)]
        Note: no delta_r parameter.
    wavelengths, L_muscle, L_vein, d_dc : as in residual_function.

    Returns
    -------
    residuals : ndarray, shape (n_wl * n_pos,)
    """
    n_wl = len(wavelengths)
    SvO2 = p_vec[0]
    Hb_mM = p_vec[1]
    mua_musc = p_vec[2: 2 + n_wl]
    lnI0 = p_vec[2 + n_wl: 2 + 2 * n_wl]

    lnI_pred = predict_dc(SvO2, Hb_mM, mua_musc, lnI0,
                          wavelengths, L_muscle, L_vein)
    return lnI_pred.ravel() - d_dc


# ─────────────────────────────────────────────────────────────────────────────
# Inverse Solver
# ─────────────────────────────────────────────────────────────────────────────

def solve_svo2(wavelengths, L_muscle, L_vein,
               I_baseline, I_perturbed=None,
               I_bulk=None, L_bulk_muscle=None,
               SvO2_init=0.70, Hb_init=2.0,
               verbose=True):
    """
    Solve the inverse problem to recover SvO2 from measured intensities.

    This is the main function of the pipeline. It takes the measured
    (or simulated) data and returns the best-fit tissue parameters.

    Parameters
    ----------
    wavelengths : ndarray, shape (n_wl,)
        Wavelengths in nm (e.g., [660, 850, 940]).

    L_muscle : ndarray, shape (n_wl, n_pos)
        Mean pathlength in muscle from MC simulation [cm].

    L_vein : ndarray, shape (n_wl, n_pos)
        Mean pathlength in vessel from MC simulation [cm].

    I_baseline : ndarray, shape (n_wl, n_pos)
        Detected intensities from the baseline vessel simulation
        (Dataset B: vessel at normal radius).

    I_perturbed : ndarray, shape (n_wl, n_pos), optional
        Detected intensities from the expanded vessel simulation
        (Dataset C: vessel at respiratory-perturbed radius).
        If provided, respiratory-AC data is included in the fit.
        If None, DC-only fit is performed.

    I_bulk : ndarray, shape (n_wl,), optional
        Detected intensities from bulk muscle simulation (Dataset A).
        Used to initialize the calibration constants lnI0.

    L_bulk_muscle : ndarray, shape (n_wl,), optional
        Mean pathlength from bulk muscle simulation.
        Used with I_bulk to initialize lnI0.

    SvO2_init : float
        Initial guess for SvO2.

    Hb_init : float
        Initial guess for [Hb] in mM.

    verbose : bool
        If True, print optimizer progress.

    Returns
    -------
    dict with keys:
        'SvO2'          : float — recovered venous oxygen saturation
        'Hb_mM'         : float — recovered hemoglobin concentration [mM]
        'mua_musc'      : ndarray (n_wl,) — recovered muscle absorption [cm^-1]
        'lnI0'          : ndarray (n_wl,) — recovered calibration constants
        'delta_r_over_r': float or None — recovered pulsation amplitude
        'cost'          : float — final cost function value
        'residuals'     : ndarray — final residual vector
        'success'       : bool — whether optimizer converged
        'message'       : str — optimizer status message
        'p_vec'         : ndarray — full parameter vector at solution
        'jacobian'      : ndarray — Jacobian at solution (for uncertainty)
        'mode'          : str — 'dc_only' or 'dc_plus_resp_ac'
    """
    n_wl = len(wavelengths)
    n_pos = L_muscle.shape[1]
    use_resp_ac = I_perturbed is not None

    mode = 'dc_plus_resp_ac' if use_resp_ac else 'dc_only'
    print(f"\n  Solver mode: {mode}")
    print(f"  Channels: {n_wl} wavelengths × {n_pos} positions = "
          f"{n_wl * n_pos} per data type")

    # ── Step A: Prepare measured data vectors ──────────────────────────────

    # DC data: natural log of baseline intensities
    d_dc = np.log(I_baseline).ravel()  # shape: (n_wl * n_pos,)

    # Respiratory-AC data: modulation ratio
    if use_resp_ac:
        M_resp_measured = (I_perturbed - I_baseline) / I_baseline
        d_resp = M_resp_measured.ravel()  # shape: (n_wl * n_pos,)

        # Compute weight to balance DC and AC residuals.
        # DC residuals are ~O(10) in magnitude (log-intensity).
        # AC residuals are ~O(0.01-0.10) (modulation ratio).
        # Without weighting, optimizer would ignore the AC data.
        mean_dc_magnitude = np.mean(np.abs(d_dc))
        mean_resp_magnitude = np.mean(np.abs(d_resp))
        if mean_resp_magnitude > 0:
            w_resp = mean_dc_magnitude / mean_resp_magnitude
        else:
            w_resp = 1.0
        print(f"  AC weight: {w_resp:.2f} "
              f"(DC magnitude: {mean_dc_magnitude:.2f}, "
              f"AC magnitude: {mean_resp_magnitude:.4f})")
        print(f"  Total data points: {2 * n_wl * n_pos}")
    else:
        print(f"  Total data points: {n_wl * n_pos}")

    # ── Step B: Set initial guesses ───────────────────────────────────────

    # Muscle absorption: start from Jacques model (literature values)
    mua_musc_init = np.array([op.muscle_mua(wl) for wl in wavelengths])

    # Calibration constants: use bulk muscle data if available
    if I_bulk is not None and L_bulk_muscle is not None:
        # In bulk muscle: ln I_bulk = lnI0 - mua_musc * L_bulk
        # Rearrange:       lnI0 = ln I_bulk + mua_musc * L_bulk
        lnI0_init = np.log(I_bulk) + mua_musc_init * L_bulk_muscle
    else:
        # Rough estimate from the vessel data at far-from-vessel positions.
        # At large |x|, L_vein ≈ 0, so ln I ≈ lnI0 - mua_musc * L_muscle.
        far_idx = -1  # last position (farthest from vessel center)
        lnI0_init = np.zeros(n_wl)
        for k in range(n_wl):
            lnI0_init[k] = (np.log(I_baseline[k, far_idx])
                            + mua_musc_init[k] * L_muscle[k, far_idx])

    # Build initial parameter vector
    if use_resp_ac:
        delta_r_init = 0.10  # initial guess for Δr/r (true = 0.15)
        p0 = np.concatenate([
            [SvO2_init, Hb_init],      # 2 params
            mua_musc_init,              # n_wl params
            lnI0_init,                  # n_wl params
            [delta_r_init],             # 1 param
        ])
    else:
        p0 = np.concatenate([
            [SvO2_init, Hb_init],
            mua_musc_init,
            lnI0_init,
        ])

    # ── Step C: Set parameter bounds ──────────────────────────────────────

    if use_resp_ac:
        lb = np.concatenate([
            [0.30, 1.0],                        # SvO2, Hb
            np.full(n_wl, 0.001),               # mua_musc lower
            np.full(n_wl, -np.inf),             # lnI0 lower (unconstrained)
            [0.0],                               # delta_r/r lower
        ])
        ub = np.concatenate([
            [1.00, 3.5],
            np.full(n_wl, 1.0),
            np.full(n_wl, np.inf),
            [0.50],
        ])
    else:
        lb = np.concatenate([
            [0.30, 1.0],
            np.full(n_wl, 0.001),
            np.full(n_wl, -np.inf),
        ])
        ub = np.concatenate([
            [1.00, 3.5],
            np.full(n_wl, 1.0),
            np.full(n_wl, np.inf),
        ])

    print(f"  Parameters: {len(p0)} unknowns")
    print(f"  Initial SvO2 = {SvO2_init:.2f}, [Hb] = {Hb_init:.1f} mM")

    # ── Step D: Run the optimizer ─────────────────────────────────────────

    if use_resp_ac:
        result = least_squares(
            residual_function,
            p0,
            args=(wavelengths, L_muscle, L_vein, d_dc, d_resp, w_resp),
            bounds=(lb, ub),
            method='trf',
            verbose=2 if verbose else 0,
            ftol=1e-12, xtol=1e-12, gtol=1e-12,
            max_nfev=500,
            jac='3-point',
        )
    else:
        result = least_squares(
            residual_dc_only,
            p0,
            args=(wavelengths, L_muscle, L_vein, d_dc),
            bounds=(lb, ub),
            method='trf',
            verbose=2 if verbose else 0,
            ftol=1e-12, xtol=1e-12, gtol=1e-12,
            max_nfev=500,
            jac='3-point',
        )

    # ── Step E: Unpack results ────────────────────────────────────────────

    p_fit = result.x

    output = {
        'SvO2': p_fit[0],
        'Hb_mM': p_fit[1],
        'mua_musc': p_fit[2: 2 + n_wl],
        'lnI0': p_fit[2 + n_wl: 2 + 2 * n_wl],
        'delta_r_over_r': p_fit[2 + 2 * n_wl] if use_resp_ac else None,
        'cost': result.cost,
        'residuals': result.fun,
        'success': result.success,
        'message': result.message,
        'p_vec': p_fit,
        'jacobian': result.jac,
        'mode': mode,
    }

    return output


# ─────────────────────────────────────────────────────────────────────────────
# Results Display
# ─────────────────────────────────────────────────────────────────────────────

def print_results(output, wavelengths, true_SvO2=None, true_Hb=None,
                  true_delta_r=None):
    """
    Print a formatted summary of the inverse solver results.

    Parameters
    ----------
    output : dict
        Output from solve_svo2().
    wavelengths : ndarray
        Wavelengths in nm.
    true_SvO2 : float, optional
        True SvO2 for comparison.
    true_Hb : float, optional
        True [Hb] for comparison.
    true_delta_r : float, optional
        True Δr/r for comparison.
    """
    print("\n" + "=" * 65)
    print("  INVERSE SOLVER RESULTS")
    print("=" * 65)
    print(f"  Mode:            {output['mode']}")
    print(f"  Converged:       {output['success']}")
    print(f"  Cost (Σr²/2):    {output['cost']:.6e}")

    print(f"\n  --- Recovered Parameters ---")
    print(f"  SvO2:            {output['SvO2']:.4f}  ({output['SvO2']*100:.1f}%)")
    if true_SvO2 is not None:
        err = (output['SvO2'] - true_SvO2) * 100
        print(f"  SvO2 (true):     {true_SvO2:.4f}  ({true_SvO2*100:.1f}%)")
        print(f"  SvO2 error:      {err:+.2f} percentage points")

    print(f"\n  [Hb]:            {output['Hb_mM']:.4f} mM")
    if true_Hb is not None:
        print(f"  [Hb] (true):     {true_Hb:.4f} mM")

    if output['delta_r_over_r'] is not None:
        print(f"\n  Δr/r (resp):     {output['delta_r_over_r']:.4f}")
        if true_delta_r is not None:
            print(f"  Δr/r (true):     {true_delta_r:.4f}")

    print(f"\n  --- Muscle Absorption ---")
    for k, wl in enumerate(wavelengths):
        lit = op.muscle_mua(wl)
        rec = output['mua_musc'][k]
        print(f"  μₐ({wl:.0f} nm):     {rec:.5f} cm⁻¹  "
              f"(literature: {lit:.5f})")

    print("=" * 65)


# ─────────────────────────────────────────────────────────────────────────────
# Noise Injection for Robustness Testing
# ─────────────────────────────────────────────────────────────────────────────

def add_noise(I, noise_fraction, seed=None):
    """
    Add multiplicative Gaussian noise to intensity measurements.

    Simulates detector noise: I_noisy = I * (1 + N(0, sigma)).

    Parameters
    ----------
    I : ndarray
        Clean intensity measurements.
    noise_fraction : float
        Standard deviation of the noise relative to the signal.
        0.01 = 1% noise, 0.05 = 5% noise.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    ndarray
        Noisy intensities (same shape as input).
    """
    rng = np.random.default_rng(seed)
    return I * (1.0 + rng.normal(0, noise_fraction, size=I.shape))


# ─────────────────────────────────────────────────────────────────────────────
# Self-Test with Synthetic Data
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("Testing forward model + inverse solver with synthetic data...\n")

    wavelengths = np.array([660.0, 850.0, 940.0])
    n_pos = 10

    # Synthetic pathlengths (roughly realistic)
    L_muscle = np.ones((3, n_pos)) * 0.5
    L_vein = np.zeros((3, n_pos))
    for i in range(n_pos):
        x = -9.6 + i * 2.4
        L_vein[:, i] = 0.03 * np.exp(-x**2 / (2 * 5**2))

    # True parameters
    true_SvO2 = 0.75
    true_Hb = 2.3
    true_mua_musc = np.array([op.muscle_mua(wl) for wl in wavelengths])
    true_lnI0 = np.array([-8.0, -8.5, -9.0])
    true_delta_r = 0.15

    # Generate synthetic DC data
    lnI_dc = predict_dc(true_SvO2, true_Hb, true_mua_musc, true_lnI0,
                        wavelengths, L_muscle, L_vein)
    I_baseline = np.exp(lnI_dc)

    # Generate synthetic respiratory-AC data
    M_resp = predict_resp_ac(true_SvO2, true_Hb, true_delta_r,
                             wavelengths, L_vein)
    I_perturbed = I_baseline * (1.0 + M_resp)

    # Solve
    output = solve_svo2(
        wavelengths, L_muscle, L_vein,
        I_baseline, I_perturbed,
        SvO2_init=0.60, Hb_init=2.0,
        verbose=False,
    )

    print_results(output, wavelengths,
                  true_SvO2=true_SvO2, true_Hb=true_Hb,
                  true_delta_r=true_delta_r)
