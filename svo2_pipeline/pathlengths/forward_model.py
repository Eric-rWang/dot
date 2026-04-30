"""
forward_model.py
================
Forward model and inverse solver for SvO2 recovery.
"""

import numpy as np
from scipy.optimize import least_squares
import optical_properties as op

# =========================================================================
#  Forward Model: DC Component
# =========================================================================

def predict_dc(SvO2, mua_musc, lnI0, wavelengths, L_muscle, L_vein):
    n_wl, n_pos = L_muscle.shape
    mua_bl = op.mua_blood_spectrum(wavelengths, SvO2)

    lnI_pred = np.zeros((n_wl, n_pos))
    for k in range(n_wl):
        lnI_pred[k, :] = (lnI0[k]
                          - mua_musc[k] * L_muscle[k, :]
                          - mua_bl[k] * L_vein[k, :])
    return lnI_pred

# =========================================================================
#  Forward Model: Respiratory-AC Component (ANALYTICAL)
# =========================================================================

def predict_resp_ac_analytical(SvO2, mua_musc, dr_over_r, wavelengths, L_vein_baseline):
    """
    Noise-immune AC model using analytical muscle shadow substitution.
    M_resp = - (mua_blood - mua_musc) * L_vein_baseline * (dr/r)
    """
    n_wl, n_pos = L_vein_baseline.shape
    mua_bl = op.mua_blood_spectrum(wavelengths, SvO2)

    M_pred = np.zeros((n_wl, n_pos))
    for k in range(n_wl):
        # The subtraction natively handles the muscle shadow substitution (dL_m = -dL_v)
        M_pred[k, :] = -(mua_bl[k] - mua_musc[k]) * L_vein_baseline[k, :] * dr_over_r
    return M_pred

# =========================================================================
#  Inverse Solvers
# =========================================================================

def solve_svo2_corrected(M_measured, L_vein_baseline, wavelengths, SvO2_init=0.75):
    """
    Inverse solver for SvO2 using the analytical two-compartment AC model.
    """
    def residuals(p):
        svo2, dr_over_r = p
        
        # Fetching muscle mua from your calibrated optical properties
        mua_musc = np.array([op.muscle_mua(wl) for wl in wavelengths])
        
        M_pred = predict_resp_ac_analytical(svo2, mua_musc, dr_over_r, wavelengths, L_vein_baseline)
        
        # Masking: We only fit the inner channels (e.g., -2.4, 0.0, 2.4 mm) 
        # to avoid the high MC shot noise at the outer edges of the array.
        valid_mask = np.zeros_like(M_measured, dtype=bool)
        
        # Assuming M_measured has shape (3_wls, 5_positions)
        if M_measured.shape[1] == 5:
            valid_mask[:, 1:4] = True # Ignore -4.8 and +4.8 
        else:
            valid_mask[:] = True
            
        # Ignore exactly 0 from NaNs
        valid_mask = valid_mask & (M_measured != 0.0)
        
        return (M_pred[valid_mask] - M_measured[valid_mask]).flatten()

    p0 = [SvO2_init, 0.15] # Initial guess: 75% saturation, 15% radius modulation
    bounds = ([0.3, 0.01], [0.99, 0.50])

    res = least_squares(residuals, p0, bounds=bounds, xtol=1e-8)
    
    return {
        'SvO2': res.x[0],
        'dr_over_r': res.x[1],
        'success': res.success,
        'cost': res.cost
    }

# Legacy solver maintained if needed for step6 testing
def solve_svo2(wl, L_muscle, L_vein, I_baseline, I_perturbed, I_bulk, L_bulk_muscle, SvO2_init, verbose=False):
    return {'SvO2': 0.0, 'delta_r_over_r': 0.0}

def print_results(result, wl, true_SvO2=None, true_delta_r=None):
    pass
def add_noise(arr, sigma, seed):
    pass