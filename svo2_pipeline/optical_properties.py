"""
optical_properties.py
=====================

Chromophore extinction spectra, water absorption, and tissue optical property
models used throughout the SvO2 recovery pipeline.

This module answers the question: "Given a wavelength and tissue composition,
what are the absorption and scattering coefficients?"

Physical background
-------------------
Light absorption in tissue is due to chromophores — molecules that absorb
photons at specific wavelengths. The dominant chromophores in the near-infrared
(600-1000 nm) are:

    - Oxyhemoglobin (HbO2): absorbs more at ~900-940 nm
    - Deoxyhemoglobin (HbR): absorbs more at ~660 nm
    - Water: absorbs increasingly above ~900 nm

The absorption coefficient of blood depends on the *ratio* of HbO2 to HbR,
which is the oxygen saturation. This wavelength-dependent ratio is what
allows us to determine SvO2 from multi-wavelength measurements.

Data sources
------------
    - HbO2/HbR extinction: Prahl (omlc.org), compiled from Gratzer and Cope
    - Water absorption: Hale & Querry (1973), Kou et al. (1993)
    - Tissue scattering model: Jacques (2013), Phys. Med. Biol. 58, R37-R61

Units convention
----------------
    - Extinction coefficients: cm^-1 / M (base-10, as conventionally tabulated)
    - Absorption coefficients: cm^-1 (natural log, as used in transport theory)
    - Concentrations: millimolar (mM)
    - The conversion factor is ln(10) = 2.302585
"""

import numpy as np
from scipy.interpolate import PchipInterpolator


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

LN10 = np.log(10.0)       # 2.302585... — converts base-10 to natural log
MW_HB = 64500.0            # g/mol — molecular weight of hemoglobin tetramer


# ─────────────────────────────────────────────────────────────────────────────
# Tabulated Spectra
# ─────────────────────────────────────────────────────────────────────────────
# Selected wavelengths from the Prahl/Gratzer compilation.
# For publication work, use the full table from:
#   https://omlc.org/spectra/hemoglobin/summary.html
#
# Values are molar extinction coefficients in cm^-1 / M (base-10).

_WL = np.array([
    540, 560, 580, 600, 620, 630, 640, 660, 680, 700,
    720, 740, 750, 760, 780, 800, 820, 840, 850, 860,
    880, 900, 920, 940, 960, 980, 1000
], dtype=np.float64)

_EPS_HBO2 = np.array([
    52880, 34360, 50104, 1824, 640, 488, 384, 320, 248, 392,
    516, 756, 876, 1016, 1192, 1620, 1740, 2072, 2296, 2212,
    1620, 1308, 1060, 1128, 872, 652, 504
], dtype=np.float64)

_EPS_HBR = np.array([
    43440, 27984, 33560, 9544, 5268, 3744, 3440, 3226, 2652, 2780,
    2012, 1620, 1378, 1212, 1044, 1798, 796, 764, 780, 716,
    800, 756, 660, 616, 420, 360, 312
], dtype=np.float64)

# Water absorption coefficient [cm^-1] (natural log)
_WATER_WL = np.array([
    540, 560, 580, 600, 620, 640, 660, 680, 700, 720,
    740, 760, 780, 800, 820, 840, 860, 880, 900, 920,
    940, 960, 980, 1000
], dtype=np.float64)

_WATER_MUA = np.array([
    0.00040, 0.00058, 0.00068, 0.0022, 0.0028, 0.0029, 0.0032, 0.0035,
    0.0060, 0.0134, 0.0210, 0.0264, 0.0198, 0.0204, 0.0300, 0.0364,
    0.0396, 0.0420, 0.0600, 0.1600, 0.2700, 0.3200, 0.3600, 0.4500
], dtype=np.float64)

# Build smooth interpolators (PCHIP = piecewise cubic Hermite)
_interp_HbO2 = PchipInterpolator(_WL, _EPS_HBO2)
_interp_HbR = PchipInterpolator(_WL, _EPS_HBR)
_interp_water = PchipInterpolator(_WATER_WL, _WATER_MUA)


# ─────────────────────────────────────────────────────────────────────────────
# Extinction Coefficient Functions
# ─────────────────────────────────────────────────────────────────────────────

def epsilon_HbO2(wavelength_nm):
    """
    Molar extinction coefficient of oxyhemoglobin.

    Parameters
    ----------
    wavelength_nm : float or array
        Wavelength in nanometers.

    Returns
    -------
    float or array
        Extinction coefficient in cm^-1 / M (base-10).

    Example
    -------
    At 660 nm, HbO2 has low extinction (~320 cm^-1/M),
    while at 940 nm it has higher extinction (~1128 cm^-1/M).
    """
    return np.float64(_interp_HbO2(wavelength_nm))


def epsilon_HbR(wavelength_nm):
    """
    Molar extinction coefficient of deoxyhemoglobin.

    Parameters
    ----------
    wavelength_nm : float or array
        Wavelength in nanometers.

    Returns
    -------
    float or array
        Extinction coefficient in cm^-1 / M (base-10).

    Example
    -------
    At 660 nm, HbR has high extinction (~3226 cm^-1/M) — about 10x
    higher than HbO2. This large contrast is why 660 nm is sensitive
    to deoxygenation.
    """
    return np.float64(_interp_HbR(wavelength_nm))


def mua_water(wavelength_nm):
    """
    Absorption coefficient of pure water.

    Parameters
    ----------
    wavelength_nm : float or array
        Wavelength in nanometers.

    Returns
    -------
    float or array
        Absorption coefficient in cm^-1.

    Notes
    -----
    Water absorption is negligible at 660 nm (~0.003 cm^-1) but becomes
    significant at 940 nm (~0.27 cm^-1). This is why 940 nm is useful
    for water content estimation but also why it complicates the
    hemoglobin measurement.
    """
    return np.float64(_interp_water(wavelength_nm))


# ─────────────────────────────────────────────────────────────────────────────
# Blood Absorption Model
# ─────────────────────────────────────────────────────────────────────────────

def mua_blood(wavelength_nm, S, Hb_mM, W_blood=0.95):
    """
    Compute the absorption coefficient of whole blood.

    This is the central equation connecting SvO2 to the optical measurement:

        mu_a = (ln10 * [Hb]) * [S * eps_HbO2 + (1-S) * eps_HbR]
               + mu_a_water * W_blood

    The first term is the hemoglobin contribution: a weighted sum of the
    two hemoglobin species, where the weights are set by the oxygen
    saturation S. The second term is the water background.

    Parameters
    ----------
    wavelength_nm : float
        Wavelength in nm.
    S : float
        Oxygen saturation, 0 to 1.
        - For venous blood: S = SvO2 (the quantity we want to recover).
        - For arterial blood: S = SaO2 (known from clinical monitoring).
    Hb_mM : float
        Total hemoglobin concentration in millimolar.
        Typical whole blood: ~2.3 mM (~15 g/dL).
    W_blood : float, optional
        Water volume fraction in blood. Default 0.95 (plasma is ~95% water).

    Returns
    -------
    float
        Absorption coefficient in cm^-1.

    Notes
    -----
    The key insight: changing S rotates the absorption spectrum (shifting
    weight between HbO2 and HbR), while changing [Hb] scales it uniformly.
    Multi-wavelength measurements can disentangle the two.
    """
    Hb_M = Hb_mM / 1000.0  # convert millimolar to molar

    eps_O2 = epsilon_HbO2(wavelength_nm)
    eps_R = epsilon_HbR(wavelength_nm)

    # Hemoglobin contribution (base-10 extinction → natural-log absorption)
    mua_hb = LN10 * Hb_M * (S * eps_O2 + (1.0 - S) * eps_R)

    # Water contribution
    mua_w = mua_water(wavelength_nm) * W_blood

    return mua_hb + mua_w


def mua_blood_spectrum(wavelengths_nm, S, Hb_mM, W_blood=0.95):
    """
    Blood absorption at multiple wavelengths (convenience wrapper).

    Parameters
    ----------
    wavelengths_nm : array-like
        Wavelengths in nm.
    S : float
        Oxygen saturation.
    Hb_mM : float
        Hemoglobin concentration in mM.
    W_blood : float, optional
        Water fraction in blood.

    Returns
    -------
    ndarray, shape (n_wavelengths,)
        Absorption coefficients in cm^-1.
    """
    return np.array([mua_blood(wl, S, Hb_mM, W_blood)
                     for wl in wavelengths_nm])


# ─────────────────────────────────────────────────────────────────────────────
# Tissue Scattering Model — Jacques (2013)
# ─────────────────────────────────────────────────────────────────────────────

def mus_prime(wavelength_nm, a_prime, f_ray, b_mie):
    """
    Reduced scattering coefficient using the Jacques (2013) model.

        mu_s' = a' * [f_Ray * (lam/500)^-4 + (1-f_Ray) * (lam/500)^-b_Mie]

    The first term is Rayleigh scattering (small structures, strong
    wavelength dependence). The second is Mie scattering (larger
    structures, weaker wavelength dependence).

    Parameters
    ----------
    wavelength_nm : float
        Wavelength in nm.
    a_prime : float
        Reduced scattering amplitude at 500 nm [cm^-1].
    f_ray : float
        Fraction of scattering due to Rayleigh scattering (0-1).
    b_mie : float
        Mie scattering power exponent.

    Returns
    -------
    float
        Reduced scattering coefficient in cm^-1.
    """
    lam_ratio = wavelength_nm / 500.0
    return a_prime * (f_ray * lam_ratio**(-4) + (1.0 - f_ray) * lam_ratio**(-b_mie))


# ─────────────────────────────────────────────────────────────────────────────
# Default Tissue Parameters (matching MCmatlab mediaPropertiesFunc)
# ─────────────────────────────────────────────────────────────────────────────

# Muscle (compartment index 5 in the simulation)
MUSCLE = {
    'B': 0.002, 'S': 0.67, 'W': 0.65, 'F': 0.0, 'M': 0.0,
    'a_prime': 42.4, 'f_ray': 0.62, 'b_mie': 1.0, 'g': 0.9,
}

# Blood scattering (same for arterial and venous)
BLOOD_SCATTER = {
    'a_prime': 10.0, 'f_ray': 0.0, 'b_mie': 1.0, 'g': 0.9,
}


def muscle_mua(wavelength_nm):
    """Muscle absorption at a given wavelength using default parameters."""
    p = MUSCLE
    # Simplified Jacques model: mu_a = B * mu_a_blood + W * mu_a_water
    mua_bl = mua_blood(wavelength_nm, S=p['S'], Hb_mM=2.3)
    return p['B'] * mua_bl + p['W'] * mua_water(wavelength_nm)


def muscle_mus_prime(wavelength_nm):
    """Muscle reduced scattering at a given wavelength."""
    p = MUSCLE
    return mus_prime(wavelength_nm, p['a_prime'], p['f_ray'], p['b_mie'])


# ─────────────────────────────────────────────────────────────────────────────
# Self-Test
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    wavelengths = [660, 850, 940]

    print("Blood and Tissue Optical Properties at Sensor Wavelengths")
    print("=" * 72)
    print(f"{'λ [nm]':>8}  {'mua_musc':>9}  {'musp_musc':>10}  "
          f"{'mua_blood':>10}  {'mua_blood':>10}  {'mua_water':>10}")
    print(f"{'':>8}  {'[cm-1]':>9}  {'[cm-1]':>10}  "
          f"{'S=0.75':>10}  {'S=0.60':>10}  {'[cm-1]':>10}")
    print("-" * 72)
    for wl in wavelengths:
        print(f"{wl:>8.0f}  {muscle_mua(wl):>9.4f}  {muscle_mus_prime(wl):>10.2f}  "
              f"{mua_blood(wl, 0.75, 2.3):>10.4f}  {mua_blood(wl, 0.60, 2.3):>10.4f}  "
              f"{mua_water(wl):>10.4f}")

    print(f"\nBlood absorption ratio 660/940:")
    for s in [0.60, 0.75, 0.90, 0.98]:
        ratio = mua_blood(660, s, 2.3) / mua_blood(940, s, 2.3)
        print(f"  S = {s:.2f}: ratio = {ratio:.3f}")
