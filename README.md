# Continuous-Wave Diffuse Optical Spectroscopy for Venous Oxygen Saturation Measurement

## Technical Description of the Proposed Pipeline

---

## 1. Overview

This document describes a method for non-invasive measurement of venous oxygen saturation (SvO₂) using a continuous-wave (CW) diffuse optical spectroscopy (DOSI) e-tattoo sensor array. The sensor is placed over superficial vasculature, and the measurement combines optical data with ultrasound-derived anatomical geometry to recover SvO₂ through a Monte Carlo–based inverse problem.

The clinical study targets three anatomical sites on anesthetized subjects:

- **Forearm:** isolated superficial vein (e.g., cephalic vein), no major artery directly beneath the sensor.
- **Wrist:** radial artery and accompanying vein (cephalic vein or radial venae comitantes).
- **Neck:** internal jugular vein (IJV) and common carotid artery.

Ground-truth SvO₂ is obtained via invasive catheter. Arterial oxygen saturation (SaO₂) is recorded via standard clinical monitoring. Subjects are under general anaesthesia with controlled mechanical ventilation.

---

## 2. Sensor Hardware

### 2.1 Array Geometry

The sensor consists of **12 photodetector (PD) sites** arranged in a linear array with uniform spacing:

$$d_{PD} = 2.4 \text{ mm}$$

Total array length:

$$L_{array} = (N_{PD} - 1) \times d_{PD} = 11 \times 2.4 = 26.4 \text{ mm}$$

where $N_{PD} = 12$ is the number of photodetectors.

Each PD site has **3 LEDs** at wavelengths:

$$\lambda_1 = 660 \text{ nm}, \quad \lambda_2 = 850 \text{ nm}, \quad \lambda_3 = 940 \text{ nm}$$

### 2.2 Source–Detector Separations

The LEDs are positioned across from (perpendicular to) each PD along the vessel axis. The source–detector separation ρ depends on wavelength:

$$\rho(\lambda) = \begin{cases} 6.3342 \text{ mm} & \text{for } \lambda = 660 \text{ nm and } 940 \text{ nm} \\ 7.9 \text{ mm} & \text{for } \lambda = 850 \text{ nm} \end{cases}$$

Each LED is pulsed independently and only the PD directly across from the active LED is sampled. This yields **one measurement per PD per wavelength**, for a total of:

$$N_{meas} = N_{PD} \times N_\lambda = 12 \times 3 = 36 \text{ measurements}$$

### 2.3 Sensor Orientation

The LED–PD axis (the line connecting each LED to its corresponding PD) is aligned **parallel** to the underlying vessel(s). The 12 PD sites are distributed along a line **perpendicular** to the vessel axis. This means each PD site samples a different lateral offset from the vessel center.

### 2.4 Coordinate System

We adopt the following coordinate convention, consistent with the MCmatlab simulation:

| Axis | Direction | Notes |
|------|-----------|-------|
| **x** | Along the PD array (perpendicular to vessel) | PD positions indexed by $x_i$ |
| **y** | Along the LED–PD axis (parallel to vessel) | Source at $+y$, detector at $-y$ |
| **z** | Depth into tissue | $z = 0$ at tissue surface |

PD positions along the array:

$$x_i, \quad i = 1, 2, \ldots, 12$$

LED positions for PD at $x_i$:

$$\text{LED}(\lambda, x_i) = (x_i, \; +\rho(\lambda)/2, \; 0)$$

PD position:

$$\text{PD}(x_i) = (x_i, \; -\rho(\lambda)/2, \; 0)$$

---

## 3. Tissue Model

### 3.1 Compartments

The tissue volume beneath the sensor is modeled as a set of discrete compartments, each with uniform optical properties. For the simplest case (single vein in muscle), the compartments are:

| Index $j$ | Compartment | Description |
|-----------|-------------|-------------|
| 1 | Muscle (background) | Homogeneous bulk tissue surrounding the vessel |
| 2 | Venous blood | Cylindrical vessel embedded in the muscle |

For the vein + artery case, an additional compartment is added:

| Index $j$ | Compartment | Description |
|-----------|-------------|-------------|
| 1 | Muscle (background) | Homogeneous bulk tissue |
| 2 | Venous blood | Cylindrical vein |
| 3 | Arterial blood | Cylindrical artery |

Additional layers (epidermis, dermis, subcutaneous fat) can be added for more anatomically realistic models, particularly for the neck site.

### 3.2 Vessel Geometry

Each vessel is modeled as an infinite cylinder oriented along the y-axis (parallel to the LED–PD axis). The geometry is parameterized by:

| Symbol | Definition | Units |
|--------|-----------|-------|
| $r_v$ | Vein radius | cm |
| $d_v$ | Depth from tissue surface to vein center | cm |
| $x_v$ | Lateral position of vein center | cm |
| $r_a$ | Artery radius (if present) | cm |
| $d_a$ | Depth from tissue surface to artery center | cm |
| $x_a$ | Lateral position of artery center | cm |

The vessel cross-section occupies voxels satisfying:

$$\text{Vein:} \quad (x - x_v)^2 + (z - d_v)^2 \leq r_v^2$$

$$\text{Artery:} \quad (x - x_a)^2 + (z - d_a)^2 \leq r_a^2$$

These parameters are obtained from **B-mode ultrasound imaging** of the measurement site prior to optical data acquisition.

### 3.3 Optical Properties

Each compartment $j$ is characterized by four wavelength-dependent optical properties:

| Symbol | Definition | Units |
|--------|-----------|-------|
| $\mu_a^{(j)}(\lambda)$ | Absorption coefficient | cm⁻¹ |
| $\mu_s^{(j)}(\lambda)$ | Scattering coefficient | cm⁻¹ |
| $g^{(j)}$ | Scattering anisotropy factor | dimensionless |
| $n^{(j)}$ | Refractive index | dimensionless |

The reduced scattering coefficient is defined as:

$$\mu_s'^{(j)}(\lambda) = \mu_s^{(j)}(\lambda) \cdot (1 - g^{(j)})$$

#### 3.3.1 Absorption Coefficient Model

Tissue absorption is modeled as a linear combination of chromophore contributions following Jacques (2013):

$$\mu_a(\lambda) = B \cdot \mu_{a,\text{blood}}(\lambda, S) + W \cdot \mu_{a,\text{water}}(\lambda) + F \cdot \mu_{a,\text{fat}}(\lambda) + M \cdot \mu_{a,\text{melanin}}(\lambda)$$

where:

| Symbol | Definition | Range |
|--------|-----------|-------|
| $B$ | Blood volume fraction | 0–1 |
| $S$ | Blood oxygen saturation (SaO₂ or SvO₂) | 0–1 |
| $W$ | Water volume fraction | 0–1 |
| $F$ | Fat volume fraction | 0–1 |
| $M$ | Melanin volume fraction | 0–0.1 |

The blood absorption term is:

$$\mu_{a,\text{blood}}(\lambda, S) = S \cdot \mu_{a,\text{HbO}_2}(\lambda) + (1 - S) \cdot \mu_{a,\text{HbR}}(\lambda)$$

where:

$$\mu_{a,\text{HbO}_2}(\lambda) = \frac{\ln(10) \cdot \varepsilon_{\text{HbO}_2}(\lambda) \cdot [\text{Hb}]}{MW_{\text{Hb}}}$$

$$\mu_{a,\text{HbR}}(\lambda) = \frac{\ln(10) \cdot \varepsilon_{\text{HbR}}(\lambda) \cdot [\text{Hb}]}{MW_{\text{Hb}}}$$

| Symbol | Definition | Typical Value |
|--------|-----------|---------------|
| $\varepsilon_{\text{HbO}_2}(\lambda)$ | Molar extinction coefficient of oxyhemoglobin | From Prahl/Gratzer spectra (cm⁻¹ M⁻¹) |
| $\varepsilon_{\text{HbR}}(\lambda)$ | Molar extinction coefficient of deoxyhemoglobin | From Prahl/Gratzer spectra (cm⁻¹ M⁻¹) |
| $[\text{Hb}]$ | Total hemoglobin concentration | ~2.3 mM (~15 g/dL) in whole blood |
| $MW_{\text{Hb}}$ | Molecular weight of hemoglobin | 64,500 g/mol |

##### Compartment-Specific Absorption Parameters

**Muscle (background):**

| Parameter | Value | Notes |
|-----------|-------|-------|
| $B_{\text{musc}}$ | 0.002 | Low capillary blood content |
| $S_{\text{musc}}$ | 0.67 | Mixed capillary saturation |
| $W_{\text{musc}}$ | 0.65 | Muscle water content |
| $F_{\text{musc}}$ | 0 | No fat in muscle compartment |
| $M_{\text{musc}}$ | 0 | No melanin at depth |

**Venous blood:**

| Parameter | Value | Notes |
|-----------|-------|-------|
| $B_{\text{ven}}$ | 1.0 | Pure blood |
| $S_{\text{ven}}$ | **SvO₂ (unknown)** | **Primary quantity to recover** |
| $W_{\text{ven}}$ | 0.95 | Plasma water fraction |
| $F_{\text{ven}}$ | 0 | — |
| $M_{\text{ven}}$ | 0 | — |

**Arterial blood (if present):**

| Parameter | Value | Notes |
|-----------|-------|-------|
| $B_{\text{art}}$ | 1.0 | Pure blood |
| $S_{\text{art}}$ | **SaO₂ (known)** | Measured via pulse oximetry or arterial line |
| $W_{\text{art}}$ | 0.95 | Plasma water fraction |
| $F_{\text{art}}$ | 0 | — |
| $M_{\text{art}}$ | 0 | — |

#### 3.3.2 Scattering Coefficient Model

The reduced scattering coefficient follows the empirical model from Jacques (2013):

$$\mu_s'(\lambda) = a' \left[ f_{\text{Ray}} \left(\frac{\lambda}{500}\right)^{-4} + (1 - f_{\text{Ray}}) \left(\frac{\lambda}{500}\right)^{-b_{\text{Mie}}} \right]$$

where:

| Symbol | Definition | Units |
|--------|-----------|-------|
| $a'$ | Reduced scattering amplitude at 500 nm | cm⁻¹ |
| $f_{\text{Ray}}$ | Fraction of scattering due to Rayleigh scattering | dimensionless |
| $b_{\text{Mie}}$ | Mie scattering power | dimensionless |

The full scattering coefficient is recovered via:

$$\mu_s(\lambda) = \frac{\mu_s'(\lambda)}{1 - g}$$

##### Compartment-Specific Scattering Parameters

| Compartment | $a'$ (cm⁻¹) | $f_{\text{Ray}}$ | $b_{\text{Mie}}$ | $g$ |
|-------------|-------------|-------------------|-------------------|-----|
| Muscle | 42.4 | 0.62 | 1.0 | 0.9 |
| Blood (arterial and venous) | 10 | 0.0 | 1.0 | 0.9 |

Blood μₛ′ is assumed from the literature and is **not recovered** from the optical data.

---

## 4. Forward Model

### 4.1 Monte Carlo Simulation

The forward model is based on Monte Carlo (MC) simulation of photon transport through the voxelized tissue geometry. The simulation is performed using **MCmatlab**, which implements the stochastic solution to the radiative transport equation (RTE):

$$\frac{dL(\mathbf{r}, \hat{s})}{ds} = -(\mu_a + \mu_s) L(\mathbf{r}, \hat{s}) + \mu_s \int_{4\pi} p(\hat{s}, \hat{s}') L(\mathbf{r}, \hat{s}') \, d\Omega'$$

where:

| Symbol | Definition |
|--------|-----------|
| $L(\mathbf{r}, \hat{s})$ | Radiance at position $\mathbf{r}$ in direction $\hat{s}$ (W cm⁻² sr⁻¹) |
| $s$ | Path length along ray direction |
| $p(\hat{s}, \hat{s}')$ | Scattering phase function (Henyey–Greenstein with parameter $g$) |
| $d\Omega'$ | Solid angle element |

MC simulation does not impose the approximations inherent in the diffusion equation (P1 approximation to the RTE) and is therefore valid for:

- Heterogeneous media with sharp boundaries (vessel walls).
- Regions of high absorption ($\mu_a \sim \mu_s'$), such as inside blood vessels.
- Short source–detector separations ($\rho \sim 1/\mu_s'$).
- Geometries with embedded absorbers that violate the local homogeneity assumption of diffusion theory.

All of these conditions apply to the present measurement geometry.

### 4.2 Simulation Configuration

For each simulation run, the MC simulation is configured with:

- **Photon count:** $N_{\text{photons}} = 10^8$ per run.
- **Source:** Narrow Gaussian beam ($\sigma \approx 0.1$ mm) incident normal to the tissue surface at position $(x_i, +\rho(\lambda)/2, 0)$.
- **Detector (light collector):** Circular aperture of diameter 1 mm at position $(x_i, -\rho(\lambda)/2, 0)$ on the tissue surface, with numerical aperture NA = 0.22.
- **Boundary condition:** Only the top boundary ($z = 0$) is an escaping boundary (boundary type 2).
- **Deposition criterion:** `onlyCollected = true` — fluence is recorded only for photon packets that eventually reach the light collector.

### 4.3 Simulation Outputs

For each (wavelength $\lambda$, PD position $x_i$) pair, the MC simulation produces:

| Output | Symbol | Definition |
|--------|--------|-----------|
| Collected intensity | $I(\lambda, x_i)$ | Total power collected by the PD, normalized per launched photon (dimensionless) |
| Collected fluence rate | $\Phi(\mathbf{r}; \lambda, x_i)$ | Volumetric fluence rate at every voxel $\mathbf{r}$, restricted to collected photons only (cm⁻²) |
| Tissue label map | $M(\mathbf{r})$ | Integer label identifying the tissue compartment at each voxel |

### 4.4 Derived Quantities from MC Outputs

#### 4.4.1 Compartment-Integrated Fluence

For each tissue compartment $j$:

$$\Phi_j^{\text{int}}(\lambda, x_i) = \sum_{\mathbf{r} \in \text{compartment } j} \Phi(\mathbf{r}; \lambda, x_i) \cdot \Delta V$$

where $\Delta V = \frac{L_x}{n_x} \cdot \frac{L_y}{n_y} \cdot \frac{L_z}{n_z}$ is the voxel volume (cm³).

This quantity represents the total photon pathlength (weighted by collected photons) spent in compartment $j$.

#### 4.4.2 Mean Photon Pathlength per Compartment

$$\langle L_j \rangle(\lambda, x_i) = \frac{\Phi_j^{\text{int}}(\lambda, x_i)}{I(\lambda, x_i)}$$

This is the mean pathlength (cm) that detected photons spend in compartment $j$. It encodes the geometric sensitivity of each measurement to each tissue compartment.

#### 4.4.3 Jacobian (Sensitivity) of Detected Intensity to Compartment Absorption

The change in detected intensity due to a small, spatially uniform perturbation $\Delta\mu_a^{(j)}$ within compartment $j$ is:

$$\Delta I(\lambda, x_i) \approx -I(\lambda, x_i) \cdot \langle L_j \rangle(\lambda, x_i) \cdot \Delta\mu_a^{(j)}(\lambda)$$

or equivalently:

$$\frac{\partial \ln I(\lambda, x_i)}{\partial \mu_a^{(j)}(\lambda)} \approx -\langle L_j \rangle(\lambda, x_i)$$

This is the **modified Beer–Lambert law** in the multi-compartment MC context. It is valid when $\Delta\mu_a^{(j)} \cdot \langle L_j \rangle \ll 1$ (perturbation is small compared to total attenuation).

---

## 5. Signal Decomposition

### 5.1 Temporal Components

The raw photodetector signal $I(t, \lambda, x_i)$ is a continuous time series sampled at rate $f_s \geq 500$ Hz. It is decomposed into three components:

$$I(t, \lambda, x_i) = I_{DC}(\lambda, x_i) + I_{AC}^{\text{card}}(\lambda, x_i) \cdot h_c(t) + I_{AC}^{\text{resp}}(\lambda, x_i) \cdot h_r(t) + \eta(t)$$

where:

| Symbol | Definition |
|--------|-----------|
| $I_{DC}$ | Time-averaged (static) detected intensity |
| $I_{AC}^{\text{card}}$ | Amplitude of the cardiac-frequency modulation |
| $h_c(t)$ | Normalized cardiac waveform (zero-mean, unit-amplitude) |
| $I_{AC}^{\text{resp}}$ | Amplitude of the respiratory-frequency modulation |
| $h_r(t)$ | Normalized respiratory waveform (zero-mean, unit-amplitude) |
| $\eta(t)$ | Noise (shot, thermal, electronic) |

### 5.2 Extraction of Temporal Components Without External ECG/PPG

Since no time-aligned ECG or PPG reference is available, the cardiac and respiratory frequencies are extracted from the optical data itself:

1. **Power spectral density (PSD):** Compute the Welch periodogram of $I(t, \lambda, x_i)$ for each channel, using ~30 s Hanning windows with 50% overlap.

2. **Cardiac frequency identification:** The cardiac peak $f_c$ is identified as the dominant spectral peak in the band 0.8–2.5 Hz (48–150 bpm). Under general anaesthesia with mechanical ventilation, heart rate is typically stable and well-defined.

3. **Respiratory frequency identification:** The respiratory peak $f_r$ is identified in the band 0.1–0.5 Hz (6–30 breaths/min). Under mechanical ventilation, $f_r$ is fixed by the ventilator and produces a sharp, stable spectral peak.

4. **Amplitude extraction:** For each channel, the AC amplitudes are extracted via narrowband filtering (e.g., Butterworth bandpass, ±0.15 Hz around $f_c$ and $f_r$) followed by envelope detection (Hilbert transform or peak-to-trough measurement):

$$I_{AC}^{\text{card}}(\lambda, x_i) = \frac{1}{2}\left[\max_{t}(I_{\text{card-filtered}}) - \min_{t}(I_{\text{card-filtered}})\right]$$

$$I_{AC}^{\text{resp}}(\lambda, x_i) = \frac{1}{2}\left[\max_{t}(I_{\text{resp-filtered}}) - \min_{t}(I_{\text{resp-filtered}})\right]$$

5. **DC extraction:**

$$I_{DC}(\lambda, x_i) = \frac{1}{T}\int_0^T I(t, \lambda, x_i) \, dt$$

where $T$ is the acquisition window duration (recommended: 60–300 s under anaesthesia).

### 5.3 Physical Origin of Each Component

#### 5.3.1 DC Component

The static intensity reflects the time-averaged optical properties of all tissue compartments:

$$I_{DC}(\lambda, x_i) = I_0(\lambda) \cdot \exp\!\left[-\sum_j \mu_a^{(j)}(\lambda) \cdot \langle L_j \rangle(\lambda, x_i)\right]$$

where $I_0(\lambda)$ is the source intensity × coupling efficiency × detector responsivity at wavelength $\lambda$ (a calibration constant), and the sum runs over all compartments $j$ (muscle, venous blood, arterial blood if present).

This component depends on the static optical properties of all compartments: the background muscle, the venous blood (via SvO₂), and the arterial blood (via SaO₂) if present.

#### 5.3.2 Cardiac-AC Component

The cardiac-frequency modulation arises from pulsatile volume changes in blood vessels synchronized with the cardiac cycle.

**Arterial contribution:** The artery expands during systole by a fractional volume change $\Delta V_a / V_a$. This modulates the effective arterial compartment absorption:

$$\Delta I_{AC,\text{art}}^{\text{card}}(\lambda, x_i) \propto -I_{DC}(\lambda, x_i) \cdot \mu_a^{\text{art}}(\lambda) \cdot \langle L_a \rangle(\lambda, x_i) \cdot \frac{\Delta V_a}{V_a}$$

**Venous contribution:** If the vein is pulsatile (particularly relevant for the IJV):

$$\Delta I_{AC,\text{ven}}^{\text{card}}(\lambda, x_i) \propto -I_{DC}(\lambda, x_i) \cdot \mu_a^{\text{ven}}(\lambda) \cdot \langle L_v \rangle(\lambda, x_i) \cdot \frac{\Delta V_v}{V_v}$$

The total cardiac-AC signal is the superposition:

$$I_{AC}^{\text{card}}(\lambda, x_i) = A_a(\lambda, x_i) \cdot \frac{\Delta r_a}{r_a} + A_v(\lambda, x_i) \cdot \frac{\Delta r_v}{r_v}$$

where we define the **compartment-specific AC sensitivity coefficients**:

$$A_j(\lambda, x_i) \equiv I_{DC}(\lambda, x_i) \cdot \mu_a^{(j)}(\lambda) \cdot \langle L_j \rangle(\lambda, x_i)$$

and $\Delta r / r$ is the fractional radius change of the vessel (related to fractional volume change by $\Delta V / V \approx 2 \Delta r / r$ for a cylinder).

**Important:** For the forearm vein site (no nearby artery), the cardiac-AC signal is expected to be weak and dominated by capillary-bed pulsation in the muscle background. For the neck (IJV + carotid), both arterial and venous components are significant.

#### 5.3.3 Respiratory-AC Component

The respiratory-frequency modulation arises primarily from venous volume changes driven by intrathoracic pressure variation during the respiratory cycle:

$$I_{AC}^{\text{resp}}(\lambda, x_i) \approx A_v(\lambda, x_i) \cdot \frac{\Delta r_v^{\text{resp}}}{r_v}$$

Respiratory modulation of arterial diameter is typically small and can be neglected to first order. Under mechanical ventilation, the respiratory venous modulation is larger and more regular than during spontaneous breathing. This makes the respiratory-AC component the **cleanest venous signal**, largely free of arterial contamination.

This is the primary channel for SvO₂ recovery.

---

## 6. Inverse Problem

### 6.1 Data Vector

The complete measurement vector for a single acquisition at one anatomical site is:

$$\mathbf{d} = \begin{bmatrix} \mathbf{d}_{DC} \\ \mathbf{d}_{AC}^{\text{card}} \\ \mathbf{d}_{AC}^{\text{resp}} \end{bmatrix} \in \mathbb{R}^{N_{\text{total}}}$$

where:

$$\mathbf{d}_{DC} = \big[\ln I_{DC}(\lambda_k, x_i)\big], \quad k = 1 \ldots 3, \quad i = 1 \ldots 12$$

$$\mathbf{d}_{AC}^{\text{card}} = \left[\frac{I_{AC}^{\text{card}}(\lambda_k, x_i)}{I_{DC}(\lambda_k, x_i)}\right]$$

$$\mathbf{d}_{AC}^{\text{resp}} = \left[\frac{I_{AC}^{\text{resp}}(\lambda_k, x_i)}{I_{DC}(\lambda_k, x_i)}\right]$$

The DC data is taken in log-space to linearize the exponential dependence. The AC data is normalized by DC to cancel the source/detector coupling constants $I_0(\lambda)$. The AC/DC ratio is often called the **modulation ratio**.

The total number of data points is:

$$N_{\text{total}} = 3 \times 36 = 108$$

(36 DC + 36 cardiac-AC + 36 respiratory-AC)

### 6.2 Parameter Vector

The unknown parameters to be recovered are:

$$\mathbf{p} = \begin{bmatrix} \mu_a^{\text{musc}}(\lambda_1) \\ \mu_a^{\text{musc}}(\lambda_2) \\ \mu_a^{\text{musc}}(\lambda_3) \\ a'_{\text{musc}} \\ b_{\text{musc}} \\ \text{SvO}_2 \\ [\text{Hb}] \\ \Delta r_a / r_a \\ \Delta r_v^{\text{card}} / r_v \\ \Delta r_v^{\text{resp}} / r_v \\ \ln I_0(\lambda_1) \\ \ln I_0(\lambda_2) \\ \ln I_0(\lambda_3) \end{bmatrix}$$

where:

| Parameter | Count | Notes |
|-----------|-------|-------|
| $\mu_a^{\text{musc}}(\lambda_k)$ | 3 | Muscle absorption at each wavelength |
| $a'_{\text{musc}}, b_{\text{musc}}$ | 2 | Muscle scattering power-law parameters |
| SvO₂ | 1 | **Primary target** |
| [Hb] | 1 | Total hemoglobin concentration in blood |
| $\Delta r_a / r_a$ | 1 | Arterial fractional radius change (cardiac); set to 0 for vein-only site |
| $\Delta r_v^{\text{card}} / r_v$ | 1 | Venous fractional radius change (cardiac); set to 0 if vein is non-pulsatile |
| $\Delta r_v^{\text{resp}} / r_v$ | 1 | Venous fractional radius change (respiratory) |
| $\ln I_0(\lambda_k)$ | 3 | Per-wavelength coupling/calibration constants |

Total: **13 parameters** for the artery + vein case. For the simple single-vein case, $\Delta r_a / r_a = 0$ and the artery compartment is absent, reducing to **12 parameters** (or fewer if cardiac venous pulsation is also neglected).

### 6.3 Fixed (Known) Inputs

| Quantity | Source |
|----------|--------|
| SaO₂ | Clinical monitoring (pulse oximetry or arterial line) |
| Vessel geometry ($r_v, d_v, x_v, r_a, d_a, x_a$) | B-mode ultrasound |
| Blood $\mu_s'(\lambda)$ | Literature values |
| Extinction coefficients $\varepsilon_{\text{HbO}_2}(\lambda), \varepsilon_{\text{HbR}}(\lambda)$ | Prahl/Gratzer tabulated spectra |
| Water absorption $\mu_{a,\text{water}}(\lambda)$ | Hale and Querry / Kou et al. tabulated spectra |

### 6.4 Forward Model Evaluation

The forward model $\mathbf{f}(\mathbf{p})$ maps parameters to predicted data:

**DC predictions:**

$$f_{DC}(\lambda_k, x_i; \mathbf{p}) = \ln I_0(\lambda_k) - \sum_j \mu_a^{(j)}(\lambda_k; \mathbf{p}) \cdot \langle L_j \rangle(\lambda_k, x_i)$$

where the compartment mean pathlengths $\langle L_j \rangle$ are pre-computed from the MC simulation at a reference tissue state (see Section 4.4.2).

**AC predictions (modulation ratios):**

$$f_{AC}^{\text{resp}}(\lambda_k, x_i; \mathbf{p}) = \mu_a^{\text{ven}}(\lambda_k; \text{SvO}_2, [\text{Hb}]) \cdot \langle L_v \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_v^{\text{resp}}}{r_v}$$

$$f_{AC}^{\text{card}}(\lambda_k, x_i; \mathbf{p}) = \mu_a^{\text{art}}(\lambda_k; \text{SaO}_2, [\text{Hb}]) \cdot \langle L_a \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_a}{r_a} + \mu_a^{\text{ven}}(\lambda_k; \text{SvO}_2, [\text{Hb}]) \cdot \langle L_v \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_v^{\text{card}}}{r_v}$$

**Venous absorption dependence on SvO₂:**

$$\mu_a^{\text{ven}}(\lambda; \text{SvO}_2, [\text{Hb}]) = \frac{\ln(10) \cdot [\text{Hb}]}{MW_{\text{Hb}}} \left[\text{SvO}_2 \cdot \varepsilon_{\text{HbO}_2}(\lambda) + (1 - \text{SvO}_2) \cdot \varepsilon_{\text{HbR}}(\lambda)\right] + \mu_{a,\text{water}}(\lambda) \cdot W_{\text{blood}}$$

**Arterial absorption (known spectral shape):**

$$\mu_a^{\text{art}}(\lambda; \text{SaO}_2, [\text{Hb}]) = \frac{\ln(10) \cdot [\text{Hb}]}{MW_{\text{Hb}}} \left[\text{SaO}_2 \cdot \varepsilon_{\text{HbO}_2}(\lambda) + (1 - \text{SaO}_2) \cdot \varepsilon_{\text{HbR}}(\lambda)\right] + \mu_{a,\text{water}}(\lambda) \cdot W_{\text{blood}}$$

### 6.5 Cost Function

The inverse problem is formulated as weighted nonlinear least squares:

$$\hat{\mathbf{p}} = \arg\min_{\mathbf{p}} \; \chi^2(\mathbf{p})$$

$$\chi^2(\mathbf{p}) = \sum_{\alpha} w_\alpha \left[d_\alpha - f_\alpha(\mathbf{p})\right]^2 + R(\mathbf{p})$$

where:

| Symbol | Definition |
|--------|-----------|
| $\alpha$ | Index over all data points (DC, cardiac-AC, respiratory-AC at each $\lambda_k, x_i$) |
| $w_\alpha$ | Weight for data point $\alpha$ (typically $1/\sigma_\alpha^2$ where $\sigma_\alpha$ is the measurement uncertainty) |
| $R(\mathbf{p})$ | Regularization term (optional, for ill-conditioned problems) |

A Tikhonov regularization term can be included:

$$R(\mathbf{p}) = \gamma \|\mathbf{p} - \mathbf{p}_{\text{prior}}\|^2$$

where $\mathbf{p}_{\text{prior}}$ is a vector of prior estimates (e.g., literature values for muscle properties) and $\gamma$ controls the regularization strength.

### 6.6 Solver

The minimization is performed using the **Levenberg–Marquardt** algorithm (or trust-region reflective for bound-constrained problems). The Jacobian $\partial f_\alpha / \partial p_k$ is computed by finite differences or analytically from the modified Beer–Lambert expressions.

Bounds on parameters:

| Parameter | Lower | Upper |
|-----------|-------|-------|
| $\mu_a^{\text{musc}}$ | 0.001 cm⁻¹ | 1.0 cm⁻¹ |
| $a'_{\text{musc}}$ | 10 cm⁻¹ | 80 cm⁻¹ |
| $b_{\text{musc}}$ | 0.5 | 2.0 |
| SvO₂ | 0.30 | 1.00 |
| [Hb] | 1.0 mM | 3.5 mM |
| $\Delta r / r$ | 0 | 0.20 |
| $\ln I_0$ | unconstrained | unconstrained |

Initial guesses:

| Parameter | Initial Value | Source |
|-----------|---------------|--------|
| $\mu_a^{\text{musc}}$ | Jacques model with $B = 0.002, S = 0.67, W = 0.65$ | Literature |
| $a', b$ | 42.4, 1.0 | Jacques (2013) |
| SvO₂ | 0.70 | Typical resting venous saturation |
| [Hb] | 2.3 mM | Normal adult |
| $\Delta r_a / r_a$ | 0.05 | Approximately 5% arterial pulsation |
| $\Delta r_v / r_v$ | 0.02 (cardiac), 0.10 (respiratory) | Literature estimates |
| $I_0(\lambda)$ | From bulk-muscle-only MC simulation | Calibration run |

---

## 7. Site-Specific Considerations

### 7.1 Forearm (Isolated Vein)

- **Geometry:** Single superficial vein (cephalic), no major artery.
- **Compartments:** Muscle + vein (2 compartments).
- **Cardiac-AC:** Weak; dominated by capillary pulsation in muscle, not a discrete artery. Set $\Delta r_a / r_a = 0$.
- **Respiratory-AC:** Expected to be present but small. This is the primary venous channel.
- **DC contrast:** The spatial pattern of DC contrast across the 12 PDs encodes the lateral position and depth of the vein.
- **This is the simplest case** and the primary validation target.

### 7.2 Wrist (Radial Artery + Vein)

- **Geometry:** Radial artery and accompanying vein (cephalic or comitantes), laterally separated by approximately 2–4 mm.
- **Compartments:** Muscle + vein + artery (3 compartments).
- **Cardiac-AC:** Strong arterial signal. SaO₂ is known, so this constrains [Hb] and validates the model.
- **Respiratory-AC:** Primarily venous.
- **Spatial information:** The 12 PDs at different lateral offsets see different mixtures of arterial and venous signal, providing geometric leverage for separating the two.

### 7.3 Neck (IJV + Carotid)

- **Geometry:** Large IJV and common carotid artery, deeper than wrist/forearm vessels.
- **Compartments:** Skin + subcutaneous tissue + muscle + IJV + carotid (5 compartments recommended).
- **Cardiac-AC:** Both arterial and venous components are significant. The IJV is pulsatile at the cardiac frequency (a, c, v waves).
- **Respiratory-AC:** Large and clean (IJV is highly respiration-sensitive). This is the dominant venous channel.
- **Posture dependence:** IJV size varies dramatically with head-of-bed elevation. Must be standardized and recorded.
- **Venous pulsatility concern:** The IJV pulsates at the cardiac frequency, meaning the assumption "cardiac-AC = arterial only" is invalid at this site.

---

## 8. Handling Venous Pulsatility at the Cardiac Frequency

### 8.1 The Problem

At sites where the vein is pulsatile (particularly the IJV), the cardiac-AC signal contains both arterial and venous contributions:

$$I_{AC}^{\text{card}}(\lambda, x_i) = A_a(\lambda, x_i) \cdot \frac{\Delta r_a}{r_a} + A_v(\lambda, x_i) \cdot \frac{\Delta r_v^{\text{card}}}{r_v}$$

The naive assumption that this is purely arterial would yield a biased estimate of SaO₂ (if used to validate) or would incorrectly absorb venous signal into the arterial model.

### 8.2 Mitigation Strategies

**Primary strategy — respiratory-AC as the venous channel:** The respiratory-frequency modulation is predominantly venous (arterial respiratory modulation is small). By using the respiratory-AC as the primary SvO₂ channel and the cardiac-AC only as a secondary constraint, the venous pulsatility at the cardiac frequency becomes a nuisance term rather than a confound.

**Secondary strategy — spatial separation across the PD array:** Since the artery and vein are at different lateral positions (known from ultrasound), PDs at different x-positions see different ratios of arterial-to-venous sensitivity:

$$\frac{A_a(\lambda, x_i)}{A_v(\lambda, x_i)} = \frac{\mu_a^{\text{art}}(\lambda) \cdot \langle L_a \rangle(\lambda, x_i)}{\mu_a^{\text{ven}}(\lambda) \cdot \langle L_v \rangle(\lambda, x_i)}$$

This ratio varies across $x_i$ because $\langle L_a \rangle$ and $\langle L_v \rangle$ have different spatial profiles (they peak at different x-positions). The 12 PDs therefore provide 12 different mixtures of arterial and venous cardiac-AC signal, enabling decomposition via the forward model.

**Tertiary strategy — self-derived waveform analysis:** The arterial and venous cardiac waveforms have different morphologies (arterial: sharp systolic upstroke; venous IJV: a/c/v wave triplet). Even without external ECG, the dominant cardiac frequency and its harmonics can be extracted from the optical data. The ratio of harmonic amplitudes differs between arterial and venous waveforms, providing additional separability. This requires sufficient sampling rate ($\geq 500$ Hz) and SNR to resolve the first 3–4 cardiac harmonics.

---

## 9. Simulation Validation Pipeline

### 9.1 Phase 1: Bulk Muscle Validation (Current)

**Objective:** Validate the MC forward model against analytical diffusion theory for a homogeneous medium.

**Procedure:**

1. Run MC simulation for bulk muscle at each wavelength with a single PD/LED position.
2. Compute the analytical CW diffusion reflectance (Farrell/Patterson/Wilson solution) at the same $\rho$:

$$R(\rho) = \frac{1}{4\pi}\left[ z_0\!\left(\mu_{\text{eff}} + \frac{1}{r_1}\right)\frac{e^{-\mu_{\text{eff}} r_1}}{r_1^2} + (z_0 + 2z_b)\!\left(\mu_{\text{eff}} + \frac{1}{r_2}\right)\frac{e^{-\mu_{\text{eff}} r_2}}{r_2^2} \right]$$

where:

| Symbol | Definition |
|--------|-----------|
| $D$ | Diffusion coefficient: $D = \frac{1}{3(\mu_a + \mu_s')}$ |
| $\mu_{\text{eff}}$ | Effective attenuation coefficient: $\mu_{\text{eff}} = \sqrt{\mu_a / D} = \sqrt{3\mu_a(\mu_a + \mu_s')}$ |
| $z_0$ | Isotropic source depth: $z_0 = 1/\mu_s'$ |
| $z_b$ | Extrapolated boundary distance: $z_b = 2D \cdot \frac{1 + R_{\text{eff}}}{1 - R_{\text{eff}}}$ |
| $R_{\text{eff}}$ | Effective reflection coefficient at tissue–air boundary (approximately 0.493 for $n = 1.4$) |
| $r_1$ | Distance from real source to detection point: $r_1 = \sqrt{\rho^2 + z_0^2}$ |
| $r_2$ | Distance from image source to detection point: $r_2 = \sqrt{\rho^2 + (z_0 + 2z_b)^2}$ |

3. Compare MC-simulated collected intensity to Farrell $R(\rho)$. Agreement within a few percent validates the MC setup.
4. Compute and store $I_0(\lambda)$ and $\langle L_{\text{musc}} \rangle(\lambda)$ as reference values.

**Note on single-ρ limitation:** With one source–detector separation per wavelength, $\mu_a$ and $\mu_s'$ cannot be independently recovered from CW intensity alone. The bulk muscle simulation serves as a forward-model consistency check (using known input properties) and as a calibration baseline, not as a tissue characterization step.

### 9.2 Phase 2: Single Vein in Muscle — Static (DC) Inverse Problem (Current)

**Objective:** Test recovery of SvO₂ from DC measurements alone, using a single vein embedded in bulk muscle.

**Procedure:**

1. Run MC simulation with vein (SvO₂ = 0.75) at all 12 PD positions × 3 wavelengths.
2. Extract collected intensity $I(\lambda_k, x_i)$ and compartment-integrated fluences $\Phi_j^{\text{int}}(\lambda_k, x_i)$ for $j \in \{\text{muscle}, \text{vein}\}$.
3. Compute vessel contrast:

$$C(\lambda_k, x_i) = \frac{I_{\text{vessel}}(\lambda_k, x_i) - I_{\text{bulk}}(\lambda_k)}{I_{\text{bulk}}(\lambda_k)}$$

4. Construct the DC-only forward model using the modified Beer–Lambert formulation (Section 6.4).
5. Run the inverse solver on the 36 DC measurements to recover SvO₂, [Hb], muscle $\mu_a(\lambda)$, and $I_0(\lambda)$.
6. Compare recovered SvO₂ to the true value (0.75).
7. Add synthetic noise ($\sigma = 0.01$ on $\ln I$, corresponding to approximately 1% intensity noise) and re-run to characterize solver variance.

### 9.3 Phase 3: Single Vein — Perturbation Runs for AC Signals (Upcoming)

**Objective:** Generate simulated cardiac and respiratory AC data.

**Procedure:**

1. Re-run MC with vein radius increased by 5% (cardiac systole equivalent): $r_v' = 1.05 \cdot r_v$.
2. Re-run MC with vein radius increased by 15% (respiratory inspiration equivalent): $r_v' = 1.15 \cdot r_v$.
3. Compute AC amplitudes as the difference between perturbed and baseline intensities.
4. Run the full inverse solver (DC + AC) on the complete 108-element data vector.

### 9.4 Phase 4: Vein + Artery Geometry (Upcoming)

**Objective:** Test SvO₂ recovery in the presence of a known artery, using spatial separation across PDs and the known SaO₂ as constraints.

### 9.5 Phase 5: Robustness Testing

**Objective:** Characterize sensitivity to model errors and measurement noise.

**Sweep parameters:**

| Parameter | Sweep Range | Tests |
|-----------|------------|-------|
| SvO₂ (true) | 0.60 – 0.85 | Recovery accuracy across anaesthesia range |
| SaO₂ (true) | 0.95 – 1.00 | Sensitivity to SaO₂ uncertainty |
| [Hb] | 1.5 – 3.0 mM | Sensitivity to hematocrit variation |
| Vessel depth | $d_v \pm 1$ mm | Sensitivity to ultrasound depth error |
| Vessel diameter | $r_v \pm 20\%$ | Sensitivity to ultrasound diameter error |
| Lateral offset | $x_v \pm 1$ mm | Sensitivity to sensor co-registration error |
| Noise level | 0.5%, 1%, 2%, 5% | Solver noise floor |

---

## 10. Linearization for SvO₂ Sweeps Without Re-Running MC

To test the inverse problem at multiple SvO₂ values without re-running MC for each, the following linearization is used:

1. Run the full MC at a reference state ($\text{SvO}_2^{\text{ref}} = 0.75$) and extract the baseline collected intensity $I^{\text{ref}}(\lambda_k, x_i)$ and compartment pathlengths $\langle L_v \rangle(\lambda_k, x_i)$.

2. For a hypothetical SvO₂ different from 0.75, predict the intensity as:

$$\ln I(\lambda_k, x_i; \text{SvO}_2) \approx \ln I^{\text{ref}}(\lambda_k, x_i) - \langle L_v \rangle(\lambda_k, x_i) \cdot \left[\mu_a^{\text{ven}}(\lambda_k; \text{SvO}_2) - \mu_a^{\text{ven}}(\lambda_k; \text{SvO}_2^{\text{ref}})\right]$$

3. This linearization is valid when $\Delta\mu_a^{\text{ven}} \cdot \langle L_v \rangle \ll 1$. For typical vessel pathlengths of $\langle L_v \rangle \sim 0.01$–$0.1$ cm and $\Delta\mu_a^{\text{ven}} \sim 1$–$5$ cm⁻¹ (for a ±15% SvO₂ swing), the product is approximately 0.01–0.5. The lower end is well within the linear regime; the upper end should be validated against a full MC re-run.

---

## 11. Experimental Protocol Summary

### 11.1 Pre-Measurement

1. Record subject demographics (age, sex, BMI).
2. Place sensor array over target site (forearm, wrist, or neck).
3. Acquire B-mode ultrasound images of the measurement site with the sensor in place (or with fiducial markers on the skin). Record vessel diameters, depths, lateral positions, and separation.
4. Establish invasive SvO₂ monitoring (catheter placement anatomically matched to sensor site).
5. Record SaO₂ from clinical monitoring.
6. Record subject posture (especially for neck site: head-of-bed angle).

### 11.2 Data Acquisition

1. Pulse each LED independently at rate $\geq 500$ Hz per complete cycle through all 36 LED/PD combinations.
2. Interleave dark frames (all LEDs off) for ambient/dark subtraction.
3. Record raw photodetector voltages with $\geq 16$ bit resolution.
4. Acquire for $\geq 60$ s (recommended: 120–300 s) per measurement epoch.
5. Record timestamps synchronized with SvO₂ catheter readings and SaO₂ readings.
6. Repeat ultrasound at end of acquisition to detect any sensor/vessel drift.

### 11.3 Offline Processing

1. Dark-subtract all channels.
2. Extract DC, cardiac-AC, and respiratory-AC per Section 5.
3. Construct the forward model using ultrasound geometry and pre-computed MC pathlengths.
4. Solve the inverse problem per Section 6.
5. Compare recovered SvO₂ to catheter ground truth.

---

## 12. Notation Summary

| Symbol | Definition | Units |
|--------|-----------|-------|
| $\lambda$ | Optical wavelength | nm |
| $\rho$ | Source–detector separation | cm |
| $x_i$ | Lateral position of the $i$-th PD along the array | cm |
| $\mu_a$ | Absorption coefficient | cm⁻¹ |
| $\mu_s$ | Scattering coefficient | cm⁻¹ |
| $\mu_s'$ | Reduced scattering coefficient: $\mu_s' = \mu_s(1-g)$ | cm⁻¹ |
| $\mu_{\text{eff}}$ | Effective attenuation coefficient: $\sqrt{3\mu_a(\mu_a + \mu_s')}$ | cm⁻¹ |
| $g$ | Scattering anisotropy | — |
| $n$ | Refractive index | — |
| $D$ | Diffusion coefficient: $1/[3(\mu_a + \mu_s')]$ | cm |
| $\Phi(\mathbf{r})$ | Fluence rate | cm⁻² (per launched photon) |
| $I$ | Detected intensity (collected power) | per launched photon |
| $I_0$ | Source/detector coupling constant | — |
| $\langle L_j \rangle$ | Mean photon pathlength in compartment $j$ | cm |
| $\varepsilon(\lambda)$ | Molar extinction coefficient | cm⁻¹ M⁻¹ |
| [Hb] | Total hemoglobin concentration | mM |
| SvO₂ | Venous oxygen saturation | — |
| SaO₂ | Arterial oxygen saturation | — |
| $B$ | Blood volume fraction | — |
| $S$ | Oxygen saturation | — |
| $W$ | Water volume fraction | — |
| $r_v, r_a$ | Vessel radius (vein, artery) | cm |
| $d_v, d_a$ | Vessel center depth | cm |
| $\Delta r / r$ | Fractional radius change (pulsatile) | — |
| $f_c$ | Cardiac frequency | Hz |
| $f_r$ | Respiratory frequency | Hz |
| $A_j$ | AC sensitivity coefficient for compartment $j$ | — |
| $\Phi_j^{\text{int}}$ | Compartment-integrated fluence | cm |
| $\Delta V$ | Voxel volume | cm³ |
| $N_{\text{photons}}$ | Number of simulated photon packets | — |
| $w_\alpha$ | Data point weight in cost function | — |
| $\gamma$ | Regularization strength | — |
| $\chi^2$ | Cost function value | — |

---

## 13. Key Assumptions and Limitations

1. **Single ρ per wavelength.** The sensor does not provide multi-distance measurements. Separation of $\mu_a$ from $\mu_s'$ relies on spectral constraints and/or assumed $\mu_s'$ from literature, not on the ρ-dependence of R(ρ).

2. **Diffusion theory is used only for analytical cross-checks.** The forward model for vessel-containing geometries is MC-based (no P1 approximation).

3. **Modified Beer–Lambert linearization.** The perturbation forward model assumes $\Delta\mu_a \cdot \langle L \rangle \ll 1$. This should be validated against full MC re-runs for the expected range of SvO₂.

4. **Blood μₛ′ is assumed from the literature**, not recovered from data. Errors in this assumption propagate into [Hb] and SvO₂ estimates.

5. **Ultrasound geometry is assumed accurate.** Errors in vessel depth, diameter, or lateral position directly affect the MC-derived pathlengths and therefore the recovered SvO₂.

6. **Wavelength set (660, 850, 940 nm) is suboptimal for SvO₂.** Maximum HbO₂/HbR discrimination occurs near 700–750 nm. The absence of a wavelength in this range limits the spectral conditioning of the SvO₂ recovery. Expected accuracy: ±5–10% SvO₂ under ideal conditions.

7. **Tissue is modeled as muscle + vessel(s) in the simplest case.** Skin, fat, fascia, and tendon layers are omitted. These layers affect short-ρ sensitivity and should be included for anatomically accurate forward models at the wrist and neck.

8. **Anaesthesia.** The study is conducted under general anaesthesia with mechanical ventilation. This provides stable heart rate, fixed respiratory rate, elevated baseline SvO₂ (reduced metabolic demand), and potentially altered vascular tone. Results may not directly translate to awake, spontaneously breathing subjects.

9. **No ECG/PPG reference.** Cardiac and respiratory frequencies are extracted from the optical data itself. Phase-locked averaging is not available; frequency-domain methods are used instead.

10. **Venous pulsatility at the cardiac frequency** is a known confound, particularly at the IJV. The respiratory-AC channel is used as the primary venous signal to mitigate this.
