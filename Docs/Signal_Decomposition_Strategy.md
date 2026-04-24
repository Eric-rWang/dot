# Signal Decomposition Strategy for SvO₂ Recovery at Three Anatomical Sites

---

## 1. Overview of the Three Signal Components

Every photodetector in our array records a continuous time series of detected light intensity. This raw signal contains information at multiple timescales, each originating from a different physiological mechanism. We decompose the signal into three components and use each to constrain different unknowns in the inverse problem.

### 1.1 DC Component

The DC component is the time-averaged detected intensity over the acquisition window (60–300 seconds):

$$I_{DC}(\lambda, x_i) = \frac{1}{T} \int_0^T I(t, \lambda, x_i) \, dt$$

**Physical origin:** The DC intensity reflects the static optical properties of all tissue compartments — muscle, venous blood, arterial blood (if present), and any other layers (skin, fat). It depends on the absorption and scattering coefficients of every compartment, the geometry (vessel sizes, depths, positions), and the source–detector coupling.

**What it encodes:** The DC signal contains the total, time-averaged attenuation. From the two-compartment MBLL:

$$\ln I_{DC}(\lambda, x_i) = \ln I_0^{\text{eff}}(\lambda) - \mu_a^{(m)}(\lambda) \cdot \langle L_m \rangle(\lambda, x_i) - \mu_a^{(v)}(\lambda) \cdot \langle L_v \rangle(\lambda, x_i)$$

For the three-compartment case (muscle + vein + artery):

$$\ln I_{DC}(\lambda, x_i) = \ln I_0^{\text{eff}}(\lambda) - \mu_a^{(m)}(\lambda) \cdot \langle L_m \rangle - \mu_a^{\text{ven}}(\lambda) \cdot \langle L_v \rangle - \mu_a^{\text{art}}(\lambda) \cdot \langle L_a \rangle$$

The DC signal mixes contributions from all compartments. Separating them requires either spatial information (different PD positions see different compartment mixtures) or spectral information (different wavelengths weight compartments differently), or both.

**Strengths:** High SNR (time-averaging reduces noise). Available at all sites. Contains geometric information through the spatial profile across the 12 PDs.

**Weaknesses:** Depends on the unknown calibration constants $I_0^{\text{eff}}(\lambda)$ (three unknowns, one per wavelength). Mixes muscle and blood contributions, making it difficult to isolate the vessel signal from the background.

### 1.2 Cardiac-AC Component

The cardiac-AC component is the amplitude of the intensity modulation at the cardiac frequency ($f_c \approx 0.8$–$2.0$ Hz):

$$I_{AC}^{\text{card}}(\lambda, x_i) = \text{amplitude of } I(t, \lambda, x_i) \text{ at frequency } f_c$$

**Physical origin:** Blood vessels expand and contract with each heartbeat. During systole, the artery's lumen expands (increasing the effective arterial blood volume in the photon path), and the venous lumen may also change (particularly in compliant veins like the IJV). This pulsatile volume change modulates the effective absorption along the photon path.

**Mathematical formulation:** When a vessel's radius changes by $\Delta r$ from its mean value $r$, the vessel's cross-sectional area changes, which changes the effective absorption experienced by photons traversing it. To first order, the fractional change in detected intensity is:

$$\frac{I_{AC}^{\text{card}}(\lambda, x_i)}{I_{DC}(\lambda, x_i)} = \mu_a^{\text{art}}(\lambda) \cdot \langle L_a \rangle(\lambda, x_i) \cdot \frac{\Delta r_a}{r_a} + \mu_a^{\text{ven}}(\lambda) \cdot \langle L_v \rangle(\lambda, x_i) \cdot \frac{\Delta r_v^{\text{card}}}{r_v}$$

The left side is the **modulation ratio** — the AC amplitude normalized by the DC level. This normalization cancels the calibration constants $I_0(\lambda)$, which is a major advantage over using DC alone.

The right side is a sum of two terms: one from arterial pulsation and one from venous pulsation. Each term is the product of three factors: the blood absorption at that wavelength, the mean pathlength in that compartment, and the fractional radius change of that vessel.

**Strengths:** The modulation ratio is self-calibrating — dividing by $I_{DC}$ cancels $I_0(\lambda)$, eliminating three unknowns. Contains spectral information about blood absorption directly. At sites with a dominant artery, provides a strong reference signal.

**Weaknesses:** If both artery and vein pulsate at the cardiac frequency (as at the neck), the two contributions are mixed and must be separated. Requires sufficient sampling rate (≥500 Hz) to resolve the cardiac waveform. Sensitive to motion artifacts (any sensor movement produces a large artifact at cardiac-like frequencies).

### 1.3 Respiratory-AC Component

The respiratory-AC component is the amplitude of the intensity modulation at the respiratory frequency ($f_r \approx 0.1$–$0.4$ Hz):

$$I_{AC}^{\text{resp}}(\lambda, x_i) = \text{amplitude of } I(t, \lambda, x_i) \text{ at frequency } f_r$$

**Physical origin:** During the respiratory cycle, intrathoracic pressure changes modulate venous return. During mechanical inspiration (positive-pressure ventilation, as used in our anesthetized subjects), intrathoracic pressure increases, impeding venous return and causing peripheral veins to distend. During expiration, pressure drops and veins partially collapse. This produces a rhythmic change in venous blood volume at the respiratory frequency.

Arterial diameter is minimally affected by respiration (the arterial wall is much stiffer than the venous wall, and arterial pressure is much higher than the respiratory pressure swing). Therefore, the respiratory modulation is predominantly venous.

**Mathematical formulation:**

$$\frac{I_{AC}^{\text{resp}}(\lambda, x_i)}{I_{DC}(\lambda, x_i)} \approx \mu_a^{\text{ven}}(\lambda) \cdot \langle L_v \rangle(\lambda, x_i) \cdot \frac{\Delta r_v^{\text{resp}}}{r_v}$$

Note: the arterial term is absent (or negligible). This is the key advantage — the respiratory-AC modulation ratio is a direct measurement of venous blood absorption, uncontaminated by arterial blood.

**Strengths:** Predominantly venous — provides the cleanest channel for SvO₂. Self-calibrating (modulation ratio cancels $I_0$). Under mechanical ventilation (as in our anesthetized subjects), the respiratory frequency is fixed and stable, producing a clean, reproducible spectral peak. The respiratory modulation of venous volume is typically larger than the cardiac modulation of venous volume (10–20% radius change vs. 2–5%).

**Weaknesses:** Lower frequency means fewer cycles per acquisition window (at 12 breaths/min, a 60-second window contains only 12 cycles — marginal for spectral estimation; 120–300 seconds is preferred). Smaller absolute signal than cardiac-AC (because the respiratory frequency is lower and the pulsation is slower). At the forearm, respiratory venous modulation may be very small.

---

## 2. Extraction of the Three Components Without ECG/PPG

Since no external ECG or PPG reference signal is available, we extract the cardiac and respiratory frequencies from the optical data itself.

### 2.1 Step 1: Identify the Frequencies

Compute the power spectral density (PSD) of each channel using Welch's method (Hanning windows, ~30 seconds, 50% overlap). Under general anaesthesia with mechanical ventilation:

**Respiratory peak ($f_r$):** The ventilator sets a fixed respiratory rate (typically 12–16 breaths/min = 0.2–0.27 Hz). This produces a sharp, narrow spectral peak that is easily identified. All channels should show the peak at the same frequency (it is a global physiological signal). Average the PSDs across channels to improve detection.

**Cardiac peak ($f_c$):** Heart rate under anaesthesia is typically 50–80 bpm (0.83–1.33 Hz). It is more variable than the respiratory rate but still produces a clear spectral peak. The cardiac peak is usually the dominant feature in the 0.5–2.5 Hz band. Channels with strong arterial sensitivity (those near the artery, at the wrist and neck sites) will show the largest cardiac peak.

### 2.2 Step 2: Extract Amplitudes

**Bandpass filtering:** For each channel, apply a narrow bandpass filter centered on the identified frequency:

- Cardiac band: $f_c \pm 0.15$ Hz (captures the fundamental and allows for slight HR variability)
- Respiratory band: $f_r \pm 0.05$ Hz (narrower, because ventilator rate is fixed)

Use a 4th-order Butterworth filter (zero-phase, applied forward and backward to avoid phase distortion).

**Amplitude estimation:** After filtering, the signal is approximately sinusoidal. The amplitude is estimated as half the peak-to-trough range, averaged over all complete cycles in the window:

$$I_{AC} = \frac{1}{2N} \sum_{n=1}^{N} \left[\max_n(I_{\text{filtered}}) - \min_n(I_{\text{filtered}})\right]$$

where $N$ is the number of complete cycles.

Alternatively, use the Hilbert transform to compute the instantaneous envelope, then take the mean envelope amplitude:

$$I_{AC} = \frac{1}{T} \int_0^T |I_{\text{analytic}}(t)| \, dt$$

where $I_{\text{analytic}}(t) = I_{\text{filtered}}(t) + i \mathcal{H}\{I_{\text{filtered}}\}(t)$ is the analytic signal.

### 2.3 Step 3: Compute Modulation Ratios

For each channel ($\lambda_k$, $x_i$), compute the modulation ratios:

$$M_{\text{card}}(\lambda_k, x_i) = \frac{I_{AC}^{\text{card}}(\lambda_k, x_i)}{I_{DC}(\lambda_k, x_i)}$$

$$M_{\text{resp}}(\lambda_k, x_i) = \frac{I_{AC}^{\text{resp}}(\lambda_k, x_i)}{I_{DC}(\lambda_k, x_i)}$$

These modulation ratios are the primary data for the AC portion of the inverse problem.

---

## 3. The Forward Model for Each Component

### 3.1 Forward Model: DC

$$\ln I_{DC}(\lambda_k, x_i) = \ln I_0^{\text{eff}}(\lambda_k) - \sum_j \mu_a^{(j)}(\lambda_k) \cdot \langle L_j \rangle(\lambda_k, x_i)$$

where $j$ runs over all compartments (muscle, vein, artery if present), and the pathlengths $\langle L_j \rangle$ are pre-computed from Monte Carlo simulation.

**Unknowns contributed by DC:** muscle $\mu_a$ at each wavelength (3), calibration constants $\ln I_0^{\text{eff}}$ (3), SvO₂ and [Hb] (which determine venous $\mu_a$), SaO₂ (known) and [Hb] (which determine arterial $\mu_a$).

### 3.2 Forward Model: Cardiac-AC Modulation Ratio

$$M_{\text{card}}(\lambda_k, x_i) = \mu_a^{\text{art}}(\lambda_k) \cdot \langle L_a \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_a}{r_a} + \mu_a^{\text{ven}}(\lambda_k) \cdot \langle L_v \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_v^{\text{card}}}{r_v}$$

**Unknowns contributed by cardiac-AC:** arterial fractional pulsation $\Delta r_a / r_a$ (1), venous cardiac fractional pulsation $\Delta r_v^{\text{card}} / r_v$ (1), plus the blood absorption terms (which depend on SvO₂, SaO₂, [Hb]).

Note: $I_0$ does not appear — it cancels in the modulation ratio.

### 3.3 Forward Model: Respiratory-AC Modulation Ratio

$$M_{\text{resp}}(\lambda_k, x_i) = \mu_a^{\text{ven}}(\lambda_k) \cdot \langle L_v \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_v^{\text{resp}}}{r_v}$$

**Unknowns contributed by respiratory-AC:** venous respiratory fractional pulsation $\Delta r_v^{\text{resp}} / r_v$ (1), plus the venous absorption terms (which depend on SvO₂ and [Hb]).

Note: no arterial term, no $I_0$. This is the cleanest equation for SvO₂.

---

## 4. Site 1: Forearm (Isolated Vein)

### 4.1 Anatomy

A single superficial vein (cephalic vein) runs beneath the sensor. No major artery is directly under the sensor array. The surrounding tissue is primarily muscle with a thin layer of skin and subcutaneous fat.

### 4.2 Compartment Model

Two compartments: muscle (background) and vein. No arterial compartment.

### 4.3 Available Signal Components

**DC (36 measurements: 12 positions × 3 wavelengths):**

$$\ln I_{DC}(\lambda_k, x_i) = \ln I_0^{\text{eff}}(\lambda_k) - \mu_a^{(m)}(\lambda_k) \cdot \langle L_m \rangle(\lambda_k, x_i) - \mu_a^{\text{ven}}(\lambda_k) \cdot \langle L_v \rangle(\lambda_k, x_i)$$

The spatial profile of $I_{DC}$ across the 12 PDs encodes the vein's lateral position: PDs directly above the vein see more attenuation (because $\langle L_v \rangle$ is largest there), while PDs far from the vein see mostly muscle.

**Cardiac-AC (36 measurements):**

At the forearm, there is no major artery under the sensor. The cardiac-frequency modulation arises from capillary-bed pulsation in the muscle and small arterioles — this is a diffuse, spatially uniform signal that carries little vessel-specific information. The cardiac-AC signal is expected to be weak (modulation ratio < 0.1%) and dominated by capillary blood at a saturation intermediate between arterial and venous.

Practically, the forearm cardiac-AC is not useful for our purposes. We set $\Delta r_a / r_a = 0$ (no artery) and $\Delta r_v^{\text{card}} / r_v \approx 0$ (minimal venous cardiac pulsation on the forearm).

**Respiratory-AC (36 measurements):**

$$M_{\text{resp}}(\lambda_k, x_i) = \mu_a^{\text{ven}}(\lambda_k) \cdot \langle L_v \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_v^{\text{resp}}}{r_v}$$

The respiratory modulation is expected to be small but present — venous distension during the inspiratory phase of mechanical ventilation affects even peripheral veins. The modulation ratio might be 0.1–1% at 660 nm.

This is the primary venous channel. The ratio of $M_{\text{resp}}$ at different wavelengths directly encodes the spectral shape of venous blood absorption, which depends on SvO₂.

### 4.4 Unknowns at the Forearm

| Unknown | Count | Constrained by |
|---------|-------|---------------|
| $\mu_a^{(m)}(\lambda_k)$ | 3 | DC: spatial profile at PDs far from vessel (where $\langle L_v \rangle \approx 0$) |
| $\ln I_0^{\text{eff}}(\lambda_k)$ | 3 | DC: overall intensity level |
| SvO₂ | 1 | Respiratory-AC: spectral ratios |
| [Hb] | 1 | Respiratory-AC: overall amplitude; DC: vessel contrast magnitude |
| $\Delta r_v^{\text{resp}} / r_v$ | 1 | Respiratory-AC: overall amplitude |

**Total: 9 unknowns from 72 measurements** (36 DC + 36 respiratory-AC). Well over-determined in count.

### 4.5 How SvO₂ is Recovered at the Forearm

**Step 1: Constrain muscle properties from DC.** PD positions far from the vein (large $|x_i - x_v|$) have $\langle L_v \rangle \approx 0$, so their DC signal is almost purely muscle:

$$\ln I_{DC}(\lambda_k, x_{\text{far}}) \approx \ln I_0^{\text{eff}}(\lambda_k) - \mu_a^{(m)}(\lambda_k) \cdot \langle L_m \rangle(\lambda_k, x_{\text{far}})$$

These measurements, combined across wavelengths, constrain muscle $\mu_a$ and $I_0^{\text{eff}}$ (or, more precisely, the combination $\ln I_0^{\text{eff}} - \mu_a^{(m)} \langle L_m \rangle$).

**Step 2: Constrain SvO₂ from spectral ratios of respiratory-AC.** The modulation ratio at wavelength $\lambda_k$ is:

$$M_{\text{resp}}(\lambda_k, x_i) = \mu_a^{\text{ven}}(\lambda_k) \cdot \langle L_v \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_v^{\text{resp}}}{r_v}$$

Take the ratio of modulation ratios at two wavelengths (say, 660 nm and 940 nm) at the same PD position:

$$\frac{M_{\text{resp}}(\lambda_1, x_i)}{M_{\text{resp}}(\lambda_3, x_i)} = \frac{\mu_a^{\text{ven}}(\lambda_1)}{\mu_a^{\text{ven}}(\lambda_3)} \cdot \frac{\langle L_v \rangle(\lambda_1, x_i)}{\langle L_v \rangle(\lambda_3, x_i)}$$

The $\Delta r_v / r_v$ term cancels (it is the same physical pulsation at both wavelengths). The pathlength ratio $\langle L_v \rangle(\lambda_1) / \langle L_v \rangle(\lambda_3)$ is known from MC (it depends on geometry and scattering, not on SvO₂ to first order). So the left side is measured, the pathlength ratio is known, and the only unknown is the absorption ratio $\mu_a^{\text{ven}}(\lambda_1) / \mu_a^{\text{ven}}(\lambda_3)$, which is a function of SvO₂ alone (Equation 25 in the derivation document).

This is the fundamental oximetry measurement — the ratio of blood absorption at two wavelengths determines oxygen saturation. The respiratory-AC signal provides this ratio cleanly because it is purely venous.

**Step 3: Determine [Hb] from the absolute magnitude.** Once SvO₂ is known (fixing the spectral shape of $\mu_a^{\text{ven}}$), the absolute magnitude of $M_{\text{resp}}$ at any one wavelength constrains the product $[\text{Hb}] \cdot \Delta r_v^{\text{resp}} / r_v$. To separate [Hb] from $\Delta r_v^{\text{resp}} / r_v$, use the DC vessel contrast (which depends on [Hb] through $\mu_a^{\text{ven}}$ but not on $\Delta r_v / r_v$).

**Step 4: Joint fit.** In practice, all unknowns are fit simultaneously using nonlinear least squares on the combined DC + respiratory-AC data vector (72 measurements, 9 unknowns). The stepwise reasoning above describes the information flow, not the actual solver procedure.

### 4.6 Challenges at the Forearm

**Small respiratory-AC signal.** Peripheral veins have lower compliance than central veins (like the IJV), so the respiratory venous modulation may be very small. If the modulation ratio is below the noise floor (say, < 0.05%), the respiratory-AC channel provides no useful information and the problem degenerates to DC-only.

**DC-only fallback.** If the respiratory-AC is undetectable, SvO₂ must be recovered from DC alone. This is harder because DC mixes muscle and blood contributions, and the calibration constants $I_0$ are unknown. However, the spatial profile across 12 PDs still provides leverage: the vessel-induced contrast pattern encodes the vessel's lateral position and, through the spectral dependence, the blood absorption. The DC-only inverse problem is less well-conditioned but may still yield SvO₂ to ±10%.

**This site is the simplest case and the first validation target.** If the method works here (where there is only one vessel and no arterial contamination), it provides confidence for the more complex sites.

---

## 5. Site 2: Wrist (Radial Artery + Vein)

### 5.1 Anatomy

The radial artery and an accompanying vein (cephalic vein or radial venae comitantes) lie beneath the sensor. The two vessels are laterally separated by approximately 2–4 mm and may be at slightly different depths (the artery is typically deeper). The surrounding tissue is muscle, tendon, and a thin layer of skin and subcutaneous fat.

### 5.2 Compartment Model

Three compartments: muscle (background), vein, and artery. The vessels are at different lateral positions ($x_v \neq x_a$), so different PDs in the array have different sensitivities to each vessel.

### 5.3 Available Signal Components

**DC (36 measurements):**

$$\ln I_{DC}(\lambda_k, x_i) = \ln I_0^{\text{eff}}(\lambda_k) - \mu_a^{(m)} \langle L_m \rangle - \mu_a^{\text{ven}} \langle L_v \rangle - \mu_a^{\text{art}} \langle L_a \rangle$$

The spatial profile now has two dips — one centered on the vein and one on the artery. The depth and width of each dip encode the vessel's position, depth, and absorption. Because the two vessels are at different x-positions, PDs near the vein see more venous absorption, PDs near the artery see more arterial absorption, and PDs between or away from both see mostly muscle.

**Cardiac-AC (36 measurements):**

$$M_{\text{card}}(\lambda_k, x_i) = \mu_a^{\text{art}}(\lambda_k) \cdot \langle L_a \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_a}{r_a} + \mu_a^{\text{ven}}(\lambda_k) \cdot \langle L_v \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_v^{\text{card}}}{r_v}$$

At the wrist, the radial artery is large and superficial, producing a strong cardiac-AC signal. The venous cardiac pulsation is typically small (forearm veins are not strongly pulsatile at the cardiac frequency). So the cardiac-AC is dominated by the arterial term:

$$M_{\text{card}}(\lambda_k, x_i) \approx \mu_a^{\text{art}}(\lambda_k) \cdot \langle L_a \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_a}{r_a}$$

Since SaO₂ is known (measured clinically), the arterial absorption $\mu_a^{\text{art}}(\lambda_k)$ is known up to the scalar [Hb]. The cardiac-AC therefore constrains [Hb] and validates the forward model.

**Respiratory-AC (36 measurements):**

$$M_{\text{resp}}(\lambda_k, x_i) \approx \mu_a^{\text{ven}}(\lambda_k) \cdot \langle L_v \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_v^{\text{resp}}}{r_v}$$

As at the forearm, the respiratory-AC is predominantly venous. It provides the cleanest SvO₂ information.

### 5.4 Unknowns at the Wrist

| Unknown | Count | Constrained by |
|---------|-------|---------------|
| $\mu_a^{(m)}(\lambda_k)$ | 3 | DC: PDs far from both vessels |
| $\ln I_0^{\text{eff}}(\lambda_k)$ | 3 | DC: overall intensity level |
| SvO₂ | 1 | Respiratory-AC: spectral ratios |
| [Hb] | 1 | Cardiac-AC: absolute amplitude (with known SaO₂) |
| $\Delta r_a / r_a$ | 1 | Cardiac-AC: absolute amplitude |
| $\Delta r_v^{\text{card}} / r_v$ | 1 | Cardiac-AC: residual after arterial model subtraction (may be negligible) |
| $\Delta r_v^{\text{resp}} / r_v$ | 1 | Respiratory-AC: absolute amplitude |

**Total: 11 unknowns from 108 measurements** (36 DC + 36 cardiac-AC + 36 respiratory-AC).

### 5.5 How SvO₂ is Recovered at the Wrist

**Step 1: Use cardiac-AC to determine [Hb].** Since SaO₂ is known, the arterial absorption spectrum is known up to the scalar [Hb]:

$$\mu_a^{\text{art}}(\lambda_k) = \frac{\ln(10) \cdot [\text{Hb}]}{MW_{\text{Hb}}} \left[\text{SaO}_2 \cdot \varepsilon_{\text{HbO}_2}(\lambda_k) + (1-\text{SaO}_2) \cdot \varepsilon_{\text{HbR}}(\lambda_k)\right] + \mu_{a,w}(\lambda_k) \cdot 0.95$$

The cardiac-AC modulation ratio at each wavelength is:

$$M_{\text{card}}(\lambda_k, x_i) \approx \mu_a^{\text{art}}(\lambda_k; [\text{Hb}]) \cdot \langle L_a \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_a}{r_a}$$

Everything on the right is known except [Hb] and $\Delta r_a / r_a$. Taking the ratio at two wavelengths at the same PD position cancels $\Delta r_a / r_a$:

$$\frac{M_{\text{card}}(\lambda_1, x_i)}{M_{\text{card}}(\lambda_2, x_i)} = \frac{\mu_a^{\text{art}}(\lambda_1) \cdot \langle L_a \rangle(\lambda_1, x_i)}{\mu_a^{\text{art}}(\lambda_2) \cdot \langle L_a \rangle(\lambda_2, x_i)}$$

The pathlength ratio is known from MC. The absorption ratio depends on SaO₂ (known) but not on [Hb] (it cancels in the ratio). So this ratio is a **consistency check** — it should match the predicted ratio based on the known SaO₂ and MC pathlengths. If it doesn't, something is wrong with the geometry model or the MC simulation.

Once validated, the absolute magnitude of $M_{\text{card}}$ at any one wavelength and PD position gives the product $[\text{Hb}] \cdot \Delta r_a / r_a$. With one more constraint (e.g., from DC vessel contrast, or from a second PD position that gives a different $\langle L_a \rangle$), [Hb] and $\Delta r_a / r_a$ can be separated.

**Step 2: Use respiratory-AC to determine SvO₂.** Exactly as at the forearm (Section 4.5, Step 2): the spectral ratio of respiratory-AC modulation ratios encodes the venous absorption ratio, which is a function of SvO₂ alone.

The key advantage over the forearm: [Hb] is now constrained by the cardiac-AC (Step 1), so it does not need to be determined from the respiratory-AC or DC. This reduces the number of unknowns the respiratory-AC channel must resolve, improving the conditioning of the SvO₂ estimate.

**Step 3: Use DC to constrain muscle properties and validate.** Same as the forearm: PDs far from both vessels constrain muscle $\mu_a$ and $I_0^{\text{eff}}$.

**Step 4: Joint fit.** All 108 measurements are fit simultaneously, but the information flow is:

- Cardiac-AC → SaO₂ validation + [Hb] + $\Delta r_a / r_a$
- Respiratory-AC → SvO₂ + $\Delta r_v^{\text{resp}} / r_v$
- DC → muscle properties + $I_0^{\text{eff}}$ + geometric consistency

### 5.6 The Spatial Leverage at the Wrist

Because the artery and vein are at different lateral positions, the 12 PDs provide **spatial separation** of the two vessels. PDs nearer the artery's lateral position ($x \approx x_a$) have large $\langle L_a \rangle$ and small $\langle L_v \rangle$; PDs nearer the vein ($x \approx x_v$) have the opposite. This means:

- The cardiac-AC signal varies across PDs in a pattern that reflects the artery's position (peak near $x_a$, falling off laterally).
- The respiratory-AC signal varies in a pattern that reflects the vein's position (peak near $x_v$).
- These patterns are different, which helps deconvolve the two contributions.

If we define the **sensitivity ratio** at each PD as:

$$R_i = \frac{\langle L_a \rangle(\lambda, x_i)}{\langle L_v \rangle(\lambda, x_i)}$$

then PDs with high $R_i$ are artery-dominated and PDs with low $R_i$ are vein-dominated. The cardiac-AC should correlate with the artery sensitivity pattern, and the respiratory-AC with the vein sensitivity pattern. This is verifiable from the MC simulation even before collecting real data.

### 5.7 Why the Wrist is Better-Conditioned than the Forearm

Three factors make the wrist a stronger site for SvO₂ recovery:

1. **The artery provides [Hb] independently.** At the forearm, [Hb] must be estimated from the venous signal or assumed from literature. At the wrist, the arterial cardiac-AC gives [Hb] directly (because SaO₂ is known), removing a confound.

2. **Spatial separation breaks degeneracies.** The two-vessel geometry means the 12 PDs see genuinely different mixtures of arterial and venous signal, providing geometric leverage that the single-vessel forearm lacks.

3. **Larger cardiac-AC signal.** The radial artery produces a strong pulsatile signal (modulation ratio ~0.5–2% at 660 nm), providing high-SNR data for the cardiac-AC channel.

---

## 6. Site 3: Neck (IJV + Carotid)

### 6.1 Anatomy

The internal jugular vein (IJV) and common carotid artery are large vessels at moderate depth beneath the sensor. The IJV is lateral to the carotid, with a center-to-center separation of approximately 5–10 mm. Overlying tissue includes skin, subcutaneous fat, the sternocleidomastoid muscle, and possibly the platysma muscle.

The IJV is the largest peripheral vein typically accessible for non-invasive optical measurement and has the strongest venous pulsatility of the three sites.

### 6.2 Compartment Model

Five compartments are recommended for anatomical accuracy: skin, subcutaneous fat, muscle (including SCM), IJV, and carotid artery. For the simplified model, three compartments (muscle + IJV + carotid) may suffice if skin and fat are thin and can be absorbed into the muscle background.

### 6.3 The Venous Pulsatility Problem

The IJV is strongly pulsatile at the cardiac frequency. Unlike the forearm vein or the wrist cephalic vein, the IJV shows prominent a, c, and v waves that are synchronized with the cardiac cycle. This means the assumption "cardiac-AC = arterial" is invalid at the neck.

Specifically, the cardiac-AC modulation ratio at the neck is:

$$M_{\text{card}}(\lambda_k, x_i) = \mu_a^{\text{art}}(\lambda_k) \cdot \langle L_a \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_a}{r_a} + \mu_a^{\text{ven}}(\lambda_k) \cdot \langle L_v \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_v^{\text{card}}}{r_v}$$

Both terms are significant. The arterial and venous contributions cannot be separated by frequency alone (both are at $f_c$). If we naively attributed the entire cardiac-AC to the artery (as pulse oximetry does), the "recovered SaO₂" would be a contaminated mixture of SaO₂ and SvO₂, biased toward SvO₂ when venous pulsation is large.

### 6.4 Available Signal Components

**DC (36 measurements):**

$$\ln I_{DC}(\lambda_k, x_i) = \ln I_0^{\text{eff}}(\lambda_k) - \mu_a^{(m)} \langle L_m \rangle - \mu_a^{\text{ven}} \langle L_v \rangle - \mu_a^{\text{art}} \langle L_a \rangle$$

Similar to the wrist, but vessels are deeper (less contrast) and the overlying tissue layers (skin, fat) may attenuate the signal more.

**Cardiac-AC (36 measurements):**

$$M_{\text{card}}(\lambda_k, x_i) = \mu_a^{\text{art}}(\lambda_k) \cdot \langle L_a \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_a}{r_a} + \mu_a^{\text{ven}}(\lambda_k) \cdot \langle L_v \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_v^{\text{card}}}{r_v}$$

Both terms are significant. However, the spatial profiles $\langle L_a \rangle(x_i)$ and $\langle L_v \rangle(x_i)$ are different (the artery and IJV are at different lateral positions), so the 12 PDs sample different mixtures. This provides the leverage to separate the two.

**Respiratory-AC (36 measurements):**

$$M_{\text{resp}}(\lambda_k, x_i) \approx \mu_a^{\text{ven}}(\lambda_k) \cdot \langle L_v \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_v^{\text{resp}}}{r_v}$$

The IJV shows large respiratory modulation — this is the strongest respiratory-AC signal of all three sites. Under positive-pressure mechanical ventilation, the IJV can distend by 15–30% in radius during inspiration. This produces modulation ratios of 1–5% at 660 nm, well above the noise floor.

This is the primary channel for SvO₂ at the neck, and it is the cleanest of the three sites for respiratory-AC signal quality.

### 6.5 Unknowns at the Neck

| Unknown | Count | Constrained by |
|---------|-------|---------------|
| $\mu_a^{(m)}(\lambda_k)$ | 3 | DC: PDs far from both vessels |
| $\ln I_0^{\text{eff}}(\lambda_k)$ | 3 | DC: overall intensity level |
| SvO₂ | 1 | Respiratory-AC (primary); cardiac-AC (secondary, after decomposition) |
| [Hb] | 1 | Cardiac-AC + known SaO₂ |
| $\Delta r_a / r_a$ | 1 | Cardiac-AC (arterial component) |
| $\Delta r_v^{\text{card}} / r_v$ | 1 | Cardiac-AC (venous component) |
| $\Delta r_v^{\text{resp}} / r_v$ | 1 | Respiratory-AC: amplitude |

**Total: 11 unknowns from 108 measurements.**

### 6.6 How SvO₂ is Recovered at the Neck

**Step 1: Use respiratory-AC as the primary SvO₂ channel.** Because the respiratory-AC is predominantly venous (arterial respiratory pulsation is negligible), and because the IJV's respiratory modulation is large and clean under mechanical ventilation, the spectral ratio of respiratory-AC modulation ratios directly encodes SvO₂:

$$\frac{M_{\text{resp}}(\lambda_1, x_i)}{M_{\text{resp}}(\lambda_3, x_i)} = \frac{\mu_a^{\text{ven}}(\lambda_1) \cdot \langle L_v \rangle(\lambda_1, x_i)}{\mu_a^{\text{ven}}(\lambda_3) \cdot \langle L_v \rangle(\lambda_3, x_i)}$$

With pathlength ratios from MC and extinction coefficients from the literature, SvO₂ is determined from this ratio alone. This is the most robust SvO₂ estimate at any of the three sites, because:

- The IJV respiratory-AC signal is the largest of the three sites.
- Mechanical ventilation produces clean, regular respiratory modulation.
- The signal is purely venous (no arterial contamination at the respiratory frequency).

**Step 2: Use cardiac-AC for [Hb] and model validation.** The cardiac-AC is a mixture of arterial and venous contributions, but we can use the forward model to decompose them:

For each PD position $x_i$ and wavelength $\lambda_k$, the cardiac-AC modulation ratio depends on four blood-related unknowns: SvO₂ (now known from Step 1), SaO₂ (known clinically), [Hb], and the two pulsation amplitudes ($\Delta r_a / r_a$ and $\Delta r_v^{\text{card}} / r_v$).

With SvO₂ and SaO₂ both known, the spectral shapes of arterial and venous absorption are fully determined. The cardiac-AC equation at each PD becomes:

$$M_{\text{card}}(\lambda_k, x_i) = [\text{Hb}] \cdot \left[\alpha_a(\lambda_k) \cdot \langle L_a \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_a}{r_a} + \alpha_v(\lambda_k) \cdot \langle L_v \rangle(\lambda_k, x_i) \cdot \frac{\Delta r_v^{\text{card}}}{r_v}\right]$$

where $\alpha_a(\lambda_k) = (\ln 10 / MW_{\text{Hb}}) [SaO_2 \cdot \varepsilon_{\text{HbO}_2}(\lambda_k) + (1-SaO_2) \cdot \varepsilon_{\text{HbR}}(\lambda_k)] + \mu_{a,w}(\lambda_k) \cdot 0.95 / [\text{Hb}]$ and similarly for $\alpha_v$.

This is linear in [Hb], $\Delta r_a / r_a$, and $\Delta r_v^{\text{card}} / r_v$ (since $\alpha_a$ and $\alpha_v$ are known after Step 1). With 36 cardiac-AC measurements and only 3 unknowns, this is heavily over-determined and can be solved robustly via linear least squares.

**Spatial decomposition:** The key to separating arterial and venous cardiac-AC contributions is that $\langle L_a \rangle(x_i)$ and $\langle L_v \rangle(x_i)$ have different spatial profiles (the IJV and carotid are at different lateral positions). PDs closer to the carotid see more arterial cardiac-AC; PDs closer to the IJV see more venous cardiac-AC. With 12 PD positions, you are effectively fitting two spatial templates (arterial and venous sensitivity profiles) to the cardiac-AC spatial pattern — this is equivalent to a two-component spatial regression.

**Step 3: Use DC for muscle properties and calibration.** Same as the other sites.

**Step 4: Cross-validate.** The cardiac-AC, after decomposition, provides a secondary estimate of SvO₂ (from the venous cardiac component's spectral ratios). This should agree with the respiratory-AC estimate from Step 1. Disagreement indicates model error — possible causes include incorrect vessel positions, incorrect tissue layering, or the presence of additional anatomical structures (e.g., the external jugular vein, thyroid, lymph nodes) not included in the model.

### 6.7 Why the Neck is the Most Clinically Important but Hardest Site

**Clinical importance:** The IJV carries blood draining from the brain. Jugular venous SvO₂ (SjvO₂) is a surrogate for cerebral oxygen extraction and is monitored in neurocritical care. Non-invasive measurement of SjvO₂ would replace invasive jugular bulb catheters.

**Difficulty:**

1. **Venous pulsatility at the cardiac frequency.** The IJV pulsates strongly, making the cardiac-AC a mixture of arterial and venous contributions. This is the central technical challenge.

2. **Deeper vessels.** The IJV and carotid are deeper than the forearm or wrist vessels (typically 5–15 mm to the center), which reduces the signal contrast and the sensitivity of the measurement.

3. **Complex overlying anatomy.** The sternocleidomastoid muscle, subcutaneous fat, and skin create a multilayered tissue model that is more complex than the forearm or wrist.

4. **Posture sensitivity.** The IJV diameter varies dramatically with head-of-bed elevation. Supine position gives the largest IJV and strongest signal; head-up position can nearly collapse the IJV.

5. **Carotid dominance.** The carotid artery is very large and produces strong cardiac-AC. The arterial signal can be 5–10× larger than the venous cardiac component, making it difficult to extract the smaller venous contribution from the cardiac-AC.

Despite these challenges, the respiratory-AC channel largely bypasses the venous pulsatility problem and provides a clean path to SvO₂. The large IJV respiratory modulation under mechanical ventilation is actually an advantage, making this the site with the strongest respiratory-AC signal.

---

## 7. Summary: Information Flow at Each Site

### 7.1 Forearm (Simplest)

```
DC ──────────────────────> muscle μₐ, I₀
                           (PDs far from vessel constrain muscle)

Respiratory-AC ──────────> SvO₂ (from spectral ratios)
                           [Hb] × Δrᵥ/rᵥ (from absolute amplitude)
                           
DC (vessel contrast) ────> [Hb] (combined with resp-AC amplitude)

Cardiac-AC ──────────────> [not used — too weak without a major artery]
```

### 7.2 Wrist (Intermediate)

```
Cardiac-AC ──────────────> [Hb] (using known SaO₂)
                           Δrₐ/rₐ (arterial pulsation amplitude)
                           SaO₂ validation (spectral ratio check)

Respiratory-AC ──────────> SvO₂ (from spectral ratios, same as forearm)
                           Δrᵥʳᵉˢᵖ/rᵥ

DC ──────────────────────> muscle μₐ, I₀
                           Geometric consistency of both vessels

Spatial separation ──────> 12 PDs see different artery/vein mixtures
                           → separates arterial and venous contributions
```

### 7.3 Neck (Most Complex)

```
Respiratory-AC ──────────> SvO₂ (PRIMARY — cleanest venous channel)
                           Δrᵥʳᵉˢᵖ/rᵥ
                           (No arterial contamination at respiratory freq)

Cardiac-AC ──────────────> [Hb] (after spatial decomposition using known
                           SvO₂ from resp-AC and known SaO₂)
                           Δrₐ/rₐ and Δrᵥᶜᵃʳᵈ/rᵥ
                           (Separated by spatial profiles from MC)
                           Secondary SvO₂ estimate (cross-validation)

DC ──────────────────────> muscle μₐ, I₀
                           Tissue layer properties

Spatial separation ──────> Critical for cardiac-AC decomposition
                           IJV and carotid at different x-positions
                           → 12 PDs provide spatial regression
```

---

## 8. Comparison Across Sites

| Feature | Forearm | Wrist | Neck |
|---------|---------|-------|------|
| **Compartments** | 2 (muscle + vein) | 3 (muscle + vein + artery) | 3–5 (skin, fat, muscle, IJV, carotid) |
| **Primary SvO₂ channel** | Respiratory-AC | Respiratory-AC | Respiratory-AC |
| **Cardiac-AC role** | Not useful | [Hb] determination, SaO₂ validation | [Hb] + spatial decomposition, secondary SvO₂ |
| **[Hb] from** | DC contrast + resp-AC amplitude | Cardiac-AC (with known SaO₂) | Cardiac-AC (with known SaO₂ + SvO₂) |
| **Venous cardiac pulsation** | Negligible | Small | Large (the central problem) |
| **Respiratory-AC amplitude** | Small | Moderate | Large |
| **Spatial leverage** | Moderate (1 vessel) | Strong (2 vessels at different x) | Strong (2 vessels at different x) |
| **Expected SvO₂ accuracy** | ±8–12% | ±5–8% | ±3–6% |
| **Clinical relevance** | Low (peripheral vein) | Moderate (local tissue oxygenation) | High (cerebral venous drainage) |

---

## 9. The Role of Known SaO₂ Across All Sites

SaO₂ is measured clinically (pulse oximeter or arterial line) and provided as a fixed input. Its impact varies by site:

**Forearm:** SaO₂ is not directly used (no artery in the model). It could be used to set the capillary bed saturation in the muscle model ($S_{\text{musc}} \approx (SaO_2 + SvO_2)/2$), but this is a minor refinement.

**Wrist:** SaO₂ is essential. It determines the spectral shape of arterial absorption, which is needed to interpret the cardiac-AC signal. Without known SaO₂, the cardiac-AC would have two more unknowns (SaO₂ and [Hb] combined), making [Hb] recovery impossible.

**Neck:** SaO₂ is essential for the same reason as the wrist, and additionally important because the cardiac-AC must be decomposed into arterial and venous components. Without known SaO₂, the spectral decomposition of cardiac-AC into "arterial-shaped" and "venous-shaped" contributions would be ambiguous.

---

## 10. The Role of Ultrasound Geometry Across All Sites

Ultrasound provides the vessel parameters ($r_v$, $d_v$, $x_v$, $r_a$, $d_a$, $x_a$) that define the geometry for the Monte Carlo simulation. The MC simulation then produces the compartment pathlengths $\langle L_j \rangle(\lambda, x_i)$ that appear in all three forward models (DC, cardiac-AC, respiratory-AC).

**Without ultrasound:** The vessel geometry becomes unknown, and the pathlengths become free parameters. With only 1 source–detector separation per wavelength, the problem becomes severely under-determined. The spatial pattern across PDs provides some geometric information, but not enough to reconstruct the full 3D vessel geometry.

**With ultrasound:** The geometry is fixed, the MC pathlengths are known constants, and the only free parameters are the tissue optical properties and pulsation amplitudes. The problem is well-determined.

Ultrasound accuracy requirements: vessel depth to ±1 mm, vessel diameter to ±0.5 mm, lateral position to ±1 mm. These are achievable with standard clinical B-mode ultrasound at frequencies of 7–15 MHz.

The robustness of the SvO₂ estimate to ultrasound measurement errors should be characterized through the sensitivity analysis described in the simulation validation plan (Section 9.5 of the companion document).
