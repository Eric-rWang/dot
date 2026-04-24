# Derivation of the Two-Compartment Modified Beer–Lambert Law for CW-DOSI Venous Oximetry

---

## 1. Background and Motivation

### 1.1 The Problem

We seek to measure venous oxygen saturation (SvO₂) non-invasively using a continuous-wave (CW) near-infrared sensor placed on the skin surface above a superficial blood vessel. The sensor consists of LED sources at three wavelengths (660, 850, 940 nm) and photodetector (PD) arrays. Each LED–PD pair measures the intensity of light that has traveled through the tissue between source and detector.

The tissue between source and detector contains two optically distinct compartments: a bulk background (muscle) and an embedded blood vessel (vein). Each compartment absorbs and scatters light differently, and the blood vessel's absorption depends on the oxygen saturation of the hemoglobin it contains. Our goal is to extract SvO₂ from the measured intensities.

The central theoretical challenge is: how does the detected intensity depend on the absorption properties of each compartment, given that photons take complex, random paths through the scattering tissue?

### 1.2 Why Not Classical Beer–Lambert?

In a transparent (non-scattering) medium — such as a cuvette of dye in a spectrophotometer — light travels in a straight line from source to detector. The classical Beer–Lambert law applies:

$$I = I_0 \cdot e^{-\mu_a \cdot d}$$

where $I_0$ is the incident intensity, $\mu_a$ is the absorption coefficient of the medium (cm⁻¹), and $d$ is the geometric path length through the medium (cm).

This law assumes:

- Light travels a single, known path of length $d$.
- No scattering occurs (or scattering is negligible).
- The medium is homogeneous along the path.

None of these assumptions hold in biological tissue. Tissue is a highly scattering medium (reduced scattering coefficient $\mu_s' \approx 5$–$15$ cm⁻¹ in the NIR), which means photons undergo many scattering events between source and detector. Each detected photon has traveled a different, random path through the tissue, and the actual distance traveled is much longer than the geometric source–detector separation $\rho$.

We need a formulation that accounts for the distribution of photon paths through a scattering, heterogeneous medium.

---

## 2. Photon Transport in a Scattering Medium

### 2.1 The Microscopic Beer–Lambert Law

While the macroscopic Beer–Lambert law fails for bulk tissue, a microscopic version still holds: along any individual photon path, absorption acts exponentially. If a photon travels a total path length $\ell$ through a medium with absorption coefficient $\mu_a$, the probability that it survives (is not absorbed) is:

$$P_{\text{survive}}(\ell) = e^{-\mu_a \cdot \ell}$$

This is exact. It follows from the fundamental definition of $\mu_a$: the probability of absorption per unit path length is $\mu_a$, so over a path of length $\ell$, the survival probability is the product of independent survival probabilities at each infinitesimal step, giving $e^{-\mu_a \ell}$.

Scattering does not violate this law — it changes the *direction* of the photon (and therefore the *path* it takes), but along whatever path the photon follows, absorption still acts as $e^{-\mu_a \ell}$.

**Key reference:** Tsuchiya (2001), "Photon path distribution and optical responses of turbid media: theoretical analysis based on the microscopic Beer–Lambert law," *Phys. Med. Biol.* 46, 2067–2084.

### 2.2 The Photon Path Length Distribution

In a scattering medium, photons launched from the source take random-walk trajectories determined by the scattering coefficient $\mu_s$, the scattering anisotropy $g$, and the geometry of the medium. Some photons are scattered back toward the surface and reach the detector; most are scattered away or absorbed.

Consider the ensemble of all possible photon paths from the source to the detector. Each path has a total length $\ell$ (the sum of all step lengths along the zigzag trajectory). The **photon path length distribution** $p(\ell)$ describes the probability density that a photon arriving at the detector has traveled a total path length $\ell$.

Critically, $p(\ell)$ is determined by the **scattering geometry alone**: it depends on $\mu_s$, $g$, the medium geometry, and the source–detector configuration, but **not** on the absorption coefficient $\mu_a$. This is because scattering determines which paths exist and how likely they are; absorption then acts as a filter that preferentially removes photons on longer paths.

This separation — geometry determines the paths, absorption filters them — is the foundation of the entire framework.

**Key reference:** Tsuchiya (2001), as above. The path length distribution is what Tsuchiya calls the "photon path distribution" (PPD).

---

## 3. Detected Intensity as a Laplace Transform

### 3.1 The Exact Expression for Detected Intensity

Consider a homogeneous scattering medium with absorption coefficient $\mu_a$. The detected intensity is the sum over all photon paths of the source intensity times the survival probability along each path:

$$I(\mu_a) = I_0 \int_0^\infty p(\ell) \cdot e^{-\mu_a \cdot \ell} \; d\ell \tag{1}$$

where:

- $I_0$ is the source intensity times all coupling factors (LED power, detector responsivity, geometric collection efficiency).
- $p(\ell)$ is the path length distribution of photons reaching the detector in the absence of absorption (determined by scattering and geometry alone).
- $e^{-\mu_a \ell}$ is the survival probability (microscopic Beer–Lambert law) for a photon traveling path length $\ell$.

The integral sums over all possible path lengths, weighting each by the probability that a photon takes that path (given by $p(\ell)$) and survives absorption (given by $e^{-\mu_a \ell}$).

**Equation (1) is exact.** No approximations have been made beyond the assumption that absorption and scattering are separable (i.e., the path length distribution is independent of absorption). This separation is exact for media where scattering is elastic and isotropic in its effect on the path geometry, which is an excellent approximation for biological tissue in the NIR.

**Mathematical note:** The integral $\int_0^\infty p(\ell) \cdot e^{-\mu_a \ell} \, d\ell$ is the **Laplace transform** of the path length distribution $p(\ell)$, evaluated at $s = \mu_a$. Denoting this as:

$$\mathcal{L}\{p\}(\mu_a) = \int_0^\infty p(\ell) \cdot e^{-\mu_a \cdot \ell} \; d\ell$$

we can write:

$$I(\mu_a) = I_0 \cdot \mathcal{L}\{p\}(\mu_a) \tag{2}$$

### 3.2 Log-Intensity

Taking the natural logarithm of Equation (2):

$$\ln I(\mu_a) = \ln I_0 + \ln \mathcal{L}\{p\}(\mu_a) \tag{3}$$

This is still exact. The log-intensity is the log of a Laplace transform of the path length distribution.

---

## 4. The Mean Pathlength of Detected Photons

### 4.1 Definition

The detected photons are a subset of all photons that reached the detector position — specifically, those that survived absorption along their path. The surviving population is biased: photons with shorter paths are more likely to survive than those with longer paths.

The **mean path length of detected photons** is the average path length over this biased population:

$$\langle L \rangle(\mu_a) = \frac{\int_0^\infty \ell \cdot p(\ell) \cdot e^{-\mu_a \cdot \ell} \; d\ell}{\int_0^\infty p(\ell) \cdot e^{-\mu_a \cdot \ell} \; d\ell} \tag{4}$$

The numerator weights each path length $\ell$ by both the path probability $p(\ell)$ and the survival probability $e^{-\mu_a \ell}$, then integrates. The denominator is the total detected intensity (the normalization). This is the standard definition of a weighted mean.

**Important:** $\langle L \rangle$ depends on $\mu_a$. As $\mu_a$ increases, long-path photons are preferentially absorbed, so the surviving population is enriched in short-path photons, and $\langle L \rangle$ decreases. This is a real physical effect, not a mathematical artifact.

**In the limit $\mu_a \to 0$:** No photons are absorbed, and the mean path length reduces to:

$$\langle L \rangle_0 = \int_0^\infty \ell \cdot p(\ell) \; d\ell$$

which is simply the first moment of the unweighted path length distribution. This is the longest possible mean path length for the given geometry.

### 4.2 Relationship to the Derivative of Log-Intensity

We now derive the central result: the mean path length of detected photons equals the negative derivative of log-intensity with respect to absorption coefficient.

Starting from Equation (3):

$$\ln I(\mu_a) = \ln I_0 + \ln \mathcal{L}\{p\}(\mu_a)$$

Differentiate with respect to $\mu_a$:

$$\frac{d \ln I}{d \mu_a} = \frac{d}{d\mu_a} \ln \mathcal{L}\{p\}(\mu_a) = \frac{\mathcal{L}\{p\}'(\mu_a)}{\mathcal{L}\{p\}(\mu_a)} \tag{5}$$

Now compute $\mathcal{L}\{p\}'(\mu_a)$, the derivative of the Laplace transform with respect to its argument:

$$\mathcal{L}\{p\}'(\mu_a) = \frac{d}{d\mu_a} \int_0^\infty p(\ell) \cdot e^{-\mu_a \ell} \; d\ell = \int_0^\infty p(\ell) \cdot (-\ell) \cdot e^{-\mu_a \ell} \; d\ell = -\int_0^\infty \ell \cdot p(\ell) \cdot e^{-\mu_a \ell} \; d\ell \tag{6}$$

Substituting Equation (6) into Equation (5):

$$\frac{d \ln I}{d \mu_a} = \frac{-\int_0^\infty \ell \cdot p(\ell) \cdot e^{-\mu_a \ell} \; d\ell}{\int_0^\infty p(\ell) \cdot e^{-\mu_a \ell} \; d\ell} = -\langle L \rangle(\mu_a) \tag{7}$$

where the last equality follows from the definition of $\langle L \rangle$ in Equation (4).

**Result:** The derivative of log-intensity with respect to absorption coefficient is the negative mean path length of detected photons:

$$\boxed{\frac{d \ln I}{d \mu_a} = -\langle L \rangle(\mu_a)} \tag{7}$$

**This result is exact.** It requires no approximations beyond the separability of absorption and scattering (which is exact for elastic scattering).

**Key references:**

- Arridge, Cope, & Delpy (1992), "The theoretical basis for the determination of optical pathlengths in tissue: temporal and frequency analysis," *Phys. Med. Biol.* 37, 1531–1560. (Derives the mean pathlength from the temporal moments of the Green's function.)
- Sassaroli & Fantini (2004), "Comment on the modified Beer–Lambert law for scattering media," *Phys. Med. Biol.* 49, N255–N257. (States the relationship $\langle L \rangle = -\partial \ln I / \partial \mu_a$ explicitly and clarifies notation.)
- Tsuchiya (2001). (Derives the framework using the photon path distribution.)

### 4.3 Connection to Time-of-Flight Measurements

The mean path length $\langle L \rangle$ is experimentally accessible through time-resolved measurements. If a short pulse of light is injected into the tissue, the detected photons arrive spread out in time. The arrival time $t$ of a photon is related to its path length by:

$$\ell = \frac{c}{n} \cdot t$$

where $c$ is the speed of light in vacuum and $n$ is the refractive index of the medium. The mean path length is therefore:

$$\langle L \rangle = \frac{c}{n} \cdot \langle t \rangle$$

where $\langle t \rangle$ is the mean arrival time (first moment of the temporal point spread function). This relationship was established experimentally and validated by Monte Carlo simulation in:

- Delpy et al. (1988), "Estimation of optical pathlength through tissue from direct time of flight measurement," *Phys. Med. Biol.* 33, 1433–1442.

### 4.4 The Differential Pathlength Factor (DPF)

For a homogeneous semi-infinite medium, the total mean path length $\langle L \rangle$ is much longer than the geometric source–detector separation $\rho$. The ratio is defined as the **differential pathlength factor**:

$$\text{DPF} = \frac{\langle L \rangle}{\rho}$$

Typical values are DPF $\approx$ 4–6 for tissue in the NIR, meaning detected photons travel 4–6 times the direct source–detector distance due to the random-walk nature of their trajectories.

The DPF can be estimated from diffusion theory (Arridge et al. 1992) or measured experimentally via time-of-flight (Delpy et al. 1988).

**Important distinction:** DPF $\times$ $\rho$ gives the *total* mean path length in a *homogeneous* medium. In a heterogeneous medium (such as tissue with an embedded vessel), there is no single DPF. The mean path length must be decomposed into contributions from each compartment, which cannot be expressed as simple multiples of $\rho$. This is addressed in Section 6.

---

## 5. The Modified Beer–Lambert Law (One Compartment)

### 5.1 Taylor Expansion of Log-Intensity

From Equation (3), $\ln I(\mu_a)$ is a smooth function of $\mu_a$. We can expand it as a Taylor series around any reference absorption value $\mu_a^{\text{ref}}$:

$$\ln I(\mu_a) = \ln I(\mu_a^{\text{ref}}) + \frac{d \ln I}{d\mu_a}\bigg|_{\mu_a^{\text{ref}}} \cdot (\mu_a - \mu_a^{\text{ref}}) + \frac{1}{2}\frac{d^2 \ln I}{d\mu_a^2}\bigg|_{\mu_a^{\text{ref}}} \cdot (\mu_a - \mu_a^{\text{ref}})^2 + \cdots \tag{8}$$

### 5.2 First-Order Approximation

Truncating the Taylor series at first order and using Equation (7):

$$\ln I(\mu_a) \approx \ln I(\mu_a^{\text{ref}}) - \langle L \rangle(\mu_a^{\text{ref}}) \cdot (\mu_a - \mu_a^{\text{ref}}) \tag{9}$$

This is the **modified Beer–Lambert law (MBLL)** for a homogeneous medium, linearized around the reference state $\mu_a^{\text{ref}}$.

**Special case: reference state is zero absorption** ($\mu_a^{\text{ref}} = 0$). Then $I(\mu_a^{\text{ref}}) = I_0 \cdot \int p(\ell) \, d\ell = I_0$ (assuming $p(\ell)$ is normalized), and the mean path length is $\langle L \rangle_0$, the unweighted mean. Equation (9) becomes:

$$\ln I \approx \ln I_0 - \langle L \rangle_0 \cdot \mu_a \tag{10}$$

or equivalently:

$$I \approx I_0 \cdot e^{-\mu_a \cdot \langle L \rangle_0} \tag{11}$$

This looks like the classical Beer–Lambert law with the geometric distance $d$ replaced by the mean optical path length $\langle L \rangle_0$. This is the form most commonly encountered in the NIRS literature.

**Key references:**

- Delpy et al. (1988) — introduced the concept of replacing the geometric path with the mean optical pathlength.
- Cope (1991), PhD thesis, University College London — worked out the MBLL formalism in detail.
- Sassaroli & Fantini (2004) — clarified that the MBLL is strictly valid only for differential changes, and that the "absolute" form (Equation 11) is an approximation.

### 5.3 Validity and Error of the Linearization

The error of the first-order approximation (Equation 9) is given by the second-order term in the Taylor expansion (Equation 8). From Equation (7), the second derivative of $\ln I$ with respect to $\mu_a$ is:

$$\frac{d^2 \ln I}{d\mu_a^2} = -\frac{d\langle L \rangle}{d\mu_a} = \text{Var}(L) \tag{12}$$

where $\text{Var}(L) = \langle L^2 \rangle - \langle L \rangle^2$ is the **variance** of the path length distribution of detected photons.

Since variance is always non-negative, $d^2 \ln I / d\mu_a^2 \geq 0$, meaning $\ln I(\mu_a)$ is a **convex function** of $\mu_a$. The first-order approximation (the tangent line) always lies below the true curve, meaning the MBLL **overestimates** the attenuation (underestimates $I$) when $\mu_a$ deviates from the reference value.

The relative error of the linearization is approximately:

$$\text{error} \approx \frac{1}{2} \text{Var}(L) \cdot (\Delta\mu_a)^2$$

This is small when $\Delta\mu_a \cdot \sigma_L \ll 1$, where $\sigma_L = \sqrt{\text{Var}(L)}$ is the standard deviation of the path length distribution.

---

## 6. Extension to Two Compartments

### 6.1 Two-Compartment Geometry

We now consider the geometry relevant to our measurement: tissue containing two compartments, muscle (index $m$) and a blood vessel (index $v$), each with a distinct absorption coefficient.

Each photon path passes through both compartments. For a single photon path, define:

- $\ell_m$ = total distance traveled in the muscle compartment along this path.
- $\ell_v$ = total distance traveled in the vessel compartment along this path.
- $\ell = \ell_m + \ell_v$ = total path length.

The fractions $\ell_m$ and $\ell_v$ are determined by the scattering geometry and the spatial arrangement of the compartments. A photon that happens to pass through the vessel will have $\ell_v > 0$; a photon whose path misses the vessel entirely will have $\ell_v = 0$.

### 6.2 Joint Path Length Distribution

Define the **joint path length distribution** $p(\ell_m, \ell_v)$: the probability density that a photon reaching the detector (in the absence of absorption) has traveled a distance $\ell_m$ in muscle and $\ell_v$ in the vessel.

This joint distribution is determined by:

- The scattering properties of each compartment ($\mu_s$, $g$).
- The geometry (vessel position, size, depth, source/detector locations).
- The refractive indices of each compartment.

It does **not** depend on the absorption coefficients $\mu_a^{(m)}$ or $\mu_a^{(v)}$.

### 6.3 Detected Intensity: Exact Expression

Extending Equation (1) to two compartments, the survival probability along a path with $\ell_m$ in muscle and $\ell_v$ in vessel is:

$$P_{\text{survive}} = e^{-\mu_a^{(m)} \cdot \ell_m - \mu_a^{(v)} \cdot \ell_v}$$

This is the product of the survival probabilities in each compartment, which follows from the microscopic Beer–Lambert law applied separately in each region.

The detected intensity is the integral over all possible path partitions:

$$I(\mu_a^{(m)}, \mu_a^{(v)}) = I_0 \int_0^\infty \int_0^\infty p(\ell_m, \ell_v) \cdot e^{-\mu_a^{(m)} \ell_m - \mu_a^{(v)} \ell_v} \; d\ell_m \; d\ell_v \tag{13}$$

**This is exact.** It is the two-dimensional Laplace transform of the joint path length distribution:

$$I = I_0 \cdot \mathcal{L}\{p\}(\mu_a^{(m)}, \mu_a^{(v)}) \tag{14}$$

### 6.4 Compartment Mean Path Lengths

By direct analogy with Section 4, define the **mean path length in each compartment** for detected photons:

$$\langle L_m \rangle = \frac{\int\int \ell_m \cdot p(\ell_m, \ell_v) \cdot e^{-\mu_a^{(m)} \ell_m - \mu_a^{(v)} \ell_v} \; d\ell_m \; d\ell_v}{\int\int p(\ell_m, \ell_v) \cdot e^{-\mu_a^{(m)} \ell_m - \mu_a^{(v)} \ell_v} \; d\ell_m \; d\ell_v} \tag{15}$$

$$\langle L_v \rangle = \frac{\int\int \ell_v \cdot p(\ell_m, \ell_v) \cdot e^{-\mu_a^{(m)} \ell_m - \mu_a^{(v)} \ell_v} \; d\ell_m \; d\ell_v}{\int\int p(\ell_m, \ell_v) \cdot e^{-\mu_a^{(m)} \ell_m - \mu_a^{(v)} \ell_v} \; d\ell_m \; d\ell_v} \tag{16}$$

These are the weighted means of $\ell_m$ and $\ell_v$ over the surviving (detected) photon population. Each depends on **both** absorption coefficients, because the survival weighting $e^{-\mu_a^{(m)} \ell_m - \mu_a^{(v)} \ell_v}$ couples the two compartments.

The total mean path length of detected photons is:

$$\langle L \rangle = \langle L_m \rangle + \langle L_v \rangle$$

### 6.5 Partial Derivatives of Log-Intensity

Taking the logarithm of Equation (14) and differentiating with respect to each absorption coefficient:

$$\frac{\partial \ln I}{\partial \mu_a^{(m)}} = -\langle L_m \rangle(\mu_a^{(m)}, \mu_a^{(v)}) \tag{17}$$

$$\frac{\partial \ln I}{\partial \mu_a^{(v)}} = -\langle L_v \rangle(\mu_a^{(m)}, \mu_a^{(v)}) \tag{18}$$

**Derivation of Equation (17):**

$$\frac{\partial \ln I}{\partial \mu_a^{(m)}} = \frac{1}{I} \frac{\partial I}{\partial \mu_a^{(m)}} = \frac{I_0 \int\int p(\ell_m, \ell_v) \cdot (-\ell_m) \cdot e^{-\mu_a^{(m)} \ell_m - \mu_a^{(v)} \ell_v} \; d\ell_m \, d\ell_v}{I_0 \int\int p(\ell_m, \ell_v) \cdot e^{-\mu_a^{(m)} \ell_m - \mu_a^{(v)} \ell_v} \; d\ell_m \, d\ell_v} = -\langle L_m \rangle$$

The derivation of Equation (18) is identical with $\ell_m$ replaced by $\ell_v$.

**Equations (17) and (18) are exact.** They are the multi-compartment generalization of Equation (7).

**Key reference:** Arridge (1999), "Optical tomography in medical imaging," *Inverse Problems* 15, R41–R93, Section 4. (Develops the Jacobian/sensitivity formalism for optical tomography, of which Equations 17–18 are special cases.)

### 6.6 The Two-Compartment Modified Beer–Lambert Law

Expand $\ln I(\mu_a^{(m)}, \mu_a^{(v)})$ as a first-order Taylor series around a reference state $(\mu_{a,0}^{(m)}, \mu_{a,0}^{(v)})$:

$$\ln I(\mu_a^{(m)}, \mu_a^{(v)}) \approx \ln I_0^{\text{ref}} - \langle L_m \rangle_0 \cdot (\mu_a^{(m)} - \mu_{a,0}^{(m)}) - \langle L_v \rangle_0 \cdot (\mu_a^{(v)} - \mu_{a,0}^{(v)}) \tag{19}$$

where the subscript "$0$" on the pathlengths indicates evaluation at the reference state, and:

$$\ln I_0^{\text{ref}} \equiv \ln I(\mu_{a,0}^{(m)}, \mu_{a,0}^{(v)})$$

is the log-intensity at the reference state (a known quantity if the reference state has been simulated or measured).

Rearranging:

$$\ln I \approx \underbrace{\ln I_0^{\text{ref}} + \langle L_m \rangle_0 \cdot \mu_{a,0}^{(m)} + \langle L_v \rangle_0 \cdot \mu_{a,0}^{(v)}}_{\equiv \; \ln I_0^{\text{eff}}} - \mu_a^{(m)} \cdot \langle L_m \rangle_0 - \mu_a^{(v)} \cdot \langle L_v \rangle_0 \tag{20}$$

Defining the **effective source term** $I_0^{\text{eff}}$ to absorb the reference-state constants:

$$\boxed{\ln I(\lambda, x_i) \approx \ln I_0^{\text{eff}}(\lambda) - \mu_a^{(m)}(\lambda) \cdot \langle L_m \rangle(\lambda, x_i) - \mu_a^{(v)}(\lambda) \cdot \langle L_v \rangle(\lambda, x_i)} \tag{21}$$

**This is the two-compartment modified Beer–Lambert law.** It states that the log-intensity is approximately a linear function of the absorption coefficients of each compartment, with the compartment mean path lengths as the coefficients.

### 6.7 What Each Term Means

| Term | Meaning |
|------|---------|
| $\ln I_0^{\text{eff}}(\lambda)$ | Effective calibration constant. Absorbs: true source power, detector responsivity, geometric coupling, and the reference-state attenuation. It is a free parameter in the inverse problem (one per wavelength). |
| $\mu_a^{(m)}(\lambda) \cdot \langle L_m \rangle(\lambda, x_i)$ | Attenuation due to the muscle compartment. The product of muscle absorption coefficient and the mean distance that detected photons travel through muscle. |
| $\mu_a^{(v)}(\lambda) \cdot \langle L_v \rangle(\lambda, x_i)$ | Attenuation due to the vessel compartment. The product of blood absorption coefficient and the mean distance that detected photons travel through the vessel. |

### 6.8 Validity Conditions

Equation (21) is a first-order approximation. It is accurate when:

1. **$|\Delta\mu_a^{(m)}| \cdot \sigma_{L_m} \ll 1$**: The change in muscle absorption from the reference state is small relative to the spread of muscle path lengths.

2. **$|\Delta\mu_a^{(v)}| \cdot \sigma_{L_v} \ll 1$**: The change in blood absorption from the reference state is small relative to the spread of vessel path lengths.

3. **Cross-terms are small**: $|\Delta\mu_a^{(m)}| \cdot |\Delta\mu_a^{(v)}| \cdot \text{Cov}(L_m, L_v)$ is negligible.

In practice, condition (1) is usually satisfied because muscle absorption is low and varies little. Condition (2) is more restrictive: blood absorption varies strongly with wavelength and SvO₂, and the vessel path length distribution may be broad relative to the mean. This is why the linearization can fail for large changes in SvO₂ (see Section 8).

### 6.9 The Compartment Path Lengths Are Not DPF × ρ

The differential pathlength factor (DPF) is defined for a homogeneous medium as:

$$\text{DPF} = \frac{\langle L \rangle}{\rho}$$

In a two-compartment medium, there is no single DPF. The quantities $\langle L_m \rangle$ and $\langle L_v \rangle$ are:

- **Not** proportional to $\rho$ in any simple way.
- **Not** derivable from diffusion theory or any analytical formula for the geometry in question (vessel + muscle).
- **Not** related to each other by a simple volume-fraction scaling.
- **Dependent** on the specific geometry: vessel position, depth, radius, and source/detector locations.

$\langle L_m \rangle$ and $\langle L_v \rangle$ must be computed numerically, either from Monte Carlo simulation or from a full numerical solution of the radiative transport equation for the specific geometry.

In particular, $\langle L_v \rangle$ depends strongly on whether the photon "banana" (the most probable photon path region between source and detector) intersects the vessel, which is a geometric question that varies with the lateral position of the source–detector pair relative to the vessel.

---

## 7. Computing Compartment Path Lengths from Monte Carlo Simulation

### 7.1 Monte Carlo Output: Collected Fluence Rate

Monte Carlo simulation tracks individual photon packets through the voxelized tissue geometry. With the "collected-only" deposition criterion (used in our MCmatlab simulations), the fluence rate $\Phi(\mathbf{r})$ recorded at each voxel $\mathbf{r}$ represents the total weight deposited by photon packets that eventually reach the detector, per unit volume, per launched photon.

Physically, $\Phi(\mathbf{r}) \cdot \Delta V$ (where $\Delta V$ is the voxel volume) is proportional to the total path length spent by detected photons in the voxel at $\mathbf{r}$.

### 7.2 Compartment-Integrated Fluence

Summing the fluence over all voxels belonging to compartment $j$:

$$\Phi_j^{\text{int}} = \sum_{\mathbf{r} \in j} \Phi(\mathbf{r}) \cdot \Delta V \tag{22}$$

This gives the total path length (weighted by detection probability) spent by all detected photons in compartment $j$.

### 7.3 Mean Path Length per Compartment

The collected intensity $I_{\text{det}}$ is the total detected power per launched photon (the `collectedIntensity` output from MCmatlab). The mean path length in compartment $j$ for detected photons is:

$$\langle L_j \rangle = \frac{\Phi_j^{\text{int}}}{I_{\text{det}}} \tag{23}$$

This is the MC estimate of the quantities defined in Equations (15) and (16).

**Important:** These path lengths are evaluated at the specific absorption coefficients used in the MC simulation (the reference state). They are the coefficients for the linearization around that reference state, not around zero absorption.

### 7.4 Practical Computation in MATLAB

From the MCmatlab simulation outputs:

```
voxel_volume = (Lx/nx) * (Ly/ny) * (Lz/nz);

% Compartment masks
muscle_mask = (model.MC.M == 5);
vessel_mask = (model.MC.M == 4);

% Compartment-integrated fluence
Phi_muscle = sum(detectedFluence(muscle_mask)) * voxel_volume;
Phi_vessel = sum(detectedFluence(vessel_mask)) * voxel_volume;

% Mean path lengths
L_muscle = Phi_muscle / collectedIntensity;
L_vessel = Phi_vessel / collectedIntensity;
```

---

## 8. The Circularity Problem and Solutions

### 8.1 The Problem

The two-compartment MBLL (Equation 21) uses $\langle L_m \rangle$ and $\langle L_v \rangle$ as known constants in the forward model. But as shown in Section 4.1, these path lengths depend on the absorption coefficients — including $\mu_a^{(v)}$, which depends on SvO₂, the quantity we are trying to determine.

Specifically: the MC simulation that computes $\langle L_v \rangle$ must assign an absorption coefficient to the vessel, which requires assuming an SvO₂ value. If the true SvO₂ differs from the assumed value, the path lengths used in the linearization are evaluated at the wrong reference state, introducing systematic error.

### 8.2 Magnitude of the Error

The error is proportional to the second-order terms in the Taylor expansion:

$$\text{Error} \approx \frac{1}{2} \text{Var}(L_v) \cdot (\Delta\mu_a^{(v)})^2$$

For a small vessel (1 mm diameter at 3 mm depth, with $\langle L_v \rangle \sim 0.01$–$0.1$ cm), the variance of $L_v$ is comparable to $\langle L_v \rangle^2$, so:

$$\text{Error} \sim \frac{1}{2} \langle L_v \rangle^2 \cdot (\Delta\mu_a^{(v)})^2$$

For a change in SvO₂ of ±15% (e.g., from the reference 0.75 to the true value 0.60), $\Delta\mu_a^{(v)}$ at 660 nm is approximately 1.5 cm⁻¹. With $\langle L_v \rangle \sim 0.05$ cm, the error term is approximately $\frac{1}{2}(0.05)^2(1.5)^2 \approx 0.003$ in log-intensity, or roughly 0.3%. This is small but not negligible compared to the vessel contrast (which might be 5–15%).

For larger vessels or larger SvO₂ deviations, the error grows quadratically and can become significant.

### 8.3 Solution: Lookup Table Approach

Instead of linearizing, run MC simulations at multiple SvO₂ values and store the detected intensity directly as a lookup table:

$$I_{\text{table}}(\lambda_k, x_i; S_n) \quad \text{for } n = 1, \ldots, N_S$$

where $S_n$ are the tabulated SvO₂ values (e.g., 0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90).

In the inverse problem, the forward model interpolates the table at the trial SvO₂:

$$I_{\text{predicted}}(\lambda_k, x_i; \text{SvO}_2) = \text{interp}(I_{\text{table}}(\lambda_k, x_i; S_1), \ldots, I_{\text{table}}(\lambda_k, x_i; S_{N_S}); \text{SvO}_2)$$

This approach captures the full nonlinear dependence of detected intensity on SvO₂, including the path-length-varies-with-absorption effect, without any linearization.

### 8.4 Hybrid Approach

The linearization error is dominated by the **vessel compartment** (high absorption, large $\Delta\mu_a$ with SvO₂ changes). The muscle compartment has low absorption that varies little, so the linearization is accurate for muscle.

A practical hybrid approach:

- Use the **lookup table** for the SvO₂ dependence (vessel compartment).
- Use the **MBLL linearization** for perturbations in muscle absorption around the reference state.

$$\ln I(\lambda_k, x_i) \approx \ln I_{\text{table}}(\lambda_k, x_i; \text{SvO}_2) - \Delta\mu_a^{(m)}(\lambda_k) \cdot \langle L_m \rangle(\lambda_k, x_i) \tag{24}$$

where $\Delta\mu_a^{(m)} = \mu_a^{(m)} - \mu_{a,\text{ref}}^{(m)}$ is the deviation of muscle absorption from the value used in the MC runs.

---

## 9. Application to SvO₂ Recovery

### 9.1 Blood Absorption as a Function of Oxygen Saturation

The absorption coefficient of blood depends on the hemoglobin composition:

$$\mu_a^{\text{blood}}(\lambda, S) = \frac{\ln(10) \cdot [\text{Hb}]}{MW_{\text{Hb}}} \left[S \cdot \varepsilon_{\text{HbO}_2}(\lambda) + (1-S) \cdot \varepsilon_{\text{HbR}}(\lambda)\right] + \mu_{a,\text{water}}(\lambda) \cdot W_{\text{blood}} \tag{25}$$

where:

- $S$ is the oxygen saturation (SvO₂ for venous blood, SaO₂ for arterial blood).
- $[\text{Hb}]$ is the total hemoglobin concentration (mM).
- $MW_{\text{Hb}} = 64{,}500$ g/mol is the molecular weight of hemoglobin.
- $\varepsilon_{\text{HbO}_2}(\lambda)$, $\varepsilon_{\text{HbR}}(\lambda)$ are the tabulated molar extinction coefficients of oxy- and deoxyhemoglobin (cm⁻¹ M⁻¹). Source: Prahl (omlc.org), compiled from Gratzer and Cope.
- $\mu_{a,\text{water}}(\lambda)$ is the absorption coefficient of water (cm⁻¹). Source: Hale & Querry (1973), Kou et al. (1993).
- $W_{\text{blood}} \approx 0.95$ is the water volume fraction in blood.

### 9.2 Spectral Sensitivity to SvO₂

Taking the derivative of blood absorption with respect to saturation:

$$\frac{\partial \mu_a^{\text{blood}}}{\partial S} = \frac{\ln(10) \cdot [\text{Hb}]}{MW_{\text{Hb}}} \left[\varepsilon_{\text{HbO}_2}(\lambda) - \varepsilon_{\text{HbR}}(\lambda)\right] \tag{26}$$

This derivative is wavelength-dependent. At wavelengths where $\varepsilon_{\text{HbO}_2} \neq \varepsilon_{\text{HbR}}$ (i.e., away from the isosbestic point near 800 nm), the blood absorption changes with saturation, providing spectral leverage for determining SvO₂.

At 660 nm: $\varepsilon_{\text{HbR}} \gg \varepsilon_{\text{HbO}_2}$, so $\partial\mu_a/\partial S < 0$ (absorption *decreases* with increasing saturation).

At 940 nm: $\varepsilon_{\text{HbO}_2} > \varepsilon_{\text{HbR}}$, so $\partial\mu_a/\partial S > 0$ (absorption *increases* with increasing saturation).

The ratio of detected signals at these two wavelengths therefore changes in opposite directions as SvO₂ varies, which is the basis for oximetry.

### 9.3 The Inverse Problem

Substituting Equation (25) into the two-compartment MBLL (Equation 21) at three wavelengths and 12 PD positions gives 36 equations. The unknowns are SvO₂, [Hb], muscle $\mu_a$ at each wavelength, and calibration constants. The system is solved by minimizing the weighted residuals between predicted and measured log-intensities using nonlinear least squares (Levenberg–Marquardt).

The mathematical details of the inverse problem are described in the companion document ("Technical Description of the Proposed Pipeline," Section 6).

---

## 10. Summary of the Derivation Chain

| Step | Result | Equation | Exact or Approximate | Key Reference |
|------|--------|----------|---------------------|---------------|
| 1. Microscopic Beer–Lambert | Survival probability: $e^{-\mu_a \ell}$ | — | Exact | Tsuchiya (2001) |
| 2. Total detected intensity | $I = I_0 \int p(\ell) e^{-\mu_a \ell} d\ell$ | (1) | Exact | Tsuchiya (2001) |
| 3. Log-intensity as Laplace transform | $\ln I = \ln I_0 + \ln \mathcal{L}\{p\}(\mu_a)$ | (3) | Exact | — |
| 4. Mean pathlength of detected photons | $\langle L \rangle = \int \ell \, p(\ell) e^{-\mu_a \ell} d\ell \; / \; \int p(\ell) e^{-\mu_a \ell} d\ell$ | (4) | Definition | Arridge et al. (1992) |
| 5. Pathlength = derivative of log-I | $d\ln I / d\mu_a = -\langle L \rangle$ | (7) | Exact | Sassaroli & Fantini (2004) |
| 6. One-compartment MBLL | $\ln I \approx \ln I_0 - \mu_a \langle L \rangle_0$ | (10) | First-order approximation | Delpy et al. (1988), Cope (1991) |
| 7. Two-compartment intensity | $I = I_0 \iint p(\ell_m, \ell_v) e^{-\mu_a^{(m)} \ell_m - \mu_a^{(v)} \ell_v} d\ell_m d\ell_v$ | (13) | Exact | Tsuchiya (2001) |
| 8. Compartment partial derivatives | $\partial \ln I / \partial \mu_a^{(j)} = -\langle L_j \rangle$ | (17–18) | Exact | Arridge (1999) |
| 9. Two-compartment MBLL | $\ln I \approx \ln I_0^{\text{eff}} - \mu_a^{(m)} \langle L_m \rangle - \mu_a^{(v)} \langle L_v \rangle$ | (21) | First-order approximation | — |
| 10. MC computation of pathlengths | $\langle L_j \rangle = \Phi_j^{\text{int}} / I_{\text{det}}$ | (23) | Numerical (MC noise) | Hiraoka et al. (1993) |
| 11. Blood absorption model | $\mu_a^{\text{blood}}(\lambda, S) = f(S, [\text{Hb}], \varepsilon)$ | (25) | Empirical model | Jacques (2013) |

The only approximation in the entire chain is the **first-order Taylor expansion** at Steps 6 and 9. Everything else is either exact or a definition. The quality of the final result depends on how well the linearization holds, which depends on how far the true tissue properties deviate from the reference state used to compute the path lengths.

---

## 11. Complete Reference List

1. Arridge, S.R. (1999). "Optical tomography in medical imaging." *Inverse Problems* 15(2), R41–R93.

2. Arridge, S.R., Cope, M., & Delpy, D.T. (1992). "The theoretical basis for the determination of optical pathlengths in tissue: temporal and frequency analysis." *Physics in Medicine & Biology* 37(7), 1531–1560.

3. Arridge, S.R., Hiraoka, M., & Schweiger, M. (1995). "Statistical basis for the determination of optical pathlength in tissue." *Physics in Medicine & Biology* 40(9), 1539–1558.

4. Bigio, I.J. & Fantini, S. (2016). *Quantitative Biomedical Optics.* Cambridge University Press.

5. Cope, M. (1991). "The development of a near infrared spectroscopy system and its application for non-invasive monitoring of cerebral blood and tissue oxygenation in the newborn infant." PhD thesis, University College London.

6. Delpy, D.T., Cope, M., van der Zee, P., Arridge, S.R., Wray, S., & Wyatt, J.S. (1988). "Estimation of optical pathlength through tissue from direct time of flight measurement." *Physics in Medicine & Biology* 33(12), 1433–1442.

7. Hiraoka, M., Firbank, M., Essenpreis, M., Cope, M., Arridge, S.R., van der Zee, P., & Delpy, D.T. (1993). "A Monte Carlo investigation of optical pathlength in inhomogeneous tissue and its application to near-infrared spectroscopy." *Physics in Medicine & Biology* 38(12), 1859–1876.

8. Jacques, S.L. (2013). "Optical properties of biological tissues: a review." *Physics in Medicine & Biology* 58(11), R37–R61.

9. Sassaroli, A. & Fantini, S. (2004). "Comment on the modified Beer–Lambert law for scattering media." *Physics in Medicine & Biology* 49(14), N255–N257.

10. Tsuchiya, Y. (2001). "Photon path distribution and optical responses of turbid media: theoretical analysis based on the microscopic Beer–Lambert law." *Physics in Medicine & Biology* 46(8), 2067–2084.
