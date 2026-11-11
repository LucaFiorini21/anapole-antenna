# The Anapole Antenna

**An antenna that radiates on paper… but not in reality**

Edoardo Fera, Luca Fiorini
Università degli Studi di Roma Tor Vergata · Wireless Electromagnetic Technologies 

---

## Summary

We try to improve a short dipole at 104.2 MHz by adding a toroidal coil around it. The toroidal dipole radiates with the same $\sin\theta$ pattern and polarization as the electric dipole, so a first analysis that adds the two contributions predicts a doubled effective length and **+6 dB** of radiated power at the same current.

A complete analysis, in which the whole structure is treated as a single current distribution and every sign is tracked, shows the opposite: with our winding the toroidal term **subtracts** from the dipole. For the design value $h = k_0^2 N A R / 2$ the far field vanishes in every direction. The structure is a **dynamic anapole**: it carries real currents but does not radiate.

The winding direction alone decides between the two outcomes. With the opposite winding the contributions add and the +6 dB result is recovered. The error of the first analysis is not in the numbers, but in assuming the relative sign instead of computing it.

This repository contains a theoretical study. The result is checked with an independent numerical integration of the wire currents and, in [`simulations/`](simulations/), with a full-wave NEC-2 model (Section 11). The antenna has not been built.

## Contents

- `README.md`: this report.
- `docs/WET_Fera_Fiorini_anapolo.pptx`: the presentation slides.
- `simulations/`: full-wave verification (NEC-2), with code, `.nec` decks, figures and its own report.
- `LICENSE`: MIT.

---

## 1. Goal and parameters

| Quantity | Value |
|---|---|
| Design frequency $f_0$ | 104.2 MHz |
| Wavelength $\lambda_0 = c/f_0$ | 2.877 m |
| Wavenumber $k_0 = 2\pi/\lambda_0$ | 2.184 rad/m |

A short dipole ($h \approx 3$ cm $\ll \lambda_0$) radiates very little: $R_\text{rad} \approx 85$ mΩ. The idea is to add a second radiator with the same pattern, in order to raise the far field. The candidate is a toroidal dipole, built as a toroidal wire coil wound around the dipole.

## 2. Geometry and model

The structure has three parts:

- **Short dipole** along $z$, length $h$, with end plates so that the current $I$ is uniform: $\int \mathbf{J}\,dV = I h\,\hat z$.
- **Toroidal coil**: $N$ circular loops of radius $a$ (area $A = \pi a^2$), lying in meridian planes, with centres on a ring of radius $R$. Each loop is a magnetic dipole $\mathbf m_i = -I A\,\hat\varphi_i$. The minus sign comes from the winding: the current flows in $+z$ on the outer side of each loop when it flows in $+z$ in the dipole.
- **Return wire** along the ring. A continuous winding also makes one full turn around the ring, which is equivalent to an extra loop of radius $R$ with $m_z = I\pi R^2$. The return wire cancels it.

Assumptions: $kR \lesssim 0.5$, $a < R$, $N \ge 3$, time convention $e^{j\omega t}$. All loops carry the same current $I$.

Design: $R = 15$ cm, $N = 12$, $a = 4.7$ cm, so $A = 6.94\cdot10^{-3}$ m² and $kR = 0.33$.

## 3. Radiation vector

All far-field quantities follow from the radiation vector

$$\mathbf N(\theta,\varphi) = \int \mathbf J(\mathbf r')\,e^{jk\,\mathbf n\cdot\mathbf r'}\,dV',$$

where $\mathbf n$ is the direction of observation. Only the transverse part $\mathbf N_\perp = \mathbf N - \mathbf n(\mathbf n\cdot\mathbf N)$ radiates:

$$\mathbf E_\perp = -\frac{j\omega\mu_0\,e^{-jkr}}{4\pi r}\,\mathbf N_\perp,\qquad U = \frac{\eta k^2}{32\pi^2}|\mathbf N_\perp|^2,\qquad R_\text{rad} = \frac{\eta k^2}{16\pi^2}\int\frac{|\mathbf N_\perp|^2}{I^2}\,d\Omega,$$

with $\eta = 120\pi$ Ω. Since $\mathbf N$ is an integral over the whole structure, the contributions of its parts add as complex vectors, with their signs and phases. The modulus is taken only at the end, to compute power.

## 4. The analysis "on paper"

### Dipole

$\hat z_\perp = -\sin\theta\,\hat\theta$, so $|\mathbf N_\perp| = I h \sin\theta$. With $\int\sin^2\theta\,d\Omega = 8\pi/3$:

$$R_\text{rad} = 20\,k^2 h^2 \approx 85\ \text{mΩ}\quad(h = 29.8\ \text{mm}),\qquad D = 3/2\ (1.76\ \text{dBi}).$$

### Coil

Each loop, displaced by $\mathbf d_i = R\hat\rho_i$ from the axis, has a toroidal moment $\mathbf T_i = \tfrac12\,\mathbf d_i\times\mathbf m_i$ (derived in Section 5). The $N$ loops add along $z$, giving $|\mathbf T| = \tfrac12 N A R I$.

A toroidal dipole radiates like an equivalent electric dipole $(ik/c)\,\mathbf T$ [2], which in the radiation vector becomes $k^2\mathbf T$. Hence

$$|\mathbf N_{T\perp}| = k^2|\mathbf T|\sin\theta \quad\Rightarrow\quad h_T = \frac{k^2}{2}NAR = 29.8\ \text{mm}.$$

These are moduli: each source is characterized on its own, and for an isolated source the sign has no meaning.

### Sum

Two first-pass hypotheses:

- **A1**: the coil is in series with the dipole, so it carries the same current $I$ with the same phase.
- **A2**: the two sources have the same $\sin\theta$ pattern and polarization.

Adding the moduli gives $|\mathbf N_\perp| = (h + h_T)\,I\sin\theta$. With the design choice $h_T = h$, the effective length doubles to 59.6 mm.

| Quantity | Bare dipole | "Paper" antenna |
|---|---|---|
| Effective length | 29.8 mm | 59.6 mm |
| Far-field amplitude at fixed $I$ | 1 | ×2 |
| $R_\text{rad}$ | 85 mΩ | 339 mΩ |
| Radiated power at fixed $I$ | 1 | ×4 (+6 dB) |
| Directivity | 1.5 | 1.5 |

The directivity is unchanged because it depends only on the shape of the pattern: scaling $\mathbf N_\perp$ by a constant multiplies $U_\text{max}$ and $P_\text{rad}$ by the same factor.

### What this does not check

A1 states that the *currents* are in phase. For the *radiated contributions* it only guarantees that their ratio is real, i.e. 0° or 180°. Which of the two depends on the sign of $T_z$, set by the winding direction, and on the coefficient of the second-order term in the expansion of $\mathbf N$. The modulus of a sum is not the sum of the moduli: writing $h + h_T$ assumes the sign.

## 5. The complete analysis

### 5.1 Expansion of the radiation vector

The phase factor is expanded in a Taylor series:

$$\mathbf N = \sum_m \frac{(jk)^m}{m!}\int(\mathbf n\cdot\mathbf r')^m\,\mathbf J\,dV' = \int\mathbf J\,dV' + jk\int(\mathbf n\cdot\mathbf r')\mathbf J\,dV' - \frac{k^2}{2}\int(\mathbf n\cdot\mathbf r')^2\mathbf J\,dV' + \dots$$

The exponential series converges uniformly on the bounded source, so sum and integral can be exchanged. The condition $kR \lesssim 0.5$ is not needed for convergence; it makes higher orders small, so that the expansion can be truncated at second order (the residual is quantified in Section 6).

Each order contains different sources: order 0 the electric dipole, order 1 the magnetic dipole and electric quadrupole, order 2 the toroidal dipole and magnetic quadrupole. The second-order coefficient is $-k^2/2$, but the sign of the contribution is not fixed by it alone: $(\mathbf n\cdot\mathbf r')^2$ is positive, while $\mathbf J$ is a vector whose direction depends on the winding.

### 5.2 Orders 0 and 1

- **Order 0.** For closed loops $\oint d\mathbf l = 0$: only the dipole contributes, $\mathbf N^{(0)} = I h\,\hat z$.
- **Order 1.** Splitting $(\mathbf n\cdot\mathbf r)\mathbf J$ into antisymmetric and symmetric parts gives the magnetic dipole term $\mathbf N^{(1)}_M = jk\,\mathbf m\times\mathbf n$ with $\mathbf m = \tfrac12\int\mathbf r\times\mathbf J\,dV$, plus the electric quadrupole. The electric quadrupole vanishes for closed currents (Lemma 1) and, for the centred dipole, because $\int z J_z\,dz = 0$.
- For the loops, $\mathbf m = -IA\sum_i\hat\varphi_i = 0$ for $N \ge 2$. The extra moment $m_z = I\pi R^2$ of a continuous winding is cancelled by the return wire. Without it the term would dominate ($|\mathbf N_M|_\text{max} = kI\pi R^2 = 0.154\,I \approx 5.2\,Ih$); being $\varphi$-polarized, it could never cancel the dipole.

### 5.3 Lemmas for closed currents

For divergence-free currents ($\nabla\cdot\mathbf J = 0$), the divergence theorem applied to suitable products gives:

1. $\int(x_bJ_a + x_aJ_b)\,dV = 0$: no electric quadrupole, and $\int x_aJ_c\,dV = \varepsilon_{acd}m_d$.
2. With $Q_{ab,c} = \int x_ax_bJ_c\,dV$: $Q_{ab,c} + Q_{bc,a} + Q_{ca,b} = 0$.
3. $\int r^2\mathbf J\,dV = -2\int(\mathbf r\cdot\mathbf J)\mathbf r\,dV$.

The lemmas apply to the coil only. The dipole carries charge on its plates; its own second-order term is of order $(kh)^2/24 \approx 2\cdot10^{-4}$ and is neglected.

With the standard definition $\mathbf T = \tfrac{1}{10}\int[(\mathbf r\cdot\mathbf J)\mathbf r - 2r^2\mathbf J]\,dV$, Lemma 3 gives $\int(\mathbf r\cdot\mathbf J)\mathbf r\,dV = 2\mathbf T$ and $\int r^2\mathbf J\,dV = -4\mathbf T$.

### 5.4 Second order

In components, $N^{(2)}_c = -\tfrac{k^2}{2}n_an_bQ_{ab,c}$. Using Lemma 2, $Q_{ab,c}$ is rewritten in terms of the antisymmetric combinations

$$Q_{ab,c} - Q_{cb,a} = 2\,\varepsilon_{acd}M_{db},\qquad M_{db} = \tfrac12\int x_b(\mathbf r\times\mathbf J)_d\,dV,$$

which gives $\mathbf N^{(2)} = \tfrac{2k^2}{3}\,\mathbf n\times(M\mathbf n)$, automatically transverse.

Splitting $M = M^S + M^A$ ($\operatorname{tr}M = 0$), the antisymmetric part is described by a single vector $\mathbf w = \int[(\mathbf r\cdot\mathbf J)\mathbf r - r^2\mathbf J]\,dV = 6\mathbf T$, and $\mathbf n\times(M^A\mathbf n) = \tfrac14\mathbf w_\perp$. Therefore

$$\mathbf N^{(2)} = k^2\,\mathbf T_\perp + \frac{2k^2}{3}\,\mathbf n\times(M^S\mathbf n):$$

a toroidal dipole term, which radiates like an electric dipole, plus a magnetic quadrupole. The term $k^2\mathbf T$ is the one used in Section 4, now derived from the whole structure together with its sign.

### 5.5 Our geometry

$M^S$ is symmetric, traceless and invariant under rotations by $2\pi/N$; for $N \ge 3$ this forces $M^S = \operatorname{diag}(\alpha,\alpha,-2\alpha)$. Its $zz$ component, $M_{zz} = \tfrac12\int z\rho J_\varphi\,dV$, vanishes: the loops lie in meridian planes ($J_\varphi = 0$), the dipole lies on the axis and the return wire lies in the plane $z = 0$. Hence $M^S = 0$ and

$$\mathbf N_\perp = (Ih + k^2T_z)\,\hat z_\perp,\qquad \hat z_\perp = -\sin\theta\,\hat\theta.$$

Dipole and coil have the same pattern, as assumed in Section 4; the result now depends only on the sign of $T_z$.

### 5.6 Toroidal moment of one loop

On a closed path $\mathbf r\cdot d\mathbf l = d(r^2/2)$, so $\oint(\mathbf r\cdot d\mathbf l)\mathbf r = -\tfrac12\oint r^2d\mathbf l$ and $\mathbf T = -\tfrac{I}{4}\oint r^2\,d\mathbf l$. For a circular loop of radius $a$ centred at $\mathbf d$, with normal $\hat n$, this gives exactly

$$\mathbf T = \tfrac12\,\mathbf d\times\mathbf m,\qquad \mathbf m = IA\,\hat n.$$

With $\mathbf m_i = -IA\,\hat\varphi_i$, $\mathbf d_i = R\hat\rho_i$ and $\hat\rho\times\hat\varphi = \hat z$:

$$\mathbf T_i = -\tfrac12 RAI\,\hat z.$$

### 5.7 The anapole condition

Summing the $N$ identical loops, $T_z = -\tfrac12 NARI$, and

$$\mathbf N_\perp = \Big(h - \frac{k^2}{2}NAR\Big)I\,\hat z_\perp.$$

The toroidal term subtracts. The far field vanishes in every direction when

$$\boxed{h = \frac{k_0^2}{2}NAR}$$

which at 104.2 MHz is exactly the design value $h = 29.8$ mm. Then $R_\text{rad} = 20k^2(h - h_T)^2 = 0$.

| | "Paper" (sign assumed) | Complete (signs tracked) |
|---|---|---|
| Toroidal term | $+k^2NAR\,I/2$ | $-k^2NAR\,I/2$ |
| $\lvert\mathbf N_\perp\rvert$ | $(h + h_T)\,I\sin\theta$ | $(h - h_T)\,I\sin\theta$ |
| At $h_T = h$ | $2h\,I\sin\theta$ | 0 |
| $R_\text{rad}$ | 339 mΩ (+6 dB) | 0 (≈ 12 µΩ with all orders) |

**The winding direction decides.** Reversing the winding (or swapping the coil terminals in the series connection) flips the sign of $T_z$: the terms add, and the numerical check gives a power ratio of 3.95 (+6.0 dB). With our winding they cancel.

## 6. Verification

### Numerical check

The far field was computed by direct numerical integration of all wire currents, without any series expansion, over 4368 directions (Python; not a NEC-2 simulation). The residual is the peak $|\mathbf N_\perp|/(Ih)$, so the bare dipole is 0 dB.

| $R$ / $N$ / $a$ | $A$ | $h = k^2NAR/2$ | $kR$ | Peak residual | Tuned $h$ | Residual (tuned) |
|---|---|---|---|---|---|---|
| 10 cm / 12 / 3.0 cm | 2.83·10⁻³ m² | 8.09 mm | 0.22 | −44 dB | 8.05 mm (−0.5%) | −57 dB |
| 15 cm / 12 / 4.7 cm | 6.94·10⁻³ m² | 29.8 mm | 0.33 | −37 dB | 29.46 mm (−1.1%) | −50 dB |
| 20 cm / 16 / 6.0 cm | 11.3·10⁻³ m² | 86.3 mm | 0.44 | −32 dB | 84.65 mm (−1.9%) | −44 dB |

Truncating the numerical expansion at second order leaves only $7\cdot10^{-5}$ (−83 dB), which confirms the derivation. $T_z = -\tfrac12NAR\,I$, $\mathbf m = 0$ and $M^S = 0$ were also checked numerically.

### Beyond second order, by hand

In the equatorial plane ($\theta = 90°$) each loop integral can be done exactly, $N_z^{(i)} = 2\pi j\,a\,I\,J_1(ka\,u_i)\,e^{jkRu_i}$ with $u_i = \cos(\varphi - \varphi_i)$. Expanding to fourth order and summing over the loops ($\sum u_i^2 = N/2$, $\sum u_i^4 = 3N/8$, odd powers vanish):

$$N_z = -\frac{k^2NAR}{2}I\Big[1 - \frac{(kR)^2}{8} - \frac{3(ka)^2}{32}\Big],\qquad \varepsilon = \frac{(kR)^2}{8} + \frac{3(ka)^2}{32} = 1.44\%\ (-36.8\ \text{dB}).$$

The numerical value is 1.43% (−36.9 dB). The exact equatorial null is at $h = h_T(1-\varepsilon) = 29.4$ mm.

## 7. Bandwidth and tolerances

### Frequency

The coil term grows as $f^2$ while the dipole length is fixed, $h_T(f) = h\,(f/f_0)^2$:

$$\frac{R_\text{rad}(f)}{R_\text{dipole}(f)} = \big[1 - (f/f_0)^2\big]^2 \approx 4\Big(\frac{\Delta f}{f_0}\Big)^2.$$

The suppression is below −20 dB within ±5% (≈ ±5 MHz), −30 dB within ±1.6% and −40 dB within ±0.5%. Including all orders, the null moves up by about 0.5–0.7% (≈ 104.8 MHz) and $R_\text{rad}(f_0) \approx 12$ µΩ instead of 0.

### Construction errors

Let $\delta$ be the relative error on the coil term $h_T$ with respect to $h$. Then $|\mathbf N_\perp| = |\delta|\,hI\sin\theta$, i.e. a residual of $20\log_{10}|\delta|$ dB and $R_\text{rad} = 85\ \text{mΩ}\cdot\delta^2$. Since $h_T \propto k^2Na^2R$, $\delta \approx 2\Delta a/a + \Delta R/R + 2\Delta f/f_0$ (plus the error on $h$).

| Error on the coil term | Ideal: $20\log\delta$ | Ideal: $R_\text{rad}$ | Numerical (all orders, $h$ tuned) |
|---|---|---|---|
| ±0.1% | −60 dB | 85 nΩ | −48 dB |
| ±1% | −40 dB | 8.5 µΩ | −40 / −38 dB |
| ±5% | −26 dB | 0.21 mΩ | −27 / −26 dB |
| ±10% | −20 dB | 0.85 mΩ | −20 dB |

The fourth-order terms set a floor of about −50 dB ($3\cdot10^{-3}$ of $Ih$ after tuning $h$): deeper suppression is not reachable, however precise the construction.

## 8. Impedance, losses and practical limits

### Input impedance

$$Z_\text{in} = R_\text{loss} + j\Big(\omega L - \frac{1}{\omega C_a}\Big),\qquad L \approx \frac{\mu_0N^2A}{2\pi R} = 1.33\ \text{µH},\quad \omega_0L = 872\ \text{Ω}.$$

$R_\text{rad} \approx 0$ does not appear. The inductance formula holds for a closely wound toroid; ours is not, so the value is an order of magnitude only. The radiation efficiency $R_\text{rad}/(R_\text{rad} + R_\text{loss}) \to 0$: at $f_0$ the input power is dissipated or stored, and the antenna behaves as a dummy load. The near field is not zero: the large stored reactive energy gives a high $Q \approx \omega_0L/R_\text{loss}$, a narrow band and strong sensitivity to nearby objects.

### What must be controlled

- **Return wire**: it must cancel $m_z = I\pi R^2$ to ≲ 0.1%; a 1% mismatch already gives −26 dB ($\varphi$-polarized).
- **Uniform current**: the model assumes all loops carry the same $I$. Twelve loops in series need about 3.5 m of wire, more than $\lambda_0 = 2.88$ m, so in a real series winding the current varies in amplitude and phase along the coil. This is the main reason why the result must be re-checked with NEC-2 using the real connections.
- **Dipole–coil coupling**: the currents follow $I = YV$, so the feed voltages must be tuned.
- **Ohmic loss**: ≈ 1.5 Ω for 3.5 m of 2 mm copper wire at 104 MHz (skin effect), far above $R_\text{rad}$.

### Feeding options

- One generator: dipole and coil in series, $I_c = I_d = I$, with the return wire.
- Two generators: $I_c/I_d = 2h/(k^2NAR)$, a real ratio (0° or 180°) whose amplitude and phase can be tuned.
- $N$ independent loops fed in parallel: no advance and no return wire, at the cost of $N$ feed points.

## 9. Reception, uses and limits

**Reception.** By reciprocity, the same effective length governs transmission and reception: the open-circuit voltage induced by a plane wave is proportional to $\mathbf h_\text{eff}\cdot\mathbf E_\text{inc}$. At the anapole condition $\mathbf h_\text{eff} = 0$, so the structure neither transmits nor receives (apart from the residual of Sections 6–7). As an FM antenna it is useless by design.

**Where anapoles are used.**
- Nanophotonics: dielectric nanodisks that scatter very little but store a lot of energy, used for third-harmonic generation, sensors and nanolasers.
- Metamaterials: resonant transparency and scattering suppression [3].
- Narrow-band stealth: reduction of the radar cross-section of protruding elements, as suggested in [1].
- Theory: non-radiating sources show that the inverse source problem is not unique. A non-radiating current can always be added to a source without changing its far field [4].

**Limits.** Narrow band, tight tolerances, residual multipoles, ohmic loss, high $Q$. The cancellation must be designed into the object and cannot be applied to an existing metal body. We are not aware of commercial products; the topic is mostly research and teaching.

## 10. Relation to the literature

Nemkov, Basharin and Fedotov [1] show that electric and toroidal dipoles interfering destructively form a point-like non-radiating source, the dynamic anapole, with condition $\mathbf D = -ik\,\mathbf T_p\,(1 + O(a^2k^2))$. In our system the coil is the toroidal dipole and the central rod is the electric dipole. Our condition $\int\mathbf J\,dV = -k^2\mathbf T$ is their $\mathbf D = -ik\,\mathbf T_p$ with $\mathbf T_p = \mathbf T/c$; the conventions differ ($e^{j\omega t}$ versus $e^{-i\omega t}$), but the real condition $h = k_0^2NAR/2$ is the same. Like theirs, our result holds up to $O(a^2k^2)$ corrections, which are the $\varepsilon$ of Section 6.

## 11. Full-wave verification (NEC-2)

The folder [`simulations/`](simulations/) checks this report with exact field integrals and with NEC-2. The numbers are radiated power relative to the bare dipole with the same current, at 104.2 MHz.

| Model | Radiated power |
|---|---|
| Ideal currents, $h = 29.8$ mm / $29.46$ mm | −38.0 dB / −50.7 dB |
| Ideal currents, winding reversed | +6.0 dB |
| NEC-2, 12 independent loops, current ratio tuned | −47.8 dB |
| NEC-2, dipole, winding and return wire in series, one generator | +22 … +28 dB |
| NEC-2, two generators, 30 series capacitors in the winding, trimmed | −33.9 dB |

The derivation of Sections 3–7 is confirmed: $h_T = 29.8$ mm, the role of the winding direction, the −37 / −50 dB residuals and the shift of the null to 104.8 MHz. Three statements of this report are corrected by the simulations:

1. **Return wire (Sections 2 and 5.2).** The continuous winding has $m_z = I\pi(R^2 + a^2/2)$, not $I\pi R^2$. The return ring needs radius $\sqrt{R^2 + a^2/2} = 153.6$ mm; with radius $R$ the suppression stops at −12.6 dB.
2. **Uniform current (Section 2).** The current must be uniform inside each loop as well. A 0.10 λ loop fed at one point has a 6.5 % current variation, which gives the coil a spurious electric dipole moment larger than that of the central dipole.
3. **Feeding (Section 8).** The single-generator series option does not produce an anapole: the 5 m of wire carry a standing wave and the structure radiates more than the bare dipole. The two-generator option works if series capacitors keep the winding current uniform, with a −20 dB bandwidth of 290 kHz and a radiation efficiency of about −40 dB once losses are included.

With real currents the condition must be written on the moments, $\int J_z\,dV + k^2T_z = 0$ over the whole structure.

## Takeaway

**A model is only as good as its assumptions.**

1. On paper: same current ⇒ same phase ⇒ the fields add: $h_\text{eff} = 2h$, $R_\text{rad} = 339$ mΩ, +6 dB.
2. Complete analysis: the sign of $T_z$ is set by the winding, and $\mathbf N_\perp = (h - k^2NAR/2)\,I\,\hat z_\perp = 0$ for $h = k_0^2NAR/2$.
3. Numerics: −37 dB (−50 dB with $h$ tuned); the null shifts with frequency. In NEC-2 the anapole survives (−34 dB) only with two generators, series capacitors and trimming (Section 11).

Moduli characterize each source; signed vectors combine them; the squared modulus gives the power. The cancellation is a design condition: $h = k_0^2NAR/2$.

## References

1. A. A. Nemkov, A. A. Basharin, V. A. Fedotov, "Non-radiating sources, dynamic anapole, and Aharonov–Bohm effect", *Phys. Rev. B* **95**, 165134 (2017). arXiv:1605.09033.
2. N. Papasimakis, V. A. Fedotov, V. Savinov, T. A. Raybould, N. I. Zheludev, "Electromagnetic toroidal excitations in matter and free space", *Nature Materials* **15**, 263 (2016).
3. V. A. Fedotov, A. V. Rogacheva, V. Savinov, D. P. Tsai, N. I. Zheludev, "Resonant transparency and non-trivial non-radiating excitations in toroidal metamaterials", *Scientific Reports* **3**, 2967 (2013).
4. A. J. Devaney, E. Wolf, "Radiating and nonradiating classical current distributions and the fields they generate", *Phys. Rev. D* **8**, 1044 (1973).

## License

Released under the MIT License. See `LICENSE`.
