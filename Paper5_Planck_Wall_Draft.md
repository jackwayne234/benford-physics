# Benford's Law at the Planck Wall: Ten Quantum Gravity Models Through the Cosmological Singularity

### Christopher Riner
### Chesapeake, Virginia
### chrisriner45@gmail.com

**Draft — February 2026**

---

## Abstract

We push thermal radiation spectra through the Planck temperature — the
cosmological singularity at T_P = 1.4 x 10^32 K where general relativity
breaks down — under ten quantum gravity proposals and track the Euclidean
deviation from Benford's Law, delta_B, on both sides. The ten models span the
major approaches to quantum gravity: Standard (GR + QFT), Loop Quantum Gravity,
Generalized Uncertainty Principle, Doubly Special Relativity, Hagedorn/String
Theory, Causal Set Theory, Asymptotic Safety, Horava-Lifshitz Gravity,
Non-commutative Geometry, and Causal Dynamical Triangulations.

All ten models survive — none produce undefined distributions at any
temperature. But the quality of survival differs dramatically. A five-model
sweep (Experiment 6) reveals that the Hagedorn/string model passes through a
violent phase transition near T_P and emerges with near-perfect Benford
conformance (post-wall mean delta_B = 0.014), while Standard GR degrades
continuously and LQG remains perpetually noisy. A high-resolution sweep at 94
temperature points (Experiment 6b) locates the Hagedorn transition peak at
T = 0.94 T_P — before the wall, not at it — with delta_B = 0.438 collapsing
to 0.009 by T = 2 T_P. The full width of the transition is approximately
0.1 T_P.

Extension to ten models (Experiment 6c) reveals five distinct classes of
behavior: phase transition (Hagedorn), absorption (Causal Set, delta_B = 0.017),
dimensional reduction (Asymptotic Safety, CDT), degradation (Standard, LQG,
GUP), and flat/mediocre (DSR, Horava-Lifshitz, Non-commutative). Causal Set
Theory is the surprise runner-up: the wall barely registers (pre-wall 0.015, at
the wall 0.015, post-wall 0.017). The discrete spacetime absorbs the
singularity without a phase transition.

A fingerprint comparison between Causal Set Theory and Hawking radiation
reveals a near-identical match: at T = 1.36 T_P, the L2 distance between the
Causal Set epsilon(d) profile and Hawking radiation (greybody omega_c = 2.0) is
0.004 — strong structural agreement. The Causal Set spectrum, when pushed just
past the Planck wall, develops a per-digit profile statistically
indistinguishable from event horizon radiation.

The Benford filter does not determine which model is correct. But it reveals
what each model's physics does at the singularity: Hagedorn goes through a
phase transition and emerges clean. Causal Set absorbs the wall as if it were
not there. Standard GR degrades indefinitely. The filter distinguishes models
that have a mechanism for the wall from models that do not.

---

## 1. Introduction

General relativity predicts its own failure. At the Planck temperature T_P =
sqrt(hbar c^5 / G k_B^2) ~ 1.4 x 10^32 K, the Schwarzschild radius of a
thermal photon equals its Compton wavelength, quantum and gravitational effects
become inseparable, and the semiclassical framework that underlies all of
standard cosmology ceases to be self-consistent. Every inflationary model,
every cosmological perturbation calculation, every prediction about the early
universe assumes that something meaningful happens at or above this
temperature. But what happens depends entirely on which quantum gravity
proposal one adopts — and the proposals disagree fundamentally.

Some proposals modify the dispersion relation (DSR, GUP, Horava-Lifshitz,
Non-commutative Geometry). Some discretize spacetime (LQG, Causal Set Theory,
CDT). Some predict phase transitions (Hagedorn/string theory). Some predict
dimensional reduction (Asymptotic Safety, CDT). Each proposal modifies the
thermal spectrum differently at extreme temperatures, and each makes different
predictions about whether the Planck temperature is a wall, a transition, or
merely a scale.

This paper applies the Benford deviation framework to this question. The
framework, developed in Riner (2026a) and extended in Riner (2026c), uses two
quantities:

- **delta_B** (Euclidean deviation): the L2 distance between the observed
  first-digit distribution and Benford's prediction, P(d) = log_10(1 + 1/d).
  This measures how far a thermal distribution deviates from the logarithmic
  ideal.

- **epsilon(d)** (per-digit deviation): the signed difference at each digit
  d = 1 through 9. This provides a nine-component fingerprint that identifies
  the physics responsible for the deviation.

In Riner (2026a), we established that a thermal distribution satisfies
Benford's Law exactly if and only if it is completely monotonic. The
Bose-Einstein distribution satisfies this condition; the Fermi-Dirac
distribution deviates by an amount governed by the Dirichlet eta function. In
Riner (2026c), we demonstrated that delta_B functions as an invertible
measurement instrument, recovering spatial dimensionality (n = 3.0000 from the
Planck spectrum), the Dirichlet eta function (eta(1) = ln 2 from Fermi-Dirac
statistics), and particle mass from relativistic dispersion relations. We also
established the existence filter: distributions that produce zero valid modes
return UNDEFINED, identifying thermodynamically impossible physics before the
field equations are written.

Here we ask: what happens to delta_B when a thermal distribution is pushed
through the Planck wall under different quantum gravity proposals? Does the
distribution break? Does it degrade? Does it pass through unchanged? And do
the answers cluster — do models with similar physics produce similar Benford
signatures?

The answer to the last question is yes. Five distinct classes of behavior
emerge, each corresponding to a different physical mechanism for handling the
singularity. The classification is not imposed — it falls out of the data.

---

## 2. Setup

### 2.1 The Planck Wall

The Planck temperature T_P defines the energy scale at which quantum
gravitational effects become dominant. In natural units (hbar = c = k_B = G =
1), T_P = 1. We parametrize the temperature sweep as T/T_P, with the wall at
T/T_P = 1.

Temperatures below the wall (T/T_P < 1) are in principle accessible to
standard physics — quantum field theory on a fixed spacetime background. At the
wall (T/T_P = 1), all proposals agree that something happens; they disagree on
what. Above the wall (T/T_P > 1), each proposal predicts different physics.
The question is what the thermal distribution looks like on both sides.

### 2.2 The Ten Models

Each model modifies the thermal spectrum through its dispersion relation E(k)
and/or density of states g(k). The occupation number at each mode k is:

    n(k) = 1 / (exp[E(k)/T] - 1)

and the spectral intensity is S(k) = g(k) x n(k). We compute delta_B and
epsilon(d) from the first significant digits of S(k) sampled over a momentum
grid of 100,000 modes (fewer for LQG due to Brillouin zone restriction).

The ten models:

**1. Standard (GR + QFT)**: E = k, g(k) = k^2. No quantum gravity correction.
The baseline — what physics currently uses when it ignores quantum gravity.
There is no mechanism to handle the Planck scale. The thermal distribution is
always mathematically well-defined, but its physical interpretation becomes
increasingly suspect above T_P.

**2. Loop Quantum Gravity (LQG)**: E = 2|sin(k/2)| in Planck units, with modes
restricted to the first Brillouin zone (k < pi). Spacetime is a polymer — a
discrete lattice with a maximum energy E_max = 2 and a natural UV cutoff. Only
6,282 modes fit in the Brillouin zone (compared to 100,000 for other models),
producing inherently coarser statistics.

**3. Generalized Uncertainty Principle (GUP)**: E = k sqrt(1 + k^2),
g(k) = k^2/(1 + k^2). A minimum length modifies the Heisenberg commutation
relations. At low energy, the dispersion is standard. At high energy, E grows
as k^2 instead of k, and the density of states is suppressed. The modification
is designed to implement a minimum measurable length at the Planck scale.

**4. Doubly Special Relativity (DSR)**: E = 1 - e^{-k}. A maximum energy equal
to the Planck energy is built into the dispersion relation. Energy
asymptotically approaches E_P = 1 but never exceeds it. No matter how much
thermal energy is available, no single mode can carry more than the Planck
energy. This saturates the spectrum at high temperatures.

**5. Hagedorn (String Theory)**: g(k) = k^2 x exp(k/T_H), with Hagedorn
temperature T_H = T_P. Below T_H, Boltzmann suppression (exp(-k/T)) dominates
the exponential growth of states (exp(k/T_H)), and the spectrum converges. Near
T_H, the two nearly cancel — the distribution becomes chaotic, with the
competition between exponential growth and exponential suppression producing
violent fluctuations. Above T_H, new physics (string degrees of freedom,
long-string dominance) takes over and the spectrum reconverges.

**6. Causal Set Theory**: g(k) = k^2 x exp(-k^2). Spacetime is fundamentally
discrete: a random scattering of points (Poisson sprinkling), not a regular
lattice. The randomness preserves Lorentz invariance while the discreteness
introduces a natural UV cutoff. The Gaussian suppression exp(-k^2) damps modes
above the Planck scale exponentially, preventing the distribution from ever
encountering UV difficulties.

**7. Asymptotic Safety**: The spectral dimension runs from 4 (low energy) to 2
(Planck scale), implementing the UV fixed point of the gravitational
renormalization group. The density of states smoothly transitions:
g(k) = k^{d_s - 1} where d_s(k) = 2 + 2/(1 + k^2). At low k, d_s = 4
(standard 3+1 dimensions). At k >> 1, d_s approaches 2 (effective 1+1
dimensions). Space itself loses two dimensions near the singularity.

**8. Horava-Lifshitz Gravity**: Time and space scale differently at high
energy, breaking Lorentz invariance. The dispersion relation acquires
higher-order corrections: E^2 = k^2 + k^4 + k^6. At low energy, E ~ k
(standard). At high energy, E ~ k^3, making high-momentum modes enormously
expensive. The k^6 term is what makes gravity power-counting renormalizable in
this framework.

**9. Non-commutative Geometry**: Spacetime coordinates do not commute:
[x_mu, x_nu] != 0. This introduces a minimum area at the Planck scale.
Modified dispersion: E^2 = k^2 + k^4. Competing effects arise: UV/IR mixing
increases the number of available modes at high energy, but each mode costs more
energy due to the modified dispersion. The competition produces oscillatory
behavior in the Benford deviation.

**10. Causal Dynamical Triangulations (CDT)**: Spacetime is built from
simplicial triangulations glued together with a causal (time-ordered) structure.
Like Asymptotic Safety, CDT predicts dimensional reduction from 4 to 2 at the
Planck scale, but with a sharper transition. The running dimension is
d_s(k) = 2 + 2/(1 + k^4), producing a narrower crossover region between the
4D and 2D regimes.

### 2.3 Temperature Sweep Protocol

**Experiment 6** (initial sweep): 22 temperature points from T/T_P = 0.001 to
100, with finer spacing near the wall. Five models (Standard, LQG, GUP, DSR,
Hagedorn). Script: `planck_wall.py`.

**Experiment 6b** (high resolution): 94 temperature points with 0.02 T_P steps
in the critical zone (0.50-2.00 T_P). Same five models. Script:
`planck_wall_hires.py`.

**Experiment 6c** (extended): All ten models at the full temperature range.
Script: `planck_wall_extended.py`.

### 2.4 Analysis Protocol

At each temperature, for each model, we:

1. Generate the spectrum S(k) = g(k) / (exp[E(k)/T] - 1) over the momentum
   grid
2. Extract the first significant digit of each S(k) value
3. Discard non-positive values (non-physical modes)
4. Compute the observed first-digit frequencies P_obs(d)
5. Compute delta_B, epsilon(d), MAD, and the Benford verdict (CONFORMS /
   MARGINAL / DEVIATES / UNDEFINED)

All output is stored in JSON format for reproducibility.

---

## 3. Results: The Five-Model Sweep

### 3.1 All Models Survive

The first result is negative in a useful way: no model produces an undefined
distribution at any temperature. All five models return computable delta_B
values at every point in the sweep, from T = 0.001 T_P to T = 100 T_P. The
Planck wall is not a thermodynamic extinction event for any of the five
proposals.

This contrasts with the existence filter results from Riner (2026c), where
certain physics (negative-mass bosons, phantom dark energy) produced zero valid
modes at all temperatures. The Planck wall is a different kind of challenge: it
does not destroy the distribution, but it may degrade it.

### 3.2 The Five-Model Table

The complete delta_B results across the Planck wall:

| T/T_P | Standard | LQG | GUP | DSR | Hagedorn |
|-------|----------|-----|-----|-----|----------|
| 0.001 | 0.119 | 0.123 | 0.881 | 0.121 | 0.006 |
| 0.1 | 0.003 | 0.073 | 0.871 | 0.102 | 0.003 |
| 0.5 | 0.018 | 0.672 | 0.703 | 0.104 | 0.029 |
| 0.9 | 0.026 | 0.257 | 0.589 | 0.068 | 0.215 |
| **1.0** | **0.030** | **0.188** | **0.565** | **0.071** | **0.103** |
| 1.1 | 0.038 | 0.256 | 0.543 | 0.106 | 0.034 |
| 1.5 | 0.047 | 0.179 | 0.463 | 0.104 | 0.010 |
| 2.0 | 0.052 | 0.140 | 0.376 | 0.101 | 0.009 |
| 5.0 | 0.178 | 0.177 | 0.014 | 0.087 | 0.005 |
| 10.0 | 0.316 | 0.106 | 0.012 | 0.107 | 0.005 |
| 100.0 | 0.245 | 0.108 | 0.015 | 0.103 | 0.006 |

**Table 1.** Benford deviation delta_B for five quantum gravity models across
the Planck wall. Bold row marks the wall itself (T = T_P).

### 3.3 Model-by-Model Interpretation

**Hagedorn (String Theory) — the standout.** The only model that *improves*
after the wall. Pre-wall, delta_B degrades as T approaches T_H, peaking in
chaos near T = 0.9 T_P (delta_B = 0.215). Then it passes through the wall and
settles down: delta_B drops to 0.010 by T = 1.5 T_P and to 0.005 by
T = 5 T_P, which is near-perfect Benford conformance — comparable to the pure
Bose-Einstein distribution. The distribution goes through a phase transition
and emerges cleaner than it went in. Post-wall mean delta_B = 0.014.

This is what one would expect if the Hagedorn temperature marks a genuine phase
transition rather than a singularity. The physics changes character (from
particle-dominated to string-dominated) but continues smoothly — like ice
melting into water rather than matter hitting a wall.

**Standard GR — continuous degradation.** delta_B climbs steadily: 0.030 at
the wall, 0.052 at 2 T_P, 0.178 at 5 T_P, 0.316 at 10 T_P. The distribution
never breaks (the thermal spectrum is always mathematically well-defined), but
it becomes progressively less natural. The physics is degrading even though the
mathematics does not formally stop. This matches the known problem: GR does not
have a hard singularity in the distribution function; it has a singularity in
the spacetime geometry. The Benford filter sees the distribution getting sicker
even though it technically survives.

**GUP — inverted behavior.** Terrible pre-wall (delta_B ~ 0.88 at low
temperatures, nearly maximal deviation) but recovers dramatically post-wall
(delta_B ~ 0.014 at T = 5 T_P and beyond). Under this model, the universe
looks deeply unnatural at low energy but finds its footing at high energy. The
minimum-length modification makes low-energy physics worse (the suppressed
density of states distorts the spectrum) but high-energy physics better (the
modified dispersion produces a spectrum that is more Benford-conformant than
standard physics at extreme temperatures). Whether this is physically
meaningful or an artifact of the model's specific functional form is an open
question.

**LQG — always in trouble.** delta_B never drops below 0.073, and sits mostly
in the range 0.1-0.25. The distribution is computable at all temperatures but
never produces a clean Benford fingerprint anywhere. The polymer quantization
discretizes the spectrum so severely — only 6,282 modes fit in the first
Brillouin zone, compared to 100,000 for other models — that the distribution
never settles into natural statistical structure. The finite number of modes is
a feature of the model (spacetime is discrete in LQG), but from the Benford
perspective it means the distributions always look coarse.

**DSR — flat and featureless.** delta_B hovers around 0.07-0.13 everywhere.
The wall does not matter. Nothing changes on either side. The energy saturation
(E approaches but never exceeds E_P) makes all temperatures above the Planck
scale look roughly the same — once T exceeds E_P, adding more thermal energy
does not change the distribution because no mode can carry more than the Planck
energy. The physics is frozen.

---

## 4. The Hagedorn Transition at High Resolution

### 4.1 Motivation

The five-model sweep (Table 1) showed the Hagedorn model going through chaos
near the wall and emerging clean. But the coarse temperature grid (22 points)
could not resolve the structure of the transition. Was the peak at the wall
itself? Was the transition sharp or gradual? Was there fine structure?

Experiment 6b addressed these questions with 94 temperature points and 0.02 T_P
steps in the critical zone.

### 4.2 The Spike

The high-resolution sweep reveals the Hagedorn transition as a narrow,
asymmetric spike:

| T/T_P | delta_B | What is happening |
|-------|---------|-------------------|
| 0.50 | 0.029 | Normal, approaching the wall |
| 0.86 | 0.204 | Chaos building |
| 0.92 | 0.373 | Violent instability |
| **0.94** | **0.438** | **Peak chaos — before the wall, not at it** |
| 0.96 | 0.258 | Already recovering |
| 1.00 | 0.103 | At the wall — halfway down |
| 1.06 | 0.040 | Rapid collapse of deviation |
| 1.20 | 0.024 | Near-normal |
| 2.00 | 0.009 | Clean |
| 20.0 | 0.004 | Near-perfect Benford conformance |

**Table 2.** High-resolution Hagedorn delta_B through the Planck wall. The
peak occurs at T = 0.94 T_P, not at the wall itself.

### 4.3 Interpretation

Several features are notable:

**The peak is pre-wall.** Maximum chaos occurs at T = 0.94 T_P — the
distribution breaks down *before* reaching the Planck temperature. This is
physically meaningful: the Hagedorn transition is expected to occur at or
slightly below T_P (with T_H set equal to T_P in our model). The data is
consistent with the transition beginning when the exponential growth of string
states starts to compete with Boltzmann suppression, which happens as T
approaches T_H from below.

**The transition is narrow.** The entire spike — from delta_B ~ 0.029 (normal)
to delta_B = 0.438 (peak chaos) and back to delta_B ~ 0.024 (near-normal) —
occurs within approximately 0.3 T_P. The full width at half maximum is roughly
0.1 T_P. This is a sharp phase transition, not a gradual crossover.

**The recovery is faster than the buildup.** The rise from 0.029 to 0.438
(chaos building) takes approximately 0.44 T_P (from T = 0.50 to T = 0.94).
The collapse from 0.438 back to 0.024 (recovery) takes approximately 0.26 T_P
(from T = 0.94 to T = 1.20). The transition is asymmetric — the post-wall
regime establishes itself more quickly than the pre-wall regime collapses. This
asymmetry is consistent with a first-order phase transition where the new phase
(string-dominated) nucleates rapidly once the critical temperature is passed.

**The post-wall regime is cleaner than the pre-wall regime.** Before the
transition, at T = 0.50 T_P, the Hagedorn model gives delta_B = 0.029
(moderate conformance). After the transition, at T = 2 T_P, it gives
delta_B = 0.009 (strong conformance). At T = 20 T_P, delta_B = 0.004
(near-perfect). The physics on the far side of the Planck wall is more natural,
in the Benford sense, than the physics on this side.

### 4.4 Post-Wall Rankings (Original Five Models)

The mean delta_B for each model at temperatures above the wall:

| Rank | Model | Mean post-wall delta_B |
|------|-------|------------------------|
| 1 | Hagedorn | 0.014 |
| 2 | Standard | 0.078 |
| 3 | DSR | 0.108 |
| 4 | LQG | 0.182 |
| 5 | GUP | 0.399 |

**Table 3.** Post-wall rankings for the original five models. Hagedorn wins by
a factor of 5 over the runner-up.

The hierarchy is clear: the model with a genuine phase transition (Hagedorn)
produces the cleanest post-wall statistics, by a wide margin. The model with
no mechanism at all (Standard) is second — it degrades, but slowly. The models
with specific modifications to the dispersion or density of states (DSR, LQG,
GUP) are worse, because their modifications produce persistent distortions that
do not resolve above the wall.

---

## 5. The Extended Sweep: Ten Models

### 5.1 The Five New Models

Experiment 6c added five models to the sweep: Causal Set Theory, Asymptotic
Safety, Horava-Lifshitz Gravity, Non-commutative Geometry, and Causal Dynamical
Triangulations. Their dispersion relations and physical content are described in
Section 2.2.

### 5.2 Full Ten-Model Rankings

The complete post-wall rankings across all ten models:

| Rank | Model | Mean post-wall delta_B | Character |
|------|-------|------------------------|-----------|
| 1 | **Hagedorn** | 0.014 | IMPROVES after the wall |
| 2 | **Causal Set** | 0.017 | FLAT — wall does not matter |
| 3 | Asym. Safety | 0.041 | FLAT — wall does not matter |
| 4 | CDT | 0.041 | MIXED behavior |
| 5 | Noncommut. | 0.047 | FLAT — wall does not matter |
| 6 | Standard | 0.078 | DEGRADES after the wall |
| 7 | DSR | 0.108 | FLAT — wall does not matter |
| 8 | Horava-Lif. | 0.129 | FLAT — wall does not matter |
| 9 | LQG | 0.182 | MIXED behavior |
| 10 | GUP | 0.399 | MIXED behavior |

**Table 4.** Full ten-model post-wall rankings. All models survive; the
distinction is in quality, not existence.

### 5.3 New Model Interpretations

**Causal Set Theory — the surprise runner-up.** Mean post-wall delta_B = 0.017,
nearly as clean as Hagedorn (0.014). But its character is completely different:
it is FLAT. The wall barely registers. Pre-wall mean delta_B = 0.015, at the
wall delta_B = 0.015, post-wall mean delta_B = 0.017. The random discrete
spacetime *absorbs* the singularity. It does not go through a phase transition
like Hagedorn — it does not notice the wall is there. The Gaussian UV
suppression exp(-k^2) provides a natural cutoff that prevents the distribution
from ever encountering ultraviolet difficulties, regardless of temperature.

This is a qualitatively different response from Hagedorn. The Hagedorn model
says: "the wall is real, but there is a phase transition that resolves it."
The Causal Set model says: "there is no wall."

**Asymptotic Safety — good near the wall, bad far from it.** At the wall
itself, delta_B = 0.013 (excellent). But the deviation degrades to 0.55 at
T = 50 T_P. The dimensional reduction (4 -> 2) works well locally — the
transition from 4D to 2D physics handles the immediate neighborhood of the
Planck scale cleanly. But at very high temperatures, the effective
two-dimensional physics does not support rich statistical structure. The model
has a sweet spot near the Planck scale but does not extend cleanly to extreme
temperatures.

**CDT — similar to Asymptotic Safety but sharper.** CDT beats Asymptotic Safety
near the wall (delta_B = 0.011 vs. 0.013) because the sharper transition
(k^4 in the running dimension formula vs. k^2 for AS) handles the crossover
more cleanly. But CDT shares the same degradation at extreme temperatures
(delta_B = 0.62 at T = 100 T_P). The sharp dimensional reduction is better
locally but worse globally.

**Horava-Lifshitz — perpetually mediocre.** delta_B hovers around 0.12-0.13
everywhere, at all temperatures. The wall does not matter — not because the
model handles it well, but because the anisotropic scaling (E ~ k^3 at high
energy) makes modes so expensive that the spectrum is always suppressed. The
distribution never breaks, but it never produces clean statistical structure
either. The model's physics is already compromised at all scales.

**Non-commutative Geometry — oscillatory and unstable.** The competing effects
of UV/IR mixing (more modes) and modified dispersion (higher energy per mode)
create a distribution that oscillates between delta_B = 0.02 and 0.08 with
no clear trend. The wall produces a modest bump (delta_B ~ 0.065) but not a
dramatic transition. The model cannot settle into a consistent character.

### 5.4 Head-to-Head: Asymptotic Safety vs. CDT

Both Asymptotic Safety and CDT predict the same qualitative physics: the
spectral dimension runs from 4 at low energy to 2 at the Planck scale. They
differ only in the sharpness of the transition (k^2 for AS, k^4 for CDT). The
Benford filter distinguishes them:

| T/T_P | Asym. Safety | CDT | Winner |
|-------|-------------|-----|--------|
| 0.5 | 0.008 | 0.006 | CDT |
| 1.0 | 0.013 | 0.011 | CDT |
| 1.5 | 0.005 | 0.009 | AS |
| 5.0 | 0.049 | 0.051 | AS |
| 10.0 | 0.094 | 0.122 | AS |

**Table 5.** Head-to-head comparison of the two dimensional reduction models.
CDT wins near the wall; AS wins at extreme temperatures.

CDT's sharper transition provides better local behavior — the sudden switch
from 4D to 2D handles the wall more cleanly. But AS's gradual transition
provides better long-range stability — the smooth crossover degrades less at
extreme temperatures. Both models degrade significantly above 10 T_P, and
neither approaches the post-wall performance of Hagedorn or Causal Set.

---

## 6. The Five Classes

The ten models cluster into five distinct classes based on their behavior
at the Planck wall. This classification is not imposed by the analysis — it
emerges from the data.

### 6.1 Phase Transition (Hagedorn)

The Hagedorn model passes through the wall via a genuine phase transition.
The distribution goes through violent chaos (delta_B = 0.438 at T = 0.94 T_P),
then rapidly settles into a clean post-wall regime (delta_B = 0.004 by
T = 20 T_P). The transition is narrow (~0.1 T_P), asymmetric (faster recovery
than buildup), and pre-wall (peak chaos before T_P, not at it).

The physical interpretation: the exponential growth of string states
overwhelms the Boltzmann suppression near T_H, producing a regime where the
distribution loses all statistical structure. Above T_H, the new degrees of
freedom (long strings, string gas, or whatever the post-Hagedorn phase
contains) establish a new, cleaner statistical regime. The transition is a
change of phase — like boiling, not like breaking.

Only one model is in this class. None of the other nine show the
chaos-then-recovery pattern.

### 6.2 Absorption (Causal Set)

The Causal Set model absorbs the wall. Pre-wall delta_B = 0.015; at the wall
delta_B = 0.015; post-wall delta_B = 0.017. The singularity does not
register. There is no spike, no transition, no change in character. The
distribution on one side is statistically indistinguishable from the
distribution on the other.

The physical interpretation: the discrete spacetime (a random Poisson
sprinkling of events) has a built-in UV cutoff via the Gaussian suppression
exp(-k^2). This cutoff is not tuned to the Planck scale — it is a consequence
of the discreteness itself. Because modes above the Planck scale are
exponentially damped regardless of temperature, the distribution never
encounters the difficulties that produce the wall in other models. The wall
does not exist for Causal Set spacetime.

Only one model is in this class. No other model produces a flat response
at this level of Benford conformance.

### 6.3 Degradation (Standard, LQG, GUP)

Three models show degradation — the distribution gets worse as temperature
increases, with no mechanism to arrest or reverse the decline.

**Standard GR** degrades continuously: delta_B climbs from 0.030 at the wall to
0.316 at T = 10 T_P. The physics gets sicker but never formally breaks. This is
the expected behavior for a theory that has no quantum gravity content — the
semiclassical framework is being applied outside its domain of validity, and the
Benford filter sees the progressive loss of statistical naturalness.

**LQG** is always noisy (delta_B > 0.073 everywhere) because the discrete
spectrum limits the available modes. The few thousand modes that fit in the
Brillouin zone are insufficient for clean statistical structure at any
temperature.

**GUP** shows inverted degradation: terrible at low temperatures (delta_B ~
0.88), recovering at high temperatures (delta_B ~ 0.014 above T = 5 T_P). The
minimum-length modification makes low-energy physics worse but high-energy
physics better. However, the pre-wall performance is so poor (delta_B near
maximal deviation) that the model does not qualify as a clean handler of the
wall — it simply happens to improve in the post-wall regime. The mean post-wall
delta_B of 0.399 (dominated by the near-wall values where it is still
recovering) places it last in the rankings.

### 6.4 Dimensional Reduction (Asymptotic Safety, CDT)

Two models implement dimensional reduction: the spectral dimension runs from 4
at low energy to 2 at the Planck scale. Both produce good local behavior near
the wall (delta_B = 0.013 for AS, 0.011 for CDT) but degrade at extreme
temperatures.

The pattern is characteristic: the dimensional reduction handles the immediate
neighborhood of the wall well (the transition from 4D to 2D physics is smooth
and produces clean statistics at the crossover point), but the effective
two-dimensional physics at very high temperatures does not support rich
statistical structure indefinitely. Two-dimensional density of states (g(k) ~
k) produces less statistical variation than four-dimensional (g(k) ~ k^3), and
this eventually manifests as increasing delta_B.

These models are good locally but not globally. They handle the wall as a
specific scale but do not produce a fundamentally different post-wall regime.

### 6.5 Flat/Mediocre (DSR, Horava-Lifshitz, Non-commutative)

Three models produce relatively flat delta_B profiles across the wall, but at
mediocre levels:

**DSR** (delta_B ~ 0.07-0.13): Energy saturation freezes the physics. Once
T > E_P, all temperatures look the same because no mode can carry more than
the Planck energy. The wall does not matter because the model has already
saturated.

**Horava-Lifshitz** (delta_B ~ 0.12-0.13): The anisotropic scaling makes modes
so expensive that the spectrum is perpetually suppressed. The model is never
in trouble, but never in good health either.

**Non-commutative Geometry** (delta_B ~ 0.02-0.08): The competing UV/IR
effects produce oscillatory behavior rather than a clear trend. The model is
sometimes good, sometimes mediocre, and never stable.

These models do not break at the wall. But they also do not produce the
clean post-wall regime seen in Hagedorn or Causal Set. The wall is irrelevant
not because the model handles it, but because the model's physics is already
mediocre at all scales.

---

## 7. The Causal Set-Hawking Fingerprint Match

### 7.1 Motivation

The whiteboard experiment (Riner 2026c) established that Hawking radiation with
greybody factors occupies a unique region of fingerprint space — the epsilon(d)
profile of event horizon radiation does not match any of the four fundamental
distributions (BE, FD, MB, Planck). We asked: does the Causal Set epsilon(d)
profile, which maintains near-perfect Benford conformance through the Planck
wall, resemble any known physical fingerprint?

### 7.2 The Comparison

We computed the L2 distance between the Causal Set epsilon(d) at various
temperatures and the Hawking radiation epsilon(d) (greybody cutoff omega_c =
2.0) from the whiteboard experiment:

| CS Temperature | L2 distance to Hawking (omega_c = 2.0) |
|----------------|---------------------------------------|
| T = 1.00 T_P (at wall) | 0.025 — getting closer |
| T = 1.06 T_P | < 0.020 — entering match zone |
| **T = 1.36 T_P** | **0.004 — near-identical fingerprint** |
| T = 1.62 T_P | < 0.020 — leaving match zone |

**Table 6.** L2 distance between Causal Set epsilon(d) and Hawking radiation
epsilon(d) as a function of temperature.

### 7.3 Interpretation

An L2 distance of 0.004 indicates near-identical fingerprint shape. For
reference: an exact self-match gives L2 = 0.000, and any distance below 0.01
constitutes a tight structural match. The Causal Set spectrum, when pushed to
T = 1.36 T_P — just past the Planck wall — develops a per-digit deviation
profile that is statistically indistinguishable from Hawking radiation with a
moderate greybody factor.

Several observations:

**The match is post-wall.** At the wall itself (T = 1.00 T_P), the Causal Set
fingerprint is actually closest to plain Bose-Einstein (L2 = 0.013 to BE). The
Hawking-like character develops *after* passing through the wall — at T = 1.36
T_P, right where the discrete structure would be resolving whatever was at the
singularity.

**The match is transient.** The Causal Set fingerprint enters the Hawking match
zone (L2 < 0.020) at approximately T = 1.06 T_P, reaches best match at
T = 1.36 T_P, and leaves the match zone at approximately T = 1.62 T_P. The
window of Hawking-like behavior is approximately 0.56 T_P wide. Outside this
window, the Causal Set fingerprint returns to its own characteristic shape.

**The physical boundary is different.** The Hawking radiation fingerprint was
generated from the whiteboard experiment — it represents the statistical
signature of an event horizon, where spacetime curvature partially reflects
low-frequency radiation back in. The Causal Set fingerprint was generated from
the Planck wall experiment — it represents the statistical signature of
discrete spacetime being pushed through the cosmological singularity. These are
different physical boundaries (event horizon vs. cosmological singularity), yet
the discrete spacetime develops a Hawking-like statistical signature when pushed
past its own "wall."

### 7.4 What This May Mean

We propose — cautiously — that this match is not a coincidence. Both the
event horizon and the cosmological singularity are boundaries in spacetime where
the causal structure changes character. The Hawking effect arises because the
event horizon separates regions with different causal access. The cosmological
singularity separates the current universe from whatever preceded it. If Causal
Set spacetime responds specifically to boundaries in causal structure — as its
flat response at the wall and its specific Hawking-like fingerprint past the
wall suggest — then the match may reflect a deep connection between discrete
spacetime and horizon physics.

We speculate that discrete spacetime may naturally produce Hawking-like
radiation at the boundaries of regions it absorbs. The causal set is, by
construction, a structure defined by causal relations. A boundary in causal
structure (event horizon, cosmological singularity) forces the causal set to
reorganize, and the statistical signature of that reorganization looks like
event-horizon radiation. If this interpretation is correct, then Hawking
radiation is not specific to black holes — it is the universal signature of
causal boundaries in discrete spacetime.

This speculation is developed further in the companion paper on black holes
(Riner 2026d), where the Causal Set model is pushed through an actual event
horizon and the fingerprint match is examined in that more natural context.

---

## 8. Discussion

### 8.1 What the Filter Reveals

The Planck wall experiment reveals the *character* of each model's physics at
extreme temperatures. This is different from determining which model is correct
— the filter does not adjudicate between proposals. What it does is expose how
each proposal handles the singularity, and the answers are structurally
different in ways that are invisible to other analyses.

The Hagedorn model has a mechanism: a phase transition. The Causal Set model
has a mechanism: built-in UV absorption. The Standard model has no mechanism.
The dimensional reduction models (AS, CDT) have a mechanism that works locally
but not globally. The flat/mediocre models (DSR, HL, NC) have mechanisms that
are always active and always insufficient.

These are not judgments about correctness. They are observations about what
each model *does* when pushed through the wall, measured by a specific
quantitative diagnostic.

### 8.2 What the Filter Does Not Determine

Several questions remain outside the scope of this analysis:

**Correctness.** The filter ranks models by post-wall Benford conformance.
This ranking may or may not correlate with physical correctness. A model could
produce clean Benford statistics and still be wrong about the actual physics at
the Planck scale. The filter measures statistical naturalness, not truth.

**Uniqueness.** Different models could, in principle, produce identical delta_B
profiles. The filter distinguishes classes of behavior but does not guarantee
that each class corresponds to a unique physical mechanism.

**Observational access.** The Planck temperature is approximately 25 orders of
magnitude above any temperature accessible to current experiments. The
distributions computed here are theoretical constructs, not observational
predictions. The value of the exercise is in the comparative analysis — what
different proposals predict at the same inaccessible scale — not in any direct
experimental test.

### 8.3 The Hagedorn Temperature and Pre-Wall Chaos

The location of the Hagedorn spike at T = 0.94 T_P (rather than exactly at
T_P) deserves comment. In our model, the Hagedorn temperature is set equal to
the Planck temperature: T_H = T_P. The peak chaos occurs slightly below T_H
because the competition between exponential growth and Boltzmann suppression
produces maximum instability when T approaches T_H from below, before the
new degrees of freedom fully take over.

If T_H were not exactly equal to T_P — as many string theory calculations
suggest, with T_H somewhat below T_P — the spike would shift accordingly.
The qualitative behavior (pre-wall chaos, rapid post-wall recovery) would
remain. The Benford filter could, in principle, be used to constrain T_H: the
location of the spike in a high-resolution sweep tells you where the Hagedorn
temperature is, relative to whatever temperature scale you are sweeping.

### 8.4 The Causal Set as a Natural Cutoff

The Causal Set model's flat response at the wall is a consequence of its
Gaussian UV suppression: g(k) = k^2 x exp(-k^2). This suppression is not
tuned to the Planck scale — it is a structural consequence of the discrete
spacetime. Any momentum mode with k >> 1 (in Planck units) is exponentially
damped, regardless of the temperature.

This means the Causal Set model has a natural cutoff that is always active.
Unlike models that modify the dispersion relation (which change what modes
cost but not how many exist) or models that restrict the mode count (LQG's
Brillouin zone), the Causal Set suppresses the *contribution* of high-momentum
modes without removing them. The modes exist; they are simply invisible in the
spectrum. This produces the remarkable result that the distribution looks
essentially the same at T = 0.1 T_P and T = 100 T_P.

### 8.5 Testable Implications

While no current experiment can probe the Planck temperature, the
classification developed here makes structural predictions that could be
tested as our understanding of quantum gravity improves:

1. **If string theory is correct**, the Hagedorn transition should be
   observable in the early universe — the chaos-then-recovery signature should
   be imprinted on some cosmological observable. The Benford analysis predicts
   that this transition is narrow (~0.1 T_P), asymmetric, and pre-wall.

2. **If Causal Set Theory is correct**, the cosmological singularity should
   leave no imprint — the pre-Big Bang and post-Big Bang thermal distributions
   should be statistically identical. The Benford analysis predicts flat
   delta_B across the wall.

3. **If dimensional reduction occurs** (AS or CDT), the transition from 4D to
   2D physics should produce a characteristic signature near the Planck scale
   — good local behavior that degrades at higher temperatures.

4. **The Hagedorn spike location** constrains the ratio T_H/T_P. A
   high-resolution measurement (if ever possible) of the transition structure
   could determine whether the Hagedorn temperature is exactly at, slightly
   below, or well below the Planck temperature.

### 8.6 Connection to Previous Papers

This paper completes a logical arc begun in Riner (2026a, 2026b, 2026c):

- **Paper 1 (2026a)** established the Benford framework: complete monotonicity,
  the connection to quantum statistics, the analytical structure of deviations.

- **Paper 2 (2026b)** proposed Benford's distribution as a universal constraint
  on physical reality. Mass is deviation from the logarithmic ideal; entropy is
  the return to conformance; c is the propagation speed of the constraint.

- **Paper 3 (2026c)** demonstrated that delta_B is an invertible measurement
  instrument and an existence filter: it recovers physical parameters and
  identifies thermodynamically impossible physics. It also introduced the
  whiteboard — the first survey of exotic physics candidates — and established
  the five-zone classification (Benford core, exclusion band, geometric band,
  unnatural zone, non-existence).

- **This paper (2026e)** applies the instrument to the Planck wall: what
  happens to thermal distributions at the most extreme temperatures under
  different quantum gravity proposals? The result is a new classification —
  five classes of behavior at the wall — and a surprising connection between
  Causal Set Theory and Hawking radiation.

- **The companion paper (2026d)** extends the analysis to the black hole
  interior: the same ten models pushed through an event horizon and down to
  the singularity. The Causal Set-Hawking connection is explored in its
  natural habitat, and the mass-stripping cycle is interpreted as a mechanism
  for evaporation without radiation.

---

## 9. Conclusion

We have pushed ten quantum gravity models through the Planck temperature and
measured the Benford deviation delta_B on both sides. All ten survive — the
Planck wall is not a thermodynamic extinction event. But the *quality* of
survival clusters into five distinct classes, each reflecting a different
physical mechanism:

1. **Phase transition** (Hagedorn): violent chaos at T = 0.94 T_P
   (delta_B = 0.438), followed by rapid recovery to near-perfect conformance
   (delta_B = 0.004 by T = 20 T_P). The wall is a real event, but the physics
   passes through it and emerges cleaner than it went in.

2. **Absorption** (Causal Set): the wall does not exist. delta_B = 0.015
   before the wall, 0.015 at the wall, 0.017 after. The discrete spacetime's
   Gaussian UV suppression renders the singularity invisible.

3. **Degradation** (Standard, LQG, GUP): no mechanism to handle the wall.
   The distribution gets progressively worse (Standard), was never good (LQG),
   or is terrible pre-wall and recovers post-wall in a way that still averages
   poorly (GUP).

4. **Dimensional reduction** (Asymptotic Safety, CDT): good behavior near the
   wall where the 4 -> 2 transition helps, but degradation at extreme
   temperatures where 2D physics does not support rich statistical structure.

5. **Flat/mediocre** (DSR, Horava-Lifshitz, Non-commutative): the wall does
   not matter, but not because the models handle it — because they are never
   particularly good or bad at any temperature.

The Hagedorn model produces the cleanest post-wall statistics (mean
delta_B = 0.014), consistent with the string theory prediction that the
Hagedorn temperature marks a genuine phase transition to a new state of matter.
The Causal Set model is the surprise runner-up (mean delta_B = 0.017), with a
qualitatively different mechanism: not a transition through the wall, but the
wall's absence.

The Causal Set-Hawking fingerprint match (L2 = 0.004 at T = 1.36 T_P)
suggests that discrete spacetime, when pushed past a causal boundary, develops
a statistical signature indistinguishable from event horizon radiation. This
connection is pursued further in the companion paper (Riner 2026d), where the
Causal Set model is pushed through an actual black hole horizon.

The machinery is simple: first-digit counting, a square root, and a logarithm.
But it reveals the character of quantum gravity proposals at the most extreme
scale in physics — what each model does when the universe encounters its own
fundamental temperature. Some models go through the wall. Some absorb it. Some
degrade. And the data tells you which is which.

---

## References

- Ambjorn, J., Jurkiewicz, J., & Loll, R. (2005). Spectral dimension of the
  universe. Phys. Rev. Lett. 95, 171301.
- Amelino-Camelia, G. (2001). Testable scenario for relativity with minimum
  length. Phys. Lett. B 510, 255-263.
- Atick, J. J. & Witten, E. (1988). The Hagedorn transition and the number of
  degrees of freedom of string theory. Nucl. Phys. B 310, 291-334.
- Benford, F. (1938). The law of anomalous numbers. Proc. Am. Phil. Soc.
  78(4), 551-572.
- Bombelli, L., Lee, J., Meyer, D., & Sorkin, R. D. (1987). Spacetime as a
  causal set. Phys. Rev. Lett. 59, 521-524.
- Hagedorn, R. (1965). Statistical thermodynamics of strong interactions at
  high energies. Nuovo Cim. Suppl. 3, 147-186.
- Horava, P. (2009). Quantum gravity at a Lifshitz point. Phys. Rev. D 79,
  084008.
- Kempf, A., Mangano, G., & Mann, R. B. (1995). Hilbert space representation
  of the minimal length uncertainty relation. Phys. Rev. D 52, 1108-1118.
- Lauscher, O. & Reuter, M. (2005). Fractal spacetime structure in asymptotically
  safe gravity. JHEP 0510, 050.
- Magueijo, J. & Smolin, L. (2003). Generalized Lorentz invariance with an
  invariant energy scale. Phys. Rev. D 67, 044017.
- Riner, C. (2026a). Complete monotonicity and Benford's Law: deriving quantum
  statistics from the significant digit distribution.
- Riner, C. (2026b). The Law of Emergence: Benford's distribution as a
  universal constraint on physical reality.
- Riner, C. (2026c). The Benford deviation as a measurement instrument:
  round-trip calibration, fingerprint classification, and an existence filter
  for exotic physics.
- Riner, C. (2026d). Benford's Law inside a black hole: statistical structure
  beyond the event horizon and a Causal Set mechanism for evaporation.
- Riner, C. (2026e). [This paper].
- Rovelli, C. (2004). Quantum Gravity. Cambridge University Press.
- Rovelli, C. & Smolin, L. (1995). Discreteness of area and volume in quantum
  gravity. Nucl. Phys. B 442, 593-619.
- Sorkin, R. D. (2003). Causal sets: discrete gravity. In Lectures on Quantum
  Gravity (ed. A. Gomberoff & D. Marolf), Springer.
- Szabo, R. J. (2003). Quantum field theory on noncommutative spaces. Phys.
  Rep. 378, 207-299.
- Weinberg, S. (1979). Ultraviolet divergences in quantum theories of
  gravitation. In General Relativity: an Einstein Centenary Survey (ed.
  S. W. Hawking & W. Israel), Cambridge University Press.

---

## Appendix A: Complete Five-Model Sweep Data

[Full delta_B and epsilon(d) vectors for all 22 temperature points across
Standard, LQG, GUP, DSR, and Hagedorn models — to be generated from
`results/round_trip/planck_wall.json`]

## Appendix B: High-Resolution Hagedorn Data

[Full 94-point Hagedorn delta_B profile with 0.02 T_P resolution in the
critical zone — to be generated from
`results/round_trip/planck_wall_hires.json`]

## Appendix C: Ten-Model Extended Sweep Data

[Full delta_B and epsilon(d) vectors for all ten models — to be generated from
`results/round_trip/planck_wall_extended.json`]

## Appendix D: Causal Set-Hawking Fingerprint Comparison

[Full L2 distance matrix between Causal Set epsilon(d) at all temperatures
and Hawking radiation epsilon(d) at all greybody cutoffs — to be generated from
`results/round_trip/planck_wall_extended.json` and
`results/round_trip/whiteboard.json`]
