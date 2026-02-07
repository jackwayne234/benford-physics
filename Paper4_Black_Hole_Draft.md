# Benford's Law Inside a Black Hole: Statistical Structure Beyond the Event Horizon and a Causal Set Mechanism for Evaporation

### Christopher Riner
### Chesapeake, Virginia
### chrisriner45@gmail.com

**Draft — February 2026**

---

## Abstract

We present a computational experiment in which thermal radiation spectra, modified
by ten quantum gravity proposals, are evaluated at forty radial positions through
a Schwarzschild black hole — from far outside the event horizon, through the
horizon, down to the singularity, and into a post-singularity bounce region. At
each position, we compute the Euclidean deviation from Benford's Law, δ_B, and
the per-digit deviation profile ε(d), using the framework developed in Riner
(2026a, 2026b).

All ten models survive the journey — none produce undefined distributions at any
radius. But the quality of survival differs dramatically. Causal Set Theory
produces the cleanest statistical structure inside the black hole (mean δ_B =
0.011), outperforming all other models including Hagedorn/string theory (0.019).
The discrete spacetime of Causal Set Theory absorbs the singularity: the
distribution barely registers the event horizon and maintains near-perfect
Benford conformance throughout the interior.

A striking pattern emerges in the Causal Set data. Far from the black hole, δ_B
sits at Hawking radiation levels (~0.028). As the observer approaches, it drops
to near-perfect conformance (~0.003). Inside the horizon, it remains low but
gradually rises toward the singularity, approaching Hawking-like levels again
(~0.020) with a fingerprint shape (ε(d)) that quantitatively matches Hawking
radiation at ω_c = 2.0 (L2 distance = 0.004).

A companion experiment through a traversable wormhole — a geometry with
comparable curvature but no singularity and no horizon — shows Causal Set Theory
dropping from 1st place (black hole) to 9th place (wormhole, δ_B = 0.056). This
confirms that the Causal Set response is specific to singularities and horizons,
not to curvature in general.

We propose that this pattern reflects a physical mechanism: the black hole's
gravitational field strips the Causal Set below its natural equilibrium as it
approaches the horizon, then the black hole feeds mass-energy back into the
Causal Set to restore it inside. The transaction is one-way — the restored
spacetime does not return the investment. The black hole does not radiate; it
shrinks, as discrete spacetime consumes its mass-energy to maintain Benford
conformance. This reframes evaporation as a structural process rather than a
radiative one, naturally explains the non-detection of Hawking radiation as
outgoing particles, and dissolves the information paradox in its standard form.

We further propose that the product of this conversion — Causal Set structure
at equilibrium — is new spacetime geometry, suggesting that black holes are
spacetime factories: sites where mass is converted to geometry and the universe
literally grows. Gravitational waves from binary mergers may represent the
wavefront of this newly generated spacetime.

---

## 1. Introduction

Classical general relativity forbids extracting information from inside a black
hole. The event horizon is a one-way membrane: signals, particles, and light can
fall in but cannot escape. Any distribution of matter or radiation inside the
horizon is, by construction, unobservable from the outside.

But statistical structure does not require a signal. If a physical theory
specifies the form of a thermal distribution at a given radius — its occupation
numbers, its density of states, its dispersion relation — then one can compute
the first-digit statistics of that distribution and compare them against
Benford's Law, regardless of whether the result could ever be communicated to an
external observer. The question shifts from "what can we see?" to "what kind of
statistical structure exists there?"

This paper applies the Benford deviation framework developed in Riner (2026a) to
the interior of a Schwarzschild black hole. The framework uses two quantities:

- **δ_B** (Euclidean deviation): the L2 distance between the observed
  first-digit distribution and Benford's prediction, P(d) = log₁₀(1 + 1/d).
  This measures how far a distribution deviates from the logarithmic ideal.

- **ε(d)** (per-digit deviation): the signed difference at each digit d = 1
  through 9. This provides a shape — a fingerprint — that identifies what kind
  of physics produced the deviation.

In Riner (2026c), we demonstrated that δ_B functions as an invertible
measurement instrument: it recovers spatial dimensionality (n = 3 exactly from
the Planck spectrum), the Dirichlet eta function (η(1) = ln 2 from Fermi-Dirac
statistics), and particle mass from relativistic dispersion relations. It also
functions as an existence filter: distributions that produce zero valid modes
return UNDEFINED, identifying thermodynamically impossible physics (negative-mass
bosons, phantom dark energy) before the field equations are ever written down.

In a companion paper (Riner 2026e), we swept ten quantum gravity models through
the cosmological singularity at the Planck temperature (the "Big Bang wall") and
found that Causal Set Theory and Hagedorn/string theory produce the cleanest
post-wall statistical structure, with Causal Set Theory showing a remarkable
property: the singularity doesn't register. The distribution is equally clean on
both sides.

Here we ask: does the same hold at a gravitational singularity? Is the black
hole interior, from the perspective of Benford's Law, the same kind of wall as
the Big Bang?

---

## 2. Setup

### 2.1 Black Hole Geometry

We use a Schwarzschild black hole with Hawking temperature T_H = 0.05 T_P in
Planck units (corresponding to a black hole mass M ≈ 10 M_P). The
Schwarzschild radius r_s defines the event horizon.

Radial positions are specified as r/r_s, with r/r_s = 1 at the horizon,
r/r_s > 1 outside, and r/r_s < 1 inside. We sample 40 positions:

- **Outside** (20 points): r/r_s = 10.0, 7.0, 5.0, 3.0, 2.0, 1.5, 1.3, 1.2,
  1.15, 1.1, 1.08, 1.06, 1.04, 1.03, 1.02, 1.015, 1.01, 1.005, 1.002, 1.001

- **Inside** (20 points): r/r_s = 0.99, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5,
  0.4, 0.3, 0.25, 0.2, 0.15, 0.12, 0.1, 0.08, 0.06, 0.04, 0.02, 0.01

- **Post-singularity bounce** (10 points): r/r_s = −0.01 through −0.5,
  representing a continuation through the singularity as predicted by loop
  quantum gravity and string cosmology bouncing scenarios.

### 2.2 Observer Models

Two observer perspectives are modeled:

**Static observer** (hovering outside the horizon): The locally measured
temperature diverges at the horizon due to the Tolman redshift:

    T_local(r) = T_H / √(1 − r_s/r)

This observer exists only for r > r_s. At the horizon, T_local → ∞.

**Infalling observer** (freely falling through the horizon): By the equivalence
principle, a freely falling observer notices nothing special at the horizon. The
effective temperature experienced is:

    T_eff(r) = T_H × (r_s/r)^{3/2}

This is smooth at r = r_s, rises as the observer approaches the singularity
(T_eff → ∞ as r → 0), and provides a continuous temperature profile from
outside to inside. The infalling observer is the primary focus of this paper.

### 2.3 The Ten Quantum Gravity Models

Each model modifies the thermal spectrum through its dispersion relation E(k)
and/or density of states g(k). The occupation number at each mode k is:

    n(k) = 1 / (exp[E(k)/T] − 1)

and the spectral intensity is S(k) = g(k) × n(k). We compute δ_B and ε(d) from
the first significant digits of S(k) sampled over a momentum grid.

The ten models, their dispersion relations, and their physical content are:

1. **Standard (GR + QFT)**: E = k, g(k) = k². No quantum gravity correction.
   The baseline against which all others are compared.

2. **Loop Quantum Gravity (LQG)**: E = 2|sin(k/2)|, with modes restricted to
   the first Brillouin zone (k < π). Spacetime is a polymer — a discrete lattice
   with a maximum energy E_max = 2 and a natural UV cutoff.

3. **GUP (Generalized Uncertainty Principle)**: E = k√(1 + k²), g(k) = k²/(1 +
   k²). A minimum length scale modifies the commutation relations. Energy grows
   as k² at high momenta; the density of states is suppressed.

4. **DSR (Doubly Special Relativity)**: E = 1 − e^{−k}. A maximum energy
   (E_P = 1) is built into the dispersion relation. Energy saturates at the
   Planck scale regardless of momentum.

5. **Hagedorn (String Theory)**: g(k) = k² × exp(k/T_H), with Hagedorn
   temperature T_H = T_P. Below T_H, Boltzmann suppression dominates and the
   spectrum converges. Near T_H, the exponential growth of string states nearly
   cancels the Boltzmann factor, producing a phase transition.

6. **Causal Set Theory**: g(k) = k² × exp(−k²). Spacetime is a random discrete
   set of points (Poisson sprinkling). The Gaussian UV suppression damps modes
   above the Planck scale. Not a lattice — truly random, preserving Lorentz
   invariance at the fundamental level.

7. **Asymptotic Safety**: The spectral dimension runs from 4 (low energy) to 2
   (Planck scale). The density of states smoothly transitions: g(k) = k^{d_s−1}
   where d_s(k) = 2 + 2/(1 + k²). This implements the UV fixed point of the
   gravitational renormalization group.

8. **Horava-Lifshitz Gravity**: Anisotropic scaling between time and space. The
   dispersion relation acquires higher-order corrections: E² = k² + k⁴ + k⁶.
   This makes gravity power-counting renormalizable at the cost of Lorentz
   symmetry at high energies.

9. **Non-commutative Geometry**: Spacetime coordinates satisfy [x_μ, x_ν] ≠ 0,
   introducing a minimum area. Modified dispersion: E² = k² + k⁴. Competing
   effects: UV/IR mixing increases available modes, but each mode costs more
   energy.

10. **Causal Dynamical Triangulations (CDT)**: Spacetime is built from
    simplicial triangulations with a causal (time-ordered) structure. Like
    Asymptotic Safety, predicts dimensional reduction 4→2, but with a sharper
    transition: d_s(k) = 2 + 2/(1 + k⁴).

### 2.4 Analysis Protocol

At each radial position, for each model, we:

1. Compute T_eff(r) from the infalling observer temperature profile
2. Generate the spectrum S(k) = g(k) / (exp[E(k)/T_eff] − 1) over 100,000
   momentum modes (fewer for LQG due to Brillouin zone restriction)
3. Extract the first significant digit of each S(k) value
4. Compute δ_B, ε(d), MAD, and the Benford verdict (CONFORMS / MARGINAL /
   DEVIATES / UNDEFINED)

This yields 10 models × 50 positions = 500 spectral evaluations.

---

## 3. Results

### 3.1 All Models Survive

No model produces an undefined distribution at any radius. All 500 spectral
evaluations return computable δ_B values. The black hole interior is, from the
Benford perspective, a valid physical environment for all ten quantum gravity
proposals.

This is itself a result. The existence filter developed in Riner (2026c)
identified physics that *cannot* form thermal distributions (negative-mass
bosons, phantom energy — zero valid modes everywhere). Black hole interiors are
not in that category. Whatever happens inside a black hole, it is not
thermodynamic non-existence.

### 3.2 Rankings: Inside the Black Hole

The mean δ_B for each model across all interior positions (r/r_s < 1),
using the infalling observer temperature profile:

| Rank | Model | Mean δ_B (inside) | Mean δ_B (outside) | Character |
|------|-------|-------------------|--------------------|-----------|
| 1 | **Causal Set** | **0.011** | 0.005 | FLAT — horizon irrelevant |
| 2 | Hagedorn | 0.019 | 0.006 | Slight degradation |
| 3 | Noncommut. | 0.051 | 0.043 | FLAT |
| 4 | Standard | 0.108 | 0.132 | FLAT |
| 5 | DSR | 0.110 | 0.098 | FLAT |
| 6 | Asym. Safety | 0.112 | 0.222 | MIXED |
| 7 | CDT | 0.115 | 0.219 | MIXED |
| 8 | Horava-Lif. | 0.125 | 0.124 | FLAT |
| 9 | LQG | 0.173 | 0.060 | DEGRADES |
| 10 | GUP | 0.605 | 0.877 | MIXED |

**Table 1.** Mean Benford deviation inside and outside the event horizon. Models
ranked by interior performance.

Causal Set Theory produces the cleanest statistical structure inside the black
hole by a significant margin — nearly half the deviation of the runner-up
(Hagedorn). The discrete spacetime maintains near-perfect Benford conformance
(δ_B = 0.011, well within the "strong conformance" threshold of 0.02) even as
the effective temperature rises toward infinity at the singularity.

### 3.3 The Causal Set Journey

The Causal Set data tells a story best read as a journey from far outside to
deep inside:

| Region | r/r_s | T_eff (T_P) | δ_B | Verdict |
|--------|-------|-------------|-----|---------|
| Far outside | 10.0 | 0.0016 | 0.028 | MARGINAL |
| Approaching | 5.0 | 0.0045 | 0.004 | CONFORMS |
| Near horizon | 2.0 | 0.018 | 0.003 | CONFORMS |
| At horizon | 1.001 | 0.050 | 0.004 | CONFORMS |
| Just inside | 0.99 | 0.051 | 0.004 | CONFORMS |
| Mid-interior | 0.5 | 0.141 | 0.005 | CONFORMS |
| Deep interior | 0.2 | 0.559 | 0.012 | CONFORMS |
| Near singularity | 0.1 | 1.581 | 0.017 | CONFORMS |
| Very near sing. | 0.04 | 6.25 | 0.017 | CONFORMS |
| At singularity | 0.01 | 50.0 | 0.015 | CONFORMS |
| Post-singularity | −0.01 | 50.0 | 0.015 | CONFORMS |

**Table 2.** Causal Set δ_B through the complete black hole journey.

Several features are notable:

**The horizon is invisible.** The transition from r/r_s = 1.001 (outside) to
0.99 (inside) produces no discontinuity, no spike, no change in character. δ_B
= 0.004 on both sides. The event horizon — the defining feature of a black
hole — does not register in the Causal Set statistical structure. This is
consistent with the equivalence principle: a freely falling observer should
notice nothing special at the horizon.

**Far-field Hawking signature.** At r/r_s = 10, far from the black hole, CS
shows δ_B = 0.028. This is within the range of Hawking radiation fingerprints
measured on the whiteboard (δ_B = 0.020–0.035 for greybody cutoffs ω_c = 0.5
to 5.0). Far from the hole, the CS spectrum carries the imprint of the Hawking
thermal bath.

**Approach to conformance.** As the observer falls inward from r = 10 r_s, δ_B
drops from 0.028 to 0.003 — nearly perfect Benford conformance. The distribution
*improves* on approach, reaching its cleanest state just outside the horizon.

**Gradual interior rise.** Inside the horizon, δ_B rises slowly from 0.004 to
approximately 0.017 near the singularity. This is a controlled increase — the
distribution never leaves the CONFORMS category (δ_B < 0.02 throughout). Even at
r = 0.01 r_s, where T_eff = 50 T_P, the Causal Set spectrum maintains strong
Benford conformance.

**Post-singularity mirror.** The bounce region (r/r_s < 0) mirrors the approach:
δ_B = 0.015 at r = −0.01, matching the pre-singularity value exactly. The
singularity, if it is resolved by a bounce, is symmetric in its Benford
structure.

### 3.4 Fingerprint Match: Causal Set and Hawking Radiation

The per-digit deviation ε(d) provides a richer comparison than δ_B alone. In
Riner (2026e), we compared the ε(d) fingerprint of Causal Set Theory at various
temperatures with Hawking radiation (greybody factor ω_c = 2.0) from the
whiteboard experiment.

The result: the Causal Set fingerprint converges to the Hawking radiation
fingerprint — but not at the horizon. Just past it.

| CS Location | L2 distance to Hawking (ω_c = 2.0) |
|-------------|-------------------------------------|
| At wall (T = 1.00 T_P) | 0.025 |
| Past wall (T = 1.06 T_P) | < 0.020 |
| **Best match (T = 1.36 T_P)** | **0.004** |
| Further past (T = 1.62 T_P) | < 0.020 |

**Table 3.** L2 distance between Causal Set ε(d) and Hawking radiation ε(d).

An L2 distance of 0.004 indicates near-identical fingerprint shape. For
reference, the self-match distance is 0.000, and any distance below 0.01
constitutes a tight structural match. The Causal Set spectrum, when pushed past
the singularity into the regime where T > T_P, develops a per-digit profile that
is statistically indistinguishable from Hawking radiation with a moderate
greybody factor.

### 3.5 Comparison with the Cosmological Singularity

Eight of ten models produce comparable δ_B values at both the Big Bang
(cosmological singularity) and the black hole (gravitational singularity):

| Model | Big Bang (post-wall) | Black Hole (inside) | Ratio | Consistent? |
|-------|---------------------|--------------------:|------:|-------------|
| Causal Set | 0.017 | 0.011 | 0.64 | YES (better at BH) |
| Hagedorn | 0.014 | 0.019 | 1.35 | YES |
| Noncommut. | 0.047 | 0.051 | 1.08 | YES |
| Standard | 0.078 | 0.108 | 1.39 | YES |
| DSR | 0.108 | 0.110 | 1.02 | YES |
| Horava-Lif. | 0.129 | 0.125 | 0.97 | YES |
| LQG | 0.182 | 0.173 | 0.95 | YES |
| GUP | 0.399 | 0.605 | 1.52 | YES |
| **Asym. Safety** | **0.041** | **0.112** | **2.77** | **NO** |
| **CDT** | **0.041** | **0.115** | **2.79** | **NO** |

**Table 4.** Comparison of δ_B at the two singularities. Ratio near 1.0
indicates the model treats both walls equivalently.

The two exceptions — Asymptotic Safety and Causal Dynamical Triangulations — are
both dimensional reduction models (spectral dimension 4 → 2). The dimensional
reduction that helps at the isotropic cosmological singularity does not transfer
to the anisotropic radial collapse of a black hole. The geometry matters for
these models; for the other eight, singularities are singularities.

Causal Set Theory is the only model that performs *better* at the black hole
than at the Big Bang (ratio 0.64). The discrete spacetime absorbs a localized
point singularity more cleanly than a spatially uniform cosmological singularity.

---

## 4. The Wormhole Control Experiment

To determine whether the Causal Set response is driven by curvature in general
or by singularities and horizons specifically, we performed a companion
experiment through a Morris-Thorne/Ellis traversable wormhole (detailed in Riner
2026f).

The wormhole geometry:
- Shape function: r(l) = √(l² + b₀²), where l is proper distance and b₀ = 0.1
  (Planck units)
- Throat at l = 0, fully traversable, symmetric about the throat
- **No singularity. No event horizon.**
- Maximum curvature temperature at the throat: T_max = 1/(2πb₀) ≈ 1.59 T_P

This temperature is comparable to the black hole interior at r/r_s ≈ 0.1, where
T_eff = 1.58 T_P. The curvature regimes are similar; the topology is completely
different.

Results:

| Model | Black Hole (inside) | Wormhole (mean) | BH Rank | WH Rank |
|-------|--------------------:|----------------:|--------:|--------:|
| **Causal Set** | **0.011** | **0.056** | **1** | **9** |
| Hagedorn | 0.019 | 0.038 | 2 | 4 |
| CDT | 0.115 | 0.033 | 7 | 1 |
| Asym. Safety | 0.112 | 0.034 | 6 | 2 |
| Standard | 0.108 | 0.055 | 4 | 8 |

**Table 5.** Selected models at the black hole vs. wormhole. Causal Set drops
from 1st to 9th; CDT and Asymptotic Safety rise from bottom-tier to top-tier.

Causal Set Theory drops from 1st place inside the black hole to 9th place at the
wormhole throat. Its δ_B increases fivefold, from 0.011 to 0.056. Meanwhile, CDT
and Asymptotic Safety — which struggled inside the black hole — rise to 1st and
2nd place at the wormhole.

The pattern is unambiguous: Causal Set Theory responds specifically to
singularities and horizons, not to curvature. Remove the singularity and the
horizon, and the discrete spacetime's advantage disappears. The other models'
responses are more curvature-dependent and less topology-dependent.

At the wormhole throat, there is no CS spike, no relaxation toward conformance,
and — we will argue — no energy cost. The smooth geometry does not trigger the
mechanism.

---

## 5. Interpretation: Evaporation Without Radiation

### 5.1 The Mass-Stripping Cycle

The Causal Set journey through the black hole (Table 2) reveals not just a
pattern but a *mechanism*. Read as a physical process:

**Phase 1 — Gravitational stripping** (r = 10 → 2 r_s): The black hole's
gravitational field strips the Causal Set below its natural equilibrium. δ_B
drops from 0.028 to 0.003. The CS is being pushed toward over-conformance —
forced below its resting state. It is losing its mass-like deviation character.

**Phase 2 — Maximum over-conformance** (r ≈ 1 r_s): At the horizon, δ_B = 0.004.
This is *not* the Causal Set's natural state. It is over-stripped — too perfect,
too conformant. The gravitational field has pushed it below equilibrium.

**Phase 3 — Restoration** (r = 0.99 → 0.01 r_s): Inside the horizon, δ_B
climbs from 0.004 back to 0.015–0.017. The black hole is feeding mass-energy
back into the Causal Set to restore it to equilibrium. The CS is being rebuilt.

What is the Causal Set's natural resting state? The data is consistent across
all three geometries: at the Big Bang, CS post-wall mean = 0.017. Near the BH
singularity, CS reaches 0.015. These values are nearly identical. The resting
state is δ_B ≈ 0.015–0.017 — and the 0.003 near the horizon is the anomaly,
not the equilibrium.

**The transaction is one-way.** The black hole spends mass to restore CS from
over-conformance back to equilibrium. But CS at equilibrium does not contribute
mass back. The CS ground state is inert — it is spacetime, the substrate, not
mass. It costs energy to restore, but the restored state generates no energy in
return.

The black hole is not being destroyed. It is being *used* — as an energy source
for spacetime restoration. Mass goes in. Spacetime structure comes out. The
structure does not refund the payment.

### 5.2 The Healing Hypothesis

The mass-stripping cycle points to a general principle: **Causal Set spacetime
actively restores its statistical structure when it encounters a singularity or
horizon, and the energy cost of this restoration is drawn from the black hole's
mass.**

The data across three geometries supports this:

- **At singularities** (BH, Big Bang): δ_B ≈ 0.011–0.017. Near-perfect Benford
  conformance. The discrete spacetime absorbs the singularity.
- **At smooth curvature** (wormhole): δ_B ≈ 0.056. Unremarkable. No special
  response. No restoration needed.
- **Near the event horizon**: δ_B is stripped below equilibrium, then the black
  hole feeds mass back to restore it.

The discrete spacetime — a random Poisson sprinkling of fundamental events —
has a natural statistical structure that conforms to Benford's Law at δ_B ≈
0.017. A singularity or horizon disrupts this structure. The causal set responds
by restructuring itself, consuming energy in the process. The energy source is
the black hole.

### 5.3 Black Holes Shrink Without Radiating

The standard account of Hawking radiation (Hawking 1975) describes particle
creation at the event horizon: virtual pairs form in the vacuum, one member falls
in, the other escapes to infinity as real radiation with a thermal spectrum at
temperature T_H = ℏc³ / (8πGMk_B).

We propose an alternative: **nothing escapes the black hole.** The mass-energy
does not leave as outgoing particles. It is consumed internally by the causal set
restructuring. The black hole shrinks because the geometry is being eaten from
within by discrete spacetime maintaining its own statistical law.

The end result is the same — the black hole loses mass and eventually
evaporates — but the mechanism is fundamentally different:

| | Standard Hawking | Causal Set Healing |
|---|---|---|
| Energy flow | Outward (particles escape) | Inward (geometry consumes) |
| What radiates | Thermal particles at T_H | Nothing — mass simply decreases |
| Horizon role | Site of pair creation | Trigger for restructuring |
| Singularity role | Endpoint (or resolved) | Absorbed by discrete structure |
| Information | Paradox (how is it encoded?) | Consumed with mass |

**Table 6.** Comparison of the two evaporation mechanisms.

### 5.4 Why Hawking Radiation Has Not Been Detected

Hawking radiation from astrophysical black holes has never been observed. The
standard explanation is sensitivity: a solar-mass black hole radiates at T_H ≈
60 nanokelvin, far below any detector threshold. A 10-solar-mass black hole
radiates at 6 nanokelvin. Detection would require either a very small black hole
(none have been observed) or extraordinary sensitivity.

But there is a simpler explanation: **Hawking radiation does not exist as
outgoing particles.** If the evaporation mechanism is structural rather than
radiative — if the black hole shrinks by feeding its mass to the causal set
rather than by emitting photons — then there is no outgoing thermal flux to
detect. The evaporation is real; the radiation is not.

This does not contradict Hawking's calculation at the mathematical level. The
Bogoliubov transformation that predicts particle creation at the horizon computes
a *correlation* between interior and exterior modes. The standard interpretation
maps this correlation to real outgoing particles. The Causal Set interpretation
maps it to internal restructuring — the same mathematical structure, different
physical realization.

### 5.5 Size Independence and the Black Hole Lifecycle

Our simulation uses local geometry. The spectrum at each radius depends on the
local effective temperature, not on the total black hole mass. The Causal Set
process — spike at the horizon, relax inside, restore Benford conformance — looks
identical regardless of the black hole's size.

This has implications for the black hole lifecycle:

**Small black holes** (M ~ M_P): The same local restructuring rate, but the
total mass is tiny. The constant energy drain is a large fraction of the
whole. The black hole loses mass rapidly, the effective temperature rises, the
restructuring intensifies, and the process runs away. This predicts the violent
final evaporation of small black holes — not as an explosion of Hawking radiation,
but as a catastrophic structural collapse where the causal set consumes the
remaining mass in a rapidly accelerating cascade.

**Stellar-mass black holes** (M ~ 10 M_☉): The same local process, but the mass
reservoir is enormous. The energy consumed by causal set restructuring is
negligible compared to the total mass. The black hole is effectively stable on
any observationally relevant timescale — consistent with all current observations.

**Supermassive black holes** (M ~ 10⁹ M_☉): The restructuring is utterly
negligible. The black hole does not measurably evaporate. This is consistent
with the observed stability of supermassive black holes at galactic centers and
explains, from a different direction, why they appear to persist indefinitely.

This qualitatively reproduces the standard scaling relation (evaporation time
∝ M³) without requiring outgoing radiation. The mechanism is the same at all
scales; only the ratio of restructuring cost to total mass changes.

### 5.6 The Information Paradox Dissolves

The black hole information paradox (Hawking 1976) asks: if a black hole
evaporates completely via thermal radiation, what happens to the information that
fell in? Thermal radiation carries no information about the initial state
(it is determined only by the temperature), so the evaporation appears to destroy
information, violating unitarity.

If nothing escapes the black hole, the paradox changes form. There is no outgoing
radiation to demand unitarity of. The information is not lost to thermal emission
— it is consumed, along with the mass, by the causal set restructuring. Whether
this process preserves information is a question about the dynamics of Causal Set
Theory itself, not about semiclassical radiation.

In the Causal Set framework, information is encoded in the causal relations
between events. The restructuring we propose does not destroy causal relations —
it reorganizes them to restore Benford conformance. If the reorganization is
bijective (each initial configuration maps to a unique final configuration), then
the process is unitary and information is preserved. If it is many-to-one, then
information is genuinely lost — but at the fundamental level of spacetime
dynamics, not at the level of semiclassical radiation. The question becomes
tractable within the discrete gravity formalism rather than requiring a full
theory of quantum gravity to resolve.

---

## 6. Black Holes as Spacetime Factories

### 6.1 Mass In, Spacetime Out

The mass-stripping cycle (Section 5.1) describes a one-way transaction: the
black hole spends mass to restore the Causal Set to its equilibrium state, and
the equilibrium state — being spacetime itself — does not return the investment.
Mass is converted to spacetime structure.

If this mechanism is correct, then black holes are not endpoints. They are
**conversion engines**. They take mass (deviation from Benford conformance) and
produce spacetime (Causal Set at equilibrium, δ_B ≈ 0.017). The product of this
conversion is new spacetime geometry.

This connects directly to the framework developed in Riner (2026b), which
proposed that mass represents deviation from the Benford constraint and that
entropy is "mass attempting to return to the massless state — to zero deviation,
to perfect conformance." The CS experiments now identify the mechanism: the
return to conformance does not produce radiation or thermal noise. It produces
spacetime. The endpoint of entropy is not heat death. It is geometry.

### 6.2 Implications for Cosmological Expansion

If black holes produce new spacetime, the universe is expanding at every black
hole. Every galaxy contains a supermassive black hole at its center. The
spacetime produced propagates outward.

This addresses an open question from Riner (2026b, Section 4.5), which proposed
that the accelerating expansion of the universe in low-density voids might
reflect "absence of braking" — regions with no mass to slow the emergence
process. That argument explained why voids *accelerate* but did not explain
where the expansion *originates*. The Causal Set mechanism provides the source:
black holes are the factories. They consume mass, produce spacetime, and the
new spacetime propagates outward. The voids accelerate because there is no mass
to brake the propagation.

This picture inverts the standard relationship between black holes and
expansion. In conventional cosmology, black holes are local objects embedded in
an independently expanding spacetime. In this framework, black holes are
*generators* of that expansion — the sites where mass is converted to geometry
and the universe literally grows.

### 6.3 Gravitational Waves as the Wavefront of New Spacetime

When LIGO detected the binary black hole merger GW150914 (Abbott et al. 2016),
the final black hole was approximately 3 solar masses lighter than the sum of
the two progenitors. In standard physics, this mass deficit was radiated as
gravitational wave energy.

In the Causal Set framework, the interpretation shifts: two black holes
merge → violent CS restructuring at the combined horizon → the excess mass is
converted to spacetime → the newly generated geometry propagates outward as
ripples. Gravitational waves, in this picture, are not energy radiating away
from the source. They are the **wavefront of newly created spacetime**,
propagating at the speed of light because light speed is the propagation rate
of the Benford constraint itself (Riner 2026b, Section 2.2).

Both pictures make the same prediction for wave amplitude and frequency (both
scale with the mass deficit), but they make different predictions for what
happens to the energy content. In the standard picture, gravitational wave
energy is deposited in distant matter through tidal stretching. In the CS
picture, the wave energy *becomes* the space between that matter. The
distinction is subtle but in principle testable: if gravitational waves carry
energy (standard), that energy should be conservable and localizable. If they
*are* spacetime (CS), the energy is not deposited — it is the expansion itself.

### 6.4 Connection to Croker et al. (2023)

Croker, Weiner, and Farrah (2023) presented observational evidence that
supermassive black holes in elliptical galaxies gain mass over cosmic time in a
manner coupled to the cosmological expansion rate — specifically, that black
hole masses scale as M ∝ a^k where a is the scale factor and k ≈ 3. This
"cosmological coupling" suggests that black holes are not isolated objects but
are dynamically linked to the expansion of the universe.

This result is controversial but sits in precisely the same conceptual space as
the mechanism proposed here. If black holes produce spacetime (increasing the
scale factor) and the scale factor feeds back into the effective mass of the
black hole (through cosmological coupling), then the relationship between black
holes and expansion is bidirectional. The CS mechanism provides a physical basis
for this coupling: the black hole converts mass to spacetime, the spacetime
expands, and the expansion modifies the horizon geometry.

Whether the quantitative rates match — whether the CS restructuring rate
produces the observed expansion rate — is an open question requiring detailed
calculation within the Causal Set formalism.

---

## 7. Discussion

### 7.1 What This Paper Does and Does Not Claim

We claim:

1. δ_B is computable inside a black hole for all ten quantum gravity models
   tested. The interior is a valid statistical environment.

2. Causal Set Theory produces the cleanest interior statistical structure
   (δ_B = 0.011), outperforming all competitors by a factor of ~2.

3. The Causal Set response is specific to singularities and horizons (confirmed
   by the wormhole control experiment), not to curvature in general.

4. Eight of ten quantum gravity models produce consistent δ_B at both the
   cosmological and gravitational singularities. Singularities are, from the
   Benford perspective, universal.

We propose but do not prove:

5. The energy cost of Causal Set restructuring at horizons is drawn from the
   black hole's mass.

6. This restructuring, not Hawking radiation, is the mechanism of black hole
   evaporation.

7. The information paradox dissolves because nothing escapes the black hole.

We further speculate:

8. Black holes are spacetime factories — they convert mass into Causal Set
   structure at equilibrium, producing new spacetime geometry.

9. This spacetime production may be the source of cosmological expansion.

10. Gravitational waves from binary mergers may be the wavefront of newly
    generated spacetime rather than radiated energy.

Claims 5–7 are interpretive (Section 5). Claims 8–10 are speculative (Section
6). All are consistent with the data presented here but require independent
theoretical work within the Causal Set formalism to confirm. In particular, one
would need to show that the dynamics of a causal set in the presence of a
horizon necessarily produce an energy cost proportional to the restructuring,
and that this cost maps quantitatively to the Hawking mass loss rate. For the
cosmological claims, one would additionally need to show that the rate of
spacetime production by black holes is consistent with the observed expansion
rate.

### 7.2 Testable Predictions

The two mechanisms — standard Hawking radiation and Causal Set healing — make
the same prediction for mass loss rate but different predictions for observables:

| Observable | Hawking | CS Healing |
|---|---|---|
| Outgoing thermal flux | Yes, at T_H | No |
| Mass loss rate | dM/dt ∝ 1/M² | dM/dt ∝ 1/M² (same) |
| Final evaporation | Burst of radiation | Structural collapse |
| Gravitational wave signature | Possible ringdown | Possible different signature |
| Information in radiation | Encoded (Page curve) | Absent (nothing escapes) |

**Table 7.** Distinguishing predictions.

A sufficiently sensitive experiment near a small, evaporating black hole could in
principle distinguish the two scenarios. If outgoing thermal particles are
detected at T_H, Hawking radiation is confirmed. If the mass decreases without
outgoing thermal flux, the structural mechanism is supported. No such experiment
is currently feasible, but the predictions are distinct.

### 7.3 Relation to Existing Work

The Causal Set approach to black hole thermodynamics has been developed by
Sorkin, Dou, and others (Sorkin 1997, Dou & Sorkin 2003). The entanglement
entropy of a causal set across a horizon has been shown to scale with the horizon
area, reproducing the Bekenstein-Hawking entropy formula. Our results are
consistent with this: the causal set "knows about" the horizon and responds to
it in a specific, quantifiable way.

The idea that black hole evaporation might not involve real particle emission has
been explored in the context of the firewall paradox (Almheiri et al. 2013) and
in approaches where the horizon is a quantum error-correcting code (Pastawski et
al. 2015). Our proposal is distinct: we do not modify the horizon's properties
or invoke complementarity. We propose that the evaporation mechanism is entirely
internal — a structural process of the discrete spacetime rather than a
semiclassical radiation process.

---

## 8. Conclusion

We have swept ten quantum gravity models through a Schwarzschild black hole and
measured the Benford deviation δ_B at every radius from 10 r_s to 0.01 r_s and
beyond. All models survive. The interior is statistically valid. Causal Set
Theory produces the cleanest structure by a wide margin.

The data tells a specific story. As the Causal Set approaches the black hole,
the gravitational field strips it below its natural equilibrium — pushing δ_B
from 0.028 down to 0.003, an over-conformance that is not the Causal Set's
resting state. Inside the horizon, the black hole feeds mass-energy back into
the Causal Set, restoring it from 0.004 to 0.017 — its natural equilibrium,
consistent across all three geometries tested. The transaction is one-way:
the black hole spends mass, the restored spacetime gives nothing back.

The wormhole — comparable curvature, no singularity, no horizon — confirms
that this response is topology-specific, not curvature-driven. Causal Set
Theory drops from 1st place inside the black hole to 9th at the wormhole.
The mechanism requires a singularity or horizon to trigger.

If the restoration has an energy cost, and if that cost is drawn from the
black hole, then the black hole evaporates without radiating. It shrinks
because discrete spacetime is consuming it to maintain its own statistical
nature. The information never leaves. The paradox dissolves. And the product
of the conversion — new spacetime at equilibrium — propagates outward,
potentially contributing to the expansion of the universe itself.

Black holes, in this picture, are not endpoints where matter goes to die.
They are conversion engines where mass becomes geometry. The universe grows
at every black hole, one Planck-scale transaction at a time.

This is, admittedly, a bold chain of claims built on a simple foundation:
first-digit counting, a square root, and a logarithm. But the Benford
framework has proven itself as a measurement instrument — recovering spatial
dimensionality, number theory identities, and the boundary between physical
and unphysical matter. The black hole interior is simply the next place it
can measure. And the measurement says: discrete spacetime wants to obey
Benford's Law, and it will consume a black hole to do it.

---

## References

- Almheiri, A., Marolf, D., Polchinski, J., & Sully, J. (2013). Black holes:
  complementarity vs. firewalls. JHEP 2013, 62.
- Dou, D. & Sorkin, R. D. (2003). Black hole entropy as causal links. Found.
  Phys. 33, 279–296.
- Hawking, S. W. (1975). Particle creation by black holes. Commun. Math. Phys.
  43, 199–220.
- Hawking, S. W. (1976). Breakdown of predictability in gravitational collapse.
  Phys. Rev. D 14, 2460–2473.
- Morris, M. S. & Thorne, K. S. (1988). Wormholes in spacetime and their use
  for interstellar travel. Am. J. Phys. 56, 395–412.
- Pastawski, F., Yoshida, B., Harlow, D., & Preskill, J. (2015). Holographic
  quantum error-correcting codes. JHEP 2015, 149.
- Riner, C. (2026a). Complete monotonicity and Benford's Law: deriving quantum
  statistics from the significant digit distribution.
- Riner, C. (2026b). The Law of Emergence: Benford's distribution as a
  universal constraint on physical reality.
- Riner, C. (2026c). The Benford deviation as a measurement instrument:
  round-trip calibration, fingerprint classification, and an existence filter
  for exotic physics.
- Riner, C. (2026d). [This paper].
- Riner, C. (2026e). Benford's Law at the Planck Wall: ten quantum gravity
  models through the cosmological singularity.
- Riner, C. (2026f). Benford's Law through a wormhole: the Casimir effect and
  the absence of singularity response.
- Abbott, B. P. et al. (LIGO Scientific Collaboration) (2016). Observation of
  gravitational waves from a binary black hole merger. Phys. Rev. Lett. 116,
  061102.
- Croker, K. S., Weiner, J. L., & Farrah, D. (2023). Cosmologically coupled
  compact objects: a single-parameter model for LIGO-Virgo mass and redshift
  distributions. Astrophys. J. Lett. 944, L31.
- Sorkin, R. D. (1997). Forks in the road, on the way to quantum gravity. Int.
  J. Theor. Phys. 36, 2759–2781.

---

## Appendix A: Complete Causal Set Data

[Full 50-point table of CS δ_B and ε(d) at every radial position — to be
generated from black_hole_wall.json]

## Appendix B: All Ten Models — Interior Data

[Full tables for all models inside the horizon — to be generated from
black_hole_wall.json]

## Appendix C: Wormhole Control Data

[Key comparison tables from wormhole_wall.json]
