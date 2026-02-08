# Benford's Law Inside a Black Hole: Statistical Structure Beyond the Event Horizon and a Causal Set Mechanism for Evaporation

### Christopher Riner
### Chesapeake, Virginia
### chrisriner45@gmail.com

**Draft — February 2026**

---

## Abstract

We present a computational experiment in which thermal radiation spectra, modified
by ten quantum gravity proposals, are evaluated at fifty-five radial positions through
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

**Abbreviations:** BH = black hole; CS = Causal Set; GR = general relativity;
LQG = Loop Quantum Gravity; GUP = Generalized Uncertainty Principle; DSR =
Doubly Special Relativity; CDT = Causal Dynamical Triangulations; δ_B =
Euclidean (L2) deviation from Benford's Law; ε(d) = per-digit deviation
profile; r_s = Schwarzschild radius; BE = Bose-Einstein; FD = Fermi-Dirac;
T_P = Planck temperature; T_H = Hawking temperature.

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

Here we ask: what does that framework reveal inside a black hole?

---

## 2. Setup

### 2.1 Black Hole Geometry

We use a Schwarzschild black hole with Hawking temperature T_H = 0.05 T_P in
Planck units (corresponding to a black hole mass M ≈ 10 M_P). The
Schwarzschild radius r_s defines the event horizon.

Radial positions are specified as r/r_s, with r/r_s = 1 at the horizon,
r/r_s > 1 outside, and r/r_s < 1 inside. We sample 55 positions:

- **Outside** (20 points): r/r_s = 10.0, 7.0, 5.0, 3.0, 2.0, 1.5, 1.3, 1.2,
  1.15, 1.1, 1.08, 1.06, 1.04, 1.03, 1.02, 1.015, 1.01, 1.005, 1.002, 1.001

- **Inside** (20 points): r/r_s = 0.99, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5,
  0.4, 0.3, 0.25, 0.2, 0.15, 0.12, 0.1, 0.08, 0.06, 0.04, 0.02, 0.01

- **Post-singularity bounce** (15 points): r/r_s = −0.01 through −1.0,
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

### 2.3 The Quantum Gravity Models

Each model modifies the thermal spectrum through its dispersion relation E(k)
and/or density of states g(k). The occupation number at each mode k is:

    n(k) = 1 / (exp[E(k)/T] − 1)

and the spectral intensity is S(k) = g(k) × n(k). We compute δ_B and ε(d) from
the first significant digits of S(k) sampled over a momentum grid.

We test ten models. The primary focus of this paper is:

**Causal Set Theory**: g(k) = k² × exp(−k²). Spacetime is a random discrete
set of points (Poisson sprinkling). The Gaussian UV suppression damps modes
above the Planck scale. Not a lattice — truly random, preserving Lorentz
invariance at the fundamental level.

The baseline is **Standard (GR + QFT)**: E = k, g(k) = k². No quantum gravity
correction.

The remaining eight models — Loop Quantum Gravity (LQG), Generalized
Uncertainty Principle (GUP), Doubly Special Relativity (DSR), Hagedorn/String
Theory, Asymptotic Safety, Horava-Lifshitz, Non-commutative Geometry, and
Causal Dynamical Triangulations (CDT) — each modify the dispersion relation
and/or density of states to implement various quantum gravity proposals. Full
dispersion relations and model descriptions are given in Appendix B. All ten
models are evaluated at every radial position; Section 3.1 reports comparative
rankings.

### 2.4 Analysis Protocol

At each radial position, for each model, we:

1. Compute T_eff(r) from the infalling observer temperature profile
2. Generate the spectrum S(k) = g(k) / (exp[E(k)/T_eff] − 1) over 100,000
   momentum modes (fewer for LQG due to Brillouin zone restriction)
3. Extract the first significant digit of each S(k) value
4. Compute δ_B, ε(d), MAD, and the Benford verdict (CONFORMS / MARGINAL /
   DEVIATES / UNDEFINED)

This yields 10 models × 55 positions = 550 spectral evaluations.

---

## 3. Results

### 3.1 All Models Survive; Causal Set Ranks First

All 550 spectral evaluations return computable δ_B values — no model produces an
undefined distribution at any radius. The black hole interior is a valid
statistical environment for all ten quantum gravity proposals.

The rankings inside the black hole (mean δ_B for r/r_s < 1):

| Rank | Model | Mean δ_B (inside) | Mean δ_B (outside) |
|------|-------|-------------------|--------------------|
| 1 | **Causal Set** | **0.011** | 0.005 |
| 2 | Hagedorn | 0.019 | 0.006 |
| 3 | Noncommut. | 0.051 | 0.043 |
| 4–10 | Standard, DSR, AS, CDT, Horava, LQG, GUP | 0.108–0.605 | 0.060–0.877 |

**Table 1.** Mean Benford deviation inside and outside the event horizon. Full
model-by-model data in Appendix B.

Causal Set Theory produces the cleanest statistical structure inside the black
hole by a significant margin — nearly half the deviation of the runner-up
(Hagedorn). The discrete spacetime maintains near-perfect Benford conformance
(δ_B = 0.011, well within the "strong conformance" threshold of 0.02) even as
the effective temperature rises toward infinity at the singularity. The
remaining sections focus on the Causal Set data.

### 3.2 The Causal Set Journey

The Causal Set data tells a story best read as a journey from far outside to
deep inside:

| Region | r/r_s | T_eff (T_P) | δ_B | Character |
|--------|-------|-------------|-----|-----------|
| Far outside | 10.0 | 0.0016 | 0.028 | Hawking-like deviation |
| Approaching | 5.0 | 0.0045 | 0.004 | Dropping toward conformance |
| Near horizon (min) | 1.1 | 0.043 | 0.002 | Near-perfect conformance |
| At horizon | 1.001 | 0.050 | 0.004 | Conformance holds |
| Just inside | 0.99 | 0.051 | 0.004 | No change at crossing |
| Mid-interior | 0.5 | 0.141 | 0.005 | Slow rise begins |
| Deep interior | 0.2 | 0.559 | 0.012 | Rising toward equilibrium |
| Near singularity | 0.1 | 1.581 | 0.017 | Approaching CS equilibrium |
| Very near sing. | 0.04 | 6.25 | 0.017 | At CS equilibrium (~0.017) |
| At singularity | 0.01 | 50.0 | 0.015 | Stable at equilibrium |
| Post-singularity | −0.01 | 50.0 | 0.015 | Mirror of approach |

**Table 2.** Causal Set δ_B through the complete black hole journey. The deviation
drops from Hawking radiation levels far outside, reaches near-perfect conformance
at the horizon, then gradually rises to the CS equilibrium (~0.017) near the
singularity. (See Figure 1 for the full visual.)

Several features are notable:

**The horizon is invisible.** The transition from r/r_s = 1.001 (outside) to
0.99 (inside) produces no discontinuity, no spike, no change in character. δ_B
= 0.004 on both sides. The event horizon — the defining feature of a black
hole — does not register in the Causal Set statistical structure. This is
consistent with the equivalence principle: a freely falling observer should
notice nothing special at the horizon.

**Far-field Hawking signature.** At r/r_s = 10, far from the black hole, CS
shows δ_B = 0.028. This is within the range of Hawking radiation fingerprints
measured in the exotic physics survey (Riner 2026c) (δ_B = 0.020–0.035 for greybody cutoffs ω_c = 0.5
to 5.0). Far from the hole, the CS spectrum carries the imprint of the Hawking
thermal bath.

**Approach to conformance.** As the observer falls inward from r = 10 r_s, δ_B
drops from 0.028 to 0.002 — nearly perfect Benford conformance. The distribution
*improves* on approach, reaching its cleanest state near the horizon (minimum
δ_B = 0.002 at r ≈ 1.1 r_s) and maintaining near-identical low values just
inside the horizon and through the inner region to about r ≈ 0.5 r_s.

**Gradual interior rise.** Inside the horizon, δ_B rises slowly from 0.004 to
approximately 0.017 near the singularity. This is a controlled increase — the
distribution never leaves the CONFORMS category (δ_B < 0.02 throughout). Even at
r = 0.01 r_s, where T_eff = 50 T_P, the Causal Set spectrum maintains strong
Benford conformance.

**Post-singularity mirror.** The bounce region (r/r_s < 0) mirrors the approach:
δ_B = 0.015 at r = −0.01, matching the pre-singularity value exactly. The
singularity, if it is resolved by a bounce, is symmetric in its Benford
structure.

### 3.3 Fingerprint Match: Causal Set and Hawking Radiation

The per-digit deviation ε(d) provides a richer comparison than δ_B alone. In
Riner (2026e), we compared the ε(d) fingerprint of Causal Set Theory at various
temperatures with Hawking radiation (greybody factor ω_c = 2.0) from the
exotic physics survey (Riner 2026c).

The result: the Causal Set fingerprint converges to the Hawking radiation
fingerprint — but not at the horizon. Just past it.

| CS Location | L2 distance to Hawking (ω_c = 2.0) |
|-------------|-------------------------------------|
| At wall (T = 1.00 T_P) | 0.025 |
| Past wall (T = 1.06 T_P) | 0.020 |
| **Best match (T = 1.36 T_P)** | **0.004** |
| Further past (T = 1.62 T_P) | 0.020 |

**Table 3.** L2 distance between Causal Set ε(d) and Hawking radiation ε(d).

An L2 distance of 0.004 indicates near-identical fingerprint shape. For
reference, the self-match distance is 0.000, and any distance below 0.01
constitutes a tight structural match. The Causal Set spectrum, when pushed past
the singularity into the regime where T > T_P, develops a per-digit profile that
is statistically indistinguishable from Hawking radiation with a moderate
greybody factor.

### 3.4 Hawking Radiation vs. Causal Set: Side-by-Side

The following table places Hawking radiation (ω_c = 2.0, constant δ_B ≈ 0.028)
beside the Causal Set at each position through the black hole. The visual
pattern is immediate: CS starts at Hawking levels, drops far below, then rises
back toward a Hawking-like fingerprint in the deep interior.

| r/r_s | Zone | Hawking δ_B | CS δ_B | CS relative to Hawking |
|-------|------|-------------|--------|------------------------|
| 10.0 | Far outside | 0.028 | 0.028 | Equal — Hawking imprint |
| 5.0 | Approaching | 0.028 | 0.004 | CS 7× cleaner |
| 2.0 | Near horizon | 0.028 | 0.003 | CS 9× cleaner |
| 1.1 | Just outside | 0.028 | 0.002 | CS 14× cleaner (minimum) |
| 1.001 | At horizon | 0.028 | 0.004 | CS 7× cleaner |
| 0.99 | Just inside | 0.028 | 0.004 | No change at crossing |
| 0.5 | Mid-interior | 0.028 | 0.005 | CS 6× cleaner |
| 0.2 | Deep interior | 0.028 | 0.012 | CS rising toward Hawking |
| 0.1 | Near sing. | 0.028 | 0.017 | CS approaching Hawking |
| 0.04 | Very near sing. | 0.028 | 0.017 | Fingerprint match (L2 = 0.004) |
| 0.01 | At singularity | 0.028 | 0.015 | CS at equilibrium |

**Table 4.** Hawking radiation (constant) vs. Causal Set (evolving) δ_B through
the black hole. The CS never reaches Hawking's level but develops the same
fingerprint shape in the deep interior (see Figure 2).

### 3.5 Comparison with the Cosmological Singularity

Eight of ten models produce comparable δ_B at both the cosmological singularity
(Riner 2026e) and the black hole interior. Causal Set Theory is the only model
that performs *better* at the black hole (0.011) than at the Big Bang (0.017),
suggesting the discrete spacetime absorbs a localized point singularity more
cleanly than a spatially uniform one. Full comparison data is in Riner (2026e).

---

## 4. Interpretation: Evaporation Without Radiation

### 4.1 The Mass-Stripping Cycle

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
climbs from 0.004 back to 0.015–0.017. The black hole is feeding mass back
into the Causal Set to restore it to equilibrium. The CS is being rebuilt.

What is the Causal Set's natural resting state? The data points toward a
consistent answer: at the Big Bang wall (Riner 2026e), CS post-wall mean =
0.017. Near the BH singularity, CS reaches 0.015. The resting state is δ_B ≈
0.015–0.017 — and the 0.002 near the horizon is the anomaly, not the
equilibrium.

**The transaction is one-way.** The black hole spends mass to restore CS from
over-conformance back to equilibrium. But CS at equilibrium does not contribute
mass back. The CS ground state is inert — it is spacetime, the substrate, not
mass. It costs mass to restore, but the restored state generates no mass in
return.

The black hole is not being destroyed. It is being *used* — as a mass source
for spacetime restoration. Mass goes in. Spacetime structure comes out. The
structure does not refund the payment.

### 4.2 The Whirlpool: Mass Goes Down the Drain, Light Does Not

The mass-stripping cycle has a natural physical image: the black hole is a
whirlpool. Everything spirals inward. But the drain — the CS restructuring
mechanism at the singularity — only processes one kind of material.

The key distinction is deviation from the CS equilibrium (~0.017):

- **Mass** corresponds to deviation *above* the CS equilibrium. In Riner
  (2026a, 2026b), we established that mass is the signature of deviation from
  Benford conformance. The greater the deviation above ~0.017, the more
  mass-like the distribution.

- **Light** sits at δ_B ≈ 0.006 (pure Bose-Einstein), which is *below* the CS
  equilibrium. Light is already more conformant than the CS resting state.

The drain processes mass-deviation into geometry. It cannot process light,
because light has nothing to give — its deviation is below the threshold. Light
reaches the center of the whirlpool but cannot go down the drain. It
accumulates, trapped by the geometry that the CS built from the mass it already
consumed.

Two classes of behavior at the singularity:

- **Mass** (deviation > ~0.017): consumed by the CS restructuring, converted to
  spacetime geometry. This is the fuel.
- **Light** (deviation < ~0.017): trapped by the closed geometry but not
  consumed. This accumulates.

The whirlpool image clarifies why E = mc² does not apply straightforwardly
inside a black hole. Einstein's mass-energy equivalence treats mass and energy
as interchangeable. The CS drain does not. It is selective: it consumes
deviation, not energy. Light carries energy but no deviation above equilibrium.
Mass carries both. Inside the event horizon, this distinction matters.

### 4.3 The Healing Hypothesis

The mass-stripping cycle points to a general principle: **Causal Set spacetime
actively restores its statistical structure when it encounters a singularity or
horizon, and the cost of this restoration is drawn from the black hole's mass.**

The data supports this:

- **At the black hole singularity**: δ_B ≈ 0.011–0.017. Near-perfect
  conformance. The discrete spacetime absorbs the singularity.
- **At smooth curvature** (wormhole, Riner 2026f): δ_B ≈ 0.056. No special
  response. No restoration needed.
- **Near the event horizon**: δ_B is stripped below equilibrium, then the black
  hole feeds mass back to restore it.

The discrete spacetime — a random Poisson sprinkling of fundamental events —
has a natural statistical structure that conforms to Benford's Law at δ_B ≈
0.017. A singularity or horizon disrupts this structure. The causal set responds
by restructuring itself, consuming mass in the process. The mass source is the
black hole.

### 4.4 Black Holes Shrink Without Radiating

The standard account of Hawking radiation (Hawking 1975) describes particle
creation at the event horizon: virtual pairs form in the vacuum, one member falls
in, the other escapes to infinity as real radiation with a thermal spectrum at
temperature T_H = ℏc³ / (8πGMk_B).

We propose an alternative: **nothing escapes the black hole.** The mass does not
leave as outgoing particles. It is consumed internally by the causal set
restructuring. The black hole shrinks because the geometry is being eaten from
within by discrete spacetime maintaining its own statistical law.

The end result is the same — the black hole loses mass and eventually
evaporates — but the mechanism is fundamentally different:

| | Standard Hawking | Causal Set Healing |
|---|---|---|
| Mass flow | Outward (particles escape) | Inward (geometry consumes mass) |
| What radiates | Thermal particles at T_H | Nothing — mass simply decreases |
| Horizon role | Site of pair creation | Trigger for restructuring |
| Singularity role | Endpoint (or resolved) | Absorbed by discrete structure |
| Information | Paradox (how is it encoded?) | Consumed with mass |

**Table 5.** Comparison of the two evaporation mechanisms.

### 4.5 Why Hawking Radiation Has Not Been Detected

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

### 4.6 Size Independence and the Black Hole Lifecycle

Our simulation uses local geometry. The spectrum at each radius depends on the
local effective temperature, not on the total black hole mass. The Causal Set
process — spike at the horizon, relax inside, restore Benford conformance — looks
identical regardless of the black hole's size.

This has implications for the black hole lifecycle:

**Small black holes** (M ~ M_P): The same local restructuring rate, but the
total mass is tiny. The constant mass drain is a large fraction of the whole.
The black hole loses mass rapidly, the effective temperature rises, the
restructuring intensifies, and the process runs away.

The whirlpool framework (Section 4.2) predicts what happens at the end. Over
the black hole's lifetime, the drain consumes mass-deviation into geometry, but
light — which cannot be processed — accumulates at the center. When the mass
runs out, the drain has nothing left to consume. The closed geometry was
maintained by CS restructuring of mass. No more mass means no more
restructuring, which means the closed geometry opens. All accumulated light
releases at once. This gives a physical mechanism for Hawking's prediction that
small black holes explode in their final moments (T ∝ 1/M, diverging as M → 0).
Hawking's math predicts the burst; the whirlpool framework explains *why*: the
drain ran out of fuel, the geometry opened, and the trapped light escaped.

**Stellar-mass black holes** (M ~ 10 M_☉): The same local process, but the mass
reservoir is enormous. The mass consumed by causal set restructuring is
negligible compared to the total mass. The black hole is effectively stable on
any observationally relevant timescale — consistent with all current observations.

**Supermassive black holes** (M ~ 10⁹ M_☉): The restructuring is utterly
negligible. The black hole does not measurably evaporate. This is consistent
with the observed stability of supermassive black holes at galactic centers and
explains, from a different direction, why they appear to persist indefinitely.

This qualitatively reproduces the standard scaling relation (evaporation time
∝ M³) without requiring outgoing radiation. The mechanism is the same at all
scales; only the ratio of restructuring cost to total mass changes.

### 4.7 The Information Paradox Dissolves

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

## 5. Black Holes as Spacetime Factories

### 5.1 Mass In, Spacetime Out

The mass-stripping cycle (Section 4.1) describes a one-way transaction: the
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

### 5.2 Implications for Cosmological Expansion

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

### 5.3 Gravitational Waves as the Wavefront of New Spacetime

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

### 5.4 Connection to Croker et al. (2023)

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

## 6. Discussion

### 6.1 What This Paper Does and Does Not Claim

We claim:

1. δ_B is computable inside a black hole for all ten quantum gravity models
   tested. The interior is a valid statistical environment.

2. Causal Set Theory produces the cleanest interior statistical structure
   (δ_B = 0.011), outperforming all competitors by a factor of ~2.

3. The Causal Set response is specific to singularities and horizons, not to
   curvature in general. A companion experiment through a traversable wormhole
   (Riner 2026f) — comparable curvature, no singularity, no horizon — shows CS
   dropping from 1st place (BH, δ_B = 0.011) to 9th place (wormhole, δ_B =
   0.056). The topology, not the curvature, drives the response.

We propose but do not prove:

5. The mass cost of Causal Set restructuring at horizons is drawn from the
   black hole itself.

6. This restructuring, not Hawking radiation, is the mechanism of black hole
   evaporation.

7. The information paradox dissolves because nothing escapes the black hole.

We further speculate:

8. Black holes are spacetime factories — they convert mass into Causal Set
   structure at equilibrium, producing new spacetime geometry.

9. This spacetime production may be the source of cosmological expansion.

10. Gravitational waves from binary mergers may be the wavefront of newly
    generated spacetime rather than radiated energy.

Claims 5–7 are interpretive (Section 4). Claims 8–10 are speculative (Section
5). All are consistent with the data presented here but require independent
theoretical work within the Causal Set formalism to confirm. In particular, one
would need to show that the dynamics of a causal set in the presence of a
horizon necessarily produce a mass cost proportional to the restructuring,
and that this cost maps quantitatively to the Hawking mass loss rate. For the
cosmological claims, one would additionally need to show that the rate of
spacetime production by black holes is consistent with the observed expansion
rate.

### 6.2 Testable Predictions

The two mechanisms — standard Hawking radiation and Causal Set healing — make
the same prediction for mass loss rate but different predictions for observables:

| Observable | Hawking | CS Healing |
|---|---|---|
| Outgoing thermal flux | Yes, at T_H | No |
| Mass loss rate | dM/dt ∝ 1/M² | dM/dt ∝ 1/M² (same) |
| Final evaporation | Burst of radiation | Structural collapse |
| Gravitational wave signature | Possible ringdown | Possible different signature |
| Information in radiation | Encoded (Page curve) | Absent (nothing escapes) |

**Table 6.** Distinguishing predictions.

A sufficiently sensitive experiment near a small, evaporating black hole could in
principle distinguish the two scenarios. If outgoing thermal particles are
detected at T_H, Hawking radiation is confirmed. If the mass decreases without
outgoing thermal flux, the structural mechanism is supported. No such experiment
is currently feasible, but the predictions are distinct.

### 6.3 Relation to Existing Work

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

### 6.4 Why General Relativity Breaks Down at Singularities

Einstein's field equations treat mass and energy identically through the
stress-energy tensor T_μν. A kilogram of matter and a kilogram-equivalent of
light produce the same spacetime curvature. This is the foundation of general
relativity.

The CS mechanism distinguishes them. Mass (deviation above ~0.017) is consumed
by the drain. Light (deviation below ~0.017) is trapped but not consumed. GR
does not know about this selectivity. It models the whirlpool — the curvature,
the infall, the closed geometry — but it cannot model what happens at the center
because it treats everything going down the drain as the same material.

This is a specific, statable reason for the breakdown. The singularity is not
where "quantum effects become important" in a vague sense. It is where the
missing variable — the distinction between mass-deviation and light-deviation —
becomes the dominant physics, and GR's equations, which lack this variable,
cannot describe what is actually happening. The equations blow up not because
the math fails, but because the model is incomplete at the point where the
drain's selectivity matters most.

---

## 7. Conclusion

We have swept ten quantum gravity models through a Schwarzschild black hole and
measured the Benford deviation δ_B at every radius from 10 r_s to 0.01 r_s and
beyond. All models survive. The interior is statistically valid. Causal Set
Theory produces the cleanest structure by a wide margin.

The data tells a specific story. As the Causal Set approaches the black hole,
the gravitational field strips it below its natural equilibrium — pushing δ_B
from 0.028 down to 0.003, an over-conformance that is not the Causal Set's
resting state. Inside the horizon, the black hole feeds mass back into
the Causal Set, restoring it from 0.004 to 0.017 — its natural equilibrium,
consistent with the Big Bang wall (Riner 2026e). The transaction is one-way:
the black hole spends mass, the restored spacetime gives nothing back. A
wormhole control experiment (Riner 2026f) confirms this response is
topology-specific — CS drops from 1st place at the black hole to 9th at the
wormhole, where there is no singularity and no horizon.

If the restoration has a mass cost, and if that cost is drawn from the
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

## Figures

**Figure 1.** δ_B for all ten quantum gravity models through the Schwarzschild
black hole, from r = 10 r_s to r = 0.01 r_s (infalling observer). Causal Set
Theory (blue) maintains the lowest deviation throughout the interior. The
horizon (r/r_s = 1) is marked; note the absence of any discontinuity in the CS
curve. (See 08_black_hole_wall.png; interactive version: 10a_black_hole_interactive.html)

**Figure 2.** Fingerprint comparison: Causal Set ε(d) at T = 1.36 T_P vs.
Hawking radiation ε(d) at ω_c = 2.0. The per-digit deviation profiles are
near-identical (L2 distance = 0.004). The CS spectrum in the deep black hole
interior develops the same statistical fingerprint as Hawking radiation.
(See 10_hawking_vs_causal_set.png)

---

## Appendix A: Complete Causal Set Data

[Full 50-point table of CS δ_B and ε(d) at every radial position — to be
generated from black_hole_wall.json]

## Appendix B: All Ten Models — Interior Data

[Full tables for all models inside the horizon — to be generated from
black_hole_wall.json]

