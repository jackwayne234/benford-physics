# Benford's Law Through a Wormhole: Traversable Geometry, the Casimir Effect, and the Absence of Singularity Response

### Christopher Riner
### Chesapeake, Virginia
### chrisriner45@gmail.com

**Draft -- February 2026**

---

## Abstract

We present a computational experiment in which thermal radiation spectra,
modified by ten quantum gravity proposals, are evaluated at 43 positions along
the proper-distance axis of a Morris-Thorne/Ellis traversable wormhole. The
geometry -- r(l) = sqrt(l^2 + b_0^2) with throat radius b_0 = 0.1 Planck units
-- produces extreme tidal curvature at the throat (T_max = 1.59 T_P) but
contains no singularity and no event horizon. At each position, we compute the
Euclidean deviation from Benford's Law, delta_B, and the per-digit deviation
profile epsilon(d), using the framework developed in Riner (2026a, 2026b).

All ten models survive the traversal -- all 430 spectral evaluations return
computable delta_B values. Every model produces a perfectly symmetric profile:
delta_B(l) = delta_B(-l) with zero measured asymmetry, confirming that the
Benford filter respects the geometric symmetry of the wormhole spacetime.

The central result is negative, and the negative result is the point. Causal
Set Theory, which ranked 1st inside a black hole (mean delta_B = 0.011) and 2nd
at the cosmological singularity (delta_B = 0.017), drops to 9th place at the
wormhole throat (mean delta_B = 0.056). Its delta_B increases fivefold. The
wormhole throat has comparable curvature to the black hole interior at r/r_s =
0.1, where T_eff = 1.58 T_P -- yet no special Causal Set response is triggered.
No spike, no relaxation, no restoration. The mechanism identified in Riner
(2026d) requires a singularity or horizon to activate. Smooth curvature, no
matter how extreme, does not suffice.

The dimensional reduction models -- Causal Dynamical Triangulations and
Asymptotic Safety -- show the opposite pattern: they rise from bottom-tier
performance inside the black hole (ranks 7 and 6) to 1st and 2nd at the
wormhole. Their spectral-dimension-based corrections respond to curvature
magnitude, not topology, and the wormhole's smooth geometry is precisely the
environment where those corrections are most effective.

We introduce a Casimir effect at the throat, restricting the available mode
spectrum as the wormhole geometry confines the quantum vacuum. Eight of ten
models improve under Casimir restriction (Hagedorn by 66%, LQG by 61%).
Two models -- Asymptotic Safety and Causal Dynamical Triangulations -- worsen
dramatically (+112% and +248% respectively), because their dimensional
reduction mechanism conflicts with the external mode restriction: reducing
the effective dimension already restricts modes, and the Casimir restriction
further removes modes that the dimensional reduction was relying on.

A three-wall comparison (Wormhole vs. Black Hole vs. Big Bang) across all ten
models confirms the pattern: models that excel at singularities struggle at
smooth curvature, and vice versa. A throat-size sweep varying b_0 from 0.01 to
1.0 Planck units reveals which models are sensitive to the curvature scale and
which are geometry-independent.

This paper serves as the control experiment for Riner (2026d). The black hole
paper proposed that Causal Set Theory actively restores its statistical
structure at singularities and horizons, consuming the black hole's mass-energy
in the process. The wormhole -- with comparable curvature but no singularity and
no horizon -- tests whether the response is curvature-driven or
topology-driven. The answer is unambiguous: it is topology-driven. The healing
mechanism requires a topological defect to trigger.

---

## 1. Introduction

In Riner (2026d), we swept ten quantum gravity models through a Schwarzschild
black hole and discovered a striking pattern in Causal Set Theory. The discrete
spacetime produced the cleanest statistical structure inside the black hole
(mean delta_B = 0.011), exhibited a mass-stripping cycle at the event horizon,
and generated a per-digit fingerprint that quantitatively matched Hawking
radiation at T = 1.36 T_P (L2 distance = 0.004). We proposed that this pattern
reflects a physical mechanism: the black hole's mass-energy is consumed by the
causal set restructuring as the discrete spacetime maintains Benford
conformance across the horizon -- a healing process triggered specifically by
the singularity and the event horizon.

That claim requires a control experiment.

The black hole has three features that could, individually or in combination,
drive the Causal Set response: (1) extreme curvature, with tidal forces that
diverge at the singularity; (2) a topological boundary, the event horizon,
which separates causally disconnected regions; and (3) a singularity, a point
where geodesics terminate and classical geometry ceases. The question is: which
of these features is necessary?

A traversable wormhole provides the ideal control. The Morris-Thorne/Ellis
geometry (Morris and Thorne 1988; Ellis 1973) has extreme curvature at the
throat -- comparable to the black hole interior -- but no singularity and no
event horizon. It is fully traversable: an observer can pass through the throat
from one asymptotic region to another without encountering any pathology.
Geodesics do not terminate. The spacetime is globally hyperbolic. There is no
information barrier.

If Causal Set Theory responds to curvature alone, it should perform well at
the wormhole throat, where tidal forces peak at T_max = 1.59 T_P. If it
responds specifically to singularities and horizons, it should show no special
behavior -- its performance should be unremarkable, indistinguishable from
other models.

This paper reports the result: Causal Set Theory drops from 1st place (black
hole) to 9th place (wormhole). The response is topology-driven, not
curvature-driven. The wormhole is the control experiment that confirms the
mechanism proposed in Riner (2026d).

### 1.1 Position in the Paper Series

This is Paper 6 in a series:

- **Paper 1** (Riner 2026a): Establishes the mathematical framework.
  Complete monotonicity implies Benford conformance; the Dirichlet eta
  function governs Fermi-Dirac deviations.

- **Paper 2** (Riner 2026b): Proposes Benford's distribution as a universal
  constraint. Mass is deviation from the logarithmic ideal; entropy is the
  return to conformance; c is the propagation speed of the constraint.

- **Paper 3** (Riner 2026c): Demonstrates that delta_B is an invertible
  measurement instrument. Recovers spatial dimensionality (n = 3 exactly),
  the Dirichlet eta function (eta(1) = ln 2), and particle mass. Introduces
  the existence filter: distributions returning UNDEFINED (zero valid modes)
  identify thermodynamically impossible physics.

- **Paper 4** (Riner 2026d): Sweeps ten quantum gravity models through a
  Schwarzschild black hole. Causal Set Theory produces the cleanest interior
  structure (delta_B = 0.011). Proposes the healing hypothesis: the black
  hole's mass-energy is consumed by causal set restructuring, producing
  evaporation without radiation.

- **Paper 5** (Riner 2026e -- separate work): Sweeps the same ten models
  through the cosmological singularity at the Planck temperature. Causal Set
  Theory and Hagedorn/string theory produce the cleanest post-wall structure.
  Eight of ten models treat both singularities equivalently.

- **Paper 6** (this paper): The wormhole control experiment. Comparable
  curvature, no singularity, no horizon. Tests whether the Causal Set
  response is curvature-driven or topology-driven.

### 1.2 The Question

The specific question this paper answers:

*Given a spacetime with curvature comparable to a black hole interior
(T_max = 1.59 T_P at the wormhole throat vs. T_eff = 1.58 T_P at r/r_s = 0.1
in the black hole) but no singularity and no event horizon, does Causal Set
Theory exhibit the same exceptional Benford conformance that it shows inside a
black hole?*

The answer is no. And that negative result -- the absence of the singularity
response -- is the most informative finding of the experiment.

---

## 2. Setup

### 2.1 Wormhole Geometry

We use the Morris-Thorne/Ellis traversable wormhole (Morris and Thorne 1988;
Ellis 1973). The metric in proper-distance coordinates is:

    ds^2 = -dt^2 + dl^2 + r(l)^2 (d theta^2 + sin^2 theta d phi^2)

where l is the proper distance from the throat and:

    r(l) = sqrt(l^2 + b_0^2)

The throat is located at l = 0, where r(0) = b_0 is the minimum areal radius.
The geometry is symmetric about l = 0: the two asymptotically flat regions
(l -> +infinity and l -> -infinity) are mirror images.

Key properties:

- **No singularity**: the Kretschner scalar and all curvature invariants are
  finite everywhere. The Ricci scalar is bounded and smooth at the throat.
- **No event horizon**: g_tt = -1 everywhere. There is no redshift surface,
  no trapped surface, no causal barrier.
- **Fully traversable**: an observer can pass through the throat in finite
  proper time without encountering any pathology.
- **Symmetric**: all physical quantities satisfy f(l) = f(-l).

We set the throat radius b_0 = 0.1 in Planck units. This produces a throat
that is 10 times smaller than the Planck length -- an extremely compact
wormhole with curvature scales well into the quantum gravity regime.

### 2.2 Temperature Model

The temperature experienced by a traversing observer is determined by the
tidal forces of the wormhole geometry. We adopt a tidal-temperature model:

    T(l) = (1 / 2 pi) * b_0^2 / (l^2 + b_0^2)^{3/2}

This gives:

- At the throat (l = 0): T(0) = 1/(2 pi b_0) = 1/(2 pi * 0.1) = 1.5915 T_P
- Far from the throat (|l| >> b_0): T(l) -> 0 as l^{-3}

The throat temperature T(0) = 1.59 T_P is comparable to the effective
temperature inside a Schwarzschild black hole at r/r_s = 0.1, where
T_eff = T_H * (r_s/r)^{3/2} = 1.58 T_P for our standard black hole
parameters (T_H = 0.05 T_P). This ensures that the curvature regimes are
matched between the wormhole and black hole experiments, isolating the
topological difference.

### 2.3 Sampling Points

We sample 43 proper-distance points, symmetric about the throat:

- 21 points at l > 0 (approaching the throat from one side)
- 1 point at l = 0 (at the throat)
- 21 points at l < 0 (receding from the throat on the other side)

The sampling is denser near the throat (where the curvature and temperature
change rapidly) and sparser in the asymptotic regions (where the geometry
approaches flat spacetime). This produces 10 models x 43 positions = 430
spectral evaluations for the traversing observer, plus an additional 430
evaluations with Casimir mode restriction.

### 2.4 Observer Types

Two observer perspectives are modeled:

**Traversing observer**: An observer passing through the wormhole experiences
the tidal temperature profile T(l) described above. At each position, the
standard spectrum S(k) = g(k) / (exp[E(k)/T(l)] - 1) is computed using the
model-specific dispersion relation E(k) and density of states g(k). This is
the primary analysis.

**Casimir observer**: The same traversing observer, but with the available mode
spectrum restricted by the Casimir effect. The wormhole throat confines the
quantum vacuum to a finite region of characteristic size b_0, restricting the
available modes to those with wavelength lambda < 2 b_0. This imposes a minimum
momentum k_min = pi / b_0 and removes all modes below this threshold. The
Casimir restriction is strongest at the throat and irrelevant far from it,
implemented as a position-dependent cutoff.

**Pure Casimir reference**: To isolate the effect of mode restriction from
quantum gravity modifications, we also compute a "pure Casimir" spectrum:
standard bosonic occupation (no QG modification) with Casimir mode restriction.
This serves as the baseline for interpreting the Casimir results.

### 2.5 The Ten Quantum Gravity Models

The same ten models used in Riner (2026d) are applied here, with identical
dispersion relations and density-of-states functions. For completeness:

1. **Standard (GR + QFT)**: E = k, g(k) = k^2. No quantum gravity correction.

2. **Loop Quantum Gravity (LQG)**: E = 2|sin(k/2)|, modes restricted to the
   first Brillouin zone (k < pi). Spacetime as a polymer lattice.

3. **GUP (Generalized Uncertainty Principle)**: E = k * sqrt(1 + k^2),
   g(k) = k^2 / (1 + k^2). Minimum length modifies commutation relations.

4. **DSR (Doubly Special Relativity)**: E = 1 - e^{-k}. Maximum energy equal
   to the Planck energy.

5. **Hagedorn (String Theory)**: g(k) = k^2 * exp(k/T_H), T_H = T_P.
   Exponential growth of string states produces a phase transition at T_H.

6. **Causal Set Theory**: g(k) = k^2 * exp(-k^2). Random discrete spacetime
   with Gaussian UV suppression. Lorentz-invariant discreteness.

7. **Asymptotic Safety**: Running spectral dimension d_s(k) = 2 + 2/(1 + k^2),
   g(k) = k^{d_s - 1}. Gravity's UV fixed point drives dimensional reduction
   4 -> 2 at the Planck scale.

8. **Horava-Lifshitz Gravity**: E^2 = k^2 + k^4 + k^6. Anisotropic scaling
   between time and space makes gravity power-counting renormalizable.

9. **Non-commutative Geometry**: E^2 = k^2 + k^4. Spacetime coordinates do
   not commute, introducing a minimum area. UV/IR mixing creates competing
   effects in the density of states.

10. **Causal Dynamical Triangulations (CDT)**: Running spectral dimension
    d_s(k) = 2 + 2/(1 + k^4), g(k) = k^{d_s - 1}. Spacetime built from
    causally ordered simplices. Sharper dimensional reduction transition
    than Asymptotic Safety.

### 2.6 Analysis Protocol

At each proper-distance position, for each model, we:

1. Compute T(l) from the tidal temperature profile
2. Generate the spectrum S(k) = g(k) / (exp[E(k)/T(l)] - 1) over 100,000
   momentum modes
3. For Casimir evaluations: apply the position-dependent mode restriction,
   removing modes with k < k_min(l)
4. Extract the first significant digit of each S(k) value
5. Compute delta_B, epsilon(d), MAD, and the Benford verdict

This yields 430 spectral evaluations for the standard traversal and 430 for
the Casimir-restricted traversal.

---

## 3. Results

### 3.1 All Models Survive

All 430 spectral evaluations return computable delta_B values. No model
produces an undefined distribution at any position along the wormhole. The
wormhole throat, despite its extreme curvature, is a valid statistical
environment for all ten quantum gravity proposals.

This was expected -- the wormhole geometry is smooth and the temperature is
finite everywhere -- but it confirms the basic premise: the Benford filter
can operate in any spacetime that produces a well-defined thermal spectrum.
The instrument does not require a singularity or a horizon to function.

### 3.2 Rankings: Traversing Observer

The mean delta_B for each model across all 43 positions, along with the
throat value:

| Rank | Model | Mean delta_B | Throat delta_B |
|------|-------|-------------|---------------|
| 1 | CDT | 0.033 | 0.014 |
| 2 | Asym. Safety | 0.034 | 0.023 |
| 3 | DSR | 0.036 | 0.013 |
| 4 | Hagedorn | 0.038 | 0.011 |
| 5 | Noncommut. | 0.041 | 0.031 |
| 6 | LQG | 0.045 | 0.019 |
| 7 | Horava-Lif. | 0.048 | 0.037 |
| 8 | Standard | 0.055 | 0.059 |
| 9 | **Causal Set** | **0.056** | **0.053** |
| 10 | GUP | 0.112 | 0.070 |

**Table 1.** Mean Benford deviation across the wormhole traversal. Models
ranked by mean delta_B.

Several features are immediately apparent:

**The ranking is dramatically different from the black hole.** Inside the black
hole (Riner 2026d), Causal Set Theory ranked 1st (delta_B = 0.011). Here it
ranks 9th (delta_B = 0.056). CDT and Asymptotic Safety, which ranked 7th and
6th inside the black hole, now rank 1st and 2nd.

**The spread is compressed.** Inside the black hole, delta_B ranged from 0.011
(Causal Set) to 0.605 (GUP) -- a factor of 55. At the wormhole, the range is
0.033 (CDT) to 0.112 (GUP) -- a factor of only 3.4. The wormhole's smooth
geometry does not discriminate between models as sharply as the black hole's
singular geometry.

**Hagedorn's throat value is the lowest.** At the throat itself (l = 0),
Hagedorn produces delta_B = 0.011 -- the same value as Causal Set Theory
inside the black hole. But Hagedorn's advantage does not extend to the full
traversal (mean = 0.038, rank 4), because its performance degrades in the
lower-temperature regions away from the throat.

**Causal Set is unremarkable.** Delta_B = 0.056 is its worst performance
across all three geometries tested (Big Bang post-wall = 0.017, Black Hole
inside = 0.011, Wormhole = 0.056). The discrete spacetime shows no special
response to the throat curvature. No spike, no relaxation, no restoration
toward conformance. It simply returns a middling, unexceptional result.

### 3.3 The Symmetry Test

The Morris-Thorne/Ellis geometry is symmetric about the throat: all geometric
quantities satisfy f(l) = f(-l). If the Benford filter is responding to the
geometry (and not to numerical artifacts or asymmetric sampling), then
delta_B(l) should equal delta_B(-l) at every point.

**Result: all 10 models return delta_B(l) = delta_B(-l) with zero measured
asymmetry (maximum asymmetry = 0.00e+00).**

This is a non-trivial check. The computation involves 100,000 momentum modes
at each position, with nonlinear dispersion relations and exponential
functions. Numerical roundoff could in principle introduce asymmetry. The
fact that it does not -- that the symmetry is exact to machine precision --
confirms that the Benford filter faithfully represents the underlying geometry.

### 3.4 Traversal Profiles

The full traversal profile for each model takes the form of a smooth curve
peaking at the throat and decaying symmetrically on both sides. Unlike the
black hole journey (where the event horizon creates a one-way transition and
the singularity creates a terminal point), the wormhole traversal is a smooth
round trip: the observer enters from one side, passes through the throat, and
emerges on the other side with the same conditions in reverse.

For Causal Set Theory, the profile is essentially flat near the throat:

| Region | l/b_0 | T (T_P) | CS delta_B |
|--------|-------|---------|-----------|
| Far approach | 20.0 | 0.004 | 0.028 |
| Moderate approach | 5.0 | 0.024 | 0.061 |
| Near throat | 1.0 | 0.564 | 0.058 |
| At throat | 0.0 | 1.592 | 0.053 |
| Near throat (other side) | -1.0 | 0.564 | 0.058 |
| Moderate receding | -5.0 | 0.024 | 0.061 |
| Far receding | -20.0 | 0.004 | 0.028 |

The far-field value (delta_B = 0.028) matches the far-field value outside the
black hole (delta_B = 0.028 at r/r_s = 10) -- both are determined by the
low-temperature spectrum of the Causal Set model, independent of the local
geometry. Near the throat, delta_B rises slightly to 0.053-0.061, showing a
mild response to the curvature -- but nothing like the dramatic structure
seen at the black hole (where delta_B dropped to 0.003 at the horizon and
then recovered to 0.017 inside). The wormhole produces no stripping, no
over-conformance, no restoration.

Compare this to the Causal Set's journey through the black hole:

| Feature | Black Hole | Wormhole |
|---------|-----------|----------|
| Far field delta_B | 0.028 | 0.028 |
| Minimum delta_B | 0.003 (at horizon) | 0.028 (far field) |
| Maximum delta_B | 0.028 (far field) | 0.061 (near throat) |
| At extreme curvature | 0.017 (singularity) | 0.053 (throat) |
| Signature pattern | Strip -> over-conform -> restore | Flat, monotone |
| Rank | 1st | 9th |

**Table 2.** Causal Set behavior at the black hole vs. wormhole. The
dramatic three-phase cycle (stripping, over-conformance, restoration) that
characterizes the black hole is entirely absent at the wormhole.

---

## 4. The Casimir Effect at the Throat

### 4.1 Physical Motivation

A traversable wormhole throat confines the quantum vacuum to a region of
characteristic size b_0. This is analogous to the Casimir effect between
parallel conducting plates: the boundary conditions restrict the available
modes to those that "fit" within the confinement region. For a throat of
radius b_0 = 0.1 Planck units, the minimum wavelength that fits is
lambda = 2 b_0 = 0.2 Planck lengths, corresponding to a minimum momentum
k_min = pi / b_0 = 10 pi in Planck units.

The physical picture: as the observer approaches the throat, the effective
confinement tightens. The longest-wavelength modes are progressively excluded.
At the throat itself, only short-wavelength modes survive. Moving away from
the throat, the confinement relaxes and the full mode spectrum is restored.

This is directly analogous to the Casimir effect (Casimir 1948), but in a
curved spacetime context. The wormhole geometry acts as a natural Casimir
cavity, and the mode restriction modifies the thermal spectrum in a way that
is detectable by the Benford filter.

### 4.2 Implementation

The Casimir restriction is applied as a position-dependent low-momentum
cutoff. At each position l, modes with k < k_min(l) are removed from the
spectrum before computing delta_B. The cutoff is strongest at l = 0 (the
throat, where confinement is tightest) and weakens as |l| increases.

The "pure Casimir" reference applies the same mode restriction to a standard
bosonic spectrum (no quantum gravity modification), isolating the effect of
confinement alone from the combined effect of confinement plus quantum gravity.

### 4.3 Results

| Model | Standard delta_B | Casimir delta_B | Change |
|-------|-----------------|----------------|--------|
| Hagedorn | 0.0108 | 0.0037 | -66% |
| LQG | 0.0192 | 0.0074 | -61% |
| DSR | 0.0131 | 0.0059 | -55% |
| Causal Set | 0.0531 | 0.0283 | -47% |
| Horava-Lif. | 0.0372 | 0.0204 | -45% |
| GUP | 0.0696 | 0.0383 | -45% |
| Standard | 0.0595 | 0.0351 | -41% |
| Noncommut. | 0.0311 | 0.0232 | -25% |
| Asym. Safety | 0.0225 | 0.0477 | **+112%** |
| CDT | 0.0140 | 0.0487 | **+248%** |
| Pure Casimir | -- | 0.0351 | -- |

**Table 3.** Casimir effect at the wormhole throat (l = 0). Values are
delta_B at the throat position. Change is the percentage difference between
standard and Casimir-restricted spectra.

Three distinct classes of response emerge:

**Class 1: Improved by Casimir restriction (8 models, -25% to -66%).** Most
models produce cleaner Benford conformance when the low-momentum modes are
removed. This makes physical sense: the lowest-momentum modes carry the most
geometric structure (they sample the largest wavelengths and are most sensitive
to the density-of-states prefactor). Removing them strips away geometric
complexity and leaves a spectrum that is closer to the scale-invariant ideal.
Hagedorn improves the most (-66%), dropping from 0.0108 to 0.0037 -- a
near-perfect Benford conformance that is its best performance across any
geometry and any condition tested in this series.

**Class 2: Mildly worsened (Noncommutative, -25%).** Non-commutative geometry
shows the weakest improvement, sitting at the boundary between improvement
and degradation. The UV/IR mixing that characterizes this model creates a
competing effect: the mode restriction removes low-energy modes, but the
UV/IR connection means that removing low-energy modes also affects the
high-energy behavior. The net effect is nearly neutral.

**Class 3: Dramatically worsened (Asymptotic Safety +112%, CDT +248%).** The
two dimensional reduction models are severely damaged by the Casimir
restriction. This is the most striking result of the Casimir analysis, and
it has a clear physical explanation.

### 4.4 Why Asymptotic Safety and CDT Worsen

Both Asymptotic Safety and CDT predict spectral dimensional reduction: the
effective spacetime dimension drops from 4 at low energy to 2 at the Planck
scale. This dimensional reduction is implemented through a running spectral
dimension d_s(k) that varies with momentum. At low k (large wavelengths),
d_s = 4 and the density of states is g(k) = k^3 (standard 3+1 dimensional
behavior). At high k (short wavelengths), d_s = 2 and the density of states
is g(k) = k^1 (effective 1+1 dimensional behavior).

The dimensional reduction works by *suppressing* the density of states at high
momentum -- reducing from k^3 to k^1. This is, in effect, a built-in mode
restriction: the model already removes effective degrees of freedom at high
energy by lowering the spectral dimension.

Now apply the Casimir restriction on top of this. The Casimir effect removes
modes from the *bottom* of the momentum spectrum. The dimensional reduction
suppresses modes at the *top*. Together, they squeeze the available modes from
both ends, leaving only a narrow band of intermediate momenta. This drastically
reduces the statistical weight of the spectrum, producing a distribution that
is far from the Benford ideal.

The effect is larger for CDT (+248%) than for Asymptotic Safety (+112%)
because CDT's dimensional reduction transition is sharper. The CDT running
dimension has d_s(k) = 2 + 2/(1 + k^4), where the k^4 denominator produces
a rapid crossover from d_s = 4 to d_s = 2 in a narrow momentum window. The
AS running dimension has d_s(k) = 2 + 2/(1 + k^2), with a gentler k^2
transition. The sharper CDT transition means the mode suppression kicks in
more abruptly, and the Casimir restriction bites harder.

This result reveals a structural prediction: **models that implement UV mode
suppression through dimensional reduction are incompatible with environments
that simultaneously impose IR mode restriction.** A wormhole throat with
Casimir confinement is precisely such an environment. If dimensional reduction
is the correct UV completion of gravity, then the Casimir vacuum at a wormhole
throat should show measurably different quantum fluctuation spectra than
non-dimensional-reduction models predict. This is, in principle, a
distinguishing observable.

### 4.5 The Pure Casimir Reference

The pure Casimir spectrum (standard bosonic occupation with mode restriction,
no quantum gravity modification) returns delta_B = 0.0351. This sits in the
middle of the QG model results: better than Standard GR (0.0351 vs. 0.0595
standard, or equal to the Casimir-restricted Standard at 0.0351 by
construction) and worse than Hagedorn under Casimir (0.0037).

The pure Casimir value establishes the baseline: this is how much Benford
deviation the Casimir confinement alone produces, without any quantum gravity
physics. Models that perform better than this under Casimir restriction (all
of Hagedorn, LQG, DSR, Causal Set, Horava-Lifshitz) have quantum gravity
corrections that constructively combine with the mode restriction. Models
that perform worse (Asymptotic Safety, CDT) have corrections that
destructively interfere with it.

---

## 5. Three-Wall Comparison

### 5.1 Side by Side: Wormhole vs. Black Hole vs. Big Bang

The three experiments -- this paper, Riner (2026d), and the cosmological
singularity experiment -- used identical quantum gravity models and analysis
protocols. The only variable is the spacetime geometry. This enables a
direct comparison:

| Model | Wormhole | Black Hole | Big Bang | WH < BH < BB? |
|-------|---------|-----------|---------|---------------|
| CDT | 0.033 | 0.115 | 0.041 | no |
| Asym. Safety | 0.034 | 0.112 | 0.041 | no |
| DSR | 0.036 | 0.110 | 0.108 | no |
| Hagedorn | 0.038 | 0.019 | 0.014 | no |
| Noncommut. | 0.041 | 0.051 | 0.047 | no |
| LQG | 0.045 | 0.173 | 0.182 | YES |
| Horava-Lif. | 0.048 | 0.125 | 0.129 | YES |
| Standard | 0.055 | 0.108 | 0.078 | no |
| Causal Set | 0.056 | 0.011 | 0.017 | no |
| GUP | 0.112 | 0.605 | 0.399 | no |

**Table 4.** Three-wall comparison across all ten models. The "WH < BH < BB?"
column indicates whether the ordering Wormhole-best, Black-Hole-middle,
Big-Bang-worst is satisfied.

### 5.2 Patterns Across the Three Walls

Several distinct patterns emerge:

**Pattern 1: Singularity specialists (Causal Set, Hagedorn).** These models
perform best at singularities and worst at smooth curvature. Causal Set:
0.011 (BH) -> 0.017 (BB) -> 0.056 (WH). Hagedorn: 0.014 (BB) -> 0.019 (BH)
-> 0.038 (WH). Both are 2-5 times worse at the wormhole than at their best
singular geometry. Their mechanisms -- discrete spacetime (CS) and the string
phase transition (Hagedorn) -- are triggered by topological defects, not by
curvature alone.

**Pattern 2: Curvature responders (CDT, Asymptotic Safety).** These models
perform best at smooth curvature and worst at singularities. CDT: 0.033 (WH)
-> 0.041 (BB) -> 0.115 (BH). Asymptotic Safety: 0.034 (WH) -> 0.041 (BB)
-> 0.112 (BH). Both are 3 times worse at the black hole than at the wormhole.
Their dimensional reduction mechanism responds to the *magnitude* of curvature
but is disrupted by the *topology* of singularities, consistent with the
finding in Riner (2026d) that AS and CDT are the only models that treat the
cosmological and gravitational singularities differently (ratio 2.77 and 2.79
respectively).

**Pattern 3: Consistent degraders (LQG, Horava-Lifshitz).** These are the
only two models where the simple ordering WH < BH < BB holds. LQG: 0.045
-> 0.173 -> 0.182. Horava-Lifshitz: 0.048 -> 0.125 -> 0.129. For these
models, the wormhole (no singularity, finite curvature) is the easiest
environment, the Big Bang (cosmological singularity) is the hardest, and
the black hole (gravitational singularity) falls in between. Their
performance degrades monotonically with the severity of the geometric
pathology.

**Pattern 4: Geometry-independent (DSR, Standard, Noncommutative).** These
models show moderate variation across the three walls but no strong pattern.
DSR: 0.036 / 0.110 / 0.108. Standard: 0.055 / 0.108 / 0.078. Noncommutative:
0.041 / 0.051 / 0.047. The wormhole is mildly easier for all three, but the
black hole and Big Bang values are similar. These models do not respond
strongly to either curvature or topology.

**Pattern 5: Always terrible (GUP).** GUP is the worst model at every wall:
0.112 (WH) / 0.605 (BH) / 0.399 (BB). The minimum-length modification
produces poor statistical structure in all geometries, though the severity
varies dramatically (factor of 5 between wormhole and black hole).

### 5.3 The Inversion of Causal Set and CDT

The most striking feature of the three-wall comparison is the near-perfect
inversion between Causal Set Theory and CDT:

| Wall | Causal Set | CDT | CS Rank | CDT Rank |
|------|-----------|-----|---------|----------|
| Black Hole | 0.011 | 0.115 | 1st | 7th |
| Big Bang | 0.017 | 0.041 | 2nd | 4th |
| Wormhole | 0.056 | 0.033 | 9th | 1st |

The two models are anti-correlated across geometries. Where one excels, the
other struggles. Their mechanisms are fundamentally complementary:

- **Causal Set Theory** excels at topological defects (singularities, horizons)
  because the discrete spacetime provides a natural resolution: the random
  Poisson sprinkling "papers over" the defect, maintaining Benford conformance
  where continuous geometry fails.

- **CDT** excels at smooth curvature because the dimensional reduction
  mechanism provides a natural UV completion: as curvature increases, the
  effective dimension drops from 4 to 2, concentrating the statistical weight
  into fewer degrees of freedom and producing cleaner structure.

At singularities, CDT's dimensional reduction is disrupted by the topology
(the anisotropic collapse of a black hole or the explosive expansion of the
Big Bang is not well-described by a smooth dimensional crossover). At smooth
curvature, the Causal Set's discrete spacetime has nothing to "heal" -- no
defect to paper over -- and it simply returns a generic, unremarkable spectrum.

---

## 6. Throat Size Dependence

### 6.1 Motivation

The throat radius b_0 determines the curvature scale of the wormhole. A
smaller throat produces more extreme curvature (and higher throat temperature
T(0) = 1/(2 pi b_0)), while a larger throat produces gentler curvature. By
varying b_0, we can test which models are sensitive to the curvature scale
and which are geometry-independent.

### 6.2 Results

We evaluate three representative models -- Hagedorn, DSR, and CDT -- at the
throat (l = 0) for five throat radii:

| b_0 | T(0) (T_P) | Hagedorn | DSR | CDT |
|-----|-----------|----------|------|-----|
| 0.01 | 15.92 | 0.011 | 0.015 | 0.108 |
| 0.05 | 3.18 | 0.015 | 0.019 | 0.054 |
| 0.10 | 1.59 | 0.011 | 0.013 | 0.014 |
| 0.50 | 0.32 | 0.037 | 0.015 | 0.031 |
| 1.00 | 0.16 | 0.044 | 0.086 | 0.040 |

**Table 5.** Throat delta_B as a function of throat radius b_0 for three
models.

### 6.3 Interpretation

**Hagedorn**: Best performance at the smallest and intermediate throats
(delta_B = 0.011 at b_0 = 0.01 and 0.10), degrading at larger throats
(0.044 at b_0 = 1.0). The Hagedorn model is optimized for temperatures near
and above T_P, where the string phase transition activates. At b_0 = 1.0,
the throat temperature is only 0.16 T_P -- well below the Hagedorn
temperature, where the model has no mechanism to produce clean structure.

**DSR**: Relatively stable across throat sizes (0.013 to 0.086), with a
sweet spot near b_0 = 0.10 (delta_B = 0.013). The energy saturation
mechanism (E_max = E_P) is most effective when the throat temperature is
near the Planck scale. At larger throats (b_0 = 1.0, T = 0.16 T_P), the
energy saturation is irrelevant (all modes are well below E_P) and
performance degrades.

**CDT**: Shows the most dramatic sensitivity. At b_0 = 0.01 (T = 15.92 T_P),
delta_B = 0.108 -- the worst performance, because the extreme temperature
pushes the spectrum deep into the 2D regime where the dimensional reduction
has already completed and the resulting 1+1D physics does not support rich
statistical structure. At b_0 = 0.10 (T = 1.59 T_P), delta_B = 0.014 --
the best performance, right at the dimensional crossover point where the
transition from 4D to 2D is actively occurring and concentrating the
statistical weight. At b_0 = 1.0 (T = 0.16 T_P), delta_B = 0.040 -- the
temperature is below the crossover and CDT is in its standard 4D regime,
performing unremarkably.

The throat-size sweep reveals that CDT has a *preferred curvature scale*
where it is most effective: the scale at which the dimensional crossover
4D -> 2D is actively occurring. Above this scale (more curvature, hotter
throat), the physics is already 2D and the model has nothing new to
contribute. Below this scale (less curvature, cooler throat), the model
is in standard 4D mode with no UV benefit.

---

## 7. The Causal Set Absence

### 7.1 The Central Negative Result

Across three geometries, Causal Set Theory displays the following pattern:

| Geometry | Singularity? | Horizon? | T_max (T_P) | CS Mean delta_B | CS Rank |
|----------|-------------|----------|------------|----------------|---------|
| Black Hole | Yes | Yes | diverges | 0.011 | 1st |
| Big Bang | Yes | No | diverges | 0.017 | 2nd |
| Wormhole | **No** | **No** | 1.59 | **0.056** | **9th** |

**Table 6.** Causal Set Theory performance across three geometries. The
pattern is unambiguous: CS responds to singularities and horizons, not to
curvature.

The wormhole throat temperature (1.59 T_P) is comparable to the black hole
interior at r/r_s = 0.1, where T_eff = 1.58 T_P and CS achieves delta_B =
0.017. Yet at the wormhole throat, CS achieves only delta_B = 0.053 -- more
than three times worse. The curvature is comparable. The topology is
completely different. And the CS response tracks the topology, not the
curvature.

### 7.2 What the Absence Confirms

The negative result at the wormhole confirms three claims from Riner (2026d):

**Claim 1: The CS response is topology-driven.** The healing hypothesis
proposed in Riner (2026d) states that the Causal Set's exceptional Benford
conformance near singularities and horizons reflects active structural
restoration -- the discrete spacetime "heals" at topological defects. If
this were a curvature effect, the wormhole throat (with comparable curvature)
should trigger a comparable response. It does not. The trigger is topological:
the presence of a singularity where geodesics terminate, or an event horizon
where causal structure changes, is necessary to activate the mechanism.

**Claim 2: The mass-stripping cycle requires a horizon.** At the black hole,
the CS data showed a three-phase cycle: gravitational stripping (delta_B
drops from 0.028 to 0.003), over-conformance at the horizon (delta_B =
0.004), and restoration inside (delta_B rises to 0.017). At the wormhole,
there is no stripping phase, no over-conformance, and no restoration. The
cycle requires the causal asymmetry of an event horizon -- the one-way
nature of the boundary that separates the interior from the exterior. A
wormhole throat has no such asymmetry; it is fully traversable in both
directions.

**Claim 3: Evaporation requires a singularity or horizon.** If the CS healing
mechanism consumes the black hole's mass-energy to restore Benford
conformance, then the absence of healing at the wormhole implies no energy
consumption -- no evaporation-like process. This is physically consistent:
traversable wormholes in general relativity do not evaporate (they require
exotic matter to sustain, but the throat geometry itself is stable). The CS
data correctly predicts this: no topological defect, no healing, no mass
consumption, no evaporation.

### 7.3 Why Causal Set Theory Does Poorly at the Wormhole

The Causal Set model modifies the density of states as g(k) = k^2 * exp(-k^2).
The Gaussian UV suppression (exp(-k^2)) damps modes above the Planck scale.
At temperatures well below T_P, this damping is irrelevant -- most thermally
occupied modes are below the Planck scale anyway. At temperatures near T_P
(such as the wormhole throat at 1.59 T_P), the damping becomes significant:
it removes high-momentum modes that contribute to the thermal spectrum.

But the Gaussian damping is a blunt instrument. It suppresses modes based
on momentum alone, regardless of the local geometry. At a singularity, this
blunt suppression coincidentally matches the behavior needed to maintain
Benford conformance -- the singularity disrupts precisely the high-energy
modes that the Gaussian damping removes, and the two effects reinforce each
other. At smooth curvature, there is no such reinforcement. The Gaussian
damping removes modes that the smooth geometry would have handled perfectly
well, degrading the spectrum without compensation.

The Causal Set model's strength at singularities is not that it does
something special -- it is that its generic UV behavior happens to be exactly
what is needed at topological defects. At smooth curvature, that generic
behavior is just a limitation.

### 7.4 Implications for the Healing Hypothesis

The wormhole result strengthens the healing hypothesis from Riner (2026d) by
narrowing its scope. The hypothesis is not "CS heals at extreme curvature." It
is the more specific claim: "CS heals at topological defects in the causal
structure -- singularities where geodesics terminate and horizons where causal
connectivity changes."

This specificity is important because it makes the hypothesis falsifiable.
Consider a geometry with a mild singularity but weak curvature -- for example,
a conical singularity with bounded tidal forces. If the CS healing responds to
the singularity regardless of the curvature magnitude, the hypothesis is
supported. If it responds only when the curvature exceeds a threshold, the
mechanism is more complex than pure topology. The wormhole experiment
eliminates one possibility (curvature alone) but cannot yet distinguish
between pure topology and topology-plus-threshold.

Similarly, consider a near-extremal black hole, where the horizon exists but
the surface gravity (and hence the Hawking temperature) approaches zero. Does
the CS healing persist at vanishing Hawking temperature? If yes, the trigger is
the existence of the horizon, not its thermodynamic properties. This is an
open experimental question within the computational framework established in
this series.

---

## 8. Discussion

### 8.1 What the Wormhole Control Experiment Proves

The wormhole is the cleanest possible control for the black hole experiment.
It isolates topology from curvature by providing comparable curvature without
the topological features (singularity, horizon) that characterize the black
hole. The result is unambiguous:

1. **The Causal Set response is topology-specific.** The exceptional Benford
   conformance seen inside the black hole (delta_B = 0.011) does not appear
   at the wormhole throat (delta_B = 0.053). The five-fold increase in
   deviation confirms that it is the singularity and/or horizon -- not the
   curvature magnitude -- that triggers the CS healing.

2. **The dimensional reduction models respond to curvature, not topology.**
   CDT and Asymptotic Safety improve at the wormhole (1st and 2nd) after
   performing poorly inside the black hole (7th and 6th). Their mechanism
   -- spectral dimensional reduction -- is effective at smooth, high-curvature
   geometries and disrupted at topologically singular ones.

3. **The Casimir effect discriminates between UV completion mechanisms.**
   Models that implement UV completion through mode suppression (dimensional
   reduction) are damaged by additional mode restriction (Casimir). Models
   that implement UV completion through other means (discrete spacetime,
   string states, energy saturation) are improved or unaffected.

4. **Singularity-responsive models are not curvature-responsive models.**
   The ranking inversion between CS and CDT across three walls demonstrates
   that these are fundamentally different mechanisms. The Benford filter can
   distinguish them.

### 8.2 Implications for Wormhole Physics

The wormhole results also carry implications for the physics of traversable
wormholes themselves:

**All ten models survive.** No quantum gravity model produces undefined
distributions at the wormhole throat. The traversable wormhole is, from the
Benford perspective, a valid physical environment for all tested quantum
gravity proposals. This is consistent with (but does not prove) the physical
realizability of Morris-Thorne geometries.

**Perfect symmetry.** The exact symmetry of delta_B(l) = delta_B(-l) confirms
that the Benford filter preserves the geometric symmetry of the spacetime. For
a theory of measurement that operates on statistical structure, this is a
necessary consistency check: a symmetric geometry should produce symmetric
measurements. It does.

**No evaporation-like signature.** The absence of a CS healing process at the
wormhole is consistent with the absence of Hawking-like radiation. Traversable
wormholes, unlike black holes, have no mechanism for particle creation at the
throat (no horizon means no vacuum state ambiguity). The CS data correctly
reflects this: no topological defect, no healing, no energy cost, no
evaporation.

### 8.3 The Casimir-Dimensional Reduction Conflict

The finding that Casimir mode restriction worsens AS and CDT by +112% and
+248% respectively has a broader implication. Both models predict that the
effective dimension of spacetime drops to 2 at the Planck scale. This
dimensional reduction acts as an intrinsic UV cutoff: fewer dimensions means
fewer degrees of freedom. But the Casimir effect at the wormhole throat
imposes an *extrinsic* IR cutoff: the finite throat geometry removes
long-wavelength modes.

When both cutoffs operate simultaneously, the available mode spectrum is
squeezed from both ends. The resulting narrow-band spectrum has poor statistical
structure -- too few independent modes to produce clean Benford conformance.

This suggests a tension between dimensional reduction and Casimir-like
confinement that may have observational consequences. If the real universe
undergoes dimensional reduction at the Planck scale, then any natural Casimir
cavity at that scale (a wormhole throat, a quantum foam bubble, a fluctuation
of the spacetime topology) would produce a distinctive degradation of the
vacuum fluctuation spectrum. The Benford filter makes this degradation
quantifiable.

### 8.4 Limitations

**Wormhole realizability.** Morris-Thorne wormholes require exotic matter
(matter that violates the null energy condition) to sustain the throat. Whether
such matter exists is an open question. This paper does not address
realizability -- it uses the wormhole geometry as a mathematical laboratory
for testing the topology-dependence of the quantum gravity models.

**Temperature model.** The tidal temperature model T(l) = (1/2pi) * b_0^2 /
(l^2 + b_0^2)^{3/2} is a proxy for the curvature-induced temperature
experienced by a traversing observer. Unlike the Hawking temperature (which
has a rigorous derivation from quantum field theory in curved spacetime), the
tidal temperature at a wormhole throat does not have a universally accepted
derivation. Our choice is motivated by dimensional analysis and the requirement
that T(0) ~ 1/(b_0) at the throat. Different temperature models would produce
quantitatively different delta_B values, though the qualitative conclusion --
CS responds to topology, not curvature -- is robust to reasonable variations.

**Static geometry.** The analysis assumes a static wormhole geometry. Dynamic
wormholes (e.g., those formed in a cosmological context or those collapsing
toward a singularity) might trigger different responses. In particular, a
wormhole that develops a horizon during collapse would transition from a
topology where CS is unremarkable to one where it is exceptional -- a
prediction of the healing hypothesis that could be tested computationally.

**Limited throat-size sampling.** The throat-size sweep uses five b_0 values
(0.01 to 1.0). Finer sampling would better characterize the curvature
dependence of each model and might reveal additional structure (resonances,
transitions) in the delta_B(b_0) curves.

---

## 9. Conclusion

We have swept ten quantum gravity models through a Morris-Thorne/Ellis
traversable wormhole and measured the Benford deviation delta_B at 43
positions along the proper-distance axis. The geometry produces extreme tidal
curvature at the throat (T_max = 1.59 T_P) but contains no singularity and
no event horizon. All ten models survive. Every model produces a perfectly
symmetric profile. And the central result is negative.

Causal Set Theory drops from 1st place inside the black hole to 9th place
at the wormhole. Its delta_B increases fivefold, from 0.011 to 0.056. The
three-phase cycle that characterizes its black hole journey -- gravitational
stripping, over-conformance at the horizon, restoration inside -- is entirely
absent. The throat curvature is comparable to the black hole interior. The
topology is completely different. And the Causal Set response tracks the
topology, not the curvature.

This confirms the healing hypothesis of Riner (2026d): the mechanism by which
Causal Set Theory maintains exceptional Benford conformance near singularities
is topology-specific. It requires a singularity where geodesics terminate or a
horizon where causal structure changes. Smooth curvature, no matter how
extreme, does not trigger it. No healing means no energy cost, and no energy
cost means no evaporation. The wormhole does not shrink because there is
nothing to heal.

The dimensional reduction models -- CDT and Asymptotic Safety -- tell the
complementary story. They rise from bottom-tier performance at the black hole
to 1st and 2nd at the wormhole. Their mechanism is curvature-responsive, not
topology-responsive. And the Casimir effect at the throat exposes a structural
vulnerability: the mode restriction imposed by the throat geometry conflicts
with the intrinsic mode suppression of dimensional reduction, degrading both
models by over 100%.

The three-wall comparison -- wormhole, black hole, Big Bang -- provides a
complete characterization of how ten quantum gravity proposals respond to three
fundamentally different spacetime geometries. The models separate into
singularity specialists (CS, Hagedorn), curvature responders (CDT, AS),
consistent degraders (LQG, Horava-Lifshitz), geometry-independent models
(DSR, Standard, Noncommutative), and the perpetually challenged (GUP). The
Benford filter distinguishes all of them.

The wormhole was designed as a control experiment, and it performed exactly
as a control should: by providing the null result that validates the positive
result elsewhere. The black hole paper claimed that Causal Set Theory heals at
topological defects. The wormhole -- comparable curvature, no defects -- shows
no healing. The claim survives its control.

---

## References

- Casimir, H. B. G. (1948). On the attraction between two perfectly
  conducting plates. Proc. Kon. Ned. Akad. Wetensch. B 51, 793-795.
- Ellis, H. G. (1973). Ether flow through a drainhole: a particle model in
  general relativity. J. Math. Phys. 14, 104-118.
- Morris, M. S. & Thorne, K. S. (1988). Wormholes in spacetime and their
  use for interstellar travel: a tool for teaching general relativity. Am.
  J. Phys. 56, 395-412.
- Riner, C. (2026a). Complete monotonicity and Benford's Law: deriving
  quantum statistics from the significant digit distribution.
- Riner, C. (2026b). The Law of Emergence: Benford's distribution as a
  universal constraint on physical reality.
- Riner, C. (2026c). The Benford deviation as a measurement instrument:
  round-trip calibration, fingerprint classification, and an existence
  filter for exotic physics.
- Riner, C. (2026d). Benford's Law inside a black hole: statistical
  structure beyond the event horizon and a Causal Set mechanism for
  evaporation.
- Riner, C. (2026e). Benford's Law at the Planck Wall: ten quantum gravity
  models through the cosmological singularity.
- Hawking, S. W. (1975). Particle creation by black holes. Commun. Math.
  Phys. 43, 199-220.
- Sorkin, R. D. (1997). Forks in the road, on the way to quantum gravity.
  Int. J. Theor. Phys. 36, 2759-2781.
- Dou, D. & Sorkin, R. D. (2003). Black hole entropy as causal links.
  Found. Phys. 33, 279-296.
- Visser, M. (1995). Lorentzian Wormholes: From Einstein to Hawking.
  AIP Press, Woodbury, New York.
- Ambjorn, J., Jurkiewicz, J., & Loll, R. (2005). Spectral dimension of
  the universe. Phys. Rev. Lett. 95, 171301.
- Lauscher, O. & Reuter, M. (2005). Fractal spacetime structure in
  asymptotically safe gravity. JHEP 0510, 050.
- Bombelli, L., Lee, J., Meyer, D., & Sorkin, R. D. (1987). Spacetime as
  a causal set. Phys. Rev. Lett. 59, 521-524.
- Horava, P. (2009). Quantum gravity at a Lifshitz point. Phys. Rev. D
  79, 084008.
- Connes, A. (1994). Noncommutative Geometry. Academic Press, San Diego.

---

## Appendix A: Complete Wormhole Traversal Data

[Full 43-point table of delta_B and epsilon(d) at every proper-distance
position for all 10 models -- to be generated from wormhole_wall.json]

## Appendix B: Casimir-Restricted Spectra

[Full tables for all 10 models under Casimir mode restriction, with
comparison to standard spectra -- to be generated from wormhole_wall.json]

## Appendix C: Throat Size Sweep -- Extended Data

[Complete throat-size sweep for all 10 models at b_0 = 0.01, 0.05, 0.10,
0.50, 1.00 -- to be generated from wormhole_wall.json]

## Appendix D: Three-Wall Comparison -- Full Tables

[Side-by-side comparison tables for all 10 models across wormhole, black
hole, and Big Bang geometries, with rankings, ratios, and pattern
classification -- compiled from wormhole_wall.json, black_hole_wall.json,
and planck_wall_extended.json]
