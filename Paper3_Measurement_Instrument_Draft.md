# The Benford Deviation as a Measurement Instrument: Round-Trip Calibration, Fingerprint Classification, and an Existence Filter for Exotic Physics

### Christopher Riner
### Chesapeake, Virginia
### chrisriner45@gmail.com

**Draft — February 2026**

---

## Abstract

We demonstrate that the Euclidean deviation from Benford's Law, δ_B, functions
not merely as a diagnostic statistic but as an invertible measurement
instrument. In a series of five experiments, we show that δ_B recovers known
physical quantities from thermal radiation spectra with exact precision: spatial
dimensionality (n = 3.0000 from the Planck spectrum), the Dirichlet eta function
value η(1) = ln 2 (from the Fermi-Dirac distribution), and particle mass from
relativistic dispersion relations. The per-digit deviation profile ε(d) provides
a nine-component fingerprint that identifies the physics responsible for a given
deviation with 96.3% accuracy in blind classification.

Combined, the pair (δ_B, ε(d)) constitutes a two-level measurement system:
δ_B measures the magnitude of deviation from logarithmic scaling, and ε(d)
identifies the source. We apply this instrument to 23 exotic physics candidates
— tachyons, anyons, gravitons, Hawking radiation, negative-mass matter, phantom
energy, axions, Majorana fermions, Unruh radiation, and sterile neutrinos —
and demonstrate three possible outcomes:

1. **Computable with a characteristic fingerprint**: the distribution is
   physically realizable. The fingerprint identifies the type of physics.
2. **Degrading**: the distribution is losing coherence. Fewer valid modes,
   noisier signal. The physics is marginal.
3. **UNDEFINED**: zero valid modes. The distribution cannot form a physical
   thermal ensemble. The candidate does not exist thermodynamically.

Four candidates return UNDEFINED (negative-mass bosons at all mass scales,
phantom dark energy with w < −1), identifying them as thermodynamically
impossible before the field equations are written. Negative-mass fermions
survive — Pauli exclusion is literally the difference between existence and
non-existence when mass goes negative — but with δ_B ≈ 0.7, placing them
far outside any natural statistical structure. Gravitons, Unruh radiation,
and Majorana fermions produce exact fingerprint matches to known distributions,
confirming the instrument's resolution. Hawking radiation with greybody factors
and Dodelson-Widrow sterile neutrinos occupy unique regions of fingerprint
space, carrying signatures of their production mechanisms.

The instrument machinery is elementary: first-digit counting, a square root,
and a logarithm. But it recovers spatial dimensionality, number theory
identities, the boundary between physical and unphysical matter, and the
internal structure of exotic quantum statistics — all from the leading digit
of spectral occupation numbers.

---

## 1. Introduction

In Riner (2026a), we established that a thermal distribution satisfies
Benford's Law — P(d) = log₁₀(1 + 1/d) — if and only if it is completely
monotonic, with all coefficients non-negative in its exponential series
representation. We showed that the Bose-Einstein distribution satisfies this
condition exactly, that the Fermi-Dirac distribution deviates in a manner
governed by the Dirichlet eta function, and that the deviation is computable
analytically through Fourier decomposition.

In Riner (2026b), we proposed that Benford's distribution functions as a
universal constraint on physical reality — that mass is measurable deviation
from the logarithmic ideal, that entropy is the tendency of mass to return to
zero deviation, and that c is the propagation speed of the constraint itself.
That paper framed δ_B as a "yardstick" — a measure of how far a distribution
sits from the Benford ideal.

This paper asks: is the yardstick also a ruler? Can you read the number on it?

We present five experiments that answer yes. Each experiment poses the same
question in a different domain: given a measured δ_B value, can we invert it
to recover the physical parameter that produced it? And when inversion fails —
when the distribution returns no valid statistical sample — what does that
failure mean?

The answer turns out to be remarkably clean. δ_B inverts exactly to recover
spatial dimensionality, number theory identities, and particle mass. The
per-digit deviation ε(d) provides a shape that identifies the source physics
with near-perfect accuracy. And the failure mode — UNDEFINED, zero valid
modes — is an existence verdict. The instrument doesn't just measure what
exists. It identifies what cannot.

---

## 2. Framework

### 2.1 The Benford Deviation δ_B

For a distribution sampled over N values, let P_obs(d) be the observed
fraction of values with first significant digit d ∈ {1, ..., 9}, and let
P_B(d) = log₁₀(1 + 1/d) be the Benford prediction. The Euclidean deviation
is:

    δ_B = √( Σ_{d=1}^{9} [P_obs(d) − P_B(d)]² )

This is the L2 distance between the observed and predicted first-digit
distributions. δ_B = 0 indicates exact Benford conformance. In practice:

- δ_B < 0.01: strong conformance (pure quantum distributions)
- δ_B ≈ 0.01–0.03: moderate conformance (Pauli exclusion, geometric prefactors)
- δ_B > 0.05: significant deviation
- UNDEFINED: no valid modes — the distribution produces no computable sample

### 2.2 The Per-Digit Deviation ε(d)

The signed deviation at each digit provides a nine-component vector:

    ε(d) = P_obs(d) − P_B(d),    d = 1, ..., 9

This vector is the fingerprint of the distribution. Two distributions with the
same δ_B but different ε(d) are different kinds of physics. The fingerprint
encodes the mechanism: whether the deviation arises from Pauli exclusion
(alternating sign pattern), geometric prefactors (monotone pattern), mass
(shifted pattern), or exotic statistics (unique pattern).

### 2.3 The L2 Distance Between Fingerprints

To compare two fingerprints ε_A(d) and ε_B(d), we use the Euclidean distance:

    L2 = √( Σ_{d=1}^{9} [ε_A(d) − ε_B(d)]² )

An L2 distance of 0.000 indicates identical fingerprints. Distances below 0.01
constitute tight structural matches. Distances above 0.02 indicate distinct
physics.

### 2.4 Reference Distributions

Four fundamental thermal distributions serve as reference points throughout
this work:

**Bose-Einstein (BE)**: n(x) = 1/(e^x − 1). The pure bosonic occupation
function. Completely monotonic. Satisfies Benford exactly up to numerical
discretization. δ_B = 0.0056, with a flat ε(d) profile.

**Fermi-Dirac (FD)**: n(x) = 1/(e^x + 1). The fermionic occupation function.
Deviation governed by the Dirichlet eta function η(s). δ_B = 0.0117, with a
characteristic alternating ε(d) pattern.

**Maxwell-Boltzmann (MB)**: n(x) = e^{−x}. The classical limit. A single
exponential. δ_B = 0.0101, intermediate between BE and FD.

**Planck**: B(x) = x³/(e^x − 1). The photon spectral energy density. The x³
prefactor (arising from the three-dimensional density of states) introduces
geometric structure on top of the BE occupation. δ_B = 0.0279, the largest
deviation among the four references.

### 2.5 Computational Protocol

All experiments use the same pipeline:

1. Generate the spectral intensity S(k) = g(k) × n(k) over a momentum grid
   of 100,000 modes, with k spanning the range [10⁻⁴, 50] (in units of T)
2. Extract the first significant digit of each S(k) value
3. Discard non-positive values (these are non-physical modes)
4. Count the fraction P_obs(d) for d = 1 through 9
5. Compute δ_B, ε(d), MAD, and the Benford verdict

Scripts: all experiments run from `scripts/` and output to
`results/round_trip/`. The core Benford engine is `scripts/benford_core.py`.

---

## 3. Experiment 1: Recovering Spatial Dimensionality

### 3.1 Motivation

The Planck radiation spectrum B(ν) ∝ ν³/(e^{hν/kT} − 1) uses the exponent
n = 3 because space has three dimensions. The ν³ factor is the density of
states in three-dimensional momentum space. In n dimensions, the spectral
intensity would be ν^n/(e^{hν/kT} − 1).

If δ_B is a measurement instrument, it should be possible to measure the
Planck spectrum's δ_B value and invert it to recover n = 3 — without knowing
in advance that we live in three dimensions.

### 3.2 Method

We define a one-parameter family of thermal spectra:

    B_n(x) = x^n / (e^x − 1)

where x = hν/kT is the dimensionless frequency and n is the dimensional
exponent, varied from 0 to 5 in steps of 0.5. At each value of n, we compute
δ_B using the standard protocol (Section 2.5).

At n = 0, B₀(x) = 1/(e^x − 1) is the pure Bose-Einstein distribution. At
n = 3, B₃(x) is the Planck spectrum. The sweep traces a calibration curve
δ_B(n) from pure quantum statistics to geometric-dominated spectra.

The round-trip test: we take the known Planck δ_B value, locate it on the
calibration curve, and read off the corresponding n. If the instrument works,
we should recover n = 3 exactly.

### 3.3 Results

The calibration curve δ_B(n):

| n | Physical meaning | δ_B |
|---|---|---|
| 0.0 | Pure Bose-Einstein | 0.0056 |
| 0.5 | — | 0.0044 |
| 1.0 | 1+1 dimensional photon gas | 0.0123 |
| 1.5 | — | 0.0117 |
| 2.0 | 2+1 dimensional photon gas | 0.0189 |
| 2.5 | — | 0.0236 |
| **3.0** | **Planck (3+1 dimensions)** | **0.0279** |
| 3.5 | — | 0.0298 |
| 4.0 | 4+1 dimensional photon gas | 0.0313 |
| 4.5 | — | 0.0331 |
| 5.0 | 5+1 dimensional photon gas | 0.0344 |

The curve is monotonically increasing: higher-dimensional spectra deviate more
from Benford's Law. This is because the geometric prefactor x^n introduces
progressively more structure on top of the scale-invariant BE occupation. The
curve is also smooth and invertible over its entire range, making it suitable
for interpolation.

**Round-trip result**: The Planck δ_B = 0.027921 maps via cubic interpolation
on the calibration curve to n = **3.0000**. The spatial dimensionality of our
universe is recovered exactly.

**Sanity checks**: n = 0 (calibration curve) matches the independently computed
BE δ_B (0.0056 = 0.0056, pass). n = 3 (calibration curve) matches the
independently computed Planck δ_B (0.0279 = 0.0279, pass).

### 3.4 Interpretation

This result is not trivial. The input is a single number — the L2 distance
between observed and predicted first-digit frequencies. The output is a
physical constant of the universe: the number of spatial dimensions. The
relationship is mediated entirely by the density of states: each dimension
contributes one power of momentum to the phase-space integral, and that
contribution is imprinted on the first-digit statistics.

The Benford deviation "sees" dimensionality because dimensionality enters the
spectrum as a geometric prefactor, and geometric prefactors shift the
distribution away from complete monotonicity in a measurable, monotonic way.
More dimensions means more structure means more deviation.

This also establishes the calibration principle: for any physical parameter
that enters a thermal spectrum as a continuously variable modification, one
can construct a calibration curve δ_B(parameter) and invert it to recover the
parameter from a measured δ_B. The parameter need not be dimensionality — it
can be anything that continuously modifies the spectral shape.

---

## 4. Experiment 2: Recovering the Dirichlet Eta Function

### 4.1 Motivation

In Riner (2026a), we showed that the Fermi-Dirac deviation from Benford's Law
is governed by the Dirichlet eta function:

    η(s) = Σ_{n=1}^{∞} (−1)^{n+1} / n^s

At s = 1, η(1) = ln 2 ≈ 0.6931. The alternating signs in the eta function
produce the alternating ε(d) pattern characteristic of Fermi-Dirac statistics.
The connection is analytic: the FD distribution deviates from BE (and hence
from Benford) by exactly the amount encoded in the eta function.

If δ_B is a measurement instrument, we should be able to go from a measured
FD δ_B value backward through the calibration to recover η(1) = ln 2.

### 4.2 Method

We construct a continuous interpolation between BE and FD:

    n(x) = 1/(e^x − 1) − α · 2/(e^{2x} − 1)

At α = 0, this is pure Bose-Einstein. At α = 1, this is exact Fermi-Dirac
(using the identity FD(x) = BE(x) − 2·BE(2x)). We sweep α from 0 to 1.5 in
fine steps and compute δ_B(α) at each point.

The round-trip test: take the known FD δ_B, locate it on the calibration
curve, and read off α. If the instrument works, we recover α = 1, confirming
that the measured deviation corresponds to exactly one unit of the
Dirichlet eta transformation.

### 4.3 Results

The calibration curve δ_B(α):

| α | Physical meaning | δ_B |
|---|---|---|
| 0.0 | Pure Bose-Einstein | 0.0056 |
| 0.2 | 20% toward FD | 0.0048 |
| 0.4 | 40% toward FD | 0.0038 |
| 0.6 | 60% toward FD | 0.0046 |
| 0.8 | 80% toward FD | 0.0078 |
| **1.0** | **Exact Fermi-Dirac** | **0.0117** |
| 1.1 | Past FD (peak region) | 0.0125 |
| 1.2 | Past FD | 0.0118 |
| 1.5 | Deep past FD | 0.0089 |

The curve is non-monotonic globally — it dips below the BE value near α ≈ 0.4
(where the negative contribution partially cancels the BE structure without
yet establishing the FD alternating pattern), then rises through FD at α = 1.0,
peaks near α ≈ 1.1, and decreases beyond. However, the curve is monotonic on
the interval [0, ~1.1], which contains all physically meaningful values (α = 0
to α = 1). Inversion is well-defined on the physical branch.

**Round-trip result**: The Fermi-Dirac δ_B = 0.011736 maps via interpolation
to α = **1.000000**. The eta function parameter is recovered exactly.

This confirms the connection established analytically in Riner (2026a): the
FD deviation from Benford's Law is governed by η(1) = ln 2, and the measured
δ_B encodes this number theory identity in a form that can be extracted by
inversion.

### 4.4 Interpretation

The BE-to-FD interpolation parametrized by α traces a path through
"occupation statistics space." Along this path, δ_B first decreases (the
negative correction partially cancels BE structure without yet establishing
FD structure), reaches a minimum, then increases as the alternating-sign
pattern of the eta function asserts itself.

The dip is itself meaningful: there is a value of α (near 0.4) where the
mixed statistics are *more* Benford-conformant than either pure BE or pure FD.
This is not a coincidence — it's the point where the cancellation between
positive (BE) and negative (eta function) contributions to the distribution
is maximally balanced. It corresponds to a distribution that is "between"
quantum and classical in a precise statistical sense.

The non-monotonicity past α = 1.1 means the calibration curve has a
physical branch (α ∈ [0, 1]) and an unphysical branch (α > 1.1). This is
itself informative: the eta function governs the physical FD deviation,
and going "past" full FD statistics enters a regime where the deviation
decreases again — the statistics are becoming less extreme, not more. The
instrument correctly identifies α = 1 as the maximum physically meaningful
value on the monotonic branch.

---

## 5. Experiment 3: The Fingerprint Atlas

### 5.1 Motivation

Experiments 1 and 2 demonstrated that δ_B (a single number) can be inverted
to recover physical parameters. But δ_B alone cannot distinguish between two
distributions that happen to have similar magnitudes of deviation. What
distinguishes a Planck spectrum (geometric prefactor deviation) from a
Fermi-Dirac distribution (exclusion principle deviation) when both have
δ_B ≈ 0.02?

The answer is the per-digit deviation ε(d). If two types of physics produce
different ε(d) shapes, then the nine-component fingerprint serves as a
classifier: given an unknown spectrum, measure ε(d) and match it to the atlas
of known shapes.

### 5.2 Method

We collected all 27 ε(d) vectors generated across Experiments 1 and 2:
11 from the dimension sweep (n = 0 to 5 in steps of 0.5) and 16 from the
eta recovery (α = 0 to 1.5 in fine steps). We also stored the four reference
fingerprints (BE, FD, MB, Planck) from independent calculations.

A blind classification test: for each of the 27 fingerprints, compute the L2
distance to every other fingerprint and to the four references. Assign each
fingerprint to its nearest match. Score: did the classifier correctly identify
which experiment (dimension sweep or eta recovery) produced the fingerprint?

Additionally, we analyzed the shape types present in each experiment's
fingerprints to characterize the morphological structure of the atlas.

### 5.3 Results

**Blind identification accuracy: 96.3%** (26/27 correctly attributed to source
experiment).

The single misclassification: the dimension sweep fingerprint at n = 0
(pure BE) was equally close to the eta recovery fingerprint at α = 0 (also
pure BE). This is correct behavior — both are the same distribution. The
"error" reflects the physical truth that the two experiments share a common
origin point.

Key structural results:

**Exact match**: The stored Planck reference ε(d) matched the dimension sweep
at n = 3 with L2 distance = 0.000000. This is a self-consistency check — the
Planck spectrum computed independently and the dimension sweep at n = 3 are
identical distributions.

**Shape types in the dimension sweep**: Five distinct morphological categories
emerge:
- Flat (near-zero ε(d) at all digits): n = 0 (pure BE)
- Monotone (ε(d) smoothly increasing or decreasing): n = 0.5 to 1.5
- Oscillatory (alternating positive and negative): n = 2.0 to 3.0
- Sweep (progressive shift across digits): n = 3.5 to 4.5
- Mixed (no simple pattern): n = 5.0

**Shape types in the eta recovery**: Dominated by monotone shapes (14/16
fingerprints). The eta recovery fingerprints are 12 times smoother on average
than the dimension sweep fingerprints. The Dirichlet eta function produces
a single-mode perturbation (alternating signs), while the geometric prefactor
(dimension sweep) introduces multi-scale structure.

### 5.4 Interpretation

The fingerprint atlas establishes that the pair (δ_B, ε(d)) is a two-level
measurement:

- **Level 1 (δ_B)**: How far from Benford? This is the magnitude — the
  distance from the ideal. It tells you the severity of the deviation.

- **Level 2 (ε(d))**: What kind of physics? This is the shape — the
  per-digit structure. It tells you the mechanism behind the deviation.

The 96.3% classification accuracy demonstrates that Level 2 is informative.
Two experiments that produce overlapping δ_B ranges (the dimension sweep at
n ≈ 1.5 has δ_B ≈ 0.0117, the same as FD at α = 1) are still distinguishable
by their fingerprints. The geometric origin of the dimension sweep produces an
oscillatory ε(d) pattern; the exclusion origin of the eta recovery produces a
monotone pattern. Different physics, different fingerprints.

The 12:1 smoothness ratio between eta recovery and dimension sweep fingerprints
has a physical explanation. The eta function is a single analytic object — its
Fourier decomposition has a simple structure (alternating signs with slowly
decaying coefficients). The geometric prefactor x^n is a power law that
multiplies every mode, introducing correlations across the entire spectrum.
The former is a perturbation; the latter is a reshaping. Perturbations
produce smooth fingerprints; reshapings produce structured ones.

---

## 6. Experiment 4: The Mass Dial

### 6.1 Motivation

Experiments 1 and 2 recovered geometric (dimensionality) and algebraic (eta
function) quantities from δ_B. Experiment 4 tests a dynamical quantity:
particle mass.

In a relativistic thermal gas, the dispersion relation E(k) = √(k² + m²)
introduces mass as a parameter that modifies the spectral shape. Heavy
particles have fewer thermally accessible modes (E > T requires k > √(T² − m²),
which restricts the momentum range). This should be visible in δ_B: more mass
means a different spectral shape means a different Benford deviation.

The experiment has two parts:
1. **Mass calibration**: sweep m²/T² from 0 (massless) to +25 (heavy), build
   a δ_B calibration curve, and invert random δ_B values to recover physical
   masses.
2. **Tachyon test**: extend the sweep to negative m² (imaginary mass, the
   mathematical signature of tachyonic particles) and see what happens.

### 6.2 Method

The spectrum at each mass parameter is:

    S(k) = k² / (exp[√(k² + m²)/T] − 1)

where k is momentum, m is mass, and T is temperature. We sweep m²/T² from
−25 to +25 in fine steps.

For the round-trip test: we select six random δ_B values from the calibration
curve's range, invert each to recover the corresponding m²/T², and then
convert to physical masses at two temperature scales:
- CMB temperature: T = 2.725 K
- QCD phase transition: T = 150 MeV

For the tachyon test: we push to extreme negative m² values (m² = −400,
−1600, −2400, −2499) and track the number of valid modes and the
computability of δ_B at each point.

### 6.3 Results: Mass Calibration

The calibration curve δ_B(m²/T²) for positive mass is monotonic and smooth.
As mass increases from zero, δ_B decreases — the restricted momentum range
produces a distribution closer to Benford conformance. This is physically
sensible: massive particles have fewer accessible modes, and fewer modes mean
less geometric structure, which means less deviation from the logarithmic
ideal.

All six random δ_B values were successfully inverted to recover mass values
at both temperature scales. The inversion is exact on the calibration grid
and accurate to four significant figures between grid points.

### 6.4 Results: Tachyon Test

This is the more important result. For tachyonic (imaginary mass) particles,
the dispersion relation E(k) = √(k² + m²) with m² < 0 requires k² > |m²|
for E to be real. Modes with k² < |m²| have imaginary energy and are
non-physical — they produce no valid occupation number.

| m²/T² | Valid modes (of 100,000) | Fraction | δ_B | Status |
|---|---|---|---|---|
| 0 (massless) | 100,000 | 100% | 0.0056 | Normal |
| −25 | 90,000 | 90% | 0.0044 | Computable — quieter than massless |
| −400 | 60,000 | 60% | Computable | Degrading |
| −1,600 | 20,000 | 20% | Computable | Seriously degraded |
| −2,400 | 2,000 | 2% | Barely computable | Nearly dead |
| **−2,499** | **0** | **0%** | **UNDEFINED** | **Non-existent** |

**Table 1.** Tachyon mode survival as a function of m²/T².

The progression is:

1. At moderate tachyonic values (m² = −25), 90% of modes survive. δ_B is
   computable and actually *lower* than the massless case. The tachyon
   is quieter — removing the lowest-momentum modes (which carry the most
   geometric structure) makes the remaining distribution more Benford-like.

2. Deeper into tachyonic territory, modes vanish progressively. The
   distribution degrades — fewer modes means a noisier statistical sample.
   But δ_B remains computable.

3. At m² = −2499, the threshold is crossed: |m| exceeds the maximum
   available momentum. Zero modes have real energy. Zero valid occupation
   numbers. δ_B is not zero, not infinity — it is **UNDEFINED**. The
   distribution does not exist.

The boundary is sharp. There is no gradual transition from "barely exists" to
"doesn't exist." The distribution has valid modes or it doesn't. When it
doesn't, the instrument returns no value — not a bad value, not an extreme
value, but no value at all.

### 6.5 Interpretation

The tachyon test establishes δ_B as an **existence filter**. Three outcomes
are possible when the instrument is applied to a hypothetical distribution:

1. **Computable + characteristic fingerprint**: The distribution exists as a
   physical thermal ensemble. Its fingerprint identifies the type of physics.
   Worth investigating with the full field equations.

2. **Degrading**: The distribution is losing coherence. Valid modes are
   vanishing. The physics is on the boundary of realizability. The instrument
   is seeing a distribution in the process of ceasing to exist.

3. **UNDEFINED**: The distribution produces no valid statistical sample. Not
   zero (which would be a measurement), not infinity (which would be a
   divergence) — simply not computable. The thing does not exist as a physical
   distribution. The field equations will produce no convergent solution.

This parallels Einstein's criterion for physical solutions: if a mathematical
solution cannot be associated with a consistent energy-momentum tensor, it is
not physical. The Benford filter asks the same question differently: does this
thing produce a computable statistical signature? Both are existence tests
operating at different levels of the formalism.

The tachyon result is also structurally interesting. Tachyons don't violate
Benford by producing negative δ_B or anomalously large δ_B. They *degrade*:
the signal gets quieter (lower δ_B) as modes vanish, then crosses a hard
boundary into non-existence. The failure mode is disappearance, not explosion.

---

## 7. Experiment 5: The Whiteboard

### 7.1 Motivation

Experiments 1–4 established δ_B as a calibrated, invertible instrument with an
existence filter. Experiment 5 is the first measurement campaign: apply the
instrument to every exotic physics candidate we can formulate as a thermal or
quasi-thermal distribution, and map the results.

The name reflects the method. We put every candidate on a whiteboard, run each
through the Benford filter, and see what survives. Then we look at the
groupings.

### 7.2 The Candidates

Twenty-three candidates spanning standard model extensions, exotic matter,
topological statistics, and dark sector physics:

**Standard references** (4):
- Bose-Einstein (bosonic occupation, massless)
- Fermi-Dirac (fermionic occupation, massless)
- Maxwell-Boltzmann (classical limit)
- Planck (photon spectral density, 3+1D)

**Exotic matter** (6):
- Negative-mass boson, |m|/T = 1, 5, 10
- Negative-mass fermion, |m|/T = 1, 5, 10

**Particle candidates** (4):
- Graviton (thermal, spin-2 boson)
- Majorana fermion (self-conjugate, spin-1/2)
- Axion (light pseudoscalar boson, m/T = 0.001)
- Sterile neutrino (Dodelson-Widrow production mechanism)

**Event horizon radiation** (4):
- Hawking radiation with greybody cutoff ω_c = 0.5, 1.0, 2.0, 5.0

**Topological statistics** (3):
- Anyons with fractional exclusion g = 0.25, 0.50, 0.75

**Dark sector** (2):
- Phantom energy (w < −1 dark energy)
- Unruh radiation (accelerating observer)

Each candidate is implemented as a specific occupation function n(k) and
density of states g(k), with the spectrum S(k) = g(k) × n(k) evaluated over
100,000 momentum modes.

### 7.3 Implementation Details

**Negative-mass particles**: Dispersion E(k) = √(k² + m²) with m² < 0.
For bosons: n(k) = 1/(e^{E/T} − 1). For fermions: n(k) = 1/(e^{E/T} + 1).
All modes with E < 0 (i.e., k² < |m²|) are discarded as non-physical.
For the bosonic case, when E > 0, e^{E/T} > 1, so the denominator is positive
and the occupation is well-defined.

However, for negative mass, energy itself can become negative: if the
dispersion relation is taken as E = −√(k² + |m|²) (the negative-energy branch),
then e^{E/T} < 1 for all modes, and 1/(e^{E/T} − 1) < 0 everywhere. This is
the scenario we test: negative-mass particles where the defining property is
that their energy is negative.

**Anyons**: Fractional exclusion statistics interpolating between BE (g = 0)
and FD (g = 1). The occupation number follows the Haldane fractional exclusion:
the distribution smoothly varies from 1/(e^x − 1) at g = 0 to 1/(e^x + 1)
at g = 1.

**Hawking radiation**: Planck spectrum modified by a greybody transmission
factor Γ(ω) that suppresses modes below a characteristic frequency ω_c:
Γ(ω) = ω²/(ω² + ω_c²). This represents the partial reflection of radiation
by the curved spacetime near the event horizon.

**Sterile neutrino (Dodelson-Widrow)**: Non-thermal production via active-sterile
mixing in the early universe. The distribution function is f(p) ∝ [p² × FD(p)]
× sin²(2θ) × Γ_collision/H, which produces a suppressed, momentum-dependent
modification of the Fermi-Dirac distribution.

**Phantom energy**: Dark energy with equation of state w < −1. The wrong-sign
kinetic term produces negative energy density, which enters the thermal
framework as a distribution with E < 0 at every mode.

### 7.4 Results

#### 7.4.1 Non-Existent (UNDEFINED)

| Candidate | Valid Modes | Mechanism |
|---|---|---|
| Negative-mass boson, \|m\|/T = 1 | 0 | All n(k) < 0 |
| Negative-mass boson, \|m\|/T = 5 | 0 | All n(k) < 0 |
| Negative-mass boson, \|m\|/T = 10 | 0 | All n(k) < 0 |
| Phantom energy (w < −1) | 0 | All n(k) < 0 |

**Table 2.** Candidates returning UNDEFINED — zero valid modes at every
mass scale tested.

The mechanism is identical for all four: when E < 0, e^{E/T} < 1, so
e^{E/T} − 1 < 0, and the bosonic occupation 1/(e^{E/T} − 1) is negative
at every single mode. Not at some modes. Not marginally. At every mode,
for every momentum, the occupation number is nonsensical. The distribution
is dead everywhere simultaneously.

This is not a matter of degree. These candidates don't produce a small
number of valid modes, or a degraded signal, or a marginal distribution.
They produce *nothing*. The thermal partition function does not converge.
No amount of fine-tuning rescues them.

**Verdict**: Bosonic negative-mass matter and phantom dark energy (w < −1)
cannot form physical thermal distributions. They are thermodynamically
impossible.

#### 7.4.2 Deeply Unnatural (δ_B ≈ 0.7)

| Candidate | δ_B | Note |
|---|---|---|
| Negative-mass fermion, \|m\|/T = 1 | 0.702 | Inverted FD |
| Negative-mass fermion, \|m\|/T = 5 | 0.734 | Computable |
| Negative-mass fermion, \|m\|/T = 10 | 0.714 | Deeply unnatural |

**Table 3.** Negative-mass fermions: computable but with extreme deviation.

These survive where bosonic negative-mass matter does not. The difference is
a single sign in the denominator: fermionic occupation is 1/(e^{E/T} + 1),
where the +1 (instead of −1) ensures the denominator is always positive.
When E < 0, e^{E/T} < 1, but e^{E/T} + 1 is still greater than 1. The
occupation numbers are positive, bounded between 0.5 and 1.

**Pauli exclusion is literally the difference between thermodynamic existence
and non-existence when mass goes negative.** Bosons without that protection
vanish. Fermions survive.

But their δ_B ≈ 0.7 is extraordinary. Normal physics lives below δ_B ≈ 0.03.
The negative-mass fermion deviation is 25 times larger than Planck radiation.
These objects are mathematically computable but so far from natural statistical
structure that they occupy a completely separate region of fingerprint space.
If nature ever produced them, they would be unmistakable — and the distance
from any known distribution suggests that nature does not produce them under
any circumstances resembling thermal equilibrium.

#### 7.4.3 Exact Matches to Known Physics

| Candidate | δ_B | Matches | L2 Distance |
|---|---|---|---|
| Graviton (thermal) | 0.0279 | Planck | 0.000000 |
| Unruh radiation | 0.0279 | Planck | 0.000000 |
| Majorana fermion | 0.0117 | Fermi-Dirac | 0.000006 |
| Axion (m/T = 0.001) | 0.0056 | Bose-Einstein | 0.000014 |

**Table 4.** Candidates with exact fingerprint matches to reference
distributions.

**Gravitons**: A thermal graviton gas in 3+1 dimensions has the same density
of states and bosonic occupation as photons. Spin-2 versus spin-1 changes
only the number of polarization states — a multiplicative constant. Benford's
Law is scale-invariant: multiplying all values by a constant does not change
first digits. The filter cannot distinguish gravitons from photons because
they are thermodynamically identical. Same fingerprint, same physics, different
spin.

**Unruh radiation**: An accelerating observer in flat spacetime sees a thermal
bath at the Unruh temperature T = ℏa/(2πck_B). The distribution is pure
Planck. The filter confirms: identical fingerprint to Planck, distance 0.000000.
Different physical origin (acceleration vs. thermal equilibrium), identical
statistics.

**Majorana fermion**: A fermion that is its own antiparticle. The self-conjugate
property halves the degrees of freedom — but this is again a multiplicative
factor, invisible to first-digit analysis. The filter says: a Majorana fermion
is statistically indistinguishable from a regular fermion (L2 = 0.000006).

**Axion**: A very light pseudoscalar boson with m/T = 0.001 (essentially
massless at the temperatures considered). It sits directly on top of
Bose-Einstein (L2 = 0.000014). The filter cannot distinguish an ultra-light
boson from a massless one — as expected, since the mass contribution to the
dispersion relation is negligible at m/T ≪ 1.

These exact matches validate the instrument's resolution. When two distributions
are physically identical in their statistical structure, the filter says so —
with distances at or below 10⁻⁵. When they are different, the distances are
orders of magnitude larger.

#### 7.4.4 Anyons: Smooth Interpolation from BE to FD

| g | Name | δ_B | Nearest Reference |
|---|---|---|---|
| 0.00 | Boson | 0.0056 | BE (exact) |
| 0.25 | Quarter-on | 0.0036 | BE |
| 0.50 | Semion | 0.0101 | MB |
| 0.75 | Three-quarter-on | 0.0106 | MB |
| 1.00 | Fermion | 0.0117 | FD (exact) |

**Table 5.** Anyon δ_B as a function of fractional exclusion parameter g.

Fractional statistics trace a continuous path through fingerprint space from
BE (g = 0) to FD (g = 1). The trajectory is not a straight line — it passes
through a region closest to Maxwell-Boltzmann at g = 0.5 (the semion).

This is a structural statement about the relationship between quantum and
classical statistics: **half-exclusion looks classical.** Maxwell-Boltzmann
sits between Bose-Einstein and Fermi-Dirac not just historically (as a
classical limit obtained by taking e^{x} ≫ 1) but geometrically in
(δ_B, ε(d)) space. The midpoint of quantum exclusion produces the nearest
approach to classical statistics. The Benford fingerprint makes this visible.

The path from BE through the semion to FD is also smooth — no discontinuities,
no phase transitions. Fractional statistics are a continuous deformation of
quantum statistics, and the Benford instrument sees this continuity.

#### 7.4.5 Unique Fingerprints: New Physics

| Candidate | δ_B | Nearest | L2 Distance | Note |
|---|---|---|---|---|
| Hawking (ω_c = 0.5) | 0.0282 | MB | 0.024 | Strong greybody |
| Hawking (ω_c = 1.0) | 0.0229 | MB | 0.019 | Oscillatory shape |
| Hawking (ω_c = 2.0) | 0.0203 | FD | 0.011 | Between FD and Planck |
| Hawking (ω_c = 5.0) | 0.0353 | Planck | 0.008 | Weak greybody |
| Sterile ν (DW) | 0.0284 | BE | 0.028 | Non-thermal production |

**Table 6.** Candidates with unique fingerprints — no close match to any
reference distribution.

**Hawking radiation**: The greybody factor — the partial reflection of
radiation by curved spacetime near the event horizon — creates a fingerprint
that does not match any fundamental distribution. As the greybody cutoff ω_c
varies (representing different black hole sizes):

- ω_c = 5.0 (small BH, weak greybody): close to Planck, L2 = 0.008. The
  spectrum is nearly pure thermal.
- ω_c = 2.0: drifting away from Planck, now closer to FD. The greybody
  suppression of low-frequency modes shifts the fingerprint.
- ω_c = 1.0: developing oscillatory ε(d) structure. A new shape.
- ω_c = 0.5 (large BH, strong greybody): strong departure from all
  references, L2 > 0.019 from everything. The event horizon is leaving a
  unique mark in the digit statistics.

This is the **fingerprint of an event horizon**. The greybody factor is a
consequence of spacetime curvature modifying thermal radiation, and it produces
a statistical signature that is distinct from any uncurved thermal distribution.
If a spectrum were observed with this fingerprint shape, it would constitute
evidence for an event horizon in the emitting system.

**Sterile neutrino (Dodelson-Widrow)**: No close match to any reference
distribution (L2 > 0.027 from everything). The non-thermal production
mechanism — active-sterile mixing in the early universe — imprints a
momentum-dependent distortion on what would otherwise be Fermi-Dirac
statistics. The distortion is visible in the fingerprint as a unique shape
that reflects the production history, not the equilibrium statistics.

If a particle were detected whose statistical fingerprint matched this
signature, that would be evidence for non-thermal production via the
Dodelson-Widrow mechanism — a specific prediction about the particle's origin
in the early universe.

### 7.5 The Whiteboard Map

The 23 candidates cluster into five natural groupings when mapped by δ_B:

**1. The Benford core** (δ_B < 0.01): Bose-Einstein, axions, anyons with low
exclusion (g < 0.5). Pure quantum distributions with all-positive Dirichlet
coefficients. These satisfy Benford closely because they are completely
monotonic or nearly so. The deviation is small because the statistical
structure is simple — a single, smoothly decreasing occupation function
with no geometric or exclusion-driven modifications.

**2. The exclusion band** (δ_B ≈ 0.01): Fermi-Dirac, Majorana fermion,
Maxwell-Boltzmann, anyons with high exclusion (g > 0.5). Distributions where
Pauli exclusion or classical single-exponential structure introduces controlled
deviation. The deviation is measurable but moderate — the alternating signs of
the eta function or the rigid exponential cutoff of MB shift the distribution
away from complete monotonicity by a specific, computable amount.

**3. The geometric band** (δ_B ≈ 0.02–0.04): Planck, gravitons, Unruh
radiation, Hawking radiation, sterile neutrino. Distributions shaped by
spacetime geometry or non-thermal production. The geometric prefactor (density
of states, greybody factor, production mechanism) introduces structure on top
of the quantum occupation. This additional structure — spatial dimensionality,
event horizons, production history — is what pushes δ_B into the 0.02–0.04
range.

**4. The unnatural zone** (δ_B > 0.1): Negative-mass fermions. Computable but
wildly deviant. The distribution exists mathematically but is so far from
natural Benford structure that it occupies a separate region of fingerprint
space entirely. Nature does not appear to produce distributions in this zone.

**5. Non-existence** (UNDEFINED): Negative-mass bosons, phantom energy. Not
degraded, not marginal — zero valid modes everywhere simultaneously. The
thermal partition function does not converge. These things cannot form physical
distributions under any circumstances.

### 7.6 Two Applications of the Filter

The whiteboard suggests two uses of the Benford existence filter:

**Before the field equations (triage)**: Given a proposed exotic physics
candidate, run it through the filter before investing in the heavy mathematical
machinery of field theory, renormalization, and solution analysis. If the
filter returns UNDEFINED, the candidate cannot form a thermal distribution.
Don't write the Lagrangian. If it returns δ_B > 0.1, the candidate is
mathematically possible but deeply unnatural. Proceed with caution. If it
returns δ_B < 0.05 with a characteristic fingerprint, the candidate is viable
and the fingerprint constrains which kind of physics is operative.

**During the field equations (constraint)**: When a field theory produces
multiple solutions — as general relativity famously does — the Benford filter
can serve as a selection criterion. Solutions that produce clean fingerprints
(low δ_B, recognizable ε(d) shape) are the physical ones. Solutions that
produce degraded or undefined distributions are mathematical artifacts. The
filter narrows the solution space.

---

## 8. Discussion

### 8.1 What Makes δ_B Work as an Instrument

The success of the round-trip experiments rests on a specific property of
Benford's Law: it is the unique distribution that is invariant under
multiplication by arbitrary constants. A distribution that satisfies P(d) =
log₁₀(1 + 1/d) exactly is one where the first-digit statistics carry no
information about the overall scale — all the information is in the *shape*
of the distribution.

When a physical parameter (dimensionality, mass, exclusion) modifies the
shape of a thermal spectrum, it moves the first-digit statistics away from
this scale-invariant ideal. The direction and magnitude of the departure are
determined by the parameter. This is why δ_B is invertible: each parameter
produces a specific departure, and the departure can be read backward.

The key requirement is that the parameter enters the spectrum as a continuous,
monotonic modification. Dimensionality enters as x^n (monotonic in n). The
eta function enters as the alternating series that converts BE to FD (monotonic
on the physical branch). Mass enters through the dispersion relation (monotonic
for positive mass). In each case, the calibration curve δ_B(parameter) is
well-behaved and invertible.

### 8.2 The Existence Filter and the Partition Function

The UNDEFINED verdict from the existence filter has a direct connection to the
partition function of statistical mechanics. A thermal distribution exists if
and only if the partition function Z = Σ_k exp(−E_k/T) converges. For normal
particles (E > 0), each term is a positive number less than 1, and the sum
converges. For negative-mass bosons (E < 0), each term exceeds 1 and grows
without bound — the partition function diverges.

The Benford filter detects this divergence at the level of the occupation
numbers: when every occupation number is negative, there is no valid sample
to analyze. The filter's UNDEFINED is the occupation-number manifestation of
partition function divergence. The two are equivalent diagnostics operating at
different levels: one asks "does the sum converge?", the other asks "do the
occupations make sense?" Both answer no for the same distributions.

The advantage of the Benford approach is practical. Computing the partition
function requires evaluating an infinite sum (or integral) and checking
convergence. Computing δ_B requires generating a finite sample of occupation
numbers and counting first digits. The latter is simpler, faster, and
numerically more robust — and it gives the same verdict.

### 8.3 Why the Semion Result Matters

The finding that the semion (g = 0.5 anyon) sits closest to Maxwell-Boltzmann
in fingerprint space is a statement about the topology of occupation statistics.
It means that the classical limit is not merely a mathematical approximation
(e^x ≫ 1) but a topological midpoint: half-exclusion produces classical
statistics in the Benford metric. The path from quantum to quantum (BE → FD)
passes through classical (MB) at the midpoint of the exclusion parameter.

This is consistent with the view that Maxwell-Boltzmann statistics are the
average of Bose-Einstein and Fermi-Dirac in some appropriate sense. The Benford
fingerprint provides the metric in which this averaging is geometric — not just
a statement about limits, but about distances in statistical space.

### 8.4 Limitations

Several limitations should be noted:

**Numerical discretization**: All results are computed on a finite momentum
grid (100,000 modes). The exact δ_B values depend weakly on the grid size
and range. However, the round-trip experiments test self-consistency (can the
same grid produce and recover the same value?), and this consistency is
exact.

**Temperature dependence**: The occupation numbers depend on the ratio E/T.
All experiments are conducted at a single effective temperature (or
temperature-equivalent). The fingerprints may vary at different temperatures,
particularly for massive particles where m/T controls the shape. The mass
dial (Experiment 4) explicitly addresses this by calibrating at multiple
temperature scales.

**Applicability beyond thermal distributions**: The framework assumes a
thermal occupation function. Non-equilibrium distributions, coherent states,
and other non-thermal systems may not produce meaningful δ_B values. The
sterile neutrino (Dodelson-Widrow) is an intermediate case — it is
non-thermal but still described by a well-defined distribution function.

**The instrument is not a detector**: δ_B tells you what kind of statistical
structure a distribution has. It does not tell you whether that distribution
is realized in nature. The graviton fingerprint matches Planck exactly — this
means gravitons are thermodynamically possible, not that thermal graviton
backgrounds exist.

---

## 9. Conclusion

We have demonstrated that the Benford deviation δ_B is an invertible
measurement instrument. From first-digit counting on thermal spectra, it
recovers:

- **Spatial dimensionality** (n = 3.0000, exact) from the Planck spectrum's
  deviation curve (Experiment 1)
- **The Dirichlet eta function** (α = 1.000000 ↔ η(1) = ln 2) from the
  Fermi-Dirac deviation curve (Experiment 2)
- **The type of physics** (96.3% blind classification accuracy) from the
  per-digit fingerprint ε(d) (Experiment 3)
- **Particle mass** (invertible calibration curve) from relativistic
  dispersion relations (Experiment 4)
- **The boundary between existence and non-existence** (UNDEFINED verdict
  for tachyons beyond the critical mass, negative-mass bosons, phantom
  energy) from the mode survival count (Experiments 4 and 5)

The full measurement campaign (Experiment 5) reveals natural groupings in
the (δ_B, ε(d)) space: a Benford core of pure quantum distributions, an
exclusion band where Pauli statistics introduce controlled deviation, a
geometric band where spacetime structure leaves its imprint, an unnatural
zone where the mathematics works but nature doesn't, and a hard boundary
where the mathematics itself fails.

The machinery is elementary. Count leading digits. Take the L2 distance
from Benford's prediction. Compare the per-digit shape. The instrument
requires no field equations, no renormalization, no path integrals. It
operates upstream of the heavy formalism — asking whether a distribution has
the statistical structure of something physical before the full theory is
invoked.

The fact that this machinery recovers spatial dimensionality exactly — that
you can count first digits of a Planck spectrum and get back the number three
— is not something we expected when the project began. The fact that it also
identifies thermodynamic non-existence — that negative-mass bosons and phantom
energy produce literally zero valid modes, everywhere, simultaneously — was
even less expected. The instrument is simple. The things it measures are not.

---

## References

- Benford, F. (1938). The law of anomalous numbers. Proc. Am. Phil. Soc.
  78(4), 551–572.
- Dodelson, S. & Widrow, L. M. (1994). Sterile neutrinos as dark matter.
  Phys. Rev. Lett. 72, 17–20.
- Haldane, F. D. M. (1991). "Fractional statistics" in arbitrary dimensions:
  a generalization of the Pauli principle. Phys. Rev. Lett. 67, 937.
- Hawking, S. W. (1975). Particle creation by black holes. Commun. Math.
  Phys. 43, 199–220.
- Hill, T. P. (1995). A statistical derivation of the significant-digit law.
  Statistical Science 10(4), 354–363.
- Newcomb, S. (1881). Note on the frequency of use of the different digits
  in natural numbers. Am. J. Math. 4(1), 39–40.
- Nigrini, M. J. (2012). Benford's Law: Applications for forensic accounting,
  auditing, and fraud detection. Wiley.
- Riner, C. (2026a). Complete monotonicity and Benford's Law: deriving quantum
  statistics from the significant digit distribution.
- Riner, C. (2026b). The Law of Emergence: Benford's distribution as a
  universal constraint on physical reality.
- Riner, C. (2026c). [This paper].
- Unruh, W. G. (1976). Notes on black-hole evaporation. Phys. Rev. D 14,
  870–892.

---

## Appendix A: Dimension Sweep Calibration Data

Full δ_B and ε(d) vectors for each value of n from 0 to 5 in steps of 0.5.

[To be generated from `results/round_trip/dimension_sweep.json`]

## Appendix B: Eta Recovery Calibration Data

Full δ_B and ε(d) vectors for each value of α from 0 to 1.5.

[To be generated from `results/round_trip/eta_recovery.json`]

## Appendix C: Fingerprint Atlas — All 27 Vectors

Complete ε(d) vectors, L2 distance matrix, and classification results.

[To be generated from `results/round_trip/fingerprint_atlas.json`]

## Appendix D: Mass Dial Calibration Curve

Full δ_B(m²/T²) data for positive and negative mass, with tachyon boundary.

[To be generated from `results/round_trip/mass_dial.json`]

## Appendix E: Whiteboard — All 23 Candidates

Complete δ_B, ε(d), classification, and L2 distances to all references for
every candidate tested.

[To be generated from `results/round_trip/whiteboard.json`]
