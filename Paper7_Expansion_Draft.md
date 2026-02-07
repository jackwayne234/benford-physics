# Black Holes as Spacetime Factories: Benford Conformance, Causal Set Dynamics, and the Source of Cosmological Expansion

### Christopher Riner
### Chesapeake, Virginia
### chrisriner45@gmail.com

**Draft --- February 2026**

---

## Abstract

We propose a chain of reasoning that begins with a single axiom --- Benford's
distribution, P(d) = log_10(1 + 1/d) --- and terminates at the accelerating
expansion of the universe.

In preceding papers, we established that mass is measurable deviation from
Benford conformance (Riner 2026b), that entropy is the tendency of mass to
return to zero deviation (Riner 2026b), and that Causal Set Theory provides
the physical substrate implementing this constraint with a natural equilibrium
at delta_B approximately equal to 0.017 (Riner 2026c, 2026d). In Riner (2026d), we
demonstrated computationally that Causal Set spacetime, when pushed through a
black hole, exhibits a mass-stripping cycle: the gravitational field strips the
Causal Set below its equilibrium near the horizon, and the black hole feeds
mass-energy back into the discrete spacetime to restore it inside. The
transaction is one-way --- mass goes in, spacetime structure comes out, and the
structure does not refund the payment.

This paper follows that mechanism to its cosmological conclusion. If black
holes convert mass (Benford deviation) into spacetime (Causal Set at
equilibrium), then black holes are spacetime factories --- sites where the
universe literally grows. Every galaxy harbors a supermassive black hole. The
spacetime produced propagates outward. In voids, where no mass brakes the
propagation, the expansion accelerates. The accelerating expansion of the
universe does not require a cosmological constant, phantom energy, or any new
field. It requires black holes producing spacetime and empty space where
nothing slows the result.

We reinterpret gravitational waves from binary black hole mergers (e.g.,
GW150914, where approximately 3 solar masses of mass-energy were unaccounted for in the
final black hole) not as radiated energy but as the wavefront of newly
generated spacetime propagating at c. We connect this framework to the
observational work of Croker et al. (2023), who reported cosmological coupling
of black hole masses to the scale factor (M proportional to a^k, k approximately equal to 3), and
propose that the Causal Set mechanism provides a physical basis for this
coupling. We present testable predictions that distinguish this framework from
standard cosmology and catalog nine open questions requiring quantitative work
within the Causal Set formalism.

The argument is speculative throughout. We are explicit about this. But the
chain is internally consistent, grounded in computational results from the
preceding papers, and --- at its core --- astonishingly simple: one axiom, one
formula, mass goes in, spacetime comes out, the universe expands, and the
constraint builds itself.

---

## 1. Introduction

This paper is the speculative capstone of a series. In Papers 1 through 6,
we built the Benford deviation framework from mathematical foundations to
physical application:

- **Paper 1** (Riner 2026a): Established the mathematical connection between
  Benford's Law and quantum statistics. A thermal distribution satisfies
  P(d) = log_10(1 + 1/d) if and only if it is completely monotonic. The
  Bose-Einstein distribution satisfies this exactly. The Fermi-Dirac
  deviation is governed by the Dirichlet eta function.

- **Paper 2** (Riner 2026b): Proposed Benford's distribution as a universal
  constraint on physical reality. Mass is deviation from the Benford ideal.
  Entropy is the tendency of mass to return to zero deviation. The speed of
  light is the propagation speed of the constraint itself. Section 2.3
  defined mass as measurable delta_B. Section 2.6 proposed that "the end state
  of entropy is not disorder --- it is the universe approaching conformance
  with the axiom as completely as mass allows." Section 4.5 proposed that
  dark energy might be "absence of braking" --- expansion accelerates in
  voids because there is no mass to slow the emergence process --- but left
  open the question of where the expansion originates.

- **Paper 3** (Riner 2026c): Demonstrated that delta_B functions as an invertible
  measurement instrument, recovering spatial dimensionality (n = 3.0000
  exactly), the Dirichlet eta function, and particle mass from first-digit
  statistics. Established the existence filter: distributions that produce
  zero valid modes return UNDEFINED, identifying thermodynamically impossible
  physics (negative-mass bosons, phantom dark energy) before the field
  equations are written.

- **Paper 4** (Riner 2026d): Swept ten quantum gravity models through a
  Schwarzschild black hole. Causal Set Theory produced the cleanest interior
  statistical structure (mean delta_B = 0.011). The mass-stripping cycle
  emerged: gravitational stripping near the horizon, restoration inside via
  mass consumption. The wormhole control experiment confirmed the response is
  topology-specific (singularities and horizons), not curvature-driven.
  Sections 6.1--6.4 proposed, speculatively, that black holes are spacetime
  factories and that gravitational waves are the wavefront of new spacetime.

- **Papers 5 and 6** (Riner 2026e, 2026f): Extended the framework through the
  cosmological singularity and wormhole geometries, confirming Causal Set
  Theory's equilibrium behavior across three distinct walls.

Paper 4, Section 6 opened a door. It proposed that if black holes convert mass
to spacetime, then every black hole is a site where the universe grows. It
proposed that gravitational waves from mergers are the expansion front of newly
created geometry. It proposed that the accelerating expansion might originate
at black holes and accelerate in voids. But those proposals occupied four pages
at the end of a data paper. They were labeled as speculation and left
undeveloped.

This paper develops them. The central question is simple: **if black holes
convert mass to spacetime, where does the spacetime go?**

The answer we propose: outward. Into the universe. As expansion.

We are explicit about the epistemic status of every claim in this paper. The
computational results from Papers 1--6 are established. The mass-stripping
cycle in Causal Set Theory is a demonstrated pattern in the data. Everything
beyond that --- the physical mechanism of mass conversion, the identification
of black holes as expansion sources, the reinterpretation of gravitational
waves --- is proposal. We use "we propose" rather than "we have shown." The
chain of reasoning is internally consistent, but each link requires
independent confirmation within the Causal Set formalism. This paper presents
the chain. The confirmation is future work.

---

## 2. The Full Argument Chain

Before developing each link in detail, we state the complete chain. This is
the argument from axiom to expansion, condensed:

1. **The axiom**: P(d) = log_10(1 + 1/d) --- Benford's distribution as the
   foundational constraint on physical distributions.

2. **Mass = deviation from the axiom**: A distribution's delta_B measures how
   far it sits from Benford conformance. Mass is this deviation made physical
   (Paper 2, Section 2.3).

3. **Entropy = return to conformance**: Mass tends toward zero deviation. The
   second law of thermodynamics, reframed, says that deviation decreases over
   time (Paper 2, Section 2.6).

4. **Causal Set Theory = the physical substrate that implements the axiom**:
   Discrete spacetime --- a random Poisson sprinkling of fundamental events ---
   has a natural equilibrium at delta_B approximately equal to 0.017. This equilibrium is
   stable across cosmological singularities, black hole interiors, and extreme
   curvature regimes (Papers 4--6).

5. **Black holes = the constraint's enforcement mechanism**: Black holes
   collect mass (deviation) and convert it to spacetime (Causal Set at
   equilibrium). The mass-stripping cycle (Paper 4, Section 5.1) is the
   mechanism: gravitational stripping near the horizon, restoration inside via
   mass consumption. The transaction is one-way.

6. **Evaporation without radiation**: The black hole shrinks not because
   particles escape (Hawking radiation) but because the Causal Set consumes
   its mass to maintain Benford conformance. Nothing leaves the horizon
   (Paper 4, Section 5.3).

7. **The expansion of the universe**: The product of mass-to-spacetime
   conversion --- Causal Set structure at equilibrium --- is new geometry. It
   propagates outward from black holes. In voids, where no mass brakes the
   propagation, it accelerates. This is the source of cosmological expansion.

8. **Gravitational waves = the wavefront of newly generated spacetime**: During
   violent binary mergers, the mass deficit (e.g., 3 solar masses in GW150914)
   is converted to spacetime. The newly created geometry propagates outward at
   c as gravitational waves.

Each link builds on the previous. Remove any link and the chain breaks. The
remainder of this paper develops links 5--8 in full, connecting them to
observational data, existing theoretical work, and testable predictions.

---

## 3. The Mechanism: Mass-Stripping and One-Way Conversion

### 3.1 Recap of the Mass-Stripping Cycle

Paper 4 (Riner 2026d, Section 5.1) presented the Causal Set journey through a
Schwarzschild black hole. The data, reproduced here for reference:

| Region | r/r_s | T_eff (T_P) | delta_B | Verdict |
|--------|-------|-------------|---------|---------|
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
| Post-singularity | -0.01 | 50.0 | 0.015 | CONFORMS |

**Table 1.** Causal Set delta_B through a Schwarzschild black hole (reproduced
from Riner 2026d, Table 2).

The cycle has three phases:

**Phase 1 --- Gravitational stripping** (r = 10 to 2 r_s): The black hole's
gravitational field strips the Causal Set below its natural equilibrium. delta_B
drops from 0.028 to 0.003. The discrete spacetime is pushed toward
over-conformance --- forced below its resting state. In the language of Paper
2, it is losing its mass-like deviation character. The gravitational field is
removing the deviation.

**Phase 2 --- Maximum over-conformance** (r approximately equal to 1 r_s): At the horizon, delta_B =
0.004. This is not the Causal Set's natural state. The CS equilibrium, measured
consistently across three geometries, is delta_B approximately equal to 0.015--0.017. The value of
0.003--0.004 near the horizon represents over-stripping --- the gravitational
field has pushed the spacetime below its ground state.

**Phase 3 --- Restoration** (r = 0.99 to 0.01 r_s): Inside the horizon, delta_B
climbs from 0.004 back to 0.015--0.017. The black hole is feeding mass-energy
back into the Causal Set to restore it to equilibrium. The discrete spacetime
is being rebuilt.

### 3.2 The Equilibrium

The Causal Set's natural resting state is remarkably consistent:

- Big Bang post-wall mean: delta_B = 0.017 (Paper 5, Riner 2026e)
- Black hole near-singularity: delta_B = 0.015--0.017 (Paper 4, Riner 2026d)
- Across all three geometries tested: delta_B approximately equal to 0.017

This value is not tuned. It emerges from the Causal Set's defining property:
a random Poisson sprinkling of fundamental spacetime events with a Gaussian UV
suppression (g(k) = k^2 times exp(-k^2)). The randomness of the sprinkling
produces a distribution that naturally conforms to Benford's Law to within
delta_B approximately equal to 0.017. Deviations below this value (over-conformance) are as
unnatural as deviations above it. The equilibrium is a property of the
discrete structure, not a parameter of the model.

### 3.3 The One-Way Transaction

The critical feature of the mass-stripping cycle is its directionality.

The black hole spends mass to restore the Causal Set from over-conformance
(delta_B approximately equal to 0.003) back to equilibrium (delta_B approximately equal to 0.017). But the Causal Set
at equilibrium does not contribute mass back to the black hole. The restored
spacetime is inert --- it is the substrate, the geometry, not a source of
mass-energy.

The analogy from Session Notes: it is like heating a room to 72 degrees. You
spend energy to reach the target temperature. The room at 72 degrees does not
generate energy for you. And the energy you spent is gone.

The black hole is not being destroyed in the sense of annihilation. It is being
*used* --- as a conversion engine. Mass goes in. Spacetime structure comes out.
The structure does not refund the payment.

This is the mechanism we propose underlies black hole evaporation (Paper 4,
Section 5.3). The black hole does not radiate. It shrinks because the Causal
Set is consuming its mass-energy to maintain its own statistical law. The
present paper asks: what happens to the spacetime structure that is produced?

---

## 4. Mass In, Spacetime Out

### 4.1 The Conversion

The mass-stripping cycle converts mass (deviation from Benford conformance)
into Causal Set structure at equilibrium (delta_B approximately equal to 0.017). We propose that
this is a literal statement about geometry: the product of the conversion is
new spacetime.

In Causal Set Theory (Bombelli et al. 1987, Sorkin 2003), spacetime is not a
smooth manifold at the fundamental level. It is a discrete partial order ---
a set of events with causal relations between them. The "number of events in a
region" is the fundamental measure of spacetime volume. More events mean more
spacetime.

The mass-stripping cycle, in this picture, converts mass-energy into new
causal set events. The black hole's mass is consumed, and the product is an
increase in the number of fundamental spacetime events --- an increase in
spacetime volume. The Causal Set at equilibrium (delta_B approximately equal to 0.017) is not
just "restored spacetime." It is *additional* spacetime that did not previously
exist.

### 4.2 The Endpoint of Entropy

Paper 2 (Riner 2026b, Section 2.6) proposed that entropy is the tendency of
mass to return to zero deviation --- to perfect Benford conformance. That
paper stated: "Black hole evaporation converts the most extreme concentration
of mass back into radiation... the most dramatic return from maximum deviation
to the massless state." And: "The end state of entropy is not disorder. It is
the universe approaching conformance with the axiom as completely as mass
allows."

The Causal Set experiments now show that this return to conformance does not
produce radiation. It produces spacetime.

We therefore revise the endpoint: **the end state of entropy is not heat
death. It is geometry.** The thermodynamic arrow of time does not point toward
a universe of maximum disorder uniformly filled with cold radiation. It points
toward a universe of maximum spacetime --- a universe where all mass has been
converted, through the black hole mechanism, into Causal Set structure at
equilibrium. The axiom enforces itself, and the product of enforcement is the
substrate that implements the axiom.

This is a self-reinforcing loop: the Benford constraint converts mass into
the spacetime that carries the Benford constraint. The constraint builds
itself.

### 4.3 Quantitative Considerations

What is the scale of the conversion?

A stellar-mass black hole of 10 solar masses (approximately 2 times 10^31 kg) has a
Schwarzschild radius of approximately 30 km. The Hawking evaporation timescale under
standard theory is of order 10^67 years --- far longer than the current age
of the universe.

Under the Causal Set mechanism, the mass loss rate is proposed to scale
identically to Hawking's prediction (dM/dt proportional to 1/M^2), because the mechanism
is local: the same restructuring process occurs at every black hole regardless
of size, and only the ratio of restructuring cost to total mass changes
(Paper 4, Section 5.5). For stellar-mass and supermassive black holes, the
restructuring rate is negligible compared to the total mass.

The conversion rate relevant to cosmological expansion is therefore not the
evaporation rate of individual black holes (which is negligible for
astrophysical black holes) but the *aggregate* rate across all black holes in
the observable universe. There are estimated to be approximately 4 times 10^19 stellar-mass
black holes (Sicilia et al. 2022) and approximately 10^11 supermassive black
holes (one per galaxy). The total mass in black holes is estimated at
approximately 1% of the total baryonic mass of the universe.

Whether the aggregate conversion rate matches the observed expansion rate
(Hubble parameter H_0 approximately equal to 70 km/s/Mpc) is an open quantitative question
that we do not resolve here. The mechanism is proposed; the rate calculation
requires detailed work within the Causal Set formalism.

---

## 5. Black Holes as Expansion Sources

### 5.1 The Proposal

We propose that **every black hole is a site where the universe expands**.
The spacetime produced by the mass-stripping cycle propagates outward from
the black hole. The universe grows at every black hole, continuously, through
the conversion of mass into geometry.

This is a structural claim about the origin of cosmological expansion. In
standard cosmology, expansion is a property of the metric itself --- the scale
factor a(t) grows according to the Friedmann equations, driven by the total
energy density of the universe. Black holes are embedded in this expanding
spacetime but do not contribute to the expansion mechanism.

In the framework proposed here, the relationship is inverted. Expansion does
not happen everywhere uniformly. It *originates* at black holes --- the sites
where mass is converted to spacetime --- and propagates outward. The Friedmann
equations describe the averaged, large-scale result of this process. But the
fundamental process is local: mass goes in at a black hole, spacetime comes
out, and the new spacetime pushes the surrounding geometry apart.

### 5.2 Galactic Structure and Expansion Sources

Every galaxy in the observable universe is believed to harbor a supermassive
black hole (SMBH) at its center. The masses of these objects range from
approximately 10^5 to 10^10 solar masses. In the framework proposed here, each SMBH is
a spacetime factory, continuously converting accreted mass into new geometry.

The spatial distribution of expansion sources is therefore the spatial
distribution of galaxies. This has a specific structure:

- **Galactic centers**: Maximum expansion rate. The SMBH is the source.
  Surrounding mass (stars, gas, dark matter) brakes the propagation. The net
  expansion within a galaxy is small or zero --- gravitational binding
  overwhelms the local spacetime production.

- **Intergalactic space**: Moderate expansion. Spacetime propagates outward
  from the nearest galaxies. Some mass is present (intergalactic medium,
  dark matter filaments) to partially brake the propagation.

- **Voids**: Maximum observed expansion. The spacetime produced by distant
  black holes propagates into regions with very little mass. Nothing brakes
  the propagation. The expansion accelerates.

### 5.3 The Answer to Paper 2's Open Question

Paper 2 (Riner 2026b, Section 4.5) proposed that the accelerating expansion
of the universe in voids might reflect "absence of braking" rather than a
driving force. In that picture, expansion is not pushed by dark energy --- it
is simply not slowed down in regions without mass. The acceleration is the
natural behavior of unbraked emergence.

That proposal explained why voids accelerate but did not explain where the
expansion comes from in the first place. The present framework provides the
answer:

- **Source**: Black holes, converting mass to spacetime via the Causal Set
  mechanism.
- **Propagation**: New spacetime moves outward from black holes at a rate
  determined by the Causal Set dynamics.
- **Acceleration in voids**: The propagation encounters no mass to slow it.
  The absence of braking, proposed in Paper 2, operates on the spacetime
  produced by the black hole mechanism proposed here.

The two proposals --- Paper 2's "absence of braking" and Paper 4's "spacetime
factories" --- are complementary. One provides the source. The other explains
the acceleration. Together, they propose a complete picture of expansion
without a cosmological constant.

---

## 6. Gravitational Waves as New Spacetime

### 6.1 The Mass Deficit in Binary Mergers

On September 14, 2015, the LIGO interferometers detected the first
gravitational wave signal, GW150914 (Abbott et al. 2016). The event was a
binary black hole merger: two black holes with masses of approximately 36 and
29 solar masses spiraled together and merged into a single black hole of
approximately 62 solar masses. The mass deficit --- approximately 3 solar masses, or
5.4 times 10^47 joules --- was carried away as gravitational wave energy.

In standard general relativity, this interpretation is straightforward. The
Einstein field equations predict that accelerating masses radiate energy as
gravitational waves, analogous to accelerating charges radiating
electromagnetic waves. The 3-solar-mass deficit was radiated into the
gravitational wave field, propagated outward at the speed of light, and was
detected as a strain pattern in the LIGO interferometers approximately 1.3
billion years later.

### 6.2 The Reinterpretation

In the framework proposed here, the same event receives a different
interpretation:

Two black holes merge. The combined horizon restructures violently. The
Causal Set at the merged horizon undergoes the mass-stripping cycle at an
extreme rate --- the over-conformance near the combined horizon is severe,
and the restoration requires consuming approximately 3 solar masses of
mass-energy. That mass is converted into Causal Set structure at equilibrium.
The newly generated spacetime propagates outward at c.

The gravitational wave, in this picture, is not energy radiating away from the
source. It is **the wavefront of newly created spacetime**, propagating at
the speed of light because c is the propagation speed of the Benford
constraint itself (Paper 2, Section 2.2).

### 6.3 Same Predictions, Different Interpretation

The two pictures --- gravitational waves as radiated energy versus
gravitational waves as new spacetime --- make the same predictions for:

- **Waveform shape**: Both scale with the mass deficit and the dynamics of
  the merger. The chirp signal (increasing frequency and amplitude as the
  black holes spiral together) is determined by the orbital mechanics in both
  pictures.

- **Amplitude**: Proportional to the mass deficit divided by the distance.
  Both pictures agree.

- **Frequency**: Determined by the orbital frequency of the binary. Both
  pictures agree.

- **Propagation speed**: c in both pictures.

The predictions diverge in what happens to the wave energy after it passes:

| | Standard GR | CS Framework |
|---|---|---|
| Wave energy | Deposited in matter via tidal stretching | Becomes the space between matter |
| Net effect on spacetime | Transient distortion, then relaxation | Permanent increase in spacetime volume |
| At the source | Mass-energy radiated away | Mass-energy converted to geometry |
| Far from source | Energy density decreases as 1/r^2 | Spacetime volume increase persists |

**Table 2.** Comparison of gravitational wave interpretations.

In standard GR, a gravitational wave passing through a region temporarily
stretches and compresses spacetime, then the spacetime returns to its
original state. The energy is eventually deposited in matter (or continues
propagating). In the CS framework, the wave represents a permanent increase
in spacetime volume. The "stretching" detected by LIGO is real --- but it is
not a temporary distortion followed by relaxation. It is the leading edge of
new spacetime that, once created, persists.

### 6.4 The GW150914 Mass Budget

The mass deficit in GW150914 was approximately 3 solar masses. In Planck units,
this represents an enormous number of fundamental spacetime events. The
spacetime volume corresponding to this mass, under Causal Set theory where
the fundamental volume per event is of order the Planck volume (approximately
4.2 times 10^-105 m^3), would be:

    V approximately equal to (3 M_sun c^2) / E_P times V_P

where E_P is the Planck energy and V_P is the Planck volume. This represents
a staggeringly large number of new spacetime events produced in a fraction of
a second --- on the order of 10^78 new fundamental events.

Whether this quantity is consistent with the observed strain amplitude is
a quantitative question we do not resolve here. The point is that the mass
budget is large enough, in principle, to represent a macroscopic increase
in spacetime volume.

---

## 7. Connection to Croker et al. (2023)

### 7.1 Cosmological Coupling of Black Hole Masses

Croker, Weiner, and Farrah (2023) presented observational evidence for a
striking result: supermassive black holes in elliptical galaxies appear to
gain mass over cosmic time in a manner coupled to the cosmological scale
factor. Specifically, they reported:

    M proportional to a^k,    k approximately equal to 3

where M is the black hole mass and a is the cosmological scale factor. A
coupling constant k = 3 means that as the universe expands by a factor of 2,
the black hole mass increases by a factor of 8. This is consistent with the
black hole mass scaling as a cosmological energy density (which goes as a^3
for matter-like coupling, but their observed k approximately equal to 3 corresponds to
vacuum-energy-like coupling).

This result is controversial. It implies that black holes are not isolated
objects --- their masses are dynamically linked to the expansion of the
universe. Standard general relativity, in which black holes are described
by the Kerr metric embedded in an expanding Friedmann background, does not
predict this coupling (though the question of how to properly embed local
solutions in cosmological backgrounds remains technically subtle).

### 7.2 The Causal Set Mechanism as Physical Basis

We propose that the Causal Set mass-stripping mechanism provides a physical
basis for the Croker et al. coupling:

**Direction 1 --- Black hole to expansion**: The black hole converts mass to
spacetime via the Causal Set mechanism. This spacetime production contributes
to the increase in the scale factor a(t). More black holes, or more massive
black holes, produce more spacetime, driving faster expansion.

**Direction 2 --- Expansion to black hole**: As the scale factor increases,
the effective geometry around the black hole changes. The expanding background
modifies the horizon structure and the rate at which the Causal Set requires
restoration. If the expansion increases the effective over-conformance near
the horizon (by stretching the geometry), the restoration rate increases,
consuming more mass --- but this increased consumption itself produces more
spacetime, which increases the expansion, creating a feedback loop.

The coupling is therefore bidirectional:

    BH mass -> spacetime production -> expansion (a increases)
    expansion -> modified horizon geometry -> modified restoration rate -> BH mass change

Whether this feedback loop produces the specific scaling M proportional to a^k with
k approximately equal to 3 is a quantitative question that requires detailed calculation within
the Causal Set formalism. The point here is structural: the CS mechanism
naturally produces a coupling between black hole mass and the scale factor,
because the black hole is not a passive object embedded in an independently
expanding spacetime --- it is an active participant in producing the
expansion.

### 7.3 Consistency Check

The Croker et al. result, if confirmed, has a specific implication: black
holes in the early universe should be significantly less massive than their
present-day equivalents, even after accounting for accretion. The ratio is:

    M(z) / M(0) = [a(z) / a(0)]^3 = (1 + z)^{-3}

At redshift z = 1 (scale factor a = 0.5), the black hole mass would be 1/8
of its present value. At z = 2 (a = 1/3), it would be 1/27.

In the CS framework, this mass gain is not mysterious --- it is the flip side
of spacetime production. The expansion produced by the black hole feeds back
into the effective mass of the black hole through the cosmological coupling.
The black hole is not gaining mass from accretion; it is gaining effective
mass because the spacetime it produced has expanded the geometry in which it
sits.

This is speculative but testable. Observations of supermassive black holes
at high redshift (e.g., with JWST) can constrain whether the mass evolution
follows the predicted (1 + z)^{-3} scaling or is better explained by
standard accretion models.

---

## 8. Dark Energy Reinterpreted

### 8.1 The Standard Picture

In standard cosmology, the accelerating expansion of the universe (Riess
et al. 1998, Perlmutter et al. 1999) is attributed to dark energy --- a
component of the cosmic energy budget with negative pressure that drives
expansion. The simplest model is the cosmological constant Lambda, corresponding
to a constant vacuum energy density. More exotic models invoke quintessence
(a slowly rolling scalar field), phantom energy (w < -1), or modifications to
general relativity at cosmological scales.

Dark energy constitutes approximately 68% of the total energy density of the
universe. Its nature is unknown. No dark energy particle or field has been
detected. The cosmological constant, while mathematically simple, faces the
fine-tuning problem: the observed value of Lambda is approximately 120 orders of
magnitude smaller than the naive quantum field theory prediction.

### 8.2 The Benford Alternative

We propose that the accelerating expansion does not require a cosmological
constant, a new field, or any modification to gravity. It requires two things:

1. **A source of new spacetime**: Black holes, converting mass to Causal Set
   structure at equilibrium via the mass-stripping cycle.

2. **Regions where nothing brakes the propagation**: Voids, where the absence
   of mass allows the newly produced spacetime to propagate freely.

The combination produces acceleration without a driving force. The black holes
provide the spacetime. The voids provide the unbraked propagation. The result
is an expansion that accelerates in low-density regions and is negligible in
high-density regions --- which is precisely what is observed.

### 8.3 Why the Expansion Appears Uniform

If expansion originates at black holes (which are localized) and accelerates
in voids (which are distributed unevenly), why does the expansion appear
nearly isotropic and homogeneous on large scales?

Because on large scales (> 100 Mpc), the distribution of galaxies --- and
therefore black holes --- is approximately uniform. The cosmic web structure
(filaments, voids, clusters) averages out. The Friedmann equations, which
assume homogeneity and isotropy, correctly describe the averaged behavior.
The individual spacetime production at each black hole is a local process;
the cosmological expansion is the large-scale statistical average of
approximately 10^11 such local processes.

This is analogous to how the temperature of a gas is a macroscopic average
of individual molecular kinetic energies. Each molecule moves in a specific
direction with a specific speed. The temperature is none of those --- it is
the statistical aggregate. Similarly, each black hole produces spacetime in
its local environment. The Hubble flow is the statistical aggregate.

### 8.4 The Phantom Energy Verdict

Paper 3 (Riner 2026c, Experiment 5) tested phantom energy (w < -1) through
the Benford existence filter. The result: UNDEFINED. Zero valid modes at
every momentum. The thermal partition function does not converge. Phantom
energy cannot form a physical thermal distribution.

If the Benford existence filter is correct, phantom energy does not exist.
This eliminates the most exotic dark energy models and the "Big Rip" scenario
in which phantom energy eventually tears apart all structure. The CS framework
replaces phantom energy with a physical mechanism (mass-to-spacetime
conversion) that produces acceleration without invoking a thermodynamically
impossible energy component.

### 8.5 Connection to Paper 2's "Absence of Braking"

Paper 2 (Riner 2026b, Section 4.5) proposed that dark energy might not be a
substance at all but rather the absence of a braking mechanism. In regions
with mass, the emergence of the Benford constraint is slowed by the presence
of deviation (mass). In regions without mass (voids), there is nothing to
slow the emergence. The constraint propagates freely, and the result looks
like accelerating expansion.

The present paper adds the source. Paper 2 described the brake (or lack
thereof). This paper describes the engine: black holes producing spacetime.
The complete picture:

- **Engine**: Black holes converting mass to spacetime at every galactic center.
- **Transmission**: New spacetime propagates outward.
- **Brake**: Mass in the propagation path slows it (gravitational binding,
  mass density).
- **Open road**: Voids, where no brake operates and the propagation
  accelerates.
- **Observation**: On large scales, the average appears as a uniformly
  accelerating expansion characterized by H_0 and Lambda --- but Lambda is an
  effective parameter, not a fundamental constant.

---

## 9. The Black Hole Lifecycle Revisited

### 9.1 Size-Independent Local Process

Paper 4 (Riner 2026d, Section 5.5) established that the Causal Set
mass-stripping cycle is a local process. The simulation uses local geometry ---
the spectrum at each radius depends on the local effective temperature, not
on the total black hole mass. The spike-strip-restore pattern looks identical
regardless of the black hole's size.

The difference between black holes of different masses is not the rate of the
local process but the ratio of the restructuring cost to the total mass. This
produces a natural lifecycle:

### 9.2 Small Black Holes (M close to M_P): Catastrophic Structural Collapse

For a black hole near the Planck mass, the total mass is comparable to the
energy consumed in a single restructuring cycle. The Causal Set's demand for
restoration represents a significant fraction of the entire black hole. The
mass drops. The effective temperature rises. The restructuring demand
intensifies. The process runs away.

We propose that the "final evaporation" of a small black hole is not an
explosion of Hawking radiation but a **catastrophic structural collapse**:
the Causal Set consumes the remaining mass in a rapidly accelerating cascade.
The endpoint is not a burst of photons at the Hawking temperature. It is the
complete conversion of the black hole's mass into spacetime --- the final,
total enforcement of the Benford constraint.

What an external observer sees during this event is an open question
(Section 11). If nothing leaves the horizon at any point, the observer sees
the black hole simply disappear --- its gravitational signature vanishes as
the mass is converted. If the final conversion produces detectable
gravitational effects (a transient in the local metric), those effects would
carry information about the conversion process.

### 9.3 Stellar-Mass Black Holes (M close to 10 M_sun): Stable Engines

For a stellar-mass black hole, the restructuring cost is negligible compared
to the total mass. The black hole loses mass continuously but immeasurably
slowly. On any observationally relevant timescale (billions of years), the
black hole is effectively stable.

These are the workhorse engines of the expansion. They continuously produce
spacetime at a rate too low to affect their own mass appreciably but,
aggregated across 10^19 stellar-mass black holes in the observable universe,
potentially sufficient to drive macroscopic expansion.

The standard Hawking evaporation time for a 10-solar-mass black hole is
approximately 10^67 years. Under the CS mechanism, the mass loss rate scales
identically (dM/dt proportional to 1/M^2), so the timescale is the same. The
difference is not the rate but the product: Hawking predicts outgoing thermal
radiation; the CS mechanism predicts spacetime production with no outgoing
radiation.

### 9.4 Supermassive Black Holes (M close to 10^9 M_sun): Negligible Drain, Maximum Output

Supermassive black holes at galactic centers have masses of millions to
billions of solar masses. The restructuring cost is utterly negligible
relative to the total mass. The black hole does not measurably evaporate ---
consistent with all observations of active galactic nuclei and quiescent
galactic centers.

But the spacetime production is not negligible. The SMBH continuously
converts accreted mass into spacetime structure. The accretion rate for active
galaxies can be substantial --- up to solar masses per year for the most
luminous quasars. This accreted mass enters the black hole, is processed
through the mass-stripping cycle, and some fraction is converted to spacetime.

We propose that the balance between accretion (mass input) and Causal Set
restoration (mass consumption) determines the net growth or shrinkage of the
SMBH. During active phases (high accretion), the SMBH grows because
accretion exceeds restoration. During quiescent phases, restoration slowly
consumes existing mass. The overall evolution couples to the cosmological
expansion through the Croker et al. mechanism (Section 7).

### 9.5 The Scaling Relation

The three regimes reproduce the standard scaling relation T_H proportional to 1/M
(Hawking temperature inversely proportional to mass) without invoking outgoing
radiation:

- **Small BH**: High effective temperature, rapid mass loss, violent endpoint.
  Not because the radiation temperature is high --- because the restructuring
  demand is a large fraction of the total mass.

- **Stellar BH**: Moderate effective temperature, negligible mass loss, long
  lifetime. Not because the radiation is faint --- because the restructuring
  demand is a tiny fraction of a large mass.

- **Supermassive BH**: Very low effective temperature, no measurable mass
  loss, effectively permanent. Not because the radiation is undetectable ---
  because the restructuring demand is infinitesimal relative to an enormous
  mass.

The physics is the same at every scale. The phenomenology differs because
of the ratio of local process to global mass.

---

## 10. Testable Predictions

The framework proposed in this paper makes specific predictions that
distinguish it from standard cosmology. Some are currently testable; others
require future observational capabilities.

### 10.1 Predictions That Distinguish from Standard Hawking Radiation

**Prediction 1: No outgoing thermal flux from black holes.** If the
evaporation mechanism is structural (Causal Set consumption) rather than
radiative (Hawking emission), there is no outgoing thermal radiation at T_H.
A sufficiently sensitive experiment near a small, evaporating black hole
would detect mass loss without corresponding outgoing particles.

**Prediction 2: The final evaporation is not a photon burst.** Standard
Hawking theory predicts that the final moments of black hole evaporation
produce a burst of high-energy radiation as the temperature diverges
(T_H proportional to 1/M, so T goes to infinity as M goes to 0). The CS mechanism predicts a
catastrophic structural collapse with no outgoing radiation --- the mass
simply converts to spacetime. The two scenarios produce different
electromagnetic and gravitational wave signatures at the endpoint.

### 10.2 Predictions That Distinguish from the Cosmological Constant

**Prediction 3: Expansion rate correlates with black hole density.** If
expansion originates at black holes, regions of the universe with higher
black hole density should show subtly higher local expansion rates. Standard
Lambda-CDM predicts a uniform expansion rate (at fixed epoch) independent of local
black hole density. The effect would be small (swamped by peculiar velocities
on small scales) but potentially detectable in precision measurements of the
Hubble flow as a function of environment.

**Prediction 4: Voids expand faster than the global average.** This is
predicted by both the CS framework (unbraked propagation of black-hole-
produced spacetime) and, to some extent, by standard structure formation
models (voids expand faster because they are underdense). The CS framework
predicts a stronger effect --- the acceleration in voids is not merely
kinematic (underdense regions expand faster in any cosmology) but dynamic
(new spacetime is being deposited with nothing to slow it).

**Prediction 5: The effective Lambda is not constant.** If the "dark energy" is
really an effective parameter arising from black hole spacetime production,
then Lambda(t) should evolve as the black hole population evolves. In the early
universe (fewer and smaller black holes), the effective Lambda should be smaller.
As black holes form and grow (after the epoch of first star formation), the
effective Lambda should increase. Standard Lambda-CDM predicts a strictly constant Lambda.
Observations of the dark energy equation of state parameter w(z) at high
redshift (e.g., from DESI, Euclid, or the Roman Space Telescope) could
constrain this prediction.

### 10.3 Predictions About Gravitational Waves

**Prediction 6: Gravitational waves permanently increase spacetime volume.**
Standard GR predicts that a gravitational wave passing through a region
causes a temporary distortion that relaxes after the wave passes. The CS
framework predicts a permanent (though extremely small) increase in spacetime
volume. This could in principle be tested by precision measurements of
spatial geometry before and after a strong gravitational wave event, though
the effect would be extraordinarily small for any currently detectable wave.

**Prediction 7: The mass deficit in binary mergers corresponds to spacetime
volume.** In the CS framework, the 3-solar-mass deficit in GW150914 became
spacetime. The volume of the gravitational wave shell at any given radius
should correspond (in appropriate units) to the mass deficit converted to
spacetime events. This is a quantitative prediction that can be checked
against the detailed waveform.

### 10.4 Predictions About Black Hole Mass Evolution

**Prediction 8: The Croker et al. scaling is physical.** If the cosmological
coupling M proportional to a^k with k approximately equal to 3 is confirmed by future observations,
the CS framework predicts that this coupling arises from the bidirectional
relationship between black hole mass and expansion (Section 7.2). The
specific value of k should be derivable from the CS restructuring rate and
the cosmological background evolution.

**Prediction 9: High-redshift SMBHs should be undermassive.** The coupling
predicts that supermassive black holes at high redshift should be
significantly less massive than their low-redshift counterparts, beyond what
accretion models predict. Specifically, M(z) proportional to (1 + z)^{-3}. JWST
observations of quasars at z > 6 can constrain this prediction.

---

## 11. Open Questions

The following questions are open. Each represents a necessary step toward
making the proposals in this paper quantitatively rigorous.

**1. Can the magnitude of the CS spike at the horizon be quantitatively
related to the Hawking temperature?** The mass-stripping cycle produces a
specific delta_B profile near the horizon. Does the depth of the over-conformance
(delta_B approximately equal to 0.003) map to a specific energy cost per unit area, and does
that cost correspond to the Bekenstein-Hawking entropy S = A/(4 l_P^2)?

**2. Does the width of the spike-and-relax region scale with black hole
mass?** The simulation uses a single black hole mass. If the radial extent
of the stripping region (r approximately equal to 10 to 1 r_s) and the restoration region
(r approximately equal to 1 to 0.01 r_s) scale differently with M, this would produce
mass-dependent signatures in the conversion rate.

**3. Would a wormhole with a mild horizon trigger a partial CS response?**
The Morris-Thorne wormhole used in the control experiment has no horizon.
Near-extremal geometries (e.g., a wormhole with a very thin shell of exotic
matter near the throat) might produce a partial horizon, and the CS response
to such a geometry would test whether the trigger is strictly topological
(horizon present/absent) or admits intermediate cases.

**4. Is there a minimum curvature threshold below which CS does not respond,
or is the trigger specifically topological?** The wormhole has extreme
curvature (T_max approximately equal to 1.59 T_P) but no CS response. Is the distinction
purely topological (horizon vs. no horizon), or is there a curvature
threshold that must be exceeded in the presence of a horizon?

**5. If black holes do not radiate, what does an external observer see during
the final evaporation?** Does the black hole disappear suddenly (mass
converts to spacetime instantaneously in the final cascade)? Does it leave
a topological remnant (a residual causal structure with zero mass)? Does the
final conversion produce detectable gravitational effects (a transient in
the local metric)?

**6. Does CS healing preserve information (unitarity), or is it a genuinely
dissipative process?** The causal set restructuring proposed here reorganizes
causal relations. If this reorganization is bijective, unitarity is
preserved. If it is many-to-one, information is lost at the fundamental
level. The question is about the dynamics of the causal set, not about
semiclassical radiation.

**7. Can the rate of spacetime production by black holes be quantitatively
matched to the observed expansion rate?** The aggregate spacetime production
across all black holes in the observable universe must produce the observed
Hubble flow. This is a calculation that requires (a) the individual
production rate from the CS formalism, (b) the black hole mass function, and
(c) the propagation dynamics of newly produced spacetime.

**8. Does the mass deficit in LIGO binary mergers correspond quantitatively
to the spacetime volume of the gravitational wave shell?** The
3-solar-mass deficit in GW150914 corresponds to approximately 10^78 Planck
volumes. Does this number match the spacetime volume of the gravitational
wave shell at any given radius, computed from the observed strain?

**9. Is the Croker et al. BH-expansion coupling consistent with CS spacetime
production rates?** If M proportional to a^k with k approximately equal to 3, the CS mechanism must
produce spacetime at a rate that, fed back through the cosmological
equations, reproduces this specific scaling. Does it?

---

## 12. Discussion

### 12.1 Three Tiers of Claims

This paper makes claims at three distinct epistemic levels. We are explicit
about which is which.

**Tier 1 --- Established results** (Papers 1--6):

- Benford's Law is satisfied exactly by completely monotonic distributions.
- delta_B is an invertible measurement instrument that recovers dimensionality,
  number theory identities, and mass.
- The existence filter identifies thermodynamically impossible physics.
- Causal Set Theory produces the cleanest statistical structure inside a black
  hole (delta_B = 0.011).
- The CS response is topology-specific (singularities and horizons), not
  curvature-driven.
- The mass-stripping cycle is a demonstrated pattern in the computational data.

These are the computational facts. They have been reproduced, tested with
control experiments, and reported in the preceding papers.

**Tier 2 --- Proposals** (Paper 4, Sections 5--6, and the present paper):

- The energy cost of CS restructuring at horizons is drawn from the black
  hole's mass.
- Black holes evaporate through mass consumption, not outgoing radiation.
- The product of the conversion is new spacetime (Causal Set at equilibrium).
- Gravitational waves from mergers are the wavefront of newly generated
  spacetime.

These are interpretive claims. They are consistent with the computational
data but require independent theoretical confirmation within the Causal Set
formalism. Specifically, one must show that the dynamics of a causal set in
the presence of a horizon produce an energy cost proportional to the
restructuring, and that this cost maps to the Hawking mass loss rate.

**Tier 3 --- Speculations** (Sections 5, 7, 8 of the present paper):

- Black holes are the source of cosmological expansion.
- The accelerating expansion does not require a cosmological constant.
- The Croker et al. coupling arises from the bidirectional relationship
  between black hole mass and expansion.
- The effective Lambda evolves with the black hole population.

These are the furthest-reaching claims. They follow logically from the Tier 2
proposals but add cosmological extrapolation. Even if the Tier 2 proposals
are correct, the Tier 3 speculations could fail if the quantitative rates
do not match observations or if additional physics intervenes at cosmological
scales.

### 12.2 What Would Falsify This Framework

The framework is falsifiable at each tier:

**Tier 1 is falsified** if the computational results do not reproduce ---
if the delta_B values change under different numerical implementations, or if
the Causal Set's advantage disappears with different momentum grids or
temperature profiles. We have tested this (Paper 4, Appendix) and found
the results robust, but independent reproduction is necessary.

**Tier 2 is falsified** if Hawking radiation is detected. A confirmed
measurement of outgoing thermal particles at T_H from a black hole would
establish that the evaporation mechanism is radiative, not structural. The
CS mass consumption proposal and the standard Hawking radiation proposal
make opposite predictions about outgoing flux. One is right.

**Tier 3 is falsified** if the expansion rate is rigorously constant (no
evolution of effective Lambda), if it shows no correlation with black hole
density, or if the Croker et al. result is ruled out observationally. Any of
these would undermine the specific cosmological claims, even if the
underlying mass-to-spacetime mechanism (Tier 2) remains viable.

### 12.3 Relation to Existing Cosmological Models

The proposal that expansion originates at black holes has precursors in the
literature, though the specific mechanism we propose (Benford-driven Causal
Set restructuring) is new.

**Cosmological natural selection** (Smolin 1997): Proposed that black holes
produce new universes inside them, with slightly varied physical constants.
Our proposal is more modest --- black holes produce new spacetime in *this*
universe, not new universes.

**Cosmological coupling** (Croker et al. 2023): Proposed that black hole
masses are coupled to the scale factor. Our proposal provides a physical
mechanism for this coupling.

**Emergent gravity** (Verlinde 2011): Proposed that gravity is an entropic
force arising from information on holographic screens. Our framework shares
the idea that gravity and thermodynamics are deeply connected but proposes a
specific mechanism (Benford conformance enforcement) rather than a general
entropic principle.

**Penrose's conformal cyclic cosmology** (Penrose 2010): Proposed that the
remote future of one universe (all mass evaporated, only radiation remaining)
becomes the Big Bang of the next. Our framework's "endpoint of entropy is
geometry" resonates with Penrose's idea that the universe evolves toward a
massless state, but we propose the endpoint is new spacetime rather than
conformal rescaling.

### 12.4 The Self-Reinforcing Loop

The most striking feature of the framework, if correct, is its
self-referential structure:

1. The Benford constraint demands conformance.
2. Mass violates the constraint.
3. Black holes collect mass and convert it to spacetime.
4. The spacetime implements the Benford constraint.
5. The constraint demands conformance. (Return to step 1.)

The constraint builds the substrate that carries the constraint. The universe
expands because the mathematical law governing its statistical structure
requires expansion to enforce itself. Mass is the violation. Black holes are
the enforcement mechanism. Spacetime is the product. And the product is the
medium in which the law operates.

This circularity is not a logical flaw --- it is the structure of a
self-consistent system. The Benford constraint is not imposed from outside.
It is a property of the spacetime itself, maintained by the spacetime itself,
at the cost of the mass embedded in the spacetime. The universe is a system
that enforces its own statistical law, and the enforcement produces the
universe.

---

## 13. Conclusion

We have followed a single axiom --- P(d) = log_10(1 + 1/d) --- from Benford's
Law through quantum statistics, past the event horizon, and out to the
expansion of the universe. The chain:

**One axiom.** Benford's distribution as the foundational constraint on
physical distributions.

**One formula.** delta_B = the square root of the sum of squared deviations from the Benford
prediction. A measurement of how far any distribution sits from the
logarithmic ideal.

**Mass goes in.** Black holes collect mass --- the physical manifestation of
deviation from the axiom. The gravitational field strips the Causal Set below
its equilibrium near the horizon.

**Spacetime comes out.** Inside the black hole, the mass is consumed by the
Causal Set to restore its equilibrium at delta_B approximately equal to 0.017. The product is
new spacetime geometry. The transaction is one-way.

**The universe expands.** The newly produced spacetime propagates outward from
black holes. In voids, nothing brakes the propagation. The result is an
accelerating expansion that, on large scales, mimics a cosmological constant.

**The constraint builds itself.** The spacetime produced by the mechanism is
the substrate that carries the Benford constraint. The constraint demands
conformance. Mass violates conformance. Black holes enforce conformance. The
enforcement produces spacetime. The spacetime carries the constraint.

Every link in this chain is either established computationally (Papers 1--6)
or proposed with explicit epistemic caveats (this paper). We have not proven
that black holes produce expansion. We have constructed a framework in which
this follows naturally from a single mathematical axiom, and we have
identified nine open questions whose resolution would confirm or refute the
proposal.

The framework makes specific, falsifiable predictions: no outgoing Hawking
radiation, expansion rate correlated with black hole density, evolving
effective cosmological constant, permanent spacetime volume increase from
gravitational waves. Each can in principle be tested. Several will be
constrained by current and near-future observations (JWST, DESI, Euclid,
LISA, the Einstein Telescope).

If the framework is wrong, it will be wrong in an instructive way ---
because the computational results (Causal Set Theory's response to
singularities and horizons, the mass-stripping cycle, the wormhole control)
are real patterns that require explanation regardless of whether the
cosmological extrapolation holds.

If it is right, then the universe is simpler than we thought. One axiom.
One formula. Mass goes in, spacetime comes out, the constraint builds itself,
and the universe expands because it must.

---

## References

- Abbott, B. P. et al. (LIGO Scientific Collaboration and Virgo
  Collaboration). (2016). Observation of gravitational waves from a binary
  black hole merger. Phys. Rev. Lett. 116, 061102.

- Bombelli, L., Lee, J., Meyer, D., & Sorkin, R. D. (1987). Space-time as a
  causal set. Phys. Rev. Lett. 59, 521--524.

- Croker, K. S., Weiner, J. L., & Farrah, D. (2023). Cosmologically coupled
  compact objects: a single-parameter model for LIGO-Virgo mass and redshift
  distributions. Astrophys. J. Lett. 944, L31.

- Hawking, S. W. (1975). Particle creation by black holes. Commun. Math.
  Phys. 43, 199--220.

- Penrose, R. (2010). Cycles of Time: An Extraordinary New View of the
  Universe. Bodley Head.

- Perlmutter, S. et al. (1999). Measurements of Omega and Lambda from 42
  high-redshift supernovae. Astrophys. J. 517, 565--586.

- Riess, A. G. et al. (1998). Observational evidence from supernovae for an
  accelerating universe and a cosmological constant. Astron. J. 116,
  1009--1038.

- Riner, C. (2026a). Complete monotonicity and Benford's Law: deriving quantum
  statistics from the significant digit distribution.

- Riner, C. (2026b). The Law of Emergence: Benford's distribution as a
  universal constraint on physical reality.

- Riner, C. (2026c). The Benford deviation as a measurement instrument:
  round-trip calibration, fingerprint classification, and an existence filter
  for exotic physics.

- Riner, C. (2026d). Benford's Law inside a black hole: statistical structure
  beyond the event horizon and a Causal Set mechanism for evaporation.

- Riner, C. (2026e). Benford's Law at the Planck Wall: ten quantum gravity
  models through the cosmological singularity.

- Riner, C. (2026f). Benford's Law through a wormhole: the Casimir effect and
  the absence of singularity response.

- Riner, C. (2026g). [This paper].

- Sicilia, A. et al. (2022). The black hole mass function across cosmic time.
  Astrophys. J. 924, 56.

- Smolin, L. (1997). The Life of the Cosmos. Oxford University Press.

- Sorkin, R. D. (2003). Causal sets: discrete gravity. In Lectures on Quantum
  Gravity, ed. A. Gomberoff & D. Marolf. Springer.

- Verlinde, E. (2011). On the origin of gravity and the laws of Newton. JHEP
  2011, 29.

---

## Appendix A: The Full Argument Chain (Condensed)

For reference, the complete chain from axiom to expansion:

    1. P(d) = log_10(1 + 1/d)                     [Axiom]
    2. Mass = deviation from (1)                    [Paper 2, Sec 2.3]
    3. Entropy = return to conformance with (1)     [Paper 2, Sec 2.6]
    4. CS = substrate implementing (1) at           [Papers 4-6]
       equilibrium delta_B ~ 0.017
    5. BH = enforcement mechanism for (1):          [Paper 4, Sec 5.1]
       mass in, CS at equilibrium out
    6. Evaporation = CS consuming BH mass,          [Paper 4, Sec 5.3]
       not outgoing radiation
    7. Expansion = new spacetime propagating         [This paper, Sec 5]
       outward from BHs
    8. GW = wavefront of (7) during mergers         [This paper, Sec 6]

## Appendix B: Comparison of Expansion Mechanisms

| Feature | Lambda-CDM | CS Framework |
|---------|-----------|--------------|
| Dark energy nature | Cosmological constant or vacuum energy | Effective parameter from BH spacetime production |
| Source of expansion | Intrinsic property of spacetime | Black holes converting mass to geometry |
| Why voids accelerate | Underdensity (kinematic) | Unbraked spacetime propagation (dynamic) |
| Lambda constant? | Yes (by assumption) | No --- evolves with BH population |
| Requires new physics | Yes (unknown dark energy) | No (CS + Benford constraint) |
| Phantom energy | Allowed (w < -1 models) | Ruled out (UNDEFINED in Benford filter) |
| BH-expansion coupling | Not predicted | Natural consequence |
| Falsifiable by | Constant w(z) across all z | Detecting outgoing Hawking radiation |

## Appendix C: Summary of Open Questions

| # | Question | What would resolve it |
|---|----------|----------------------|
| 1 | CS spike magnitude related to Hawking temperature? | Analytic calculation in CS formalism |
| 2 | Spike width scales with BH mass? | Multi-mass CS simulations |
| 3 | Partial horizon triggers partial CS response? | Near-extremal geometry simulations |
| 4 | Minimum curvature threshold vs topological trigger? | Curvature sweep at fixed topology |
| 5 | External observer sees what at final evaporation? | CS dynamics at M approaching 0 |
| 6 | CS healing preserves unitarity? | Bijectivity analysis of CS restructuring |
| 7 | Aggregate BH production matches H_0? | Rate calculation + BH mass function |
| 8 | GW150914 mass deficit matches spacetime volume? | Quantitative comparison |
| 9 | Croker coupling consistent with CS rates? | Combined CS + cosmological calculation |
