# The Euler Prime Substrate: Replacing the Causal Set Dimension with Number Theory in a Five-Dimensional Metric Tensor

### Christopher Riner
### Chesapeake, Virginia
### chrisriner45@gmail.com

**Draft — February 2026**

---

## Abstract

We present a reformulation of the five-dimensional Benford metric (Riner 2026a, 2026b, 2026c) in which the Causal Set dimension is replaced by a single equation from number theory: Euler's product formula for the Riemann zeta function. The substitution is motivated by a new empirical finding — that prime numbers are perfectly anti-Benford across all integer bases 2 through 10 — and by a structural correspondence between the five axioms of Causal Set theory and the properties of the prime numbers.

The Euler product ζ(s) = Π_p 1/(1 − p^(−s)) encodes discreteness, partial ordering, local finiteness, countability, and the derivation of continuous structure from ordering relations — all five Causal Set axioms — in a single formula. Its logarithm, ln ζ(s), serves as the metric component in the fifth dimension of the tensor, replacing the previous Causal Set metric g_δ̂ = log₁₀(1 + 1/δ_B).

The reformulated tensor is validated against GPS time dilation measurements. Using entropy rate β = √(r_s/r) derived from spatial geometry — with no time dimension — the framework reproduces the measured GPS satellite clock correction of approximately 38.6 μs/day to the same precision as standard general relativity. The Euler prime substrate is inert in weak gravitational fields (the determinant at Earth's surface exceeds the Benford floor by a factor of 10³⁵) and activates only in the strong-field regime inside black holes.

A multi-base heatmap analysis of prime numbers provides the empirical foundation: for each number in a test range, leading digits in bases 2 through 10 are compared against Benford's predicted distribution. Primes show zero deviation from composite numbers across all bases — they carry no Benford signature whatsoever. This is the only dataset, out of all domains tested in this research program (quantum spectra, black hole metrics, gravitational curvature, financial and demographic data), that exhibits this property. We interpret this as evidence that primes occupy a layer beneath Benford-conformant emergence: they are the anti-Benford substrate on which Benford-like physical structure is built.

---

## 1. Introduction

The previous papers in this series established three results. First, that Benford's Law deviation δ_B provides a meaningful diagnostic for quantum gravity proposals at extreme boundaries — the Planck temperature wall, the black hole singularity, and traversable wormhole throats (Riner 2026a, 2026b). Second, that a four-dimensional spatial metric with a Benford floor constraint and a Causal Set (CS) dimension prevents the Schwarzschild singularity while leaving general relativity unchanged outside the event horizon (Riner 2026c). Third, that this framework makes testable predictions for binary black hole mergers, wormhole traversal mechanics, and the nature of black hole evaporation.

A persistent limitation of the framework was that the CS dimension, while empirically motivated, required five separate axioms for its definition — discreteness, partial ordering, local finiteness, countability, and the derivation of geometry from causal ordering. These axioms could not be condensed into a single equation suitable for direct insertion into the metric tensor. The CS metric component g_δ̂ = log₁₀(1 + 1/δ_B) was derived from a fitted curve (δ_B(r) = 0.003389|ln(r/r_s)| + 0.002508), which, while accurate, lacked a fundamental derivation.

This paper resolves both limitations. We show that prime numbers satisfy all five Causal Set axioms, that they are empirically anti-Benford (carrying zero Benford signature across nine integer bases), and that Euler's product formula — a single equation from 1737 — encodes the complete axiomatic content of the Causal Set in a form suitable for the metric tensor.

### 1.1 Structure of this paper

Section 2 presents the empirical discovery: the multi-base Benford analysis of prime numbers and the anti-Benford finding. Section 3 establishes the structural correspondence between prime numbers and the five Causal Set axioms. Section 4 introduces the Euler product as the replacement for the CS dimension and constructs the reformulated tensor. Section 5 validates the new tensor against GPS time dilation using entropy rate in place of a time dimension. Section 6 discusses the three-layer emergence hierarchy implied by the results. Section 7 summarizes predictions and next steps.

---

## 2. Primes Are Anti-Benford: An Empirical Finding

### 2.1 Method

Benford's Law states that in many naturally occurring datasets, the leading digit d appears with probability:

> P(d) = log₁₀(1 + 1/d),  d = 1, 2, ..., 9

This generalizes to arbitrary base B:

> P_B(d) = log_B(1 + 1/d),  d = 1, 2, ..., B−1

We constructed a multi-base heatmap analysis to test whether prime numbers exhibit Benford-like structure that could be exploited for prime prediction. For each integer n in the range [2, 10000], and for each base B in {2, 3, 4, 5, 6, 7, 8, 9, 10}:

1. Compute the leading digit d of n in base B.
2. Using the full dataset, compute the observed frequency of leading digit d among primes and among composites in that base.
3. Assign a log-likelihood ratio score: log₂(P(d | prime) / P(d | composite)).
4. Sum scores across all nine bases to produce a composite score for each number.

If primes deviate from Benford's prediction differently than composites, this composite score should separate the two populations.

### 2.2 Result

No separation was observed. The score distributions for primes and composites overlap completely. No single base, and no combination of bases, produced a statistically meaningful signal. Primes are indistinguishable from composites under multi-base Benford analysis.

This result is stronger than "noisy" or "weak signal." The heatmap visualization (Figure 1) shows that primes resist pattern formation more thoroughly than random noise — they do not even produce the accidental clumping that random distributions exhibit. The absence of structure is itself perfectly consistent across all nine bases.

### 2.3 Context

This is the only anti-Benford result in the research program. All other datasets tested — thermal radiation spectra, quantum gravity models, black hole metric components, gravitational curvature invariants, financial data, population data — exhibit Benford-like structure to varying degrees. Prime numbers are the singular exception.

We interpret this as follows: Benford's Law is a signature of systems that arise from multiplicative processes (Pinkham 1961, Hill 1995). Primes are defined by the absence of multiplicative structure — they are the irreducible atoms of multiplication. The sieve of Eratosthenes systematically removes every multiplicative pattern (multiples of 2, multiples of 3, multiples of 5, ...). What remains is everything that survived the deletion of all patterns. An anti-Benford result is therefore not surprising — it is predicted by the nature of primality itself.

---

## 3. Primes Satisfy All Five Causal Set Axioms

Causal Set theory (Bombelli, Lee, Meyer, Sorkin 1987) proposes that spacetime at the fundamental level is a locally finite partially ordered set. The theory is defined by five properties:

| # | Causal Set Axiom | Prime Numbers |
|---|---|---|
| 1 | **Discrete**: Elements are individual, separated points | Primes are integers — discrete, no continuum |
| 2 | **Partial order**: Reflexive, antisymmetric, transitive relation | Divisibility ordering on integers is a partial order; primes are its atoms |
| 3 | **Locally finite**: Finitely many elements between any two bounds | Finitely many primes between any two integers |
| 4 | **Countable**: The set is countably infinite | The primes are countably infinite (Euclid, ~300 BCE) |
| 5 | **Order produces geometry**: Sorkin's "Order + Number = Geometry" | The divisibility ordering of primes produces the structure of all integers via the fundamental theorem of arithmetic |

The correspondence is exact. Both primes and causal set elements are discrete, partially ordered, locally finite, countable substrates whose ordering relations generate continuous structure.

Both are also anti-Benford. The Causal Set spectrum (Experiment 6c, Riner 2026a) showed only one of five characteristic equations with Benford-like behavior. Primes show zero Benford-like behavior across nine bases.

### 3.1 Myrheim-Meyer Dimension Test

To test Axiom 5 quantitatively, we apply the Myrheim-Meyer dimension estimator to the prime divisibility poset and compare against a standard causal set constructed by random sprinkling into 2D Minkowski space.

**Method.** For a set of N elements with a partial ordering, the *ordering fraction* f is the ratio of ordered pairs (pairs where one element precedes the other) to total pairs. The Myrheim-Meyer estimator relates this fraction to the effective dimension d of the underlying space:

> f(d) = Γ(d+1) Γ(d/2) / (4 Γ(3d/2))

For reference: f(2) ≈ 0.333, f(3) ≈ 0.147, f(4) ≈ 0.070, f(5) ≈ 0.035, f(6) ≈ 0.018.

**Setup.** Two structures with 499 elements each:
1. **Causal Set (control)**: 499 points sprinkled uniformly into a 2D Minkowski region [0,10]×[0,10]. Ordering: a ≺ b iff a is in the causal past of b (timelike separation).
2. **Prime Poset**: integers 2 through 500 with divisibility ordering: a ≺ b iff a | b and a ≠ b.

**Results.**

| Property | Causal Set (2D Minkowski) | Prime Poset |
|---|---|---|
| Elements | 499 | 499 |
| Ordered pairs | 60,520 | 2,191 |
| Ordering fraction | 0.4871 | 0.0176 |
| MM dimension | 1.04 | **5.22** |
| Mean chain length | 3.46 | 1.59 |
| Mean antichain size | 18.8 | 215.0 |
| Mean interval size | 45.50 | 1.90 |

The causal set correctly recovers d ≈ 1 (reduced from the expected 2 by finite-region edge effects). The prime poset yields d ≈ 5.22.

**Interpretation.** The effective dimension of the prime divisibility poset — computed with no physical input — independently recovers approximately five dimensions. This aligns with the five-dimensional metric tensor constructed in Section 4 on entirely separate grounds. The sparse ordering (only 1.76% of pairs are related), large antichains (mean size 215 out of 499), and small intervals (mean 1.9 elements between related pairs) are characteristic of a high-dimensional causal structure where "light cones" are narrow and most elements are spacelike-separated.

**Scaling behavior.** The ordering fraction does not stabilize with N:

| N | Ordering fraction |
|---|---|
| 20 | 0.1579 |
| 100 | 0.0583 |
| 500 | 0.0176 |
| 1,000 | 0.0102 |
| 2,000 | 0.0058 |
| 5,000 | 0.0027 |

A causal set sprinkled into a fixed-dimension space would show a stable ordering fraction as N grows. The prime poset's fraction continues to decrease, implying an effective dimension that grows without bound. This distinguishes the prime substrate from any fixed-dimensional spacetime and supports the interpretation that Axiom 5 applies only on the emergence side of the Euler product — where continuous geometry with a stable dimension arises from the prime substrate. The substrate itself has no fixed dimension because it exists beneath the emergence boundary.

---

## 4. The Euler Product as a Metric Component

### 4.1 The formula

Euler's product formula (1737) states:

> ζ(s) = Π_p  1/(1 − p^(−s)) = Σ_{n=1}^{∞}  1/n^s

The left side of the product is a function over discrete atoms (primes). The right side of the sum is a continuous, analytic function over all positive integers. The equality states that the discrete substrate generates the continuous structure.

This single equation encodes all five Causal Set axioms:

1. **Discrete**: The product is indexed over discrete elements (primes).
2. **Partial order**: The product implicitly uses the divisibility ordering (fundamental theorem of arithmetic).
3. **Locally finite**: Any bounded range of the product contains finitely many terms.
4. **Countable**: The index set is countably infinite.
5. **Order = geometry**: The ordering relation (divisibility) produces the zeta function — a continuous geometric object encoding the distribution and density of primes.

### 4.2 The metric component

The previous Causal Set metric component was:

> g_δ̂ = log₁₀(1 + 1/δ_B)

where δ_B(r) = 0.003389|ln(r/r_s)| + 0.002508 was fitted to computational data.

The replacement is:

> g_ζ̂ = ln ζ(s) = −Σ_p  ln(1 − p^(−s))

Note the structural parallel: the CS metric used log(1 + 1/x), and the Euler product logged uses −log(1 − p^(−s)). Both are logarithms of ratios near unity. The CS metric derived continuous geometry from a fitted deviation parameter. The Euler metric derives continuous geometry from the prime numbers themselves.

### 4.3 The five-dimensional tensor

The reformulated metric tensor is:

```
g_μν^(5) = diag( Ṡ(r),  g_rr,  r²,  r²sin²θ,  ln ζ(s) )
```

Coordinates: {Ṡ, r, θ, φ, ζ̂}

Where:
- **Ṡ(r) = β = √(r_s/r)**: entropy rate (replaces time dimension). Derived from the Painlevé-Gullstrand cross-term — the velocity at which space flows inward.
- **g_rr = 1 + r_s/r**: radial component (PG form with river velocity absorbed). Clamped by the Benford floor inside the horizon.
- **r², r²sin²θ**: standard spherical angular components.
- **ln ζ(s)**: Euler prime substrate (5th dimension).

The Benford floor constraint acts on the spatial determinant:

> det(g_spatial) = g_rr × r⁴ sin²θ × ln ζ(s) ≥ 0.4068

This is not a component of the tensor. It is a constraint the tensor must satisfy — derived from the L2 norm of Benford's digit probability vector.

---

## 5. GPS Validation

### 5.1 The test

GPS satellites orbit at 26,571 km from Earth's center (altitude ~20,200 km). Their atomic clocks are corrected daily for a measured time dilation of approximately 38.6 μs/day relative to ground clocks. This correction has two components:

- Gravitational (general relativistic): satellite clocks run faster at higher altitude. Contributes +45.7 μs/day.
- Velocity (special relativistic): satellite motion slows clocks. Contributes −7.2 μs/day.
- Net: +38.5 μs/day.

This is the most precisely measured time dilation effect in continuous operation. We compute it three ways.

### 5.2 Method 1: Standard general relativity

Clock rate at radius r: √(1 − r_s/r), where r_s = 2GM/c² = 8.87 mm for Earth.

Gravitational dilation: (r_s/2)(1/R_earth − 1/R_gps) = 5.284 × 10⁻¹⁰

Velocity correction: −v²/(2c²) = −8.346 × 10⁻¹¹

Total: **+38.45 μs/day**

### 5.3 Method 2: Entropy rate (no time dimension)

The river velocity β = √(r_s/r) is derived entirely from spatial geometry — it is the Painlevé-Gullstrand cross-term absorbed into the effective radial metric component. No time dimension is used.

Clock rate at radius r: √(1 − β²) = √(1 − r_s/r)

This is algebraically identical to Method 1. The mechanism differs: Method 1 says "time runs slower near mass." Method 2 says "the rate of geometric change is lower near mass — less spatial activity means slower clock processing."

β at Earth surface: 3.729 × 10⁻⁵
β at GPS orbit: 1.827 × 10⁻⁵

Total: **+38.45 μs/day**

### 5.4 Method 3: Full Euler prime tensor

The Euler prime substrate contributes ln ζ(s) to the fifth dimension. At GPS distances:

- r/r_s at Earth surface: 7.19 × 10⁸
- r/r_s at GPS orbit: 3.00 × 10⁹
- det(spatial) at Earth: 3.14 × 10³⁵ (floor = 0.4068)

The floor is not active. The prime substrate contributes a nearly identical value at both radii. It factors out of the dilation ratio.

Total: **+38.45 μs/day**

### 5.5 Result

All three methods agree to within 2 × 10⁻⁶ μs/day (floating point precision). The reformulated tensor reproduces GPS time dilation without a time dimension and with the Euler prime substrate in the fifth dimension.

The prime substrate is completely inert in the weak-field regime. The spatial determinant exceeds the Benford floor by a factor of 10³⁵ at Earth's surface. The substrate activates only when det(spatial) approaches the floor — deep inside black holes at r/r_s < ~0.65.

---

## 6. The Three-Layer Emergence Hierarchy

The combined results — anti-Benford primes, Benford-conformant physics, and the Euler product bridging them — suggest a three-layer hierarchy:

### Layer 1: Below emergence (anti-Benford)

Prime numbers and Causal Set elements occupy this layer. Both are discrete, irreducible, and carry no Benford signature. The Benford diagnostic lens cannot resolve them because there is no multiplicative structure to detect. This layer is the substrate — the scaffolding on which emergent structure is built.

The substrate must be anti-Benford. If it were Benford-like, there would be no gradient between substrate and emergent structure — no "voltage difference" to drive emergence. The anti-Benford character of the substrate is not an absence of structure. It is the necessary precondition for structure to arise.

### Layer 2: At emergence (the boundary)

The Euler product formula sits at this boundary. Its right side is the discrete, anti-Benford substrate (a product over primes). Its left side is the continuous, Benford-like structure (the zeta function). The equals sign is the emergence event.

In the Causal Set analysis (Riner 2026a), one of five characteristic equations exhibited Benford-like behavior. That equation may correspond to this boundary — the single point where discrete structure transitions to continuous geometry.

### Layer 3: Above emergence (Benford-like)

Everything physical occupies this layer. Quantum spectra, black hole metrics, gravitational curvature, thermal radiation, financial data, population distributions — all carry the Benford signature. This is the domain where Benford's Law acts as a yardstick (Riner 2026, working notes): a universal invariant that all measurable systems conform to.

---

## 7. Discussion

### 7.1 What the Euler product provides

The substitution of the Euler product for the Causal Set dimension resolves two limitations of the previous framework:

1. **Axiomatic compression**: Five axioms are encoded in one formula. The metric tensor now contains a derived quantity (ln ζ(s)) rather than a fitted parameter (δ_B(r)).

2. **Fundamental derivation**: The CS equation δ_B(r) = 0.003389|ln(r/r_s)| + 0.002508 was fitted to computational data. The Euler product is a theorem of number theory — proven, exact, and independent of any dataset.

### 7.2 What remains to be determined

The parameter s in ζ(s) has not been mapped to a physical coordinate. In the GPS test, this mapping is unnecessary because the substrate is inert. Inside a black hole, where the substrate activates, the mapping of s to the radial coordinate (or to entropy rate, or to some other geometric quantity) becomes essential. This is the primary open question of the reformulation.

The Riemann Hypothesis — that all non-trivial zeros of ζ(s) lie on the critical line Re(s) = 1/2 — may have a physical interpretation in this framework. If s maps to a geometric parameter, the critical line may correspond to the emergence boundary between anti-Benford substrate and Benford-conformant physics. This is speculative and is noted only as a direction for future investigation.

### 7.3 Predictions

The reformulated tensor inherits all predictions from the previous framework:

1. **No singularity**: The Benford floor prevents det(spatial) from reaching zero. Geometry redistributes rather than collapsing to a point.
2. **Scale-free interior**: The floor activates at the same r/r_s (~0.65) for black holes of any mass.
3. **Binary merger signature**: Tidal compression during inspiral may briefly activate the floor between merging black holes (Riner 2026c, Phase 1b).
4. **GPS compatibility**: Demonstrated in Section 5.

New predictions specific to this reformulation:

5. **The substrate is number-theoretic**: The fifth dimension of the metric tensor is derived from prime numbers. If this is physical, the structure of spacetime at the Planck scale may be arithmetic rather than geometric.
6. **Anti-Benford as diagnostic**: Any physical system that exhibits anti-Benford statistics may be probing the substrate layer directly. This provides a new experimental criterion.

---

## 8. Conclusion

Prime numbers are the only dataset in this research program that exhibits zero Benford conformance across all integer bases tested. They satisfy all five axioms of Causal Set theory. Euler's product formula encodes these five axioms in a single equation that serves as a metric component in a five-dimensional tensor.

The reformulated tensor reproduces GPS time dilation without a time dimension, using entropy rate derived from spatial geometry. The Euler prime substrate is inert in weak fields and activates only inside black holes, where the Benford floor prevents singularity formation.

The result suggests a three-layer hierarchy: an anti-Benford substrate (primes, Causal Set elements) beneath an emergence boundary (the Euler product) beneath Benford-conformant physical reality. The equals sign in ζ(s) = Π_p 1/(1 − p^(−s)) may be the mathematical expression of that boundary.

---

## References

- Bombelli, L., Lee, J., Meyer, D., Sorkin, R.D. (1987). "Space-time as a causal set." *Physical Review Letters* 59(5), 521–524.
- Euler, L. (1737). "Variae observationes circa series infinitas." *Commentarii academiae scientiarum Petropolitanae* 9, 160–188.
- Hill, T.P. (1995). "A statistical derivation of the significant-digit law." *Statistical Science* 10(4), 354–363.
- Newcomb, S. (1881). "Note on the frequency of use of the different digits in natural numbers." *American Journal of Mathematics* 4(1), 39–40.
- Pinkham, R.S. (1961). "On the distribution of first significant digits." *Annals of Mathematical Statistics* 32(4), 1223–1230.
- Riner, C. (2026a). "Benford's Law as a Diagnostic for Quantum Gravity at the Planck Scale." Zenodo. DOI: 10.5281/zenodo.18553466.
- Riner, C. (2026b). "Benford's Law Inside a Black Hole." Draft, February 2026.
- Riner, C. (2026c). "Four-Dimensional Spatial Metric with Benford Floor and Causal Set Dimension." Zenodo. DOI: 10.5281/zenodo.18553466.
- Sorkin, R.D. (2003). "Causal sets: Discrete gravity." In *Lectures on Quantum Gravity*, ed. A. Gomberoff, D. Marolf. Springer.
- Verlinde, E. (2011). "On the Origin of Gravity and the Laws of Newton." *Journal of High Energy Physics* 2011(4), 29.

---

## Figures

**Figure 1**: Multi-base Benford heatmap. Columns: numbers 2–1000. Rows: bases 2–10 plus composite score. Gold cells indicate leading digits over-represented among primes; blue cells indicate under-representation. Top row shows ground truth (gold = prime). No correlation between prime location and composite score is visible. Tool: `benford_primes.html`.

**Figure 2**: Score distribution histogram. Normalized frequency of composite scores for primes (gold) and composites (blue). Complete overlap confirms zero discriminative power.

**Figure 3**: Side-by-side tensor comparison. Left: Einstein's 4×4 tensor with g_tt. Right: Riner's 5×5 tensor with entropy rate and Euler prime substrate. Bottom: clock rate vs. distance graph showing both methods produce identical curves. GPS prediction: +38.45 μs/day for both. File: `einstein_vs_riner.png`.

**Figure 4**: Five-dimensional metric tensor displayed as 5×5 matrix with component definitions. The Euler prime substrate ln ζ(s) occupies the (5,5) position. File: `benford_5d_euler_prime_tensor.png`.

---

## Data and Code Availability

All code, data, and interactive visualizations are available at:

- Heatmap tool: `benford_primes.html` (interactive, browser-based)
- GPS test: `scripts/gps_euler_prime_test.py`
- Tensor visualization: `scripts/make_tensor_matrix.py`
- Comparison figure: `scripts/einstein_vs_riner.py`
- Metric formula: `scripts/make_formula_euler_prime.py`
- Research log: `research_log.txt`
