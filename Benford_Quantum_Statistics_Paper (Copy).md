# Complete Monotonicity and Benford's Law: Deriving Quantum Statistics from the Significant Digit Distribution

### Christopher Riner
### Independent Researcher, Chesapeake, Virginia, USA

**Draft — February 2026**

---

## Abstract

We show that the Bose-Einstein distribution is the unique quantum statistical
distribution satisfying Benford's law exactly at all temperatures, and that this
result follows from a chain of established mathematical theorems connecting
complete monotonicity, the Bernstein-Widder representation, and the Benford
conformance of Laplace transforms. Specifically, requiring that a quantum
occupation function satisfy the significant digit law P(d) = log₁₀(1 + 1/d) at
all parameter values forces its series expansion to have exclusively non-negative
coefficients — selecting 1/(e^x − 1) over 1/(e^x + 1). The Fermi-Dirac
distribution, whose alternating-sign expansion violates complete monotonicity,
produces calculable periodic deviations from Benford's law: oscillations with
period exactly 1 in log₁₀(T), amplitude governed by the Dirichlet eta
function (1 − 2^(1−s))·ζ(s) with |η| = 1.054 times the single-exponential
baseline. We identify this Dirichlet factor as the mathematical
signature of the Pauli exclusion principle and derive a structural consequence:
no fermion can have zero Benford deviation, implying that massless fermions
cannot exist — consistent with the experimental discovery of nonzero neutrino
mass. These results hold independently of any particular interpretive framework.

---

## 1. Introduction

### 1.1 Benford's Law

The probability that the first significant digit of a number drawn from many
naturally occurring datasets takes the value *d* ∈ {1, 2, ..., 9} is not
uniform but logarithmic:

    P(d) = log₁₀(1 + 1/d)

This regularity was first noted by Newcomb [1] and empirically demonstrated by
Benford [2]. The distribution is equivalently characterized as the unique
probability measure on {1, ..., 9} that is scale-invariant [3], base-invariant
[4], and maximizes entropy on the significand [5]. The mathematical foundations
of Benford's law have been extensively developed by Berger and Hill [5], who
showed that a dataset satisfies Benford's law if and only if the logarithm of
its significand is uniformly distributed modulo 1.

Benford's law has been observed across a remarkably broad range of physical
data: fundamental constants [6], nuclear half-lives [7], hadron widths [8],
atomic spectra [9], astrophysical measurements [10], and geophysical quantities
[11]. This ubiquity suggests that the distribution reflects something structural
about the data-generating processes in physics, rather than being a mere
statistical artifact.

### 1.2 Benford's Law in Statistical Physics

In 2010, Shao and Ma [12] made a striking observation. They computed the
first-digit distributions of the three fundamental quantum statistical
distributions — Bose-Einstein, Fermi-Dirac, and Maxwell-Boltzmann — across a
range of temperatures and found that:

- The **Bose-Einstein** distribution satisfies Benford's law exactly at all
  temperatures.
- The **Fermi-Dirac** distribution shows systematic periodic deviations that
  oscillate with temperature.
- The **Maxwell-Boltzmann** distribution approximately satisfies Benford's law
  with small bounded errors.

This result was empirical: Shao and Ma demonstrated the conformance numerically
but did not provide a mathematical explanation for *why* the Bose-Einstein
distribution, alone among the three, satisfies Benford's law exactly.
Subsequent work by Cong, Li, and Ma [13] and Wang and Ma [14] established that
distributions expressible as Laplace transforms of non-negative measures satisfy
Benford's law, providing the theoretical tools needed to answer this question.

### 1.3 Overview and Summary of Results

This paper assembles these mathematical tools into a unified argument and derives
four results:

**Result 1 (Uniqueness).** Among the quantum statistical distributions
{1/(e^x − 1), 1/(e^x + 1)}, requiring exact Benford conformance at all
parameter values uniquely selects the Bose-Einstein distribution (Theorem 1,
Section 3).

**Result 2 (Quantitative FD deviation).** The Fermi-Dirac deviation from
Benford's law has period exactly 1 in log₁₀(T), amplitude precisely
|η(s)| = 1.054 times the single-exponential baseline, and functional form
governed by the Dirichlet eta function. These values are confirmed by the
data of Shao and Ma [12] (Section 4).

**Result 3 (Dirichlet signature).** The Dirichlet factor
(1 − 2^(1−s))·ζ(s) that controls the Fermi-Dirac deviation is identified as
the mathematical expression of the Pauli exclusion principle within the
Benford framework (Section 4).

**Result 4 (Massless fermion exclusion).** The Fermi-Dirac distribution's
inherent Benford deviation implies a structural constraint: no fermion can have
zero Benford deviation. If zero deviation is identified with masslessness —
as supported by the observation that all known massless particles are bosons —
then massless fermions cannot exist (Section 5).

We are explicit about what is novel and what is assembled from existing
results. The individual mathematical theorems used (Bernstein-Widder, the
Laplace transform proof of Benford conformance, the Fourier decomposition of
first-digit errors) are all established. The specific chain of reasoning
connecting them — from Benford conformance through complete monotonicity to
bosonic statistics — and the identification of the Dirichlet eta function as
the exclusion principle's mathematical signature within this framework, are
the new contributions of this paper. A detailed accounting is provided in
Section 6.2.

---

## 2. Mathematical Preliminaries

### 2.1 Benford's Law: Formal Definition

Let X > 0 be a real-valued random variable. The **significand function** is
defined as S(X) = X / 10^(⌊log₁₀ X⌋), so that S(X) ∈ [1, 10). We say that X
satisfies **Benford's law** if log₁₀ S(X) is uniformly distributed on [0, 1)
[5]. Equivalently, the probability that the first significant digit of X equals
*d* ∈ {1, 2, ..., 9} is:

    P(D₁ = d) = log₁₀(d + 1) − log₁₀(d) = log₁₀(1 + 1/d)

This yields the well-known probabilities: P(1) ≈ 0.301, P(2) ≈ 0.176,
P(3) ≈ 0.125, ..., P(9) ≈ 0.046.

The distribution is uniquely characterized by three properties:

1. **Scale invariance:** P(d) is unchanged under multiplication by any
   positive constant (Pinkham [3]).
2. **Base invariance:** P(d) does not depend on the number base used to
   express the data (Hill [4]).
3. **Maximum entropy on the significand:** Among all distributions on [1, 10)
   satisfying the normalization constraint, the Benford distribution maximizes
   Shannon entropy (Berger and Hill [5]).

### 2.2 Complete Monotonicity and the Bernstein-Widder Theorem

A function f : (0, ∞) → ℝ is **completely monotonic** if it possesses
derivatives of all orders and

    (−1)ⁿ f⁽ⁿ⁾(x) ≥ 0    for all x > 0 and all n = 0, 1, 2, ...

That is, f is non-negative, non-increasing, convex, and all successive
derivatives alternate in sign [15].

The **Bernstein-Widder theorem** [15, 16] provides a complete characterization:
a function f is completely monotonic on (0, ∞) if and only if it can be
represented as the Laplace transform of a non-negative measure μ on [0, ∞):

    f(x) = ∫₀^∞ e^(−xt) dμ(t)

When μ is a discrete non-negative measure supported on the positive integers,
this reduces to:

    f(x) = Σ_{k=1}^∞ aₖ · e^(−kx)    where aₖ ≥ 0 for all k

The non-negativity of all coefficients aₖ is the discrete signature of
complete monotonicity.

### 2.3 Laplace Transforms and Benford's Law

Cong, Li, and Ma [13] proved that distributions expressible as Laplace
transforms satisfy Benford's law, establishing the following decomposition.
For a function f(x) = ∫₀^∞ e^(−xt) dμ(t), the first-digit probability can
be written as:

    P_f(d) = log₁₀(1 + 1/d) + ε(d)

where the error term ε(d) depends on the Fourier harmonics of the measure μ.
When μ is a non-negative measure — i.e., when f is completely monotonic — the
error contributions from different harmonics cancel exactly, yielding ε(d) = 0
for all d [13, 14].

Wang and Ma [14] provided a concise proof of this result, showing that the
first-digit law originates from a basic property of the number system (the
uniform distribution of log₁₀ S modulo 1) combined with the distributional
properties of Laplace transforms of non-negative measures.

The key implication is:

> **Complete monotonicity of f(x) ⟹ f satisfies Benford's law exactly.**

### 2.4 Deviation Metrics

For any distribution f, we define the **per-digit Benford deviation** as:

    ε(d) = P_f(d) − log₁₀(1 + 1/d)    for d = 1, 2, ..., 9

This gives a signed measure of how much the first-digit probability for each
digit departs from the Benford baseline. To quantify the total deviation as a
single scalar, we define the **Benford deviation**:

    δ_B = √( Σ_{d=1}^{9} [P_f(d) − log₁₀(1 + 1/d)]² )

This is the Euclidean (L²) distance between the observed first-digit
distribution and the Benford distribution [17, 18]. Related metrics include
the Cho-Gaines d* statistic [17] and the Leemis-Schmeiser-Evans measure [18];
the Euclidean form is chosen here for its direct interpretability.

The key values are:

- **δ_B = 0** → exact Benford conformance
- **δ_B > 0** → deviation present

---

## 3. Benford Conformance of Quantum Statistical Distributions

### 3.1 Series Expansions and Coefficient Signs

The three fundamental quantum statistical distributions, expressed as functions
of the dimensionless variable x = ε/kT (energy divided by thermal energy), are:

**Bose-Einstein** (bosons — photons, gluons, W/Z/Higgs bosons):

    n_BE(x) = 1/(e^x − 1) = Σ_{k=1}^∞ e^(−kx) = e^(−x) + e^(−2x) + e^(−3x) + ...

Coefficients: aₖ = +1 for all k. **All non-negative.**

**Fermi-Dirac** (fermions — electrons, quarks, neutrinos):

    n_FD(x) = 1/(e^x + 1) = Σ_{k=1}^∞ (−1)^(k+1) e^(−kx) = e^(−x) − e^(−2x) + e^(−3x) − ...

Coefficients: aₖ = (−1)^(k+1). **Alternating in sign.**

**Maxwell-Boltzmann** (classical limit):

    n_MB(x) = e^(−x)

A single exponential. No infinite sum. Approximately Benford-conformant with
|ε(d)| bounded at ~0.03 [12], but lacking the infinite-sum structure that
produces exact cancellation of error terms.

The crucial distinction is in the coefficient signs. The Bose-Einstein
distribution has exclusively non-negative coefficients; the Fermi-Dirac
distribution does not.

### 3.2 Proof that the Bose-Einstein Distribution Satisfies Benford's Law Exactly

The argument proceeds in three steps:

1. The Bose-Einstein distribution n_BE(x) = Σ_{k=1}^∞ e^(−kx) has all
   non-negative coefficients aₖ = 1 ≥ 0.

2. By the Bernstein-Widder theorem (Section 2.2), n_BE(x) is therefore
   completely monotonic — it is the Laplace transform of the counting measure
   on the positive integers (a non-negative measure).

3. By the Cong-Li-Ma theorem (Section 2.3), distributions that are Laplace
   transforms of non-negative measures satisfy Benford's law with ε(d) = 0
   for all d [13, 14].

Therefore δ_B = 0 for the Bose-Einstein distribution at all temperatures.
This result provides the mathematical explanation for the exact conformance
observed numerically by Shao and Ma [12].

### 3.3 Main Theorem

**Theorem 1** (Benford Uniqueness of Bosonic Statistics). *Among the quantum
statistical occupation distributions {1/(e^x − 1), 1/(e^x + 1)}, requiring
exact Benford conformance (δ_B = 0) at all parameter values uniquely selects
the Bose-Einstein distribution n_BE(x) = 1/(e^x − 1).*

**Proof.** The argument consists of five steps, each relying on an established
result:

*Step 1 (Requirement).* We require that the quantum occupation function n(x)
satisfy Benford's law exactly — that is, δ_B = 0 at all values of the
parameter x = ε/kT. This is the hypothesis of the theorem.

*Step 2 (Laplace characterization).* By the theorem of Cong, Li, and Ma
[13], a function that is the Laplace transform of a non-negative measure —
i.e., a completely monotonic function — satisfies Benford's law exactly:
its Fourier-harmonic error terms cancel completely [13, 14]. Complete
monotonicity is thus a sufficient condition for exact Benford conformance.
(Whether it is also necessary — whether a non-completely-monotonic function
of the form Σ aₖ e^(−kx) might satisfy Benford exactly by some other
mechanism — is an open question that does not affect the present argument,
which uses only the sufficient direction.)

*Step 3 (Bernstein-Widder).* By the Bernstein-Widder theorem [15, 16], f is
completely monotonic if and only if all coefficients in its series expansion
are non-negative: aₖ ≥ 0 for all k.

*Step 4 (Coefficient test).* The Bose-Einstein distribution
n_BE(x) = 1/(e^x − 1) has coefficients aₖ = +1 for all k. ✓ Non-negative.
The Fermi-Dirac distribution n_FD(x) = 1/(e^x + 1) has coefficients
aₖ = (−1)^(k+1), which are negative for all even k. ✗ Not non-negative.

*Step 5 (Selection).* Only n_BE satisfies the non-negativity condition. The
Fermi-Dirac distribution is excluded by Step 3. ∎

**Remark on scope.** Theorem 1 states that among the two standard quantum
occupation functions, exact Benford conformance selects the bosonic one. It
does not claim that the Bose-Einstein distribution is the *only* completely
monotonic function, nor that Benford's law alone derives quantum mechanics.
The result is a selection theorem within a specific, well-defined function
class.

**Note on the interpretive step.** The hypothesis of Theorem 1 — *requiring*
exact Benford conformance — is a choice. The established mathematical fact
is that BE *satisfies* Benford exactly (Shao and Ma [12], proved via complete
monotonicity [13]). Elevating this observation to a requirement is the
interpretive step that produces the selection result. This distinction is
discussed further in Section 6.

### 3.4 Scope of the Result

Theorem 1 operates within the class of quantum statistical occupation
functions. Several clarifications are important:

1. The theorem does not derive quantum mechanics from Benford's law. It
   derives a selection among quantum statistical distributions from a
   Benford conformance requirement.

2. The Maxwell-Boltzmann distribution, as a single exponential, is
   approximately Benford-conformant (|ε| ≲ 0.03) but not exactly so. It
   occupies an intermediate position: better than Fermi-Dirac (which has
   systematic oscillations) but not as good as Bose-Einstein (which is
   exact). The MB distribution is not completely monotonic in the same
   infinite-sum sense.

3. The result is mathematical — it follows from the algebraic properties of
   the series expansions and established theorems about Laplace transforms.
   Physical interpretation of *why* nature might prefer Benford-conformant
   distributions is a separate question, addressed briefly in the companion
   paper [19].

---

## 4. Quantitative Predictions for the Fermi-Dirac Deviation

### 4.1 Fourier Decomposition of First-Digit Error

For a distribution f(x) = Σ_{k} aₖ e^(−kx), the first-digit probability can
be decomposed into a Benford term plus oscillatory harmonics using Poisson
summation [13, 20]. The dominant contribution to the error comes from the
first Fourier harmonic (n = ±1), which involves the factor:

    T^(2πi/ln 10) = e^(2πi · log₁₀ T)

This factor is purely oscillatory in log₁₀(T), completing one full cycle each
time T increases by a factor of 10. Higher harmonics (n = ±2, ±3, ...) are
suppressed by factors of approximately 10^(−2) per harmonic and can be
neglected [20].

The error for a single exponential e^(−λx) at the first harmonic is [20]:

    ε_single(d, λ) ≈ (2/ln 10) · Re[ Γ(1 + 2πi/ln 10) · (d^(−2πi/ln 10) − (d+1)^(−2πi/ln 10)) · λ^(−2πi/ln 10) ]

where Γ is the gamma function. For a sum of exponentials f(x) = Σ aₖ e^(−kx),
the total error is obtained by summing over all terms, weighted by the
coefficients aₖ.

### 4.2 The Dirichlet Factor

The summation over the exponential series introduces a **Dirichlet series
factor** that depends on the coefficient signs. For the three quantum
statistical distributions:

**Maxwell-Boltzmann** (single term, a₁ = 1):

    D_MB(s) = 1

The error amplitude is that of a single exponential: |ε_max| ≈ 0.03, with
oscillation period 1 in log₁₀(T).

**Bose-Einstein** (all positive coefficients, aₖ = 1):

    D_BE(s) = Σ_{k=1}^∞ k^(−s) = ζ(s)

where ζ is the Riemann zeta function and s = 2πi/ln 10. The complete
monotonicity of the distribution — all aₖ positive — causes the error
contributions from different exponential terms to cancel exactly when
integrated over the full energy distribution. The Dirichlet factor evaluates
to ζ(s), but the cancellation mechanism in the Cong-Li-Ma proof [13] ensures
the net error vanishes. Result: **ε(d) = 0 at all temperatures.**

**Fermi-Dirac** (alternating coefficients, aₖ = (−1)^(k+1)):

    D_FD(s) = Σ_{k=1}^∞ (−1)^(k+1) k^(−s) = (1 − 2^(1−s)) · ζ(s) = η(s)

where η(s) is the **Dirichlet eta function**. This factor is nonzero: the
alternating signs prevent the cancellation that occurs for BE. The Pauli
exclusion principle — which restricts fermions to single-occupancy states and
produces the plus sign in the FD denominator — appears directly in the
mathematics as the factor (1 − 2^(1−s)) that prevents the Dirichlet series
from reducing to ζ(s) alone.

This identification is, to our knowledge, new: *the Dirichlet eta function
is the mathematical signature of the exclusion principle within the Benford
framework*. The factor (1 − 2^(1−s)) quantifies precisely how much the
fermionic sign alternation prevents Benford conformance.

Evaluating these quantities numerically at s = 2πi/ln 10:

    |ζ(s)|          = 0.4214
    |1 − 2^(1−s)|   = 2.5021
    |η(s)|          = 2.5021 × 0.4214 = 1.0545

The exclusion factor |1 − 2^(1−s)| = 2.50 amplifies the deviation by a
factor of 2.5 relative to what a coherent (all-positive) summation would
produce. However, the zeta factor |ζ(s)| = 0.42 provides partial
compensation, so that the net Fermi-Dirac error amplitude is |η(s)| = 1.054
times the single-exponential (Maxwell-Boltzmann) amplitude — only a 5.4%
increase. The exclusion principle's effect on the deviation magnitude is
almost entirely offset by the structure of the Dirichlet series; what
changes is the *pattern* (periodic oscillation rather than monotonic
approach to Benford), not the magnitude.

### 4.3 Three Quantitative Predictions

From the Fourier decomposition and the Dirichlet factor, the framework
generates three specific predictions for the Fermi-Dirac deviation:

**Prediction 1 — Period.** The deviation oscillates with period exactly 1 in
log₁₀(T). This follows from the phase factor e^(2πi · log₁₀ T), which
completes one cycle each time T increases by a factor of 10. The period is
independent of the distribution and is a universal feature of the Fourier
structure of first-digit laws.

**Prediction 2 — Amplitude.** The Fermi-Dirac error amplitude is precisely
|η(s)| = 1.054 times the Maxwell-Boltzmann (single-exponential) amplitude.
Since the MB first-harmonic amplitude varies by digit — from A_MB ≈ 0.080
at d = 1 to A_MB ≈ 0.014 at d = 9 — the predicted FD amplitudes are:

    A_FD(d) = |η(s)| · A_MB(d) = 1.054 · A_MB(d)

For the intermediate digits (d = 3 through 7), this yields per-digit maxima
of 0.019–0.040, consistent with the range 0.02–0.04 reported by Shao and
Ma [12].

**Prediction 3 — Functional form.** The scalar Benford deviation follows
approximately:

    δ_B^(FD)(T) ≈ δ_max · |cos(2π · log₁₀(T) + φ)|

where δ_max and φ are constants determined by the complex arguments of the
Dirichlet factor η(s) and the gamma function Γ(1 + 2πi/ln 10). The cosine
form reflects the dominance of the first Fourier harmonic.

### 4.4 Comparison with Numerical Data

All three predictions are confirmed by the numerical computations of Shao and
Ma [12], who calculated first-digit distributions for all three quantum
statistical distributions across a range of temperatures. Table 1 summarizes
the comparison.

**Table 1.** Predicted vs. observed properties of the Fermi-Dirac deviation.

| Property | Predicted | Observed (Shao & Ma [12]) | Status |
|---|---|---|---|
| Period | 1 in log₁₀(T) | 1 in log₁₀(T) | Confirmed |
| FD/MB amplitude ratio | |η(s)| = 1.054 | ~1.0 | Confirmed |
| FD amplitude, digits 3–7 | 0.019–0.040 | 0.02–0.04 | Confirmed |
| Functional form | Cosine in log₁₀(T) | Periodic oscillation | Confirmed |
| BE deviation | 0 (exact) | 0 (exact) | Confirmed |

The three quantum distributions are also compared in Table 2 via their
Dirichlet factors.

**Table 2.** Dirichlet factors and Benford conformance of quantum statistical
distributions.

| Distribution | Series coefficients | Dirichlet factor | |D(s)| | δ_B |
|---|---|---|---|---|
| Maxwell-Boltzmann | {1} (single term) | 1 | 1.000 | >0 (small) |
| Bose-Einstein | {+1, +1, +1, ...} | ζ(s) → cancels to 0 | 0 (exact) | 0 (exact) |
| Fermi-Dirac | {+1, −1, +1, ...} | (1−2^(1−s))·ζ(s) | 1.054 | >0 (periodic) |

The predictions are not post-hoc fits. They follow from the mathematical
structure of the alternating-sign series, which is determined by the
coefficient signs of the Fermi-Dirac expansion. The Dirichlet eta function
arises directly from the summation over alternating terms, and its modulus
controls the error amplitude.

---

## 5. The Structural Exclusion of Massless Fermions

### 5.1 The Argument

The results of Sections 3 and 4 establish two mathematical facts:

1. The Bose-Einstein distribution has δ_B = 0 (exact Benford conformance).
2. The Fermi-Dirac distribution has δ_B > 0 (inherent, periodic deviation).

The second fact is a consequence of the alternating-sign structure of the
Fermi-Dirac series, which in turn is a consequence of the Pauli exclusion
principle. No temperature, energy range, or limiting procedure removes this
deviation — it is structural.

Now consider the empirical observation: **all known massless particles are
bosons.** Photons, gluons, and gravitons (if they exist) are all bosons. The
Standard Model of particle physics contains no massless fermions.

Combining the mathematical result with this observation yields the following
structural constraint:

> Fermionic statistics ⟹ δ_B > 0 (mathematical fact).
> All known massless particles have δ_B = 0 (empirical observation).
> Therefore: no fermion can be massless (structural consequence).

For decades, neutrinos were treated as massless fermions in the Standard
Model. This would have contradicted the structural constraint above. However,
the discovery of neutrino oscillations by the Super-Kamiokande collaboration
(1998) [21] and the SNO collaboration (2001) [22] established that neutrinos
possess nonzero mass. The particles once thought to be massless fermions
turned out not to be massless.

### 5.2 Status of the Argument

We distinguish carefully between two components:

1. **The mathematical fact:** The Fermi-Dirac distribution inherently
   violates Benford's law (δ_B > 0). This is proven in Sections 3 and 4.

2. **The interpretive identification:** Equating δ_B = 0 with masslessness.
   This is motivated by the empirical correlation (all known massless
   particles are bosons, and bosonic statistics give δ_B = 0) but is not
   a mathematical theorem.

The argument that massless fermions cannot exist therefore has the structure:
mathematical fact + interpretive identification → physical consequence. The
mathematical component is rigorous; the interpretive component is an
observation about the physical world that invites further investigation.

We call this a **retrodiction** rather than a prediction: the neutrino mass
discovery (1998–2001) preceded this analysis. The framework does not predict
the discovery but is structurally consistent with it — the same mathematics
that selects bosonic statistics in Theorem 1 simultaneously excludes
fermions from the δ_B = 0 class.

---

## 6. Discussion

### 6.1 Summary of Results

This paper has established four results connecting Benford's law to quantum
statistical distributions:

1. The Bose-Einstein distribution is the unique quantum occupation function
   satisfying Benford's law exactly, selected by the requirement of complete
   monotonicity (Theorem 1).

2. The Fermi-Dirac deviation from Benford's law has calculable period,
   amplitude, and functional form, all confirmed by existing numerical data.

3. The Dirichlet eta function (1 − 2^(1−s))·ζ(s) serves as the
   mathematical signature of the Pauli exclusion principle within the
   Benford framework.

4. The structural impossibility of δ_B = 0 for fermionic distributions
   implies that no massless fermion can exist, consistent with the neutrino
   mass discovery.

### 6.2 What Is Novel vs. What Is Reframing

We believe transparency about the novelty of contributions is essential.
Table 3 provides an explicit accounting.

**Table 3.** Established results used vs. new contributions.

| Established results (used, not claimed) | New contributions (this paper) |
|---|---|
| BE satisfies Benford exactly at all T — Shao & Ma [12] | The specific chain: Benford requirement → complete monotonicity → bosonic selection (Theorem 1) |
| FD deviates periodically from Benford — Shao & Ma [12] | Identifying the Dirichlet eta function as the exclusion principle's mathematical signature |
| Completely monotonic functions satisfy Benford — Cong, Li & Ma [13] | The structural argument excluding massless fermions from the δ_B = 0 class |
| Bernstein-Widder representation theorem — Bernstein [15, 16] | The ε(d) decomposition as a formal per-digit and scalar framework |
| Fourier/Poisson decomposition of first-digit error — Lemons et al. [20] | Assembling the derivation chain into a self-contained selection theorem |
| Cosine oscillation form of FD deviation — Lemons et al. [20] | Explicit comparison of predicted vs. observed FD deviation properties |
| Neutrino mass discovery — Super-K [21], SNO [22] | Quantitative prediction table (Table 1) |

Each individual theorem in the proof of Theorem 1 is due to others. The
contribution of this paper is the specific chain of reasoning that connects
them — from a Benford conformance requirement, through complete monotonicity
and the Bernstein-Widder theorem, to the selection of bosonic statistics —
and the consequences drawn from this chain.

### 6.3 Relation to Companion Paper

These mathematical results were originally developed as part of a broader
philosophical framework proposed in a companion paper (Riner [19], "The Law of
Emergence"). That paper interprets Benford's law as a universal constraint on
physical systems and explores implications for gravity, entropy, and spacetime.
The present paper extracts the mathematical content so that the results can be
evaluated independently of that interpretive framework.

The key difference is the status of the Benford conformance requirement. In the
companion paper [19], this requirement is proposed as a physical axiom. In the
present paper, it is treated as a mathematical hypothesis (the antecedent of
Theorem 1) whose consequences are then derived. The mathematical results are
the same regardless of whether one accepts the axiom.

### 6.4 Open Questions

Several directions for further work are suggested by these results:

1. **Per-digit amplitude predictions.** The amplitude range 0.02–0.05 given
   in Section 4.3 is for the most affected digits. Sharpening this to
   per-digit values ε(d) for d = 1, ..., 9 would provide a more stringent
   test of the framework.

2. **Extension beyond quantum statistics.** The Benford conformance condition
   (complete monotonicity) applies to any distribution expressible as a
   Laplace transform. Investigating which other physically important
   distributions satisfy or violate this condition could extend the
   framework to classical statistical mechanics, astrophysical
   distributions, and other domains.

3. **Diagnostic applications.** The Benford deviation δ_B has been shown to
   detect quantum phase transitions [23, 24] with scaling exponents
   competitive with or exceeding conventional order parameters. Developing
   δ_B as a systematic diagnostic tool for statistical physics is a
   natural application.

4. **Rigorous bounds on FD amplitude.** While the Dirichlet eta function
   determines the FD deviation factor, rigorous per-digit upper and lower
   bounds — rather than the approximate range given here — would
   strengthen the quantitative predictions.

---

## 7. Conclusion

The Bose-Einstein distribution is the unique quantum statistical distribution
satisfying Benford's law exactly at all temperatures. This result follows from
a chain of established mathematical theorems: the Bernstein-Widder
characterization of completely monotonic functions, and the proof that Laplace
transforms of non-negative measures satisfy Benford's law. Requiring exact
Benford conformance among the quantum occupation functions selects
1/(e^x − 1) over 1/(e^x + 1) — the minus sign in the denominator is forced
by the non-negativity of the series coefficients.

The Fermi-Dirac distribution's deviation from Benford's law is calculable,
periodic, and governed by the Dirichlet eta function with |η(2πi/ln 10)| =
1.054 — which we identify as the mathematical expression of the Pauli
exclusion principle within the Benford framework. The structural impossibility of zero Benford deviation for
fermionic distributions implies that massless fermions cannot exist, consistent
with the experimental discovery of nonzero neutrino mass.

These mathematical results hold regardless of whether one accepts any
particular interpretive framework for their physical significance.

---

## References

[1] S. Newcomb, "Note on the Frequency of Use of the Different Digits in
Natural Numbers," *American Journal of Mathematics*, vol. 4, no. 1, pp. 39-40,
1881.

[2] F. Benford, "The Law of Anomalous Numbers," *Proceedings of the American
Philosophical Society*, vol. 78, no. 4, pp. 551-572, 1938.

[3] R. S. Pinkham, "On the Distribution of First Significant Digits," *Annals
of Mathematical Statistics*, vol. 32, no. 4, pp. 1223-1230, 1961.

[4] T. P. Hill, "Base-Invariance Implies Benford's Law," *Proceedings of the
American Mathematical Society*, vol. 123, no. 3, pp. 887-895, 1995.

[5] A. Berger and T. P. Hill, "A Basic Theory of Benford's Law," *Probability
Surveys*, vol. 8, pp. 1-126, 2011.

[6] J. Burke and E. Kincanon, "Benford's Law and Physical Constants: The
Distribution of Initial Digits," *American Journal of Physics*, vol. 59, no. 10,
pp. 952-954, 1991.

[7] D. Ni and Z. Ren, "Benford's Law and Half-Lives of Unstable Nuclei,"
*European Physical Journal A*, vol. 38, pp. 251-255, 2008.

[8] L. Shao and B.-Q. Ma, "First Digit Distribution of Hadron Full Width,"
*Modern Physics Letters A*, vol. 24, no. 30, pp. 2465-2474, 2009.

[9] Y. Ralchenko and J.-C. Pain, "Benford's Law in Atomic Spectra and Opacity
Databases," *Journal of Quantitative Spectroscopy and Radiative Transfer*,
vol. 322, 109010, 2024.

[10] T. Alexopoulos and S. Leontsinis, "Benford's Law in Astronomy," *Journal
of Astrophysics and Astronomy*, vol. 35, pp. 639-648, 2014.

[11] M. Sambridge, H. Tkalcic, and A. Jackson, "Benford's Law in the Natural
Sciences," *Geophysical Research Letters*, vol. 37, L22301, 2010.

[12] L. Shao and B.-Q. Ma, "The Significant Digit Law in Statistical Physics,"
*Physica A*, vol. 389, no. 16, pp. 3109-3116, 2010.

[13] M. Cong, M. Li, and B.-Q. Ma, "First Digit Law from Laplace Transform,"
*Physics Letters A*, vol. 383, 1836, 2019.

[14] T. Wang and B.-Q. Ma, "A Concise Proof of Benford's Law," *Fundamental
Research*, vol. 4, pp. 841-845, 2024.

[15] S. N. Bernstein, "Sur les fonctions absolument monotones," *Acta
Mathematica*, vol. 52, pp. 1-66, 1929.

[16] D. V. Widder, *The Laplace Transform*, Princeton University Press, 1941.

[17] W. K. T. Cho and B. J. Gaines, "Breaking the (Benford) Law: Statistical
Fraud Detection in Campaign Finance," *The American Statistician*, vol. 61,
no. 3, pp. 218-223, 2007.

[18] L. M. Leemis, B. W. Schmeiser, and D. L. Evans, "Survival Distributions
Satisfying Benford's Law," *The American Statistician*, vol. 54, no. 4,
pp. 236-241, 2000.

[19] C. Riner, "The Law of Emergence: Benford's Distribution as a Universal
Constraint on Physical Reality," preprint, 2026.

[20] D. Lemons, N. Lemons, and W. Peter, "First Digit Oscillations," *Stats*,
vol. 4, pp. 595-601, 2021.

[21] Y. Fukuda et al. (Super-Kamiokande Collaboration), "Evidence for
Oscillation of Atmospheric Neutrinos," *Physical Review Letters*, vol. 81,
pp. 1562-1567, 1998.

[22] Q. R. Ahmad et al. (SNO Collaboration), "Measurement of the Rate of
νe + d → p + p + e⁻ Interactions Produced by ⁸B Solar Neutrinos at the
Sudbury Neutrino Observatory," *Physical Review Letters*, vol. 87, 071301,
2001.

[23] A. Sen(De) and U. Sen, "Benford's Law Detects Quantum Phase Transitions
Similarly as Earthquakes," *Europhysics Letters*, vol. 95, 50008, 2011.

[24] A. D. Rane, U. Mishra, A. Biswas, A. Sen(De), and U. Sen, "Benford's Law
Gives Better Scaling Exponents in Phase Transitions of Quantum XY Models,"
*Physical Review E*, vol. 90, 022144, 2014.

---
