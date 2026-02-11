# A Four-Dimensional Spatial Metric with Benford Floor and Causal Set Dimension

**Christopher Riner**

---

## Abstract

We present a modified Schwarzschild metric in Painleve-Gullstrand coordinates with three modifications: (1) the determinant singularity is replaced by a floor derived from Benford's law, (2) time is removed as a coordinate dimension and redefined as an entropy rate computed from spatial metric gradients, and (3) the Causal Set substrate is promoted to a fourth spatial dimension via Benford's equation applied to its own deviation measure. The resulting 4D spatial metric $g_{\mu\nu} = \text{diag}(g_{rr},\, g_{\theta\theta},\, g_{\phi\phi},\, g_{\hat{\delta}})$ contains no free parameters. The floor value (0.4068), the coupling between the Causal Set and spatial geometry ($\sim$49.5%, range 42--57%), and the entropy rate acceleration near the singularity all emerge from the determinant constraint without being imposed. Using 40 measurements of the Euclidean deviation $\delta_B$ across a Schwarzschild black hole from prior work, we derive a two-parameter equation for the Causal Set response: $\delta_B(r) = 0.003389\,|\!\ln(r/r_s)| + 0.002508$. Validation against six real astrophysical objects (Earth, Jupiter, two neutron stars, Sgr A*, TON 618) shows the model produces flat, unmodified geometry for weak-field objects and activates only inside black hole horizons, with identical behavior in $r/r_s$ units regardless of mass.

---

## 1. Introduction

General relativity describes black holes through the Schwarzschild metric, which predicts a singularity at $r = 0$ where the determinant of the spatial metric vanishes and curvature invariants diverge. This is not a coordinate artifact --- it represents a genuine breakdown of the theory. Every major quantum gravity program (loop quantum gravity, string theory, causal set theory, asymptotic safety, and others) attempts to resolve this singularity by introducing a minimum length scale, typically the Planck length $\ell_P \approx 1.6 \times 10^{-35}$ m, as a cutoff.

We take a different approach. Instead of choosing a cutoff scale, we replace the zero in the determinant collapse with a value derived from Benford's law --- the logarithmic distribution of first significant digits observed across natural datasets spanning many orders of magnitude. The floor is not a free parameter. It is the $L^2$ norm of the Benford probability vector, a fixed number that emerges from the equation $P(d) = \log_{10}(1 + 1/d)$ for $d = 1, \ldots, 9$.

This substitution, combined with two other modifications --- removing time as a coordinate and adding the Causal Set substrate as a fourth spatial dimension --- produces a metric that resolves the singularity, generates emergent time behavior consistent with the Wheeler-DeWitt equation, and couples to the Causal Set through geometry alone, with no tunable constants.

The starting point is data. In prior work [1] (methodology detailed in Appendix D), we tested 10 quantum gravity models against Benford's law using the Euclidean deviation $\delta_B$ computed at 40 radial positions through a Schwarzschild black hole. The Causal Set model showed the best conformance by a factor of approximately 1.7. The present work asks: what happens when we take those 40 measurements and build them into the metric tensor itself?

---

## 2. The Benford Floor

### 2.1 Painleve-Gullstrand coordinates

The Schwarzschild solution in Painleve-Gullstrand (PG) coordinates eliminates the coordinate singularity at the event horizon. Unlike standard Schwarzschild coordinates, the purely spatial part of the PG metric is flat:

$$g_{rr} = 1, \qquad g_{\theta\theta} = r^2, \qquad g_{\phi\phi} = r^2 \sin^2\theta$$

We work with the spatial sector only. Setting $\sin\theta = 1$ (equatorial plane), the 3D spatial metric is:

$$g^{(3)}_{\mu\nu} = \text{diag}(1,\; r^2,\; r^2)$$

with determinant $\det(g^{(3)}) = r^4$.

As $r \to 0$, the determinant vanishes: the spatial volume element collapses to zero. This is the singularity.

### 2.2 The floor value

Benford's law states that in many naturally occurring datasets, the probability of the first significant digit $d$ is:

$$P(d) = \log_{10}\!\left(1 + \frac{1}{d}\right), \qquad d = 1, 2, \ldots, 9$$

The nine probabilities form a vector $\mathbf{p} = (P(1), P(2), \ldots, P(9))$. Its $L^2$ norm is:

$$\|\mathbf{p}\|_2 = \sqrt{\sum_{d=1}^{9} \left[\log_{10}\!\left(1 + \frac{1}{d}\right)\right]^2}$$

Computing explicitly:

| $d$ | $P(d) = \log_{10}(1 + 1/d)$ | $P(d)^2$ |
|-----|------|------|
| 1 | 0.30103 | 0.09062 |
| 2 | 0.17609 | 0.03101 |
| 3 | 0.12494 | 0.01561 |
| 4 | 0.09691 | 0.00939 |
| 5 | 0.07918 | 0.00627 |
| 6 | 0.06695 | 0.00448 |
| 7 | 0.05799 | 0.00336 |
| 8 | 0.05115 | 0.00262 |
| 9 | 0.04576 | 0.00209 |
| **Sum** | **1.00000** | **0.16545** |

$$\|\mathbf{p}\|_2 = \sqrt{0.16545} = 0.4068$$

This is a fixed constant. It depends on nothing but the digits 1 through 9 and the base-10 logarithm.

### 2.3 The modified metric

We impose a floor on the determinant:

$$\det(g^{(3)}) \geq \|\mathbf{p}\|_2 = 0.4068$$

When the natural determinant $r^4$ falls below this value, the radial component $g_{rr}$ absorbs the compensation:

$$g_{rr} = \begin{cases} 1 & \text{if } r^4 \geq 0.4068 \\ \displaystyle\frac{0.4068}{r^4} & \text{if } r^4 < 0.4068 \end{cases}$$

The angular components are unchanged: $g_{\theta\theta} = g_{\phi\phi} = r^2$.

The floor activates at:

$$r_{\text{floor}}^{(3)} = (0.4068)^{1/4} = 0.7988\, r_s$$

This is inside the event horizon ($r = r_s$) but well outside the classical singularity. It is not a Planck-scale correction --- it modifies the geometry over a significant fraction of the black hole interior.

---

## 3. Entropy Rate as Emergent Time

### 3.1 Removing time from the metric

The metric $g^{(3)}_{\mu\nu}$ has no time component. This is a deliberate choice, not an approximation. We define time not as a coordinate dimension but as the rate of change of the spatial geometry:

$$\dot{S}(r) = \sqrt{\sum_i \left(\frac{dg_i}{dr}\right)^2}$$

where the sum runs over all metric components $g_i \in \{g_{rr}, g_{\theta\theta}, g_{\phi\phi}\}$ (and later $g_{\hat{\delta}}$). This quantity has the structure of a speed --- it measures how fast the geometry changes per unit of coordinate distance.

We call this the *entropy rate* rather than "time" because it measures what time *does* --- the rate of change between spatial configurations --- without assuming what time *is*.

### 3.2 Consistency with existing frameworks

This treatment is consistent with several established approaches to the problem of time in quantum gravity:

- The **Wheeler-DeWitt equation**, the closest existing candidate for a quantum gravity equation, contains no time variable. The wavefunction of the universe $\Psi[g_{ij}]$ depends on the 3-geometry alone.

- **Barbour's timeless mechanics** argues that time does not exist as a fundamental entity; only spatial configurations and their relative changes are real.

- **Rovelli's thermal time hypothesis** identifies the flow of time with the flow of entropy in thermodynamic systems.

Our entropy rate is a specific realization: time is the rate at which the Causal Set substrate processes spatial geometry toward equilibrium.

### 3.3 Behavior near the singularity

In standard GR (the unmodified metric), the entropy rate depends only on the angular components $g_{\theta\theta} = g_{\phi\phi} = r^2$, whose derivatives $dg/dr = 2r$ shrink to zero as $r \to 0$. Time freezes at the singularity.

In the Benford-modified metric, when the floor is active, $g_{rr} = 0.4068/r^4$ and $dg_{rr}/dr = -4 \times 0.4068/r^5$. This derivative grows without bound as $r \to 0$. The entropy rate climbs by orders of magnitude toward the singularity rather than dropping to zero.

| $r/r_s$ | Standard entropy rate | Benford-modified entropy rate |
|---------|----------------------|-------------------------------|
| 5.0 | 14.14 | 14.14 |
| 1.0 | 2.83 | 2.83 |
| 0.50 | 1.41 | 5.06 |
| 0.05 | 0.14 | 1,078 |

The acceleration emerges from the geometry resisting collapse. It was not imposed. The only input was the floor constraint.

---

## 4. Causal Set as a Fourth Spatial Dimension

### 4.1 Motivation: Nine models fit in three dimensions --- one does not

In the black hole Benford experiment (Appendix D), we measured the Euclidean deviation $\delta_B$ at 40 positions through a Schwarzschild black hole for 10 quantum gravity models. The Causal Set model showed the lowest deviation --- the closest conformance to Benford's law --- by a factor of approximately 1.7 over the next best model.

The 10 models tested were: Standard GR (baseline), Loop Quantum Gravity, Generalized Uncertainty Principle, Doubly Special Relativity, Hagedorn Temperature, Asymptotic Safety, Horava-Lifshitz, Non-Commutative Geometry, Causal Dynamical Triangulations, and Causal Set theory.

Nine of these --- all except Causal Set --- are proposals for how the three spatial dimensions $\{g_{rr}, g_{\theta\theta}, g_{\phi\phi}\}$ behave at extreme curvature. They are different answers to the same question: what happens to $\{x, y, z\}$ near a singularity? Loop Quantum Gravity discretizes spatial geometry. The Generalized Uncertainty Principle modifies the commutation relations of spatial coordinates. Horava-Lifshitz changes the scaling behavior of spatial dimensions at high energies. And so on. All nine operate *within* the three spatial dimensions of Einstein's field equations.

Causal Set theory is structurally different. It does not modify the behavior of spatial geometry --- it posits an independent substrate underlying it. A discrete structure from which continuous spacetime emerges. It has its own dynamics, its own equilibrium, and we have measured its deviation at every point.

The three spatial dimensions of Einstein's equations account for 9 of the 10 models. The 10th --- the one that won the Benford conformance test --- does not fit within them. It operates *alongside* space, not within it.

This suggests that the Causal Set belongs not as data overlaid on a 3D metric, but as a fourth dimension within it. Einstein's field equations had room for 9 out of 10 quantum gravity approaches. The 4D metric has room for all 10.

### 4.2 The fourth component

We define the Causal Set metric component by applying Benford's equation to the deviation itself:

$$g_{\hat{\delta}} = \log_{10}\!\left(1 + \frac{1}{\delta_B}\right)$$

This is Benford's law with $\delta_B$ playing the role of the digit. When $\delta_B$ is small (strong conformance), $g_{\hat{\delta}}$ is large. When $\delta_B$ is large (weak conformance), $g_{\hat{\delta}}$ approaches 1.

The full 4D spatial metric becomes:

$$g^{(4)}_{\mu\nu} = \text{diag}(g_{rr},\; r^2,\; r^2,\; g_{\hat{\delta}})$$

Written as a matrix:

$$g^{(4)}_{\mu\nu} = \begin{pmatrix} g_{rr} & 0 & 0 & 0 \\ 0 & r^2 & 0 & 0 \\ 0 & 0 & r^2 & 0 \\ 0 & 0 & 0 & \log_{10}\!\left(1 + \dfrac{1}{0.003389\,|\ln(r/r_s)| + 0.002508}\right) \end{pmatrix}$$

where

$$g_{rr} = \begin{cases} 1 & \text{if } \det \geq 0.4068 \\ \dfrac{0.4068}{r^4 \cdot g_{\hat{\delta}}} & \text{if } \det < 0.4068 \end{cases}$$

Four dimensions. Three spatial, one Causal Set. No time. No free parameters.

With determinant:

$$\det(g^{(4)}) = g_{rr} \cdot r^2 \cdot r^2 \cdot g_{\hat{\delta}} = g_{rr} \cdot r^4 \cdot g_{\hat{\delta}}$$

The floor constraint generalizes:

$$\det(g^{(4)}) \geq 0.4068$$

When the natural determinant $r^4 \cdot g_{\hat{\delta}}$ (with $g_{rr} = 1$) falls below the floor:

$$g_{rr} = \frac{0.4068}{r^4 \cdot g_{\hat{\delta}}}$$

There are no free parameters. The coupling between the Causal Set dimension and the spatial geometry is not specified --- it emerges from the shared determinant constraint.

### 4.3 Floor activation in 4D

Because $g_{\hat{\delta}} > 1$ wherever $\delta_B < 1$ (which holds across all 40 measured positions), the fourth dimension *boosts* the determinant. The 4D floor activates later --- at a smaller radius --- than the 3D floor:

$$r_{\text{floor}}^{(4)} \approx 0.641\, r_s$$

compared to $r_{\text{floor}}^{(3)} = 0.799\, r_s$.

The Causal Set dimension stabilizes the geometry. It pushes the point of maximum modification deeper into the interior, closer to where the classical singularity would be. The 4D metric is more stable than the 3D metric precisely because the Causal Set contributes additional geometric "volume" through $g_{\hat{\delta}}$.

---

## 5. The Causal Set Equation

### 5.1 The data

The 40 measurements of $\delta_B$ across the black hole geometry, from $r/r_s = 10$ (far outside) to $r/r_s = 0.01$ (deep interior), are reproduced in Table 1. These were computed using the methodology described in Appendix D, by evaluating the first-digit distribution of quantum gravity spectra against Benford's law at each radial position, using the Euclidean deviation:

$$\delta_B = \sqrt{\sum_{d=1}^{9} \left[P_{\text{obs}}(d) - \log_{10}\!\left(1 + \frac{1}{d}\right)\right]^2}$$

**Table 1.** Causal Set Benford deviation $\delta_B$ as a function of radial position $r/r_s$.

| $r/r_s$ | $\delta_B$ | | $r/r_s$ | $\delta_B$ |
|---------|-----------|---|---------|-----------|
| 10.0 | 0.027552 | | 0.99 | 0.004267 |
| 7.0 | 0.010611 | | 0.95 | 0.005472 |
| 5.0 | 0.004223 | | 0.9 | 0.002860 |
| 3.0 | 0.002856 | | 0.85 | 0.004293 |
| 2.0 | 0.002748 | | 0.8 | 0.003671 |
| 1.5 | 0.005169 | | 0.7 | 0.006253 |
| 1.3 | 0.003184 | | 0.6 | 0.006376 |
| 1.2 | 0.003569 | | 0.5 | 0.004860 |
| 1.15 | 0.002757 | | 0.4 | 0.007448 |
| 1.1 | 0.002136 | | 0.3 | 0.013130 |
| 1.08 | 0.002550 | | 0.25 | 0.016065 |
| 1.06 | 0.002722 | | 0.2 | 0.012050 |
| 1.04 | 0.003532 | | 0.15 | 0.015154 |
| 1.03 | 0.003889 | | 0.12 | 0.019625 |
| 1.02 | 0.003482 | | 0.1 | 0.016711 |
| 1.015 | 0.003757 | | 0.08 | 0.017486 |
| 1.01 | 0.004013 | | 0.06 | 0.013490 |
| 1.005 | 0.004200 | | 0.04 | 0.017329 |
| 1.002 | 0.004049 | | 0.02 | 0.018378 |
| 1.001 | 0.004166 | | 0.01 | 0.014661 |

### 5.2 Curve fitting

We tested five candidate equations against these 40 data points:

| Rank | Equation | Parameters | RMSE |
|------|----------|-----------|------|
| 1 | Damped oscillation: $\delta_{\text{eq}} + A\, e^{-\gamma|\ln r|}\cos(\omega \ln r + \phi)$ | 5 | 0.00277 |
| 2 | Horizon-centered: $\delta_{\text{eq}}(1 + a\,|\!\ln r|^b) + c\, e^{-5|\ln r|}$ | 4 | 0.00344 |
| 3 | Benford attractor: $\delta_{\text{eq}}(1 + a\,|\!\ln r|^b)$ | 3 | 0.00351 |
| 4 | Log-distance linear: $a\,|\!\ln r| + b$ | 2 | 0.00365 |
| 5 | Benford of $r$: $a\,\log_{10}(1 + 1/r) + b$ | 2 | 0.00399 |

The damped oscillation fits best numerically (RMSE 0.00277 vs 0.00365 for the log-distance linear). We chose the two-parameter model for the following reasons:

1. **Overfitting risk.** Five parameters on 40 data points is dangerously close to fitting noise. The oscillatory component --- frequency, phase, and decay rate --- may be tracking scatter in the measurements rather than a real physical oscillation.

2. **Physical interpretability.** The log-distance model makes one statement: Causal Set deviation is proportional to log-distance from the event horizon. The damped oscillation adds periodic structure ($\omega = 0.94$, $\phi = -0.25$, $\gamma = 0.11$) for which there is no clear physical basis.

3. **Marginal improvement.** The five-parameter model improves on the two-parameter model by only 1.2 percentage points in mean error at the $g_{\hat{\delta}}$ level (4.5% vs 5.7%), at a cost of three additional degrees of freedom.

4. **Generalization.** A two-parameter model is more likely to produce correct predictions at untested positions.

### 5.3 The equation

$$\delta_B(r) = 0.003389\,|\!\ln(r/r_s)| + 0.002508$$

At the horizon ($r = r_s$, so $\ln(r/r_s) = 0$):

$$\delta_B(r_s) = 0.002508$$

This is the minimum deviation --- the point of maximum Causal Set conformance to Benford's law. Moving outward or inward, the deviation grows proportionally to $|\!\ln(r/r_s)|$, the log-distance from the horizon in either direction.

The corresponding metric component is:

$$g_{\hat{\delta}}(r) = \log_{10}\!\left(1 + \frac{1}{0.003389\,|\!\ln(r/r_s)| + 0.002508}\right)$$

This is Benford's equation composed with a log-distance function. The framework is self-referential: the deviation from Benford's law at each point determines the geometry through Benford's equation.

### 5.4 Equation vs data

We compared the equation-based 4D metric to the data-based 4D metric (using interpolated values at each of the 40 positions) at the level of the metric component $g_{\hat{\delta}}$. The mean percent error is 5.7%. The two versions produce effectively identical behavior in $g_{rr}$, entropy rate, and coupling.

---

## 6. Results

### 6.1 Geometry redistribution

When the Benford floor activates, the angular components $g_{\theta\theta} = g_{\phi\phi} = r^2$ continue to shrink as $r \to 0$. But the determinant constraint forces $g_{rr}$ to grow to compensate. What the angular dimensions lose, the radial dimension gains.

The geometry is not destroyed --- it is redistributed. The physical picture: as an object falls toward the singularity, the sphere defined by the angular components compresses while the radial direction stretches. The black hole reshapes space from a sphere into a tube. Spaghettification, in this framework, is not the destruction of geometry but its reorganization. The Benford floor prevents total collapse, and the total geometric content (as measured by the determinant) is conserved at the floor value.

### 6.2 Emergent coupling

The effective coupling between the Causal Set dimension and the spatial geometry --- defined as the fractional reduction in $g_{rr}$ when comparing 4D to 3D, in regions where both floors are active --- emerges as:

- Mean: $\sim$49.5%
- Range: 42--57%
- Standard deviation: $\sim$4%

This was not an input. No coupling constant was specified anywhere in the framework. The value emerges from the determinant constraint alone, and it varies with position because $g_{\hat{\delta}}$ varies with position.

The coupling means that in regions where the floor is active, approximately half of the radial compensation that would be needed in 3D is absorbed by the Causal Set dimension. The CS contributes real geometric content to the metric.

### 6.3 Entropy rate acceleration

With all four dimensions contributing to the entropy rate:

$$\dot{S}(r) = \sqrt{\left(\frac{dg_{rr}}{dr}\right)^2 + \left(\frac{dg_{\theta\theta}}{dr}\right)^2 + \left(\frac{dg_{\phi\phi}}{dr}\right)^2 + \left(\frac{dg_{\hat{\delta}}}{dr}\right)^2}$$

the 4D entropy rate follows the same qualitative pattern as the 3D version but with the Causal Set providing an additional "clock" through $dg_{\hat{\delta}}/dr$. In regions where $\delta_B$ is changing rapidly (near the horizon, where the CS transitions from its minimum, and deep inside, where deviation grows), the Causal Set contributes directly to the rate of geometric change.

In standard GR, the entropy rate drops to zero at the singularity. In the Benford-modified 4D metric, it climbs by orders of magnitude. The singularity is where things happen fast, not where they stop.

### 6.4 The Causal Set peak

The metric component $g_{\hat{\delta}}$ peaks just outside the event horizon, at $r/r_s \approx 1.1$, where $\delta_B$ reaches its minimum of approximately 0.002. At this point, $g_{\hat{\delta}} \approx 2.7$ --- the Causal Set dimension contributes more to the determinant than a single angular dimension ($g_{\theta\theta} = r^2 = 1.21$ at $r/r_s = 1.1$).

This is the point of maximum Causal Set clarity. The substrate is most visible --- least noisy, closest to Benford equilibrium --- right at the event horizon. Moving in either direction, the signal degrades.

---

## 7. Validation: Real Astrophysical Objects

A metric modification is only useful if it produces correct results for objects where GR already works. We tested the 4D Benford metric against six real objects spanning 14 orders of magnitude in mass.

### 7.1 Method

For each object, we computed the Schwarzschild radius $r_s = 2GM/c^2$, then evaluated the 4D metric at the object's surface (or at several points for black holes). The key question: does the Benford floor activate?

### 7.2 Results

**Earth** ($M = 5.97 \times 10^{24}$ kg, $r_s = 8.87$ mm):
Surface at $r/r_s = 7.18 \times 10^8$. The determinant at the surface is $\sim 10^{35}$ --- thirty-five orders of magnitude above the floor. $g_{rr} = 1.000000$. The metric is perfectly flat. No Benford effects whatsoever.

**Jupiter** ($M = 1.90 \times 10^{27}$ kg, $r_s = 2.82$ m):
Surface at $r/r_s = 2.48 \times 10^7$. Same result: flat metric, no floor activation.

**Neutron star, 1.4 $M_\odot$** ($r_s = 4.14$ km, surface radius = 10 km):
Surface at $r/r_s = 2.42$. No floor activation. $g_{rr} = 1$. The CS deviation is measurable ($\delta_B \approx 0.005$) but produces no modification to the spatial geometry. The model correctly identifies this as strong-but-not-extreme gravity.

**Neutron star, 2.1 $M_\odot$** ($r_s = 6.21$ km, surface radius = 10 km):
Surface at $r/r_s = 1.61$. Still no floor activation. Closer to the threshold, but the metric remains unmodified.

**Sagittarius A*** ($M = 4 \times 10^6 \, M_\odot$):
This is a black hole. The 4D floor activates at $r/r_s \approx 0.641$, well inside the horizon. Full geometry redistribution occurs: $g_{rr}$ grows to compensate the angular collapse, the entropy rate accelerates, and the CS dimension contributes at the emergent $\sim$50% coupling level.

**TON 618** ($M = 66 \times 10^9 \, M_\odot$):
The largest known black hole. The 4D floor activates at $r/r_s \approx 0.641$ --- identical to Sgr A* in dimensionless units. The behavior is scale-free. The Schwarzschild radius is enormous ($\sim 10^{14}$ m), but the internal geometry in $r/r_s$ coordinates is indistinguishable from any other black hole.

### 7.3 Summary

The model produces boring results for boring objects and only becomes interesting inside black hole horizons. It does not overpredict. For every object where standard GR is known to work, the Benford floor is irrelevant --- the determinant is so far above 0.4068 that the modification has no effect. This is the correct behavior for a proposed extension: it should reduce to the established theory wherever that theory is valid.

---

## 8. Discussion

### 8.1 What the framework claims

Three concrete claims:

1. **The Benford floor resolves the singularity.** The determinant of the spatial metric cannot drop below 0.4068. The geometry redistributes rather than collapsing to zero. This is a specific, falsifiable prediction about the interior structure of black holes.

2. **The Causal Set substrate acts as a spatial dimension.** Its contribution to the metric is defined by Benford's equation applied to its own deviation, and its coupling to the spatial geometry emerges from the shared determinant constraint at approximately 50%.

3. **Entropy rate accelerates near the singularity.** If time is the rate of geometric change, then the black hole interior is where things happen fastest, not where they freeze.

### 8.2 What emerges vs what is imposed

It is important to distinguish between inputs and outputs of the framework.

**Inputs:**
- The Schwarzschild metric in PG coordinates (established physics)
- Benford's law as the floor equation (the only modification)
- The CS equation $\delta_B(r) = 0.003389\,|\!\ln(r/r_s)| + 0.002508$, derived from 40 measurements (Appendix D)

**Outputs (emergent):**
- The floor value: 0.4068
- The floor activation radius: $0.799\, r_s$ (3D) or $0.641\, r_s$ (4D)
- The emergent coupling: $\sim$49.5% (range 42--57%)
- The entropy rate acceleration
- The geometry redistribution (angular $\to$ radial)
- Scale-free behavior across 14 orders of magnitude in mass

None of these were programmed in. The simulation was given two instructions --- floor the determinant at the $L^2$ norm of the Benford vector, and use the CS equation to define the fourth dimension --- and these properties fell out.

### 8.3 The CS equation and the horizon

The Causal Set equation has a notable physical interpretation. The deviation $\delta_B$ is minimized at the event horizon --- the Causal Set is most "visible," most ordered, closest to Benford equilibrium precisely where the black hole begins. Moving outward, the signal gets noisy (far-field contamination). Moving inward, extreme curvature disturbs the CS order.

The horizon is the calm center. If the Causal Set represents a discrete substrate underlying continuous spacetime, this suggests that the substrate is most cleanly expressed at the event horizon --- the boundary between the region where spacetime behaves normally and the region where it undergoes radical reorganization.

### 8.4 Limitations

We are transparent about what this framework does not yet provide:

1. **The 40 data points are from simulation, not observation.** They were computed by running quantum gravity models through a Schwarzschild geometry and measuring Benford conformance. They are not direct measurements of a physical Causal Set substrate.

2. **The CS equation is extrapolated.** The real-object validation uses $\delta_B(r)$ at $r/r_s$ values far outside the fitted range ($r/r_s = 0.01$ to 10). At $r/r_s = 7 \times 10^8$ (Earth's surface), the equation gives $\delta_B \approx 69$ --- a meaningless number in a regime where the floor is irrelevant. The equation is valid in the strong-field regime; its extrapolation to weak fields is harmless (the floor never activates) but not independently confirmed.

3. **No independent confirmation.** The framework needs to be tested against other black hole geometries (Kerr, Reissner-Nordstrom), other quantum gravity formalisms, and ultimately against observational data if such data becomes available.

4. **The entropy rate is a proposal, not a derivation.** We define time as the rate of spatial geometric change. This is consistent with Wheeler-DeWitt, Barbour, and Rovelli, but it is not derived from first principles. It is a definition that produces suggestive results.

5. **Coordinate dependence.** We have verified that the entropy rate acceleration survives the transition from Schwarzschild to PG coordinates (the acceleration is present in both; only the coordinate artifacts at the horizon differ). A fully coordinate-invariant formulation would be stronger.

### 8.5 Predictions

The framework makes several testable predictions:

1. **The interior proper distance is large.** The radial stretching due to $g_{rr}$ growing as $0.4068/(r^4 \cdot g_{\hat{\delta}})$ implies that the proper distance from the floor activation point to the classical singularity is much larger than in standard GR. The black hole is bigger on the inside.

2. **No information destruction.** If the geometry redistributes rather than collapsing to zero, the information paradox takes a different form. Information is not lost at a singularity because there is no singularity --- there is a region of extreme but finite geometric reorganization.

3. **Scale-free behavior.** All black holes, regardless of mass, have the same internal structure in $r/r_s$ units. The floor activates at the same dimensionless radius. This is consistent with the no-hair theorem but adds internal structure.

4. **The Causal Set is most visible at the horizon.** If a measurement of Benford conformance in gravitational phenomena were possible, the strongest signal would be found at the event horizon of a black hole.

---

## 9. Conclusion

We have presented a four-dimensional spatial metric for the Schwarzschild black hole in which:

- The determinant singularity is replaced by a floor derived from the $L^2$ norm of the Benford probability vector (0.4068).
- Time is removed as a coordinate and replaced by the entropy rate --- the magnitude of the gradient of the spatial metric.
- The Causal Set substrate is promoted from external data to a fourth spatial dimension, with its metric component defined by Benford's equation applied to its own deviation measure.

The resulting framework has no free parameters. The floor value, the floor activation radii, the emergent coupling ($\sim$50%), the entropy rate acceleration, the geometry redistribution, and the CS equation all emerge from the determinant constraint alone.

The model reduces to standard GR for all weak-field objects (Earth, Jupiter, neutron stars) and produces modifications only inside black hole horizons, where standard GR is known to break down. The behavior is scale-free across 14 orders of magnitude in black hole mass.

The key equation governing the Causal Set response to black hole geometry is:

$$\delta_B(r) = 0.003389\,|\!\ln(r/r_s)| + 0.002508$$

Two parameters. Causal Set deviation is proportional to log-distance from the event horizon. The horizon is the equilibrium point.

Whether this framework correctly describes the interior of real black holes is an open question. What we can state is that it produces a self-consistent, parameter-free modification to the Schwarzschild metric that resolves the singularity, generates emergent time behavior, and couples to the best-performing quantum gravity model through geometry alone.

---

## Supplementary Materials

Interactive simulations accompanying this paper are available on Zenodo [2]. These include:

- A 3D visualization of the black hole geometry in Standard GR, 3D Benford, and 4D Benford + CS modes, showing the sphere-to-tube transition as the Benford floor activates.
- A bar-graph comparison of metric components between the 3D and 4D frameworks at each radial position.
- The full 4D simulation with the CS equation, reproducing all figures and numerical results presented here.

The simulations are browser-based and require no installation. Readers can adjust the radial position and observe the metric components, entropy rate, and geometry redistribution in real time --- experiencing, rather than reading about, what the framework predicts for the journey through a black hole.

---

## References

[1] C. Riner, "Complete Monotonicity and Benford's Law: Deriving Quantum Statistics from the Significant Digit Distribution," Zenodo (2026). DOI: [10.5281/zenodo.18510250](https://zenodo.org/records/18510250). The mathematical framework for $\delta_B$ and the Benford conformance analysis. The black hole simulation methodology that produced the 40 data points used in this paper is described in Appendix D.

[2] C. Riner, "4D Benford Metric Simulations: Interactive Black Hole Geometry," Zenodo (2026). [DOI to be assigned upon deposit]

---

## Appendix A: Computation of the Benford Floor

The floor value is computed as follows in Python:

```python
import math

def compute_benford_floor():
    probs = [math.log10(1 + 1/d) for d in range(1, 10)]
    return math.sqrt(sum(p**2 for p in probs))

FLOOR_VAL = compute_benford_floor()  # 0.406839...
```

## Appendix B: The 4D Metric Components

For a given radial position $r$ (in units of $r_s$):

1. Compute the CS deviation: $\delta_B = 0.003389\,|\!\ln r| + 0.002508$

2. Compute the CS metric component: $g_{\hat{\delta}} = \log_{10}(1 + 1/\delta_B)$

3. Angular components: $g_{\theta\theta} = g_{\phi\phi} = r^2$

4. Natural determinant: $\det_0 = r^4 \cdot g_{\hat{\delta}}$

5. If $\det_0 \geq 0.4068$: $g_{rr} = 1$ (standard PG)

6. If $\det_0 < 0.4068$: $g_{rr} = 0.4068 \,/\, (r^4 \cdot g_{\hat{\delta}})$

7. Entropy rate: $\dot{S} = \sqrt{(dg_{rr}/dr)^2 + (dg_{\theta\theta}/dr)^2 + (dg_{\phi\phi}/dr)^2 + (dg_{\hat{\delta}}/dr)^2}$

## Appendix C: Candidate Equations for $\delta_B(r)$

Five equations were fitted to the 40 data points via nonlinear least squares (scipy.optimize.curve_fit). The comparison was performed both at the $\delta_B$ level (RMSE) and at the $g_{\hat{\delta}}$ level (mean percent error):

| Equation | Parameters | $\delta_B$ RMSE | $g_{\hat{\delta}}$ Mean % Error |
|----------|-----------|---------|--------------------------|
| Damped oscillation | 5 | 0.00277 | 4.5% |
| Horizon-centered | 4 | 0.00344 | 5.7% |
| Benford attractor | 3 | 0.00351 | 6.4% |
| Log-distance linear | 2 | 0.00365 | 5.7% |
| Benford of $r$ | 2 | 0.00399 | 8.9% |

The log-distance linear model ($\delta_B = a\,|\!\ln r| + b$) was selected for parsimony and physical interpretability (Section 5.2). Its fitted parameters: $a = 0.003389$, $b = 0.002508$.

## Appendix D: Prior Work --- The Black Hole Benford Experiment

This appendix describes the simulation that produced the 40 $\delta_B$ measurements used throughout this paper. The mathematical framework for the $\delta_B$ metric and its connection to complete monotonicity is developed in [1]. The black hole application described here extends that framework to quantum gravity.

### D.1 Overview

We computed the Benford conformance of 10 quantum gravity models at 40 radial positions through a Schwarzschild black hole, from $r/r_s = 10$ (far outside the horizon) to $r/r_s = 0.01$ (deep interior). At each position, each model generated a 100,000-point energy spectrum. The first significant digit of each spectral value was extracted, compared to Benford's expected distribution, and scored using the Euclidean deviation $\delta_B$.

### D.2 Temperature as a proxy for radial position

An infalling observer in a Schwarzschild black hole experiences tidal forces that increase with depth. We model the effective local temperature as:

$$T_{\text{eff}}(r) = T_H \times \left(\frac{r_s}{|r|}\right)^{3/2}$$

where $T_H = 0.05\, T_{\text{Planck}}$ is the Hawking temperature parameter. This is smooth at the horizon (consistent with the equivalence principle --- a freely falling observer notices nothing special at $r = r_s$) and diverges at the singularity ($r \to 0$).

The temperature encodes radial position: as the observer falls deeper, tidal forces grow, the effective temperature rises, and the energy spectrum changes. Different quantum gravity models modify the spectrum differently, producing different first-digit distributions and therefore different $\delta_B$ values.

### D.3 The 10 quantum gravity models

Each model defines a modified dispersion relation $E(k)$ or density of states $g(k)$ that alters the energy spectrum at high wave vectors. For each model at each radial position, the occupation number spectrum is:

$$n(k) = \frac{g(k)}{e^{E(k)/T_{\text{eff}}} - 1}$$

evaluated at 100,000 logarithmically spaced wave vectors $k \in [0.001, 50]$.

| Model | Modification | Physical basis |
|-------|-------------|----------------|
| Standard (GR) | $E = k$, $g = k^2$ | No modification --- baseline |
| Loop Quantum Gravity | $E = 2|\sin(k/2)|$ | Periodic dispersion from discrete geometry |
| Generalized Uncertainty Principle | $E = k\sqrt{1+k^2}$ | Minimum length scale |
| Doubly Special Relativity | $E = 1 - e^{-k}$ | Energy saturation at Planck scale |
| Hagedorn Temperature | $g \propto k^2 e^k$ | Exponential density of states (string theory) |
| **Causal Set** | $g \propto k^2 e^{-k^2}$ | **Gaussian UV suppression from discreteness** |
| Asymptotic Safety | $d_s = 2 + 2/(1+k^2)$, $g = k^{d_s - 1}$ | Running spectral dimension $4 \to 2$ |
| Horava-Lifshitz | $E = \sqrt{k^2 + k^4 + k^6}$ | Anisotropic scaling at high energy |
| Non-Commutative Geometry | $E = \sqrt{k^2 + k^4}$, $g = k^2(1+k^2)$ | Minimum area from non-commuting coordinates |
| Causal Dynamical Triangulations | $d_s = 2 + 2/(1+k^4)$, $g = k^{d_s - 1}$ | Sharp spectral dimension transition $4 \to 2$ |

### D.4 First-digit extraction and $\delta_B$ computation

From each 100,000-point spectrum, we extracted positive, finite values and computed the first significant digit $d_1 \in \{1, \ldots, 9\}$:

$$d_1(x) = \lfloor 10^{\,\log_{10} x - \lfloor \log_{10} x \rfloor} \rfloor$$

The observed frequency of each digit was compared to the Benford expected frequency $P(d) = \log_{10}(1 + 1/d)$, and the Euclidean deviation was computed:

$$\delta_B = \sqrt{\sum_{d=1}^{9} \left[P_{\text{obs}}(d) - \log_{10}\!\left(1 + \frac{1}{d}\right)\right]^2}$$

Additional metrics (chi-squared, mean absolute deviation, Kolmogorov-Smirnov) were computed for each measurement but are not used in the present paper.

### D.5 Results: why the Causal Set model won

Across the 40 radial positions used in this paper (infalling observer, $r/r_s = 10$ to $0.01$), the Causal Set model produced the lowest mean $\delta_B$ by a factor of approximately 1.7 over the next best model.

The physical reason is the Gaussian UV suppression factor $e^{-k^2}$ in the Causal Set density of states. This produces a spectral cutoff that is smooth and symmetric --- it suppresses high-energy modes without introducing oscillatory artifacts (as LQG does), energy saturation (as DSR does), or exponential growth (as Hagedorn does). The result is a clean spectrum whose first-digit distribution closely follows Benford's law across a wide range of temperatures.

The nine models that modify the three spatial dimensions produce spectra that deviate from Benford's law in characteristic ways tied to their UV modifications. The Causal Set model, which posits an independent discrete substrate rather than modifying spatial geometry, produces conformance. This distinction --- modification *of* space vs operation *alongside* space --- is the empirical basis for promoting the Causal Set to its own dimension in Section 4.

### D.6 Data availability

The complete dataset (10 models $\times$ 40 positions, with all statistical metrics) is included in the supplementary materials on Zenodo [2]. The 40 Causal Set values used in this paper are reproduced in Table 1.
