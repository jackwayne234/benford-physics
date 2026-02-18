#!/usr/bin/env python3
"""Generate a clean, readable PNG of the 5D Benford metric with Euler Prime Substrate replacing CS."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(14, 24))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
fig.patch.set_facecolor('#1a1a2e')

# Title
ax.text(0.5, 0.97, "Five-Dimensional Spacetime Metric\nwith Benford Floor and Euler Prime Substrate",
        fontsize=22, fontweight='bold', color='white',
        ha='center', va='top', family='serif')

ax.text(0.5, 0.915, "Christopher Riner  |  Building on: Zenodo DOI: 10.5281/zenodo.18553466",
        fontsize=11, color='#888888', ha='center', va='top', family='serif')

ax.plot([0.1, 0.9], [0.905, 0.905], color='#444466', linewidth=1)

y = 0.89

# --- 1. Benford's Law ---
ax.text(0.05, y, "1. Benford's Law", fontsize=15, fontweight='bold', color='#66ccff',
        va='top', family='serif')
y -= 0.032
ax.text(0.5, y, r"$P(d) = \log_{10}\!\left(1 + \frac{1}{d}\right), \quad d = 1, 2, \ldots, 9$",
        fontsize=19, color='white', ha='center', va='top')
y -= 0.04

# --- 2. Floor Value ---
ax.text(0.05, y, "2. The Floor  (L\u00B2 norm of Benford probability vector)", fontsize=15,
        fontweight='bold', color='#66ccff', va='top', family='serif')
y -= 0.032
ax.text(0.5, y,
        r"$\Vert \mathbf{p} \Vert_2 \;=\; \sqrt{\sum_{d=1}^{9}\left[\log_{10}\!\left(1+\frac{1}{d}\right)\right]^2} \;=\; \sqrt{0.16545} \;=\; 0.4068$",
        fontsize=17, color='white', ha='center', va='top')
y -= 0.045

# --- 3. 5D Spacetime Metric ---
ax.text(0.05, y, "3. The 5D Spacetime Metric  (Painlev\u00E9\u2013Gullstrand coordinates)",
        fontsize=15, fontweight='bold', color='#66ccff', va='top', family='serif')
y -= 0.032
ax.text(0.5, y,
        r"$g_{\mu\nu}^{(5)} = \mathrm{diag}\!\left(\,g_{tt},\;\; g_{rr},\;\; r^2,\;\; r^2,\;\; g_{\hat{\zeta}}\,\right)$",
        fontsize=21, color='white', ha='center', va='top')
y -= 0.035
ax.text(0.5, y, r"Coordinates: $\{t,\; r,\; \theta,\; \phi,\; \hat{\zeta}\}$"
        r"$\quad$\u2014  Einstein's 4D spacetime  +  Euler Prime Substrate as 5th dimension",
        fontsize=11, color='#999999', ha='center', va='top', family='serif', style='italic')
y -= 0.035

# --- 4. Time Component ---
ax.text(0.05, y, "4. Time Component", fontsize=15, fontweight='bold',
        color='#66ccff', va='top', family='serif')
y -= 0.032
ax.text(0.5, y,
        r"$g_{tt} \;=\; -\!\left(1 \;-\; \frac{r_s}{r}\right)$",
        fontsize=21, color='white', ha='center', va='top')
y -= 0.035
ax.text(0.5, y, r"Timelike ($g_{tt} < 0$) outside horizon.   Spacelike ($g_{tt} > 0$) inside."
        "   Independent of the floor.",
        fontsize=11, color='#999999', ha='center', va='top', family='serif', style='italic')
y -= 0.04

# --- 5. Euler Product (The Prime Substrate) ---
ax.text(0.05, y, "5. The Euler Prime Substrate", fontsize=15, fontweight='bold',
        color='#66ccff', va='top', family='serif')
y -= 0.035
ax.text(0.5, y,
        r"$\zeta(s) \;=\; \prod_{p \;\mathrm{prime}} \frac{1}{1 \;-\; p^{-s}}$",
        fontsize=23, color='#ffcc44', ha='center', va='top')
y -= 0.045
ax.text(0.5, y, "Euler's product formula:  discrete prime atoms (right)  =  continuous structure (left)",
        fontsize=11, color='#999999', ha='center', va='top', family='serif', style='italic')
y -= 0.025
ax.text(0.5, y, "Encodes all 5 Causal Set axioms:  discrete, partially ordered, locally finite, countable, order = geometry",
        fontsize=10, color='#777777', ha='center', va='top', family='serif', style='italic')
y -= 0.035

# --- 6. 5th Dimension Metric Component ---
ax.text(0.05, y, "6. The 5th Dimension Metric Component", fontsize=15, fontweight='bold',
        color='#66ccff', va='top', family='serif')
y -= 0.035
ax.text(0.5, y,
        r"$g_{\hat{\zeta}} \;=\; \ln\,\zeta(s) \;=\; -\sum_{p \;\mathrm{prime}} \ln\!\left(1 \;-\; p^{-s}\right)$",
        fontsize=21, color='#ffcc44', ha='center', va='top')
y -= 0.04
ax.text(0.5, y, r"The metric component is the log of the Euler product \u2014 a sum over all primes.",
        fontsize=11, color='#999999', ha='center', va='top', family='serif', style='italic')
y -= 0.025
ax.text(0.5, y, r"Structural parallel:  $g_{\hat{\delta}} = \log_{10}(1 + 1/\delta_B)$  \u2192  $g_{\hat{\zeta}} = -\sum_p \log(1 - p^{-s})$",
        fontsize=11, color='#777777', ha='center', va='top', family='serif', style='italic')
y -= 0.04

# --- 7. Floor Constraint ---
ax.text(0.05, y, "7. Spatial Determinant Floor Constraint", fontsize=15, fontweight='bold',
        color='#66ccff', va='top', family='serif')
y -= 0.032
ax.text(0.5, y,
        r"$\det\!\left(g_{\mathrm{spatial}}\right) \;=\; g_{rr} \cdot r^4 \cdot g_{\hat{\zeta}} \;\;\geq\;\; 0.4068$",
        fontsize=19, color='white', ha='center', va='top')
y -= 0.035
ax.text(0.5, y, r"Applies to spatial dimensions only.   $g_{tt}$ is not part of this product.",
        fontsize=11, color='#999999', ha='center', va='top', family='serif', style='italic')
y -= 0.035

# g_rr piecewise
ax.text(0.5, y,
        r"$g_{rr} = 1 \qquad\qquad\qquad\;\;\; \mathrm{if}\;\; r^4 \cdot g_{\hat{\zeta}} \geq 0.4068$",
        fontsize=17, color='white', ha='center', va='top')
y -= 0.03
ax.text(0.5, y,
        r"$g_{rr} = \frac{0.4068}{r^4 \cdot g_{\hat{\zeta}}} \qquad\quad\;\; \mathrm{if}\;\; r^4 \cdot g_{\hat{\zeta}} < 0.4068$",
        fontsize=17, color='white', ha='center', va='top')
y -= 0.04

# --- 8. Entropy Rate ---
ax.text(0.05, y, "8. Entropy Rate", fontsize=15, fontweight='bold',
        color='#66ccff', va='top', family='serif')
y -= 0.035
ax.text(0.5, y,
        r"$\dot{S}(r) \;=\; \sqrt{\left(\frac{dg_{tt}}{dr}\right)^{\!2} + \left(\frac{dg_{rr}}{dr}\right)^{\!2} + \left(\frac{dg_{\theta\theta}}{dr}\right)^{\!2} + \left(\frac{dg_{\phi\phi}}{dr}\right)^{\!2} + \left(\frac{dg_{\hat{\zeta}}}{dr}\right)^{\!2}}$",
        fontsize=16, color='white', ha='center', va='top')
y -= 0.035
ax.text(0.5, y, "Rate of geometric change across all 5 metric components.   A measure of activity, not time.",
        fontsize=11, color='#999999', ha='center', va='top', family='serif', style='italic')
y -= 0.035

# Divider
ax.plot([0.1, 0.9], [y + 0.005, y + 0.005], color='#444466', linewidth=1)
y -= 0.02

# --- 9. The Substitution ---
ax.text(0.05, y, "9. The Substitution:  Causal Set \u2192 Euler Prime Substrate", fontsize=15,
        fontweight='bold', color='#66ccff', va='top', family='serif')
y -= 0.03

subs = [
    (r"Old CS dimension:   $g_{\hat{\delta}} = \log_{10}(1 + 1/\delta_B)$",
     r"  \u2014  Benford's equation applied to a fitted deviation parameter"),
    (r"New prime substrate:   $g_{\hat{\zeta}} = \ln\,\zeta(s) = -\sum_p \ln(1 - p^{-s})$",
     r"  \u2014  Euler's product, encoding all 5 CS axioms in one formula"),
]
for main, sub in subs:
    ax.text(0.08, y, main, fontsize=13, color='#ffcc44', va='top')
    y -= 0.023
    ax.text(0.08, y, sub, fontsize=10, color='#999999', va='top', family='serif', style='italic')
    y -= 0.025

y -= 0.01

# Divider
ax.plot([0.1, 0.9], [y + 0.005, y + 0.005], color='#444466', linewidth=1)
y -= 0.02

# Key properties
ax.text(0.05, y, "Properties", fontsize=15, fontweight='bold',
        color='#66ccff', va='top', family='serif')
y -= 0.025
properties = [
    "Outside the horizon:  floor inactive, prime substrate inert \u2192 reduces to standard Schwarzschild",
    "GPS time dilation:  38.62 \u03BCs/day  (matches measured value)",
    "Inside the horizon:  floor clamps spatial determinant, prime substrate holds the floor",
    "Primes are anti-Benford:  confirmed via multi-base heatmap (bases 2\u201310, zero signal)",
    "The substrate MUST be anti-Benford for emergence to occur (no gradient \u2192 no emergence)",
    "Euler product: right side = discrete anti-Benford atoms,  left side = continuous Benford-like structure",
    "Scale-free: identical behavior for any black hole mass",
    "No singularity:  geometry is redistributed, not destroyed",
]
for line in properties:
    ax.text(0.08, y, "\u2022  " + line, fontsize=11, color='#cccccc',
            va='top', family='serif')
    y -= 0.023

out = '/home/jackwayne/Desktop/Projects/Benford_Fun/results/benford_5d_euler_prime_metric.png'
plt.savefig(out, dpi=180, bbox_inches='tight', facecolor=fig.get_facecolor(), pad_inches=0.4)
plt.close()
print(f"Done: {out}")
