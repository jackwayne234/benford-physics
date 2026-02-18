#!/usr/bin/env python3
"""Generate a clean, readable PNG of the 5D Benford metric formulas."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(14, 22))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
fig.patch.set_facecolor('#1a1a2e')

# Title
ax.text(0.5, 0.97, "Five-Dimensional Spacetime Metric\nwith Benford Floor and Causal Set Dimension",
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
        r"$g_{\mu\nu}^{(5)} = \mathrm{diag}\!\left(\,g_{tt},\;\; g_{rr},\;\; r^2,\;\; r^2,\;\; g_{\hat{\delta}}\,\right)$",
        fontsize=21, color='white', ha='center', va='top')
y -= 0.035
ax.text(0.5, y, r"Coordinates: $\{t,\; r,\; \theta,\; \phi,\; \hat{\delta}\}$"
        r"$\quad$\u2014  Einstein's 4D spacetime  +  Causal Set as 5th dimension",
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
y -= 0.035

# --- 5. CS Metric Component ---
ax.text(0.05, y, "5. The Causal Set Dimension", fontsize=15, fontweight='bold',
        color='#66ccff', va='top', family='serif')
y -= 0.032
ax.text(0.5, y,
        r"$g_{\hat{\delta}} \;=\; \log_{10}\!\left(1 \;+\; \frac{1}{\delta_B}\right)$",
        fontsize=21, color='white', ha='center', va='top')
y -= 0.03
ax.text(0.5, y, r"Benford's equation applied to the deviation itself:  $\delta_B$ plays the role of the digit $d$",
        fontsize=11, color='#999999', ha='center', va='top', family='serif', style='italic')
y -= 0.035

# --- 6. CS Equation ---
ax.text(0.05, y, "6. The Causal Set Equation", fontsize=15, fontweight='bold',
        color='#66ccff', va='top', family='serif')
y -= 0.035
ax.text(0.5, y,
        r"$\delta_B(r) \;=\; 0.003389 \;\left|\,\ln\!\left(\frac{r}{r_s}\right)\right| \;+\; 0.002508$",
        fontsize=21, color='#ffcc44', ha='center', va='top')
y -= 0.04
ax.text(0.5, y, "Minimum deviation at the horizon.   Grows with log-distance in either direction.",
        fontsize=11, color='#999999', ha='center', va='top', family='serif', style='italic')
y -= 0.03

# Combined g_delta with the equation substituted
ax.text(0.5, y,
        r"$g_{\hat{\delta}}(r) \;=\; \log_{10}\!\left(1 + \frac{1}{0.003389\,\left|\ln(r/r_s)\right| + 0.002508}\right)$",
        fontsize=17, color='#ffcc44', ha='center', va='top')
y -= 0.045

# --- 7. Floor Constraint ---
ax.text(0.05, y, "7. Spatial Determinant Floor Constraint", fontsize=15, fontweight='bold',
        color='#66ccff', va='top', family='serif')
y -= 0.032
ax.text(0.5, y,
        r"$\det\!\left(g_{\mathrm{spatial}}\right) \;=\; g_{rr} \cdot r^4 \cdot g_{\hat{\delta}} \;\;\geq\;\; 0.4068$",
        fontsize=19, color='white', ha='center', va='top')
y -= 0.035
ax.text(0.5, y, r"Applies to spatial dimensions only.   $g_{tt}$ is not part of this product.",
        fontsize=11, color='#999999', ha='center', va='top', family='serif', style='italic')
y -= 0.035

# g_rr piecewise
ax.text(0.5, y,
        r"$g_{rr} = 1 \qquad\qquad\qquad\;\;\; \mathrm{if}\;\; r^4 \cdot g_{\hat{\delta}} \geq 0.4068$",
        fontsize=17, color='white', ha='center', va='top')
y -= 0.03
ax.text(0.5, y,
        r"$g_{rr} = \frac{0.4068}{r^4 \cdot g_{\hat{\delta}}} \qquad\quad\;\; \mathrm{if}\;\; r^4 \cdot g_{\hat{\delta}} < 0.4068$",
        fontsize=17, color='white', ha='center', va='top')
y -= 0.04

# --- 8. Entropy Rate ---
ax.text(0.05, y, "8. Entropy Rate", fontsize=15, fontweight='bold',
        color='#66ccff', va='top', family='serif')
y -= 0.035
ax.text(0.5, y,
        r"$\dot{S}(r) \;=\; \sqrt{\left(\frac{dg_{tt}}{dr}\right)^{\!2} + \left(\frac{dg_{rr}}{dr}\right)^{\!2} + \left(\frac{dg_{\theta\theta}}{dr}\right)^{\!2} + \left(\frac{dg_{\phi\phi}}{dr}\right)^{\!2} + \left(\frac{dg_{\hat{\delta}}}{dr}\right)^{\!2}}$",
        fontsize=16, color='white', ha='center', va='top')
y -= 0.035
ax.text(0.5, y, "Rate of geometric change across all 5 metric components.   A measure of activity, not time.",
        fontsize=11, color='#999999', ha='center', va='top', family='serif', style='italic')
y -= 0.035

# Divider
ax.plot([0.1, 0.9], [y + 0.005, y + 0.005], color='#444466', linewidth=1)
y -= 0.02

# Key properties
ax.text(0.05, y, "Properties", fontsize=15, fontweight='bold',
        color='#66ccff', va='top', family='serif')
y -= 0.025
properties = [
    "Outside the horizon:  floor inactive, CS inert \u2192 reduces to standard Schwarzschild",
    "GPS time dilation:  38.62 \u03BCs/day  (matches measured value)",
    "Inside the horizon:  floor clamps spatial determinant, CS absorbs radial growth",
    "Scale-free: identical behavior for any black hole mass",
    "No singularity:  geometry is redistributed, not destroyed",
]
for line in properties:
    ax.text(0.08, y, "\u2022  " + line, fontsize=11.5, color='#cccccc',
            va='top', family='serif')
    y -= 0.023

plt.savefig('/home/jackwayne/Desktop/benford_5d_metric_formula.png',
            dpi=180, bbox_inches='tight', facecolor=fig.get_facecolor(),
            pad_inches=0.4)
plt.close()
print("Done: ~/Desktop/benford_5d_metric_formula.png")
