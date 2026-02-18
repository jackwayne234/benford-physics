#!/usr/bin/env python3
"""Side-by-side comparison: Einstein's time vs Entropy Rate."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import math

fig = plt.figure(figsize=(20, 24))
fig.patch.set_facecolor('#0a0a14')

# ═══════════════════════════════════════════════════════════
# TOP HALF: TWO TENSORS SIDE BY SIDE
# ═══════════════════════════════════════════════════════════

ax_top = fig.add_axes([0.0, 0.42, 1.0, 0.58])
ax_top.set_xlim(0, 1)
ax_top.set_ylim(0, 1)
ax_top.axis('off')

# Main title
ax_top.text(0.5, 0.97, "Einstein  vs  Entropy Rate",
            fontsize=32, fontweight='bold', color='white',
            ha='center', va='top', family='serif')
ax_top.text(0.5, 0.925, "Same measurement.  Different mechanism.  No time dimension required.",
            fontsize=13, color='#888888', ha='center', va='top', family='serif', style='italic')

# Dividing line
ax_top.plot([0.5, 0.5], [0.05, 0.90], color='#333355', linewidth=1.5, linestyle='--', alpha=0.5)

# ── LEFT SIDE: EINSTEIN ──

ax_top.text(0.25, 0.88, "EINSTEIN", fontsize=22, fontweight='bold',
            color='#ff6b6b', ha='center', va='top', family='serif')
ax_top.text(0.25, 0.845, "General Relativity  (1915)", fontsize=11,
            color='#666666', ha='center', va='top', family='serif')

# 4×4 matrix
labels_e = [r'$t$', r'$r$', r'$\theta$', r'$\phi$']
colors_e = ['#ff6b6b', '#4ecdc4', '#ffe66d', '#ffe66d']
entries_e = [
    r'$-\!(1-\frac{r_s}{r})$',
    r'$\frac{1}{1-r_s/r}$',
    r'$r^2$',
    r'$r^2\!\sin^2\!\theta$'
]

mx, my = 0.06, 0.40
mw, mh = 0.38, 0.38
cw, ch = mw/4, mh/4

for j in range(4):
    ax_top.text(mx + (j+0.5)*cw, my + mh + 0.012, labels_e[j],
                fontsize=14, color=colors_e[j], ha='center', va='bottom', fontweight='bold')
for i in range(4):
    ax_top.text(mx - 0.015, my + mh - (i+0.5)*ch, labels_e[i],
                fontsize=14, color=colors_e[i], ha='right', va='center', fontweight='bold')

for i in range(4):
    for j in range(4):
        x = mx + j*cw
        yc = my + mh - (i+1)*ch
        if i == j:
            if i == 0:
                fc, ec, lw = '#2a0a0a', '#ff6b6b', 2.5
            else:
                fc, ec, lw = '#0a1a1a', '#4ecdc4', 1.5
            ax_top.add_patch(patches.FancyBboxPatch(
                (x+0.003, yc+0.003), cw-0.006, ch-0.006,
                boxstyle="round,pad=0.003", facecolor=fc, edgecolor=ec, linewidth=lw))
            ax_top.text(x+cw/2, yc+ch/2, entries_e[i],
                        fontsize=13, color='white', ha='center', va='center')
        else:
            ax_top.add_patch(patches.FancyBboxPatch(
                (x+0.003, yc+0.003), cw-0.006, ch-0.006,
                boxstyle="round,pad=0.003", facecolor='#08080f', edgecolor='#1a1a2a', linewidth=0.5))
            ax_top.text(x+cw/2, yc+ch/2, r'$0$',
                        fontsize=10, color='#222233', ha='center', va='center')

# Brackets
bx0, bx1 = mx - 0.005, mx + mw + 0.005
by0, by1 = my - 0.005, my + mh + 0.005
ax_top.plot([bx0+0.01, bx0, bx0, bx0+0.01], [by1, by1, by0, by0], color='white', lw=2)
ax_top.plot([bx1-0.01, bx1, bx1, bx1-0.01], [by1, by1, by0, by0], color='white', lw=2)

# Einstein's mechanism
ax_top.text(0.25, 0.37, "Time is a dimension in the matrix",
            fontsize=11, color='#ff6b6b', ha='center', va='top', fontweight='bold')
ax_top.text(0.25, 0.34, r"Clock rate  $= \sqrt{-g_{tt}} = \sqrt{1 - r_s/r}$",
            fontsize=14, color='white', ha='center', va='top')
ax_top.text(0.25, 0.30, '"Clocks tick slower because time itself slows\nnear mass."',
            fontsize=10, color='#888888', ha='center', va='top', family='serif', style='italic')

# Arrow pointing to g_tt
ax_top.annotate('', xy=(mx + cw/2, my + mh - ch/2 - 0.01),
                xytext=(mx + cw/2, my + mh - ch/2 + 0.04),
                arrowprops=dict(arrowstyle='->', color='#ff6b6b', lw=1.5))

# ── RIGHT SIDE: RINER ──

ax_top.text(0.75, 0.88, "RINER", fontsize=22, fontweight='bold',
            color='#ffcc44', ha='center', va='top', family='serif')
ax_top.text(0.75, 0.845, "Entropy Rate + Euler Prime Substrate  (2026)", fontsize=11,
            color='#666666', ha='center', va='top', family='serif')

# 5×5 matrix
labels_r = [r'$\dot{S}$', r'$r$', r'$\theta$', r'$\phi$', r'$\hat{\zeta}$']
colors_r = ['#ff9944', '#4ecdc4', '#ffe66d', '#ffe66d', '#ffcc44']
entries_r = [
    r'$\beta\!=\!\sqrt{r_s/r}$',
    r'$1+r_s/r$',
    r'$r^2$',
    r'$r^2\!\sin^2\!\theta$',
    r'$\ln\zeta(s)$'
]

mx2, my2 = 0.545, 0.35
mw2, mh2 = 0.40, 0.44
cw2, ch2 = mw2/5, mh2/5

for j in range(5):
    ax_top.text(mx2 + (j+0.5)*cw2, my2 + mh2 + 0.012, labels_r[j],
                fontsize=13, color=colors_r[j], ha='center', va='bottom', fontweight='bold')
for i in range(5):
    ax_top.text(mx2 - 0.012, my2 + mh2 - (i+0.5)*ch2, labels_r[i],
                fontsize=13, color=colors_r[i], ha='right', va='center', fontweight='bold')

for i in range(5):
    for j in range(5):
        x = mx2 + j*cw2
        yc = my2 + mh2 - (i+1)*ch2
        if i == j:
            if i == 0:
                fc, ec, lw = '#2a1a00', '#ff9944', 2.5
            elif i == 4:
                fc, ec, lw = '#2a1f00', '#ffcc44', 2.5
            else:
                fc, ec, lw = '#0a1a1a', '#4ecdc4', 1.5
            ax_top.add_patch(patches.FancyBboxPatch(
                (x+0.002, yc+0.002), cw2-0.004, ch2-0.004,
                boxstyle="round,pad=0.002", facecolor=fc, edgecolor=ec, linewidth=lw))
            ax_top.text(x+cw2/2, yc+ch2/2, entries_r[i],
                        fontsize=11, color='white', ha='center', va='center')
        else:
            ax_top.add_patch(patches.FancyBboxPatch(
                (x+0.002, yc+0.002), cw2-0.004, ch2-0.004,
                boxstyle="round,pad=0.002", facecolor='#08080f', edgecolor='#1a1a2a', linewidth=0.5))
            ax_top.text(x+cw2/2, yc+ch2/2, r'$0$',
                        fontsize=9, color='#222233', ha='center', va='center')

# Brackets
bx0, bx1 = mx2 - 0.004, mx2 + mw2 + 0.004
by0, by1 = my2 - 0.004, my2 + mh2 + 0.004
ax_top.plot([bx0+0.01, bx0, bx0, bx0+0.01], [by1, by1, by0, by0], color='white', lw=2)
ax_top.plot([bx1-0.01, bx1, bx1, bx1-0.01], [by1, by1, by0, by0], color='white', lw=2)

# Riner's mechanism
ax_top.text(0.75, 0.32, "Time is NOT in the matrix — it's computed FROM it",
            fontsize=11, color='#ff9944', ha='center', va='top', fontweight='bold')
ax_top.text(0.75, 0.29, r"Clock rate  $= \sqrt{1 - \beta^2} = \sqrt{1 - r_s/r}$",
            fontsize=14, color='white', ha='center', va='top')
ax_top.text(0.75, 0.25, '"Clocks tick slower because local geometric activity\nis lower near mass. Less processing = slower clock."',
            fontsize=10, color='#888888', ha='center', va='top', family='serif', style='italic')

# Euler product note
ax_top.add_patch(patches.FancyBboxPatch(
    (0.56, 0.175), 0.38, 0.055,
    boxstyle="round,pad=0.008", facecolor='#1a1500', edgecolor='#ffcc44',
    linewidth=1, alpha=0.7))
ax_top.text(0.75, 0.215,
            r"5th dimension:  $g_{\hat{\zeta}} = \ln\zeta(s) = -\sum_p \ln(1-p^{-s})$",
            fontsize=11, color='#ffcc44', ha='center', va='top')
ax_top.text(0.75, 0.19, "Euler's prime substrate — silent at GPS, active inside black holes",
            fontsize=9, color='#999966', ha='center', va='top', family='serif', style='italic')

# ── RESULT BOXES AT BOTTOM OF TOP SECTION ──

# Einstein result
ax_top.add_patch(patches.FancyBboxPatch(
    (0.06, 0.05), 0.38, 0.08,
    boxstyle="round,pad=0.01", facecolor='#1a0808', edgecolor='#ff6b6b',
    linewidth=2))
ax_top.text(0.25, 0.115, "GPS Prediction", fontsize=10, color='#ff6b6b',
            ha='center', va='top', fontweight='bold')
ax_top.text(0.25, 0.085, "+38.45 μs/day", fontsize=22, color='white',
            ha='center', va='top', fontweight='bold', family='serif')

# Riner result
ax_top.add_patch(patches.FancyBboxPatch(
    (0.56, 0.05), 0.38, 0.08,
    boxstyle="round,pad=0.01", facecolor='#1a1500', edgecolor='#ffcc44',
    linewidth=2))
ax_top.text(0.75, 0.115, "GPS Prediction", fontsize=10, color='#ffcc44',
            ha='center', va='top', fontweight='bold')
ax_top.text(0.75, 0.085, "+38.45 μs/day", fontsize=22, color='white',
            ha='center', va='top', fontweight='bold', family='serif')

# Equals sign between them
ax_top.text(0.5, 0.085, "=", fontsize=36, color='#44ff44',
            ha='center', va='center', fontweight='bold')


# ═══════════════════════════════════════════════════════════
# BOTTOM HALF: CLOCK RATE COMPARISON GRAPH
# ═══════════════════════════════════════════════════════════

ax = fig.add_axes([0.10, 0.08, 0.82, 0.30])
ax.set_facecolor('#0c0c18')

# Physical constants
G = 6.67430e-11
c = 2.99792458e8
M_earth = 5.9722e24
r_s = 2 * G * M_earth / c**2
R_earth = 6.3781e6
R_gps = 2.6571e7

# Compute over range of distances (in units of R_earth)
r_range = np.linspace(1.0, 8.0, 1000)  # 1 to 8 Earth radii
r_meters = r_range * R_earth

# Einstein: clock rate = sqrt(1 - r_s/r)
clock_einstein = np.sqrt(1 - r_s / r_meters)

# Entropy rate: clock rate = sqrt(1 - β²) where β = sqrt(r_s/r)
beta = np.sqrt(r_s / r_meters)
clock_entropy = np.sqrt(1 - beta**2)

# Normalize: show the fractional offset from 1 in parts per billion
offset_einstein = (clock_einstein - 1) * 1e9  # ppb
offset_entropy = (clock_entropy - 1) * 1e9    # ppb

# Plot both (they'll overlap perfectly)
ax.plot(r_range, offset_einstein, color='#ff6b6b', linewidth=3,
        label='Einstein  (from $g_{tt}$)', zorder=3)
ax.plot(r_range, offset_entropy, color='#ffcc44', linewidth=1.5,
        linestyle='--', label='Entropy Rate  (from $\\beta$)', zorder=4)

# Mark Earth surface and GPS orbit
r_gps_earth_radii = R_gps / R_earth

ax.axvline(1.0, color='#4ecdc4', alpha=0.4, linestyle=':', linewidth=1)
ax.axvline(r_gps_earth_radii, color='#4ecdc4', alpha=0.4, linestyle=':', linewidth=1)

# Earth surface annotation
offset_at_earth = (np.sqrt(1 - r_s/R_earth) - 1) * 1e9
offset_at_gps = (np.sqrt(1 - r_s/R_gps) - 1) * 1e9

ax.plot(1.0, offset_at_earth, 'o', color='#4ecdc4', markersize=10, zorder=5)
ax.annotate('Earth\nsurface', xy=(1.0, offset_at_earth),
            xytext=(1.4, offset_at_earth - 0.06),
            fontsize=10, color='#4ecdc4', ha='center',
            arrowprops=dict(arrowstyle='->', color='#4ecdc4', lw=1))

ax.plot(r_gps_earth_radii, offset_at_gps, 'o', color='#4ecdc4', markersize=10, zorder=5)
ax.annotate('GPS\norbit', xy=(r_gps_earth_radii, offset_at_gps),
            xytext=(r_gps_earth_radii - 0.6, offset_at_gps + 0.08),
            fontsize=10, color='#4ecdc4', ha='center',
            arrowprops=dict(arrowstyle='->', color='#4ecdc4', lw=1))

# The dilation between them
ax.annotate('', xy=(5.8, offset_at_earth), xytext=(5.8, offset_at_gps),
            arrowprops=dict(arrowstyle='<->', color='#44ff44', lw=2))
ax.text(6.1, (offset_at_earth + offset_at_gps) / 2,
        '38.45 μs/day',
        fontsize=12, color='#44ff44', fontweight='bold', va='center')

ax.set_xlabel('Distance from Earth center  (Earth radii)', fontsize=12, color='#aaaaaa')
ax.set_ylabel('Clock offset from infinity  (parts per billion)', fontsize=12, color='#aaaaaa')
ax.set_title('Clock Rate vs Distance — Both Methods Identical',
             fontsize=16, color='white', fontweight='bold', family='serif', pad=15)

ax.legend(loc='lower right', fontsize=11, facecolor='#0c0c18', edgecolor='#333355',
          labelcolor='white', framealpha=0.9)

ax.tick_params(colors='#888888')
ax.spines['bottom'].set_color('#333355')
ax.spines['left'].set_color('#333355')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, alpha=0.15, color='#333355')

out = '/home/jackwayne/Desktop/einstein_vs_riner.png'
plt.savefig(out, dpi=180, bbox_inches='tight', facecolor=fig.get_facecolor(), pad_inches=0.3)
plt.close()
print(f"Done: {out}")
