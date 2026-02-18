#!/usr/bin/env python3
"""Generate a 5x5 metric tensor matrix with Euler Prime Substrate."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, ax = plt.subplots(figsize=(18, 20))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
fig.patch.set_facecolor('#0e0e1a')

# Title
ax.text(0.5, 0.97, r"$g_{\mu\nu}^{(5)}$" + "  —  Five-Dimensional Metric Tensor",
        fontsize=26, fontweight='bold', color='white',
        ha='center', va='top', family='serif')
ax.text(0.5, 0.945, "with Benford Floor and Euler Prime Substrate",
        fontsize=16, color='#66ccff',
        ha='center', va='top', family='serif')
ax.text(0.5, 0.925, "Christopher Riner  |  Zenodo DOI: 10.5281/zenodo.18553466",
        fontsize=10, color='#666666', ha='center', va='top', family='serif')

# ── THE 5x5 MATRIX ──────────────────────────────────────

# Matrix layout parameters
mx = 0.12       # left edge of matrix
my = 0.52       # bottom edge of matrix
mw = 0.76       # total width
mh = 0.36       # total height
cw = mw / 5     # cell width
ch = mh / 5     # cell height

# Column/Row labels
coords = [r'$t$', r'$r$', r'$\theta$', r'$\phi$', r'$\hat{\zeta}$']
coord_colors = ['#ff6b6b', '#4ecdc4', '#ffe66d', '#ffe66d', '#ffcc44']

# Column headers
for j, (label, col) in enumerate(zip(coords, coord_colors)):
    ax.text(mx + (j + 0.5) * cw, my + mh + 0.015, label,
            fontsize=18, color=col, ha='center', va='bottom', fontweight='bold')

# Row headers
for i, (label, col) in enumerate(zip(coords, coord_colors)):
    ax.text(mx - 0.02, my + mh - (i + 0.5) * ch, label,
            fontsize=18, color=col, ha='right', va='center', fontweight='bold')

# Diagonal entries (formulas)
diag_formulas = [
    r"$-\!\left(1 - \dfrac{r_s}{r}\right)$",
    r"$1 + \dfrac{r_s}{r}$",
    r"$r^2$",
    r"$r^2 \sin^2\!\theta$",
    r"$\ln\,\zeta(s)$",
]
diag_subtexts = [
    r"$g_{tt}$",
    r"$g_{rr}$",
    r"$g_{\theta\theta}$",
    r"$g_{\phi\phi}$",
    r"$g_{\hat{\zeta}}$",
]

# Draw the matrix grid
for i in range(5):
    for j in range(5):
        x = mx + j * cw
        y_cell = my + mh - (i + 1) * ch

        # Cell background
        if i == j:
            # Diagonal — active cells
            if i == 4:
                # Euler prime substrate — special highlight
                facecolor = '#2a1f00'
                edgecolor = '#ffcc44'
                lw = 2.5
            elif i == 0:
                # Time component
                facecolor = '#1a0a0a'
                edgecolor = '#ff6b6b'
                lw = 2
            else:
                # Spatial components
                facecolor = '#0a1a1a'
                edgecolor = '#4ecdc4'
                lw = 2
        else:
            # Off-diagonal — zeros
            facecolor = '#0a0a14'
            edgecolor = '#1e2233'
            lw = 1

        rect = patches.FancyBboxPatch(
            (x + 0.003, y_cell + 0.003), cw - 0.006, ch - 0.006,
            boxstyle="round,pad=0.003",
            facecolor=facecolor, edgecolor=edgecolor, linewidth=lw)
        ax.add_patch(rect)

        if i == j:
            # Diagonal — show formula
            ax.text(x + cw/2, y_cell + ch/2 + 0.008, diag_formulas[i],
                    fontsize=16, color='white', ha='center', va='center')
            ax.text(x + cw/2, y_cell + 0.012, diag_subtexts[i],
                    fontsize=10, color='#888888', ha='center', va='bottom')
        else:
            # Off-diagonal — zero
            ax.text(x + cw/2, y_cell + ch/2, r"$0$",
                    fontsize=14, color='#333344', ha='center', va='center')

# Big brackets around the matrix
bracket_x_left = mx - 0.008
bracket_x_right = mx + mw + 0.008
bracket_top = my + mh + 0.005
bracket_bot = my - 0.005

# Left bracket
ax.plot([bracket_x_left + 0.015, bracket_x_left, bracket_x_left, bracket_x_left + 0.015],
        [bracket_top, bracket_top, bracket_bot, bracket_bot],
        color='white', linewidth=2.5, solid_capstyle='round')

# Right bracket
ax.plot([bracket_x_right - 0.015, bracket_x_right, bracket_x_right, bracket_x_right - 0.015],
        [bracket_top, bracket_top, bracket_bot, bracket_bot],
        color='white', linewidth=2.5, solid_capstyle='round')

# g_μν label
ax.text(mx - 0.06, my + mh/2, r"$g_{\mu\nu}^{(5)} \;=$",
        fontsize=22, color='white', ha='right', va='center')


# ── COMPONENT DEFINITIONS BELOW THE MATRIX ──────────────

y = my - 0.045

ax.plot([0.05, 0.95], [y + 0.015, y + 0.015], color='#333355', linewidth=1)

# Time
ax.text(0.05, y, "Time:", fontsize=13, fontweight='bold', color='#ff6b6b',
        va='top', family='serif')
ax.text(0.14, y,
        r"$g_{tt} = -(1 - r_s/r)$" + "  —  Standard Schwarzschild.  Timelike outside horizon, spacelike inside.",
        fontsize=12, color='#cccccc', va='top', family='serif')
y -= 0.032

# Radial
ax.text(0.05, y, "Radial:", fontsize=13, fontweight='bold', color='#4ecdc4',
        va='top', family='serif')
ax.text(0.14, y,
        r"$g_{rr} = 1 + r_s/r$" + "  —  Painlev\u00E9-Gullstrand form (river velocity absorbed).  Clamped by floor inside horizon.",
        fontsize=12, color='#cccccc', va='top', family='serif')
y -= 0.032

# Angular
ax.text(0.05, y, "Angular:", fontsize=13, fontweight='bold', color='#ffe66d',
        va='top', family='serif')
ax.text(0.14, y,
        r"$g_{\theta\theta} = r^2, \;\; g_{\phi\phi} = r^2\sin^2\theta$"
        + "  —  Standard spherical geometry.  Shrinks toward singularity.",
        fontsize=12, color='#cccccc', va='top', family='serif')
y -= 0.04

# Prime Substrate — the star of the show
ax.plot([0.05, 0.95], [y + 0.012, y + 0.012], color='#333355', linewidth=1)

ax.text(0.05, y, "Euler Prime Substrate  (5th dimension):", fontsize=14, fontweight='bold',
        color='#ffcc44', va='top', family='serif')
y -= 0.04

ax.text(0.5, y,
        r"$g_{\hat{\zeta}} \;=\; \ln\,\zeta(s) \;=\; -\sum_{p \;\mathrm{prime}} \ln\!\left(1 - p^{-s}\right)$",
        fontsize=22, color='#ffcc44', ha='center', va='top')
y -= 0.045

ax.text(0.5, y, "Replaces the Causal Set dimension.  One formula encoding all 5 CS axioms:",
        fontsize=11, color='#999999', ha='center', va='top', family='serif', style='italic')
y -= 0.022
ax.text(0.5, y, "discrete  \u00B7  partially ordered  \u00B7  locally finite  \u00B7  countable  \u00B7  order = geometry",
        fontsize=11, color='#777777', ha='center', va='top', family='serif', style='italic')
y -= 0.035

# The Euler product itself
ax.text(0.5, y,
        r"$\zeta(s) \;=\; \prod_{p \;\mathrm{prime}} \frac{1}{1 - p^{-s}} \;=\; \sum_{n=1}^{\infty} \frac{1}{n^s}$",
        fontsize=20, color='white', ha='center', va='top')
y -= 0.035
ax.text(0.5, y, "Right: discrete anti-Benford atoms (primes)    =    Left: continuous Benford-like structure",
        fontsize=11, color='#999999', ha='center', va='top', family='serif', style='italic')
y -= 0.04

# Divider
ax.plot([0.05, 0.95], [y + 0.012, y + 0.012], color='#333355', linewidth=1)

# Determinant floor
ax.text(0.05, y, "Determinant Floor:", fontsize=13, fontweight='bold', color='#66ccff',
        va='top', family='serif')
y -= 0.035
ax.text(0.5, y,
        r"$\det\!\left(g_{\mathrm{spatial}}\right) = g_{rr} \cdot r^2 \cdot r^2\sin^2\theta \cdot g_{\hat{\zeta}} \;\geq\; 0.4068$",
        fontsize=18, color='white', ha='center', va='top')
y -= 0.04

# Key insight box
box = patches.FancyBboxPatch(
    (0.08, y - 0.05), 0.84, 0.06,
    boxstyle="round,pad=0.008",
    facecolor='#1a1500', edgecolor='#ffcc44', linewidth=1.5, alpha=0.8)
ax.add_patch(box)

ax.text(0.5, y - 0.02,
        "The anti-Benford substrate (primes) sits in the 5th dimension.\n"
        "Everything physical (Benford-like) emerges on top of it.  The equals sign in Euler's product is the emergence.",
        fontsize=12, color='#ffcc44', ha='center', va='top', family='serif', style='italic')

out = '/home/jackwayne/Desktop/benford_5d_euler_prime_tensor.png'
plt.savefig(out, dpi=180, bbox_inches='tight', facecolor=fig.get_facecolor(), pad_inches=0.3)
plt.close()
print(f"Done: {out}")
