#!/usr/bin/env python3
"""
3D Wormhole Floor Sphere — Morris-Thorne Benford Metric Visualization

Renders the wormhole mouth as it appears in 3D space: a sphere defined by
the Benford floor boundary, where the metric determinant drops below the
CS scale. The sphere is color-mapped by g_delta-hat (CS dimension strength).

Shows:
  - Translucent floor boundary sphere
  - Cross-section slice revealing internal metric structure
  - The throat at the center
  - CS coloring (Plasma colormap)

Output: results/figures/morris_thorne_3d_sphere.png + ~/Desktop/
"""

import math
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import Normalize

try:
    from scipy.interpolate import interp1d
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# ---------------------------------------------------------------------------
# Wormhole Parameters (from morris_thorne_benford_4d.py)
# ---------------------------------------------------------------------------
B0 = 0.1

# Benford floor
def compute_benford_floor():
    probs = [math.log10(1 + 1 / d) for d in range(1, 10)]
    return math.sqrt(sum(p ** 2 for p in probs))

FLOOR_VAL = compute_benford_floor()

# CS Wormhole data — positive-l side (including throat)
CS_POSITIVE = [
    (0.0, 0.053134), (0.001, 0.053162), (0.002, 0.053086),
    (0.005, 0.052646), (0.01, 0.051016), (0.02, 0.049459),
    (0.03, 0.055167), (0.05, 0.054867), (0.07, 0.060177),
    (0.1, 0.051888), (0.15, 0.059727), (0.2, 0.032515),
    (0.3, 0.078334), (0.5, 0.057310), (0.7, 0.055346),
    (1.0, 0.058440), (1.5, 0.065096), (2.0, 0.105850),
    (3.0, 0.089252), (5.0, 0.038157), (7.0, 0.031026),
    (10.0, 0.034009),
]

def build_interpolator():
    ls = np.array([p[0] for p in CS_POSITIVE])
    dbs = np.array([p[1] for p in CS_POSITIVE])
    if HAS_SCIPY:
        interp_func = interp1d(ls, dbs, kind='linear', fill_value='extrapolate')
        def interpolator(l_arr):
            return interp_func(np.abs(np.asarray(l_arr, dtype=float)))
    else:
        def interpolator(l_arr):
            al = np.abs(np.asarray(l_arr, dtype=float))
            result = np.zeros_like(al)
            for i, a in enumerate(al.flat):
                idx = np.searchsorted(ls, a)
                if idx <= 0:
                    result.flat[i] = dbs[0]
                elif idx >= len(ls):
                    result.flat[i] = dbs[-1]
                else:
                    t = (a - ls[idx-1]) / (ls[idx] - ls[idx-1])
                    result.flat[i] = dbs[idx-1] + t * (dbs[idx] - dbs[idx-1])
            return result
    return interpolator

wh_interp = build_interpolator()

# ---------------------------------------------------------------------------
# Metric functions
# ---------------------------------------------------------------------------
def g_delta_scalar(db):
    return math.log10(1 + 1 / max(db, 1e-6))

def compute_4d_metric_at_l(l_val):
    """Compute 4D metric at a single l value. Returns dict."""
    r2 = l_val**2 + B0**2
    db = float(wh_interp(np.array([l_val]))[0])
    gd = g_delta_scalar(db)
    gtt = r2
    gpp = r2
    det_natural = r2**2 * gd
    if det_natural >= FLOOR_VAL:
        gll = 1.0
        floor_active = False
    else:
        gll = FLOOR_VAL / (r2**2 * gd)
        floor_active = True
    return {
        'gll': gll, 'gtt': gtt, 'gpp': gpp, 'gd': gd,
        'db': db, 'det': gll * gtt * gpp * gd, 'floor': floor_active
    }


# ---------------------------------------------------------------------------
# Main Visualization
# ---------------------------------------------------------------------------
def make_figure():
    print("=" * 60)
    print("3D Wormhole Floor Sphere — Morris-Thorne")
    print(f"Throat radius b₀ = {B0}, Floor = {FLOOR_VAL:.6f}")
    print("=" * 60)

    dark_bg = '#1a1a2e'
    tc = '#e0e0e0'

    fig = plt.figure(figsize=(16, 14), facecolor=dark_bg)

    # ── Panel layout: 2x2 ──
    # Top-left: Full 3D sphere at throat
    # Top-right: Full 3D sphere at l=2 (floor zone edge)
    # Bottom-left: Cross-section (theta=pi/2 slice) at throat
    # Bottom-right: l-scan showing sphere shape evolution

    # ===================================================================
    # Helper: generate sphere mesh at given l
    # ===================================================================
    def make_sphere(l_val, n_theta=50, n_phi=60):
        m = compute_4d_metric_at_l(l_val)
        eq_r = math.sqrt(m['gtt'])
        pol_r = math.sqrt(m['gll'])
        # Power-law compression: x^0.2 brings extreme ratios to viewable range
        # At throat: eq=0.1 → 0.63, pol=56 → 2.24 (ratio ~3.6:1 — a nice tall egg)
        # Far away:  eq=10 → 1.58, pol=1 → 1.0  (ratio ~0.63:1 — oblate, correct)
        disp_eq = eq_r ** 0.2
        disp_pol = pol_r ** 0.2
        disp_eq = max(disp_eq, 0.05)

        theta = np.linspace(0, np.pi, n_theta)
        phi = np.linspace(0, 2 * np.pi, n_phi)
        THETA, PHI = np.meshgrid(theta, phi, indexing='ij')

        X = disp_eq * np.sin(THETA) * np.cos(PHI)
        Y = disp_eq * np.sin(THETA) * np.sin(PHI)
        Z = disp_pol * np.cos(THETA)

        return X, Y, Z, m

    # ===================================================================
    # Panel 1: Sphere at THROAT (l=0) — maximally squeezed
    # ===================================================================
    ax1 = fig.add_subplot(2, 2, 1, projection='3d', facecolor=dark_bg)
    X0, Y0, Z0, m0 = make_sphere(0.0)

    # Color by g_delta
    norm = Normalize(vmin=1.0, vmax=1.5)
    colors0 = cm.plasma(norm(np.full_like(X0, m0['gd'])))
    ax1.plot_surface(X0, Y0, Z0, facecolors=colors0, alpha=0.7,
                     edgecolor='none', shade=True, antialiased=True)

    # Wireframe for structure
    ax1.plot_wireframe(X0, Y0, Z0, color='#ffffff', alpha=0.06,
                       rstride=5, cstride=5, linewidth=0.3)

    # Axes — use 3.0 to fit the power-law compressed shapes
    lim = 3.0
    ax1.set_xlim(-lim, lim)
    ax1.set_ylim(-lim, lim)
    ax1.set_zlim(-lim, lim)
    ax1.set_box_aspect([1, 1, 1])
    ax1.set_title(f'Throat (l = 0)\n'
                  f'g_ll = {m0["gll"]:.1f}, g_θθ = {m0["gtt"]:.4f}',
                  color=tc, fontsize=11, pad=10)
    ax1.set_xlabel('x', color='#666', fontsize=8)
    ax1.set_ylabel('y', color='#666', fontsize=8)
    ax1.set_zlabel('z', color='#666', fontsize=8)
    ax1.tick_params(colors='#444', labelsize=6)
    ax1.xaxis.pane.fill = False
    ax1.yaxis.pane.fill = False
    ax1.zaxis.pane.fill = False
    ax1.xaxis.pane.set_edgecolor('#333')
    ax1.yaxis.pane.set_edgecolor('#333')
    ax1.zaxis.pane.set_edgecolor('#333')
    ax1.grid(True, alpha=0.1)
    ax1.view_init(elev=15, azim=30)

    # ===================================================================
    # Panel 2: Sphere at l=5 — far from throat, nearly round
    # ===================================================================
    ax2 = fig.add_subplot(2, 2, 2, projection='3d', facecolor=dark_bg)
    X5, Y5, Z5, m5 = make_sphere(5.0)

    colors5 = cm.plasma(norm(np.full_like(X5, m5['gd'])))
    ax2.plot_surface(X5, Y5, Z5, facecolors=colors5, alpha=0.7,
                     edgecolor='none', shade=True, antialiased=True)
    ax2.plot_wireframe(X5, Y5, Z5, color='#ffffff', alpha=0.06,
                       rstride=5, cstride=5, linewidth=0.3)

    ax2.set_xlim(-lim, lim)
    ax2.set_ylim(-lim, lim)
    ax2.set_zlim(-lim, lim)
    ax2.set_box_aspect([1, 1, 1])
    ax2.view_init(elev=20, azim=30)
    ax2.set_title(f'Far from throat (l/b₀ = 50)\n'
                  f'g_ll = {m5["gll"]:.4f}, g_θθ = {m5["gtt"]:.2f}',
                  color=tc, fontsize=11, pad=10)
    ax2.set_xlabel('x', color='#666', fontsize=8)
    ax2.set_ylabel('y', color='#666', fontsize=8)
    ax2.set_zlabel('z', color='#666', fontsize=8)
    ax2.tick_params(colors='#444', labelsize=6)
    ax2.xaxis.pane.fill = False
    ax2.yaxis.pane.fill = False
    ax2.zaxis.pane.fill = False
    ax2.xaxis.pane.set_edgecolor('#333')
    ax2.yaxis.pane.set_edgecolor('#333')
    ax2.zaxis.pane.set_edgecolor('#333')
    ax2.grid(True, alpha=0.1)

    # ===================================================================
    # Panel 3: Equatorial cross-section at multiple l positions
    # ===================================================================
    ax3 = fig.add_subplot(2, 2, 3, facecolor=dark_bg)

    l_samples = [0.0, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0]
    phi = np.linspace(0, 2 * np.pi, 200)
    cmap = cm.plasma
    l_norm = Normalize(vmin=0, vmax=10)

    for l_val in l_samples:
        m = compute_4d_metric_at_l(l_val)
        eq_r = math.log10(1 + math.sqrt(m['gtt']))
        eq_r = max(eq_r, 0.01)
        x_circ = eq_r * np.cos(phi)
        y_circ = eq_r * np.sin(phi)
        c = cmap(l_norm(l_val))
        lw = 2.5 if l_val == 0.0 else 1.2
        alpha = 1.0 if l_val == 0.0 else 0.6
        label = f'l/b₀={l_val/B0:.0f}' if l_val > 0 else 'Throat (l=0)'
        ax3.plot(x_circ, y_circ, color=c, linewidth=lw, alpha=alpha, label=label)

    ax3.set_xlim(-2.0, 2.0)
    ax3.set_ylim(-2.0, 2.0)
    ax3.set_aspect('equal')
    ax3.set_title('Equatorial Cross-Section (θ = π/2)\nCircle radius = log₁₀(1 + √g_θθ)',
                  color=tc, fontsize=11, pad=10)
    ax3.set_xlabel('x (log scale)', color='#888', fontsize=9)
    ax3.set_ylabel('y (log scale)', color='#888', fontsize=9)
    ax3.tick_params(colors='#666', labelsize=7)
    for s in ax3.spines.values():
        s.set_color('#333')
    ax3.grid(True, alpha=0.15, color='#333')
    ax3.legend(fontsize=7, facecolor='#2a2a4e', edgecolor='#3a3a5e',
               labelcolor=tc, loc='upper right', ncol=2)

    # ===================================================================
    # Panel 4: Aspect ratio (pol/eq) vs l — shows shape evolution
    # ===================================================================
    ax4 = fig.add_subplot(2, 2, 4, facecolor=dark_bg)

    l_scan = np.linspace(0, 10, 500)
    aspect_std = []
    aspect_3d = []
    aspect_4d = []
    gll_4d_arr = []
    gtt_4d_arr = []

    for lv in l_scan:
        # Standard: gll=1, gtt=l^2+b0^2 → aspect=1/sqrt(l^2+b0^2)
        r2 = lv**2 + B0**2
        aspect_std.append(1.0 / math.sqrt(r2))

        # 3D Benford
        det3 = r2**2
        if det3 < FLOOR_VAL:
            gll_3 = FLOOR_VAL / r2**2
        else:
            gll_3 = 1.0
        aspect_3d.append(math.sqrt(gll_3) / math.sqrt(r2))

        # 4D
        m4 = compute_4d_metric_at_l(lv)
        aspect_4d.append(math.sqrt(m4['gll']) / math.sqrt(m4['gtt']))
        gll_4d_arr.append(m4['gll'])
        gtt_4d_arr.append(m4['gtt'])

    ax4.plot(l_scan / B0, aspect_std, color='#888888', linewidth=1.5,
             label='Standard GR', alpha=0.7)
    ax4.plot(l_scan / B0, aspect_3d, color='#ff4444', linewidth=2.0,
             label='3D Benford', alpha=0.9)
    ax4.plot(l_scan / B0, aspect_4d, color='#ffdd44', linewidth=2.0,
             label='4D Benford + CS', alpha=0.9)

    ax4.axhline(y=1.0, color='#44ff88', linestyle=':', linewidth=0.8,
                alpha=0.5, label='Sphere (aspect=1)')
    ax4.set_yscale('log')
    ax4.set_xlabel('l / b₀', color='#888', fontsize=10)
    ax4.set_ylabel('Polar/Equatorial aspect ratio', color=tc, fontsize=10)
    ax4.set_title('Shape Evolution: Sphere → Tall Ellipsoid at Throat',
                  color=tc, fontsize=11, pad=10)
    ax4.tick_params(colors='#666', labelsize=8)
    for s in ax4.spines.values():
        s.set_color('#333')
    ax4.grid(True, alpha=0.15, color='#333')
    ax4.legend(fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e',
               labelcolor=tc)

    # --- Colorbar ---
    sm = cm.ScalarMappable(cmap='plasma', norm=Normalize(vmin=1.0, vmax=1.5))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=[ax1, ax2], shrink=0.5, pad=0.02, aspect=20)
    cbar.set_label(r'$g_{\hat\delta}$ (CS dimension)', color=tc, fontsize=10)
    cbar.ax.tick_params(colors='#888', labelsize=8)

    fig.suptitle('Morris-Thorne Wormhole — 3D Floor Sphere Visualization\n'
                 '4D Benford Metric with Causal Set Dimension',
                 color='#ffffff', fontsize=15, fontweight='bold', y=0.98)

    fig.subplots_adjust(left=0.06, right=0.94, top=0.91, bottom=0.06, wspace=0.25, hspace=0.30)

    # ── Save ──
    fig_dir = os.path.join(os.path.dirname(__file__), '..', 'results', 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    fig_path = os.path.join(fig_dir, 'morris_thorne_3d_sphere.png')
    fig.savefig(fig_path, dpi=200, facecolor=dark_bg)
    print(f"\nFigure saved to: {fig_path}")

    desktop_path = os.path.expanduser('~/Desktop/morris_thorne_3d_sphere.png')
    fig.savefig(desktop_path, dpi=200, facecolor=dark_bg)
    print(f"Figure copied to: {desktop_path}")

    plt.close(fig)

    # ── Summary ──
    m_throat = compute_4d_metric_at_l(0.0)
    m_far = compute_4d_metric_at_l(10.0)
    print(f"\nAt throat (l=0):")
    print(f"  g_ll = {m_throat['gll']:.4f}, g_θθ = {m_throat['gtt']:.4f}")
    print(f"  g_δ̂  = {m_throat['gd']:.4f}, δ_B = {m_throat['db']:.6f}")
    print(f"  Aspect ratio (pol/eq) = {math.sqrt(m_throat['gll'])/math.sqrt(m_throat['gtt']):.2f}")
    print(f"  Floor active: {m_throat['floor']}")
    print(f"\nAt l/b₀=100 (far away):")
    print(f"  g_ll = {m_far['gll']:.4f}, g_θθ = {m_far['gtt']:.2f}")
    print(f"  Aspect ratio = {math.sqrt(m_far['gll'])/math.sqrt(m_far['gtt']):.4f}")
    print(f"  Floor active: {m_far['floor']}")


if __name__ == '__main__':
    make_figure()
