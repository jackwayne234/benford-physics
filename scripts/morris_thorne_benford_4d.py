#!/usr/bin/env python3
"""
4D Spatial Metric — Morris-Thorne Wormhole with CS as a Proper Dimension (delta-hat)

Applies the same 4D Benford floor framework used for Schwarzschild to the
Morris-Thorne/Ellis traversable wormhole geometry.

Key differences from the black hole:
  - Coordinate: proper distance l (linear, symmetric about throat)
  - No singularity, no horizon
  - r(l) = sqrt(l^2 + b0^2), throat at l=0 where r = b0
  - g_ll = 1 (radially flat), g_theta = g_phi = l^2 + b0^2
  - CS conformance is weak (delta_B ~ 0.05) vs BH (delta_B ~ 0.003)
  - Floor zone is wide, entropy rate minimises at throat

Compares:
  - 3D Benford floor (wormhole geometry only)
  - 4D with measured CS wormhole data (43 points)
  - 4D with BH equation (wrong geometry — that's the point)

Reference: "Complete Monotonicity and Benford's Law" — C. Riner (2026)
"""

import math
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

try:
    from scipy.interpolate import interp1d
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# ---------------------------------------------------------------------------
# Wormhole Parameters
# ---------------------------------------------------------------------------
B0 = 0.1  # Throat radius

# ---------------------------------------------------------------------------
# CS Wormhole Data — 43 measurements (traversing observer, Causal Set)
# From results/round_trip/wormhole_wall.json, "Causal Set" key
# Format: (l, delta_b) — proper distance from throat
# ---------------------------------------------------------------------------
CS_WORMHOLE_DATA = [
    (-10.0, 0.034009), (-7.0, 0.031026), (-5.0, 0.038157),
    (-3.0, 0.089252), (-2.0, 0.105850), (-1.5, 0.065096),
    (-1.0, 0.058440), (-0.7, 0.055346), (-0.5, 0.057310),
    (-0.3, 0.078334), (-0.2, 0.032515), (-0.15, 0.059727),
    (-0.1, 0.051888), (-0.07, 0.060177), (-0.05, 0.054867),
    (-0.03, 0.055167), (-0.02, 0.049459), (-0.01, 0.051016),
    (-0.005, 0.052646), (-0.002, 0.053086), (-0.001, 0.053162),
    (0.0, 0.053134),
    (0.001, 0.053162), (0.002, 0.053086), (0.005, 0.052646),
    (0.01, 0.051016), (0.02, 0.049459), (0.03, 0.055167),
    (0.05, 0.054867), (0.07, 0.060177), (0.1, 0.051888),
    (0.15, 0.059727), (0.2, 0.032515), (0.3, 0.078334),
    (0.5, 0.057310), (0.7, 0.055346), (1.0, 0.058440),
    (1.5, 0.065096), (2.0, 0.105850), (3.0, 0.089252),
    (5.0, 0.038157), (7.0, 0.031026), (10.0, 0.034009),
]

# Use only positive-l side (including throat) for interpolation on |l|
# The data is symmetric, so we interpolate on |l| and mirror
CS_POSITIVE = [(abs(l), db) for l, db in CS_WORMHOLE_DATA if l >= 0]
CS_POSITIVE.sort(key=lambda x: x[0])  # ascending in |l|


def build_wormhole_interpolator():
    """Build interpolator for CS delta_B as a function of |l|."""
    ls = np.array([p[0] for p in CS_POSITIVE])
    dbs = np.array([p[1] for p in CS_POSITIVE])

    if HAS_SCIPY:
        # Linear interpolation in (|l|, delta_b) space
        # Use linear (not log) since l includes 0
        interp_func = interp1d(ls, dbs, kind='linear', fill_value='extrapolate')

        def interpolator(l_arr):
            return interp_func(np.abs(l_arr))
    else:
        def interpolator(l_arr):
            result = np.zeros_like(l_arr, dtype=float)
            al = np.abs(l_arr)
            for i, ali in enumerate(al):
                idx = np.searchsorted(ls, ali)
                if idx <= 0:
                    result[i] = dbs[0]
                elif idx >= len(ls):
                    result[i] = dbs[-1]
                else:
                    t = (ali - ls[idx - 1]) / (ls[idx] - ls[idx - 1])
                    result[i] = dbs[idx - 1] + t * (dbs[idx] - dbs[idx - 1])
            return result

    return interpolator


# ---------------------------------------------------------------------------
# BH equation mapped to wormhole areal radius (for comparison — wrong geometry)
# ---------------------------------------------------------------------------
BH_CS_A = 0.003389
BH_CS_B = 0.002508


def delta_b_bh_equation(l_arr):
    """BH CS equation mapped to wormhole: delta_B = a*|ln(r(l)/b0)| + b
    where r(l) = sqrt(l^2 + b0^2). This is WRONG for the wormhole — that's the point."""
    r_areal = np.sqrt(l_arr**2 + B0**2)
    return BH_CS_A * np.abs(np.log(r_areal / B0)) + BH_CS_B


# ---------------------------------------------------------------------------
# Benford Floor
# ---------------------------------------------------------------------------
def compute_benford_floor():
    probs = [math.log10(1 + 1 / d) for d in range(1, 10)]
    return math.sqrt(sum(p ** 2 for p in probs))


FLOOR_VAL = compute_benford_floor()


# ---------------------------------------------------------------------------
# Wormhole Metrics
# ---------------------------------------------------------------------------
def compute_3d_metric(l):
    """3D Morris-Thorne spatial metric with Benford floor.
    g_ll = 1, g_theta = g_phi = l^2 + b0^2.
    det_natural = 1 * (l^2+b0^2) * (l^2+b0^2) = (l^2+b0^2)^2."""
    g_theta = l**2 + B0**2
    g_phi = l**2 + B0**2
    det_natural = (l**2 + B0**2)**2  # with g_ll=1

    g_ll = np.ones_like(l)
    floor_active = det_natural < FLOOR_VAL
    g_ll[floor_active] = FLOOR_VAL / ((l[floor_active]**2 + B0**2)**2)

    return g_ll, g_theta, g_phi, floor_active


def compute_4d_metric(l, delta_b_func):
    """4D Morris-Thorne metric: 3 spatial + g_delta = log10(1 + 1/delta_B).
    det_natural = (l^2+b0^2)^2 * g_delta."""
    delta_b = delta_b_func(l)
    delta_b = np.clip(delta_b, 1e-6, None)

    g_delta = np.log10(1.0 + 1.0 / delta_b)
    g_theta = l**2 + B0**2
    g_phi = l**2 + B0**2

    det_natural = (l**2 + B0**2)**2 * g_delta

    g_ll = np.ones_like(l)
    floor_active = det_natural < FLOOR_VAL
    g_ll[floor_active] = FLOOR_VAL / (
        (l[floor_active]**2 + B0**2)**2 * g_delta[floor_active]
    )

    return g_ll, g_theta, g_phi, g_delta, floor_active, delta_b


def entropy_rate(l, *components):
    """Entropy rate = sqrt(sum of (dg_i/dl)^2)."""
    total_sq = np.zeros_like(l)
    for g in components:
        dg_dl = np.gradient(g, l)
        total_sq += dg_dl**2
    return np.sqrt(total_sq)


# ---------------------------------------------------------------------------
# Main Simulation
# ---------------------------------------------------------------------------
def run_simulation():
    print("=" * 70)
    print("4D Spatial Metric — Morris-Thorne Wormhole")
    print(f"Throat radius b₀ = {B0}")
    print("=" * 70)

    # Linear grid, symmetric about throat
    l = np.linspace(-10, 10, 2001)

    # Build interpolator from wormhole CS data
    wh_interp = build_wormhole_interpolator()

    # --- Three metric versions ---
    g_ll_3d, g_th_3d, g_ph_3d, floor_3d = compute_3d_metric(l)

    g_ll_wh, g_th_wh, g_ph_wh, g_d_wh, floor_wh, db_wh = \
        compute_4d_metric(l, wh_interp)

    g_ll_bh, g_th_bh, g_ph_bh, g_d_bh, floor_bh, db_bh = \
        compute_4d_metric(l, delta_b_bh_equation)

    # --- Entropy rates ---
    erate_3d = entropy_rate(l, g_ll_3d, g_th_3d, g_ph_3d)
    erate_wh = entropy_rate(l, g_ll_wh, g_th_wh, g_ph_wh, g_d_wh)
    erate_bh = entropy_rate(l, g_ll_bh, g_th_bh, g_ph_bh, g_d_bh)

    # --- Coupling (where both 3D and 4D floors active) ---
    both_floor_wh = floor_3d & floor_wh
    coupling_wh = np.full_like(l, np.nan)
    coupling_wh[both_floor_wh] = 1.0 - (g_ll_wh[both_floor_wh] / g_ll_3d[both_floor_wh])

    both_floor_bh = floor_3d & floor_bh
    coupling_bh = np.full_like(l, np.nan)
    coupling_bh[both_floor_bh] = 1.0 - (g_ll_bh[both_floor_bh] / g_ll_3d[both_floor_bh])

    # --- Analysis ---
    # 3D floor boundary: (l^2 + b0^2)^2 < FLOOR_VAL => l^2 + b0^2 < sqrt(FLOOR_VAL)
    # => l^2 < sqrt(FLOOR_VAL) - b0^2
    l_floor_3d_sq = math.sqrt(FLOOR_VAL) - B0**2
    l_floor_3d = math.sqrt(l_floor_3d_sq) if l_floor_3d_sq > 0 else 0.0

    # 4D floor boundaries (from grid)
    fi_wh = np.where(floor_wh)[0]
    fi_bh = np.where(floor_bh)[0]

    # Find extent of floor zone for wormhole data
    if len(fi_wh) > 0:
        l_floor_wh_min = l[fi_wh[0]]
        l_floor_wh_max = l[fi_wh[-1]]
    else:
        l_floor_wh_min = l_floor_wh_max = None

    if len(fi_bh) > 0:
        l_floor_bh_min = l[fi_bh[0]]
        l_floor_bh_max = l[fi_bh[-1]]
    else:
        l_floor_bh_min = l_floor_bh_max = None

    # Throat values (l=0)
    idx_throat = np.argmin(np.abs(l))
    g_ll_throat_3d = g_ll_3d[idx_throat]
    g_ll_throat_wh = g_ll_wh[idx_throat]
    g_ll_throat_bh = g_ll_bh[idx_throat]
    g_d_throat_wh = g_d_wh[idx_throat]
    g_d_throat_bh = g_d_bh[idx_throat]
    db_throat_wh = db_wh[idx_throat]
    db_throat_bh = db_bh[idx_throat]

    # Critical throat radius: b_crit = FLOOR_VAL^(1/4)
    b_crit = FLOOR_VAL**0.25

    # Traversability: proper distance through floor zone
    dl = l[1] - l[0]
    sqrt_gll_wh = np.sqrt(g_ll_wh)
    _trapz = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
    proper_dist_total = _trapz(sqrt_gll_wh, l)
    if len(fi_wh) > 0:
        proper_dist_floor = _trapz(sqrt_gll_wh[fi_wh[0]:fi_wh[-1]+1],
                                   l[fi_wh[0]:fi_wh[-1]+1])
    else:
        proper_dist_floor = 0.0

    print(f"\n1. Benford floor value: {FLOOR_VAL:.6f}")
    print(f"   Critical throat radius: b_crit = FLOOR^(1/4) = {b_crit:.6f}")
    print(f"   Our b₀ = {B0} < b_crit → floor IS active")

    print(f"\n2. Floor activation zones:")
    print(f"   3D floor:       |l/b₀| < {l_floor_3d/B0:.2f}  (|l| < {l_floor_3d:.4f})")
    if l_floor_wh_min is not None:
        print(f"   4D(wh data):    l ∈ [{l_floor_wh_min:.4f}, {l_floor_wh_max:.4f}]"
              f"  (|l/b₀| < {abs(l_floor_wh_max)/B0:.2f})")
    if l_floor_bh_min is not None:
        print(f"   4D(BH eq):      l ∈ [{l_floor_bh_min:.4f}, {l_floor_bh_max:.4f}]"
              f"  (|l/b₀| < {abs(l_floor_bh_max)/B0:.2f})")

    print(f"\n3. Metric at throat (l=0):")
    print(f"   g_ll (3D):      {g_ll_throat_3d:.4f}")
    print(f"   g_ll (4D wh):   {g_ll_throat_wh:.4f}")
    print(f"   g_ll (4D BH):   {g_ll_throat_bh:.4f}")
    print(f"   g_δ̂ (wh data):  {g_d_throat_wh:.4f}  (δ_B = {db_throat_wh:.6f})")
    print(f"   g_δ̂ (BH eq):    {g_d_throat_bh:.4f}  (δ_B = {db_throat_bh:.6f})")
    print(f"   g_θθ = g_φφ:    {B0**2:.4f}")
    print(f"   det(3D) at throat: {B0**4:.6e}")
    print(f"   det(4D wh) at throat: {B0**4 * g_d_throat_wh:.6e}")

    # Coupling
    vc_wh = coupling_wh[~np.isnan(coupling_wh)]
    vc_bh = coupling_bh[~np.isnan(coupling_bh)]

    print(f"\n4. Emergent coupling:")
    if len(vc_wh) > 0:
        print(f"   Wormhole data: mean = {np.mean(vc_wh)*100:.2f}%,"
              f" range [{np.min(vc_wh)*100:.2f}%, {np.max(vc_wh)*100:.2f}%]")
    if len(vc_bh) > 0:
        print(f"   BH equation:   mean = {np.mean(vc_bh)*100:.2f}%,"
              f" range [{np.min(vc_bh)*100:.2f}%, {np.max(vc_bh)*100:.2f}%]")
    print(f"   (Black hole coupling was ~49.5%)")

    # Entropy rate at throat
    print(f"\n5. Entropy rate at throat:")
    print(f"   3D:  {erate_3d[idx_throat]:.6f}")
    print(f"   4D (wh data): {erate_wh[idx_throat]:.6f}")
    print(f"   4D (BH eq):   {erate_bh[idx_throat]:.6f}")
    print(f"   (Entropy rate should be at a MINIMUM at throat — opposite of BH singularity)")

    # Verify symmetry
    idx_plus5 = np.argmin(np.abs(l - 5.0))
    idx_minus5 = np.argmin(np.abs(l + 5.0))
    sym_diff = abs(g_ll_wh[idx_plus5] - g_ll_wh[idx_minus5])
    print(f"\n6. Symmetry check:")
    print(f"   g_ll(l=+5) - g_ll(l=-5) = {sym_diff:.6e}  {'✓' if sym_diff < 1e-6 else '⚠'}")

    # Traversability
    print(f"\n7. Traversability:")
    print(f"   g_ll is finite everywhere: "
          f"max = {np.max(g_ll_wh):.4f} ({'✓ FINITE' if np.isfinite(np.max(g_ll_wh)) else '✗ DIVERGES'})")
    print(f"   Proper distance through full path (l=-10 to +10): {proper_dist_total:.4f}")
    print(f"   Proper distance through floor zone: {proper_dist_floor:.4f}")
    print(f"   Wormhole remains TRAVERSABLE — finite proper distance everywhere")

    # g_delta range
    print(f"\n8. g_δ̂ range:")
    print(f"   Wormhole data: [{np.min(g_d_wh):.4f}, {np.max(g_d_wh):.4f}]")
    print(f"   BH equation:   [{np.min(g_d_bh):.4f}, {np.max(g_d_bh):.4f}]")
    print(f"   (BH had range [1.5, 2.6] — wormhole CS is much weaker)")

    print("\n" + "=" * 70)

    # ===================================================================
    # VISUALIZATION — 4 panels (dark theme)
    # ===================================================================
    fig = plt.figure(figsize=(14, 20), facecolor='#1a1a2e')
    gs = GridSpec(4, 1, hspace=0.25, top=0.94, bottom=0.04, left=0.10, right=0.95)

    dark_bg = '#1a1a2e'
    tc = '#e0e0e0'
    gc = '#2a2a4e'

    def style(ax, ylabel):
        ax.set_facecolor(dark_bg)
        ax.tick_params(colors=tc, which='both')
        ax.set_ylabel(ylabel, color=tc, fontsize=12)
        for s in ax.spines.values():
            s.set_color(gc)
        ax.grid(True, alpha=0.2, color=gc)
        ax.set_xlim(-10, 10)

    # --- Panel 1: All 4 metric components (wormhole data) ---
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(l, g_ll_wh, color='#ff4444', linewidth=2.0, label=r'$g_{ll}$ (4D+CS, wormhole data)')
    ax1.plot(l, g_ll_3d, color='#ff4444', linewidth=1.2, alpha=0.4, linestyle='--',
             label=r'$g_{ll}$ (3D only)')
    ax1.plot(l, g_th_wh, color='#44ff44', linewidth=2.0,
             label=r'$g_{\theta\theta}$ = $g_{\phi\phi}$ = $l^2 + b_0^2$')
    ax1.plot(l, g_d_wh, color='#ffff00', linewidth=2.0, label=r'$g_{\hat\delta}$ (CS)')
    ax1.axvline(x=0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.4, label='Throat')
    ax1.axhline(y=FLOOR_VAL, color='#ff8800', linestyle=':', linewidth=1, alpha=0.5,
                label=f'Floor = {FLOOR_VAL:.4f}')
    ax1.set_yscale('log')
    ax1.set_ylim(1e-4, 1e5)
    style(ax1, 'Metric Component')
    ax1.set_title('Panel 1: All 4 Dimensions — The Full Wormhole Metric',
                  color=tc, fontsize=13, pad=10)
    ax1.legend(fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e', labelcolor=tc,
               loc='upper right')

    # --- Panel 2: g_delta comparison (wormhole data vs BH equation) ---
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(l, g_d_wh, color='#ffdd44', linewidth=1.8,
             label=r'$g_{\hat\delta}$ (wormhole CS data)')
    ax2.plot(l, g_d_bh, color='#00ffaa', linewidth=2.0, linestyle='--',
             label=r'$g_{\hat\delta}$ (BH equation — wrong geometry)')
    # Scatter the 43 measured points
    l_pts = np.array([p[0] for p in CS_WORMHOLE_DATA])
    db_pts = np.array([p[1] for p in CS_WORMHOLE_DATA])
    gd_pts = np.log10(1.0 + 1.0 / db_pts)
    ax2.scatter(l_pts, gd_pts, color='#ffdd44', s=25, zorder=5, edgecolors='#aa9900',
                linewidths=0.5, label='43 measured points')
    ax2.axvline(x=0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.4)
    style(ax2, r'$g_{\hat\delta}$')
    ax2.set_ylim(0.8, 2.0)
    ax2.set_title('Panel 2: CS Dimension — Wormhole Data vs BH Equation',
                  color=tc, fontsize=13, pad=10)
    ax2.legend(fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e', labelcolor=tc)

    # --- Panel 3: Entropy rate ---
    ax3 = fig.add_subplot(gs[2])
    ax3.plot(l, erate_3d, color='#888888', linewidth=1.5, label='3D entropy rate')
    ax3.plot(l, erate_wh, color='#ff4444', linewidth=1.8,
             label='4D entropy rate (wormhole data)', alpha=0.8)
    ax3.plot(l, erate_bh, color='#00ffaa', linewidth=1.8,
             label='4D entropy rate (BH equation)', linestyle='--')
    ax3.axvline(x=0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.4)
    ax3.set_yscale('log')
    style(ax3, 'Entropy Rate')
    ax3.set_title('Panel 3: Entropy Rate — Minimum at Throat (opposite of BH!)',
                  color=tc, fontsize=13, pad=10)
    ax3.legend(fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e', labelcolor=tc)

    # --- Panel 4: Emergent coupling ---
    ax4 = fig.add_subplot(gs[3])
    m_wh = ~np.isnan(coupling_wh)
    m_bh = ~np.isnan(coupling_bh)
    if np.any(m_wh):
        ax4.plot(l[m_wh], coupling_wh[m_wh] * 100, color='#ff4444', linewidth=1.8,
                 label='Coupling (wormhole data)', alpha=0.8)
    if np.any(m_bh):
        ax4.plot(l[m_bh], coupling_bh[m_bh] * 100, color='#00ffaa', linewidth=1.8,
                 label='Coupling (BH equation)', linestyle='--')
    ax4.axvline(x=0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.4)
    style(ax4, 'Coupling (%)')
    ax4.set_xlabel(r'$l / b_0$ (proper distance from throat, units of $b_0$)', color=tc, fontsize=12)
    ax4.set_title('Panel 4: Emergent Coupling — Wormhole vs BH Equation',
                  color=tc, fontsize=13, pad=10)
    ax4.legend(fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e', labelcolor=tc)

    fig.suptitle('4D Benford Metric: Morris-Thorne Traversable Wormhole',
                 color='#ffffff', fontsize=16, fontweight='bold')

    # Save
    fig_dir = os.path.join(os.path.dirname(__file__), '..', 'results', 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    fig_path = os.path.join(fig_dir, 'morris_thorne_benford_4d.png')
    fig.savefig(fig_path, dpi=200, facecolor=dark_bg)
    print(f"\nFigure saved to: {fig_path}")

    desktop_path = os.path.expanduser('~/Desktop/morris_thorne_benford_4d.png')
    fig.savefig(desktop_path, dpi=200, facecolor=dark_bg)
    print(f"Figure copied to: {desktop_path}")

    plt.close(fig)

    return {
        'floor_val': FLOOR_VAL,
        'b_crit': b_crit,
        'l_floor_3d': l_floor_3d,
        'l_floor_4d_wh': (l_floor_wh_min, l_floor_wh_max) if l_floor_wh_min is not None else None,
        'g_ll_throat_3d': float(g_ll_throat_3d),
        'g_ll_throat_4d_wh': float(g_ll_throat_wh),
        'g_delta_throat_wh': float(g_d_throat_wh),
        'db_throat_wh': float(db_throat_wh),
        'coupling_wh_mean': float(np.mean(vc_wh)) if len(vc_wh) > 0 else None,
        'coupling_bh_mean': float(np.mean(vc_bh)) if len(vc_bh) > 0 else None,
        'erate_throat_3d': float(erate_3d[idx_throat]),
        'erate_throat_wh': float(erate_wh[idx_throat]),
        'proper_dist_total': float(proper_dist_total),
        'proper_dist_floor': float(proper_dist_floor),
        'traversable': bool(np.all(np.isfinite(g_ll_wh))),
    }


if __name__ == '__main__':
    results = run_simulation()
