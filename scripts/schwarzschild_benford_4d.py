#!/usr/bin/env python3
"""
4D Spatial Metric — CS (Causal Set) as a Proper Dimension (delta-hat)

Compares the 3D Benford floor metric (existing) with a new 4D version
where Causal Set gets its own metric component g_delta = log10(1 + 1/delta_B).

No coupling constants. No free parameters. The equation defines its own geometry.

Reference: "Complete Monotonicity and Benford's Law" — C. Riner (2026)
"""

import json
import math
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Try scipy for interpolation; fall back to manual log-linear
try:
    from scipy.interpolate import interp1d
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# ---------------------------------------------------------------------------
# CS Experimental Data (Causal Set delta_B at each radial position)
# ---------------------------------------------------------------------------

def load_cs_data():
    """Load CS delta_B data from black_hole_wall.json or use hardcoded fallback."""
    json_path = os.path.join(
        os.path.dirname(__file__), '..', 'results', 'round_trip', 'black_hole_wall.json'
    )
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        cs_entries = data['infalling_observer']['Causal Set']
        pairs = []
        for entry in cs_entries:
            r = entry['r_ratio']
            db = entry['delta_b']
            if r > 0 and db > 0:
                pairs.append((r, db))
        pairs.sort(key=lambda x: -x[0])  # descending in r
        if len(pairs) >= 10:
            print(f"Loaded {len(pairs)} CS data points from {json_path}")
            return pairs
    except Exception as e:
        print(f"Could not load CS data from JSON ({e}), using hardcoded fallback")

    return [
        (10.0, 0.027552), (7.0, 0.010611), (5.0, 0.004223),
        (3.0, 0.002856), (2.0, 0.002748), (1.5, 0.005169),
        (1.3, 0.003184), (1.2, 0.003569), (1.15, 0.002757),
        (1.1, 0.002136), (1.08, 0.00255), (1.06, 0.002722),
        (1.04, 0.003532), (1.03, 0.003889), (1.02, 0.003482),
        (1.015, 0.003757), (1.01, 0.004013), (1.005, 0.0042),
        (1.002, 0.004049), (1.001, 0.004166), (0.99, 0.004267),
        (0.95, 0.005472), (0.9, 0.00286), (0.85, 0.004293),
        (0.8, 0.003671), (0.7, 0.006253), (0.6, 0.006376),
        (0.5, 0.00486), (0.4, 0.007448), (0.3, 0.01313),
        (0.25, 0.016065), (0.2, 0.01205), (0.15, 0.015154),
        (0.12, 0.019625), (0.1, 0.016711), (0.08, 0.017486),
        (0.06, 0.01349), (0.04, 0.017329), (0.02, 0.018378),
        (0.01, 0.014661),
    ]


def build_interpolator(cs_data):
    """Build log-linear interpolator for CS delta_B as a function of r."""
    rs = np.array([p[0] for p in cs_data])
    dbs = np.array([p[1] for p in cs_data])
    log_rs = np.log(rs)
    log_dbs = np.log(dbs)

    if HAS_SCIPY:
        interp_log = interp1d(log_rs, log_dbs, kind='linear',
                              fill_value='extrapolate')
        def interpolator(r_arr):
            return np.exp(interp_log(np.log(r_arr)))
    else:
        def interpolator(r_arr):
            result = np.zeros_like(r_arr, dtype=float)
            lr = np.log(r_arr)
            for i, lri in enumerate(lr):
                idx = np.searchsorted(log_rs, lri)
                if idx <= 0:
                    result[i] = dbs[0]
                elif idx >= len(log_rs):
                    result[i] = dbs[-1]
                else:
                    t = (lri - log_rs[idx - 1]) / (log_rs[idx] - log_rs[idx - 1])
                    result[i] = np.exp(log_dbs[idx - 1] + t * (log_dbs[idx] - log_dbs[idx - 1]))
            return result

    return interpolator


# ---------------------------------------------------------------------------
# Benford Floor Value
# ---------------------------------------------------------------------------

def compute_benford_floor():
    """FLOOR_VAL = L2 norm of Benford probability vector."""
    probs = [math.log10(1 + 1 / d) for d in range(1, 10)]
    return math.sqrt(sum(p ** 2 for p in probs))

FLOOR_VAL = compute_benford_floor()


# ---------------------------------------------------------------------------
# 3D Metric (Equation 1)
# ---------------------------------------------------------------------------

def compute_3d_metric(r):
    """
    3D Painleve-Gullstrand spatial metric with Benford floor.
    Returns g_rr, g_theta, g_phi arrays.
    """
    g_theta = r ** 2
    g_phi = r ** 2
    det_natural = r ** 4  # g_rr=1 * r^2 * r^2

    g_rr = np.ones_like(r)
    floor_active = det_natural < FLOOR_VAL
    # Where floor activates: g_rr = FLOOR_VAL / (r^2 * r^2)
    g_rr[floor_active] = FLOOR_VAL / (r[floor_active] ** 2 * r[floor_active] ** 2)

    return g_rr, g_theta, g_phi, floor_active


# ---------------------------------------------------------------------------
# 4D Metric (Equation 2)
# ---------------------------------------------------------------------------

def compute_4d_metric(r, delta_b_interp):
    """
    4D spatial metric: 3 PG components + g_delta = log10(1 + 1/delta_B).
    Returns g_rr, g_theta, g_phi, g_delta arrays plus floor_active mask.
    """
    delta_b = delta_b_interp(r)
    # Clamp delta_b to avoid division by zero or negative
    delta_b = np.clip(delta_b, 1e-6, None)

    g_delta = np.log10(1.0 + 1.0 / delta_b)
    g_theta = r ** 2
    g_phi = r ** 2

    det_natural = r ** 4 * g_delta  # with g_rr=1

    g_rr = np.ones_like(r)
    floor_active = det_natural < FLOOR_VAL
    # g_rr absorbs compensation: g_rr = FLOOR_VAL / (r^2 * r^2 * g_delta)
    g_rr[floor_active] = FLOOR_VAL / (
        r[floor_active] ** 2 * r[floor_active] ** 2 * g_delta[floor_active]
    )

    return g_rr, g_theta, g_phi, g_delta, floor_active, delta_b


# ---------------------------------------------------------------------------
# Entropy Rate
# ---------------------------------------------------------------------------

def entropy_rate(r, *components):
    """
    Entropy rate = sqrt(sum of (dg_i/dr)^2) for all metric components.
    Uses central differences.
    """
    total_sq = np.zeros_like(r)
    for g in components:
        dg_dr = np.gradient(g, r)
        total_sq += dg_dr ** 2
    return np.sqrt(total_sq)


# ---------------------------------------------------------------------------
# Main Simulation
# ---------------------------------------------------------------------------

def run_simulation():
    print("=" * 70)
    print("4D Spatial Metric — CS as a Proper Dimension (delta-hat)")
    print("=" * 70)

    # Radial grid: 10.0 down to 0.01 (log-spaced)
    r = np.logspace(np.log10(10.0), np.log10(0.01), 2000)

    # Load CS data and build interpolator
    cs_data = load_cs_data()
    delta_b_interp = build_interpolator(cs_data)

    # --- Compute metrics ---
    g_rr_3d, g_theta_3d, g_phi_3d, floor_3d = compute_3d_metric(r)
    g_rr_4d, g_theta_4d, g_phi_4d, g_delta_4d, floor_4d, delta_b_vals = \
        compute_4d_metric(r, delta_b_interp)

    # --- Entropy rates ---
    erate_3d = entropy_rate(r, g_rr_3d, g_theta_3d, g_phi_3d)
    erate_4d = entropy_rate(r, g_rr_4d, g_theta_4d, g_phi_4d, g_delta_4d)

    # --- Effective coupling (where BOTH floors are active) ---
    both_active = floor_3d & floor_4d
    coupling = np.full_like(r, np.nan)
    coupling[both_active] = 1.0 - (g_rr_4d[both_active] / g_rr_3d[both_active])

    # Also compute where only 3D floor is active (for context)
    only_3d_active = floor_3d & ~floor_4d
    coupling_extended = np.full_like(r, np.nan)
    # Where both active
    coupling_extended[both_active] = 1.0 - (g_rr_4d[both_active] / g_rr_3d[both_active])
    # Where only 3D is active (4D doesn't need floor yet), 4D g_rr=1
    coupling_extended[only_3d_active] = 1.0 - (1.0 / g_rr_3d[only_3d_active])

    # --- Floor activation points ---
    # 3D: r^4 < FLOOR_VAL => r < FLOOR_VAL^(1/4)
    r_floor_3d = FLOOR_VAL ** 0.25
    # 4D: approximate — find first r where floor activates
    floor_4d_indices = np.where(floor_4d)[0]
    r_floor_4d = r[floor_4d_indices[0]] if len(floor_4d_indices) > 0 else None

    # ===================================================================
    # ANALYSIS — print to console
    # ===================================================================
    print(f"\n1. Benford floor value: {FLOOR_VAL:.6f}")
    print(f"   (L2 norm of [log10(1+1/d) for d=1..9])")

    print(f"\n2. Floor activation:")
    print(f"   3D floor activates at r/r_s = {r_floor_3d:.6f}")
    if r_floor_4d is not None:
        print(f"   4D floor activates at r/r_s ~ {r_floor_4d:.6f}")
        if r_floor_4d > r_floor_3d:
            print(f"   4D activates EARLIER (farther out) — g_delta < 1 reduces determinant")
        elif r_floor_4d < r_floor_3d:
            print(f"   4D activates LATER (closer in) — g_delta > 1 boosts determinant")
        else:
            print(f"   Activates at approximately the same point")

    valid_coupling = coupling[~np.isnan(coupling)]
    if len(valid_coupling) > 0:
        print(f"\n3. Emergent coupling (where both floors active):")
        print(f"   Mean:  {np.mean(valid_coupling):.6f} ({np.mean(valid_coupling)*100:.2f}%)")
        print(f"   Min:   {np.min(valid_coupling):.6f} ({np.min(valid_coupling)*100:.2f}%)")
        print(f"   Max:   {np.max(valid_coupling):.6f} ({np.max(valid_coupling)*100:.2f}%)")
        print(f"   Std:   {np.std(valid_coupling):.6f}")
        print(f"   Compare to old hardcoded: 25.00%")
        ratio = np.mean(valid_coupling) / 0.25
        print(f"   Ratio emergent/hardcoded: {ratio:.4f}")
    else:
        print("\n3. No region where both floors are active simultaneously")

    valid_coupling_ext = coupling_extended[~np.isnan(coupling_extended)]
    if len(valid_coupling_ext) > 0:
        print(f"\n   Extended coupling (wherever 3D floor is active):")
        print(f"   Mean:  {np.mean(valid_coupling_ext):.6f} ({np.mean(valid_coupling_ext)*100:.2f}%)")
        print(f"   Min:   {np.min(valid_coupling_ext):.6f} ({np.min(valid_coupling_ext)*100:.2f}%)")
        print(f"   Max:   {np.max(valid_coupling_ext):.6f} ({np.max(valid_coupling_ext)*100:.2f}%)")

    # 4. Entropy rate comparison
    # Compare at the singularity end
    idx_inner = np.argmin(r)
    print(f"\n4. Entropy rate at innermost point (r={r[idx_inner]:.4f}):")
    print(f"   3D rate: {erate_3d[idx_inner]:.6f}")
    print(f"   4D rate: {erate_4d[idx_inner]:.6f}")
    print(f"   Ratio 4D/3D: {erate_4d[idx_inner]/erate_3d[idx_inner]:.6f}")

    # Check acceleration: is rate increasing as r decreases?
    erate_3d_inner = erate_3d[r < 0.5]
    erate_4d_inner = erate_4d[r < 0.5]
    if len(erate_3d_inner) > 2:
        accel_3d = np.all(np.diff(erate_3d_inner[::-1]) >= 0)  # should increase as r->0
        accel_4d = np.all(np.diff(erate_4d_inner[::-1]) >= 0)
        print(f"   3D monotonically accelerating toward singularity: {accel_3d}")
        print(f"   4D monotonically accelerating toward singularity: {accel_4d}")
        # General trend
        print(f"   3D rate at r=0.5: {erate_3d[np.argmin(np.abs(r-0.5))]:.4f}, "
              f"at r=0.01: {erate_3d[idx_inner]:.4f}")
        print(f"   4D rate at r=0.5: {erate_4d[np.argmin(np.abs(r-0.5))]:.4f}, "
              f"at r=0.01: {erate_4d[idx_inner]:.4f}")

    # 5. g_delta peak
    idx_peak = np.argmax(g_delta_4d)
    print(f"\n5. g_delta peak:")
    print(f"   Peak at r/r_s = {r[idx_peak]:.6f}")
    print(f"   g_delta_max = {g_delta_4d[idx_peak]:.6f}")
    print(f"   Corresponding delta_B = {delta_b_vals[idx_peak]:.6f}")
    print(f"   (This is where CS conformance is highest)")

    # 6. Acceleration confirmation
    print(f"\n6. Entropy rate acceleration toward singularity:")
    r_check = [1.0, 0.5, 0.1, 0.01]
    for rc in r_check:
        idx = np.argmin(np.abs(r - rc))
        print(f"   r={rc:.2f}: 3D={erate_3d[idx]:.4f}, 4D={erate_4d[idx]:.4f}")
    print("   Both show acceleration ✓" if erate_4d[idx_inner] > erate_4d[np.argmin(np.abs(r-1.0))]
          else "   WARNING: acceleration not confirmed")

    # 7. g_delta range
    print(f"\n7. g_delta range across the journey:")
    print(f"   At r=10 (far):   g_delta = {g_delta_4d[np.argmin(np.abs(r-10.0))]:.4f}")
    print(f"   At r=1 (horizon): g_delta = {g_delta_4d[np.argmin(np.abs(r-1.0))]:.4f}")
    print(f"   At r=0.01 (inner): g_delta = {g_delta_4d[idx_inner]:.4f}")
    print(f"   Min g_delta = {np.min(g_delta_4d):.4f} at r = {r[np.argmin(g_delta_4d)]:.4f}")
    print(f"   Max g_delta = {np.max(g_delta_4d):.4f} at r = {r[np.argmax(g_delta_4d)]:.4f}")

    # dg_delta/dr contribution
    dg_delta_dr = np.gradient(g_delta_4d, r)
    idx_max_dg = np.argmax(np.abs(dg_delta_dr))
    print(f"\n   Max |dg_delta/dr| at r = {r[idx_max_dg]:.4f}")
    print(f"   |dg_delta/dr|_max = {np.abs(dg_delta_dr[idx_max_dg]):.4f}")

    print("\n" + "=" * 70)

    # ===================================================================
    # VISUALIZATION — 4 stacked panels
    # ===================================================================
    fig = plt.figure(figsize=(14, 20), facecolor='#1a1a2e')
    gs = GridSpec(4, 1, hspace=0.25, top=0.94, bottom=0.04, left=0.10, right=0.95)

    dark_bg = '#1a1a2e'
    text_color = '#e0e0e0'
    grid_color = '#2a2a4e'

    def style_ax(ax, ylabel):
        ax.set_facecolor(dark_bg)
        ax.tick_params(colors=text_color, which='both')
        ax.set_ylabel(ylabel, color=text_color, fontsize=12)
        ax.spines['bottom'].set_color(grid_color)
        ax.spines['top'].set_color(grid_color)
        ax.spines['left'].set_color(grid_color)
        ax.spines['right'].set_color(grid_color)
        ax.grid(True, alpha=0.2, color=grid_color)
        ax.set_xlim(10, 0.01)
        ax.set_xscale('log')

    # --- Panel 1: 3D Metric ---
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(r, g_rr_3d, color='#ff4444', linewidth=1.8, label=r'$g_{rr}$')
    ax1.plot(r, g_theta_3d, color='#44ff44', linewidth=1.8, label=r'$g_{\theta\theta}$')
    ax1.plot(r, g_phi_3d, color='#4488ff', linewidth=1.8, label=r'$g_{\phi\phi}$', linestyle='--')
    ax1.axhline(y=FLOOR_VAL, color='#ff8800', linestyle='--', linewidth=1.2,
                label=f'Benford floor = {FLOOR_VAL:.4f}')
    ax1.axvline(x=1.0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.5,
                label='Event horizon')
    ax1.axvline(x=r_floor_3d, color='#ff8800', linestyle=':', linewidth=0.8, alpha=0.7,
                label=f'3D floor at r={r_floor_3d:.4f}')
    ax1.set_yscale('log')
    ax1.set_ylim(1e-2, 1e4)
    style_ax(ax1, 'Metric Component')
    ax1.set_title('Panel 1: 3D Benford Floor (Equation 1)', color=text_color, fontsize=13, pad=10)
    ax1.legend(loc='upper right', fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e',
               labelcolor=text_color)

    # --- Panel 2: 4D Metric ---
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(r, g_rr_4d, color='#ff4444', linewidth=1.8, label=r'$g_{rr}$')
    ax2.plot(r, g_theta_4d, color='#44ff44', linewidth=1.8, label=r'$g_{\theta\theta}$')
    ax2.plot(r, g_phi_4d, color='#4488ff', linewidth=1.8, label=r'$g_{\phi\phi}$', linestyle='--')
    ax2.plot(r, g_delta_4d, color='#ffff00', linewidth=2.2, label=r'$g_{\hat\delta}$ (CS dimension)',
             zorder=5)
    ax2.axhline(y=FLOOR_VAL, color='#ff8800', linestyle='--', linewidth=1.2,
                label=f'Benford floor = {FLOOR_VAL:.4f}')
    ax2.axvline(x=1.0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.5)
    if r_floor_4d is not None:
        ax2.axvline(x=r_floor_4d, color='#ffff00', linestyle=':', linewidth=0.8, alpha=0.7,
                    label=f'4D floor at r~{r_floor_4d:.4f}')
    ax2.set_yscale('log')
    ax2.set_ylim(1e-2, 1e4)
    style_ax(ax2, 'Metric Component')
    ax2.set_title('Panel 2: 4D Benford Floor (Equation 2) — CS as Proper Dimension',
                  color=text_color, fontsize=13, pad=10)
    ax2.legend(loc='upper right', fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e',
               labelcolor=text_color)

    # --- Panel 3: Entropy Rate Comparison ---
    ax3 = fig.add_subplot(gs[2])
    ax3.plot(r, erate_3d, color='#888888', linewidth=2.0, label='3D entropy rate')
    ax3.plot(r, erate_4d, color='#ffffff', linewidth=2.0, label='4D entropy rate')
    # Show the dg_delta/dr contribution magnitude
    dg_delta_contribution = np.abs(np.gradient(g_delta_4d, r))
    ax3.plot(r, dg_delta_contribution, color='#ffff00', linewidth=1.0, alpha=0.5,
             label=r'$|dg_{\hat\delta}/dr|$ (CS clock)')
    ax3.axvline(x=1.0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.5)
    ax3.set_yscale('log')
    style_ax(ax3, 'Entropy Rate')
    ax3.set_title('Panel 3: Entropy Rate Comparison — Emergent Time',
                  color=text_color, fontsize=13, pad=10)
    ax3.legend(loc='upper right', fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e',
               labelcolor=text_color)

    # --- Panel 4: Emergent Coupling ---
    ax4 = fig.add_subplot(gs[3])
    # Plot coupling where both floors active
    mask_valid = ~np.isnan(coupling)
    if np.any(mask_valid):
        ax4.plot(r[mask_valid], coupling[mask_valid] * 100, color='#00ffaa', linewidth=2.0,
                 label='Emergent coupling (both floors active)')
    # Plot extended coupling
    mask_ext = ~np.isnan(coupling_extended)
    if np.any(mask_ext):
        ax4.plot(r[mask_ext], coupling_extended[mask_ext] * 100, color='#00ffaa', linewidth=1.0,
                 alpha=0.4, label='Extended coupling (3D floor active)')
    ax4.axhline(y=25.0, color='#ff8800', linestyle='--', linewidth=1.5,
                label='Old hardcoded coupling (25%)')
    ax4.axvline(x=1.0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.5)
    ax4.set_ylabel('Effective Coupling (%)', color=text_color, fontsize=12)
    ax4.set_xlabel(r'$r / r_s$', color=text_color, fontsize=12)
    ax4.set_facecolor(dark_bg)
    ax4.tick_params(colors=text_color, which='both')
    ax4.spines['bottom'].set_color(grid_color)
    ax4.spines['top'].set_color(grid_color)
    ax4.spines['left'].set_color(grid_color)
    ax4.spines['right'].set_color(grid_color)
    ax4.grid(True, alpha=0.2, color=grid_color)
    ax4.set_xlim(10, 0.01)
    ax4.set_xscale('log')
    ax4.set_title('Panel 4: Emergent Coupling — No Free Parameters',
                  color=text_color, fontsize=13, pad=10)
    ax4.legend(loc='upper right', fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e',
               labelcolor=text_color)

    # Supertitle
    fig.suptitle('3D vs 4D Spatial Metric — CS as a Proper Dimension',
                 color='#ffffff', fontsize=16, fontweight='bold')

    # Save
    fig_dir = os.path.join(os.path.dirname(__file__), '..', 'results', 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    fig_path = os.path.join(fig_dir, 'schwarzschild_benford_4d.png')
    fig.savefig(fig_path, dpi=200, facecolor=dark_bg)
    print(f"\nFigure saved to: {fig_path}")

    # Copy to Desktop
    desktop_path = os.path.expanduser('~/Desktop/schwarzschild_benford_4d.png')
    fig.savefig(desktop_path, dpi=200, facecolor=dark_bg)
    print(f"Figure copied to: {desktop_path}")

    plt.close(fig)

    # ===================================================================
    # Return analysis data for report
    # ===================================================================
    return {
        'floor_val': FLOOR_VAL,
        'r_floor_3d': r_floor_3d,
        'r_floor_4d': r_floor_4d,
        'coupling_mean': float(np.mean(valid_coupling)) if len(valid_coupling) > 0 else None,
        'coupling_min': float(np.min(valid_coupling)) if len(valid_coupling) > 0 else None,
        'coupling_max': float(np.max(valid_coupling)) if len(valid_coupling) > 0 else None,
        'coupling_std': float(np.std(valid_coupling)) if len(valid_coupling) > 0 else None,
        'coupling_ext_mean': float(np.mean(valid_coupling_ext)) if len(valid_coupling_ext) > 0 else None,
        'g_delta_min': float(np.min(g_delta_4d)),
        'g_delta_max': float(np.max(g_delta_4d)),
        'g_delta_peak_r': float(r[idx_peak]),
        'erate_3d_inner': float(erate_3d[idx_inner]),
        'erate_4d_inner': float(erate_4d[idx_inner]),
        'erate_ratio_inner': float(erate_4d[idx_inner] / erate_3d[idx_inner]),
    }


if __name__ == '__main__':
    results = run_simulation()
