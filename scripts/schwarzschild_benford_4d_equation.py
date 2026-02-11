#!/usr/bin/env python3
"""
4D Spatial Metric — CS (Causal Set) as a Proper Dimension (delta-hat)
VERSION 2: Using the |ln(r)| equation instead of 40 data points.

The CS equation:  δ_B(r) = 0.003389 × |ln(r/r_s)| + 0.002508

Compares:
  - 3D Benford floor (Equation 1)
  - 4D with data points (original)
  - 4D with equation (new)

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
# CS Data (original 40 points, for comparison)
# ---------------------------------------------------------------------------
CS_DATA = [
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

# ---------------------------------------------------------------------------
# The CS Equation: δ_B(r) = a × |ln(r)| + b
# ---------------------------------------------------------------------------
CS_A = 0.003389
CS_B = 0.002508

def delta_b_equation(r):
    """CS equation: δ_B = a × |ln(r)| + b"""
    return CS_A * np.abs(np.log(r)) + CS_B

def delta_b_data_interp(r):
    """Interpolator from the 40 data points (for comparison)."""
    rs = np.array([p[0] for p in CS_DATA])
    dbs = np.array([p[1] for p in CS_DATA])
    log_rs = np.log(rs)
    log_dbs = np.log(dbs)
    if HAS_SCIPY:
        interp_log = interp1d(log_rs, log_dbs, kind='linear', fill_value='extrapolate')
        return np.exp(interp_log(np.log(r)))
    else:
        # Simple nearest-neighbor fallback
        result = np.zeros_like(r, dtype=float)
        for i, ri in enumerate(r):
            idx = np.argmin(np.abs(rs - ri))
            result[i] = dbs[idx]
        return result

# ---------------------------------------------------------------------------
# Benford Floor
# ---------------------------------------------------------------------------
def compute_benford_floor():
    probs = [math.log10(1 + 1/d) for d in range(1, 10)]
    return math.sqrt(sum(p**2 for p in probs))

FLOOR_VAL = compute_benford_floor()

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def compute_3d_metric(r):
    g_theta = r**2
    g_phi = r**2
    det_natural = r**4
    g_rr = np.ones_like(r)
    floor_active = det_natural < FLOOR_VAL
    g_rr[floor_active] = FLOOR_VAL / (r[floor_active]**4)
    return g_rr, g_theta, g_phi, floor_active

def compute_4d_metric(r, delta_b_func):
    delta_b = delta_b_func(r)
    delta_b = np.clip(delta_b, 1e-6, None)
    g_delta = np.log10(1.0 + 1.0 / delta_b)
    g_theta = r**2
    g_phi = r**2
    det_natural = r**4 * g_delta
    g_rr = np.ones_like(r)
    floor_active = det_natural < FLOOR_VAL
    g_rr[floor_active] = FLOOR_VAL / (r[floor_active]**4 * g_delta[floor_active])
    return g_rr, g_theta, g_phi, g_delta, floor_active, delta_b

def entropy_rate(r, *components):
    total_sq = np.zeros_like(r)
    for g in components:
        dg_dr = np.gradient(g, r)
        total_sq += dg_dr**2
    return np.sqrt(total_sq)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run():
    print("=" * 70)
    print("4D Spatial Metric — Equation vs Data")
    print(f"CS Equation: δ_B(r) = {CS_A} × |ln(r)| + {CS_B}")
    print("=" * 70)

    r = np.logspace(np.log10(10.0), np.log10(0.01), 2000)

    # --- Three versions ---
    g_rr_3d, g_th_3d, g_ph_3d, floor_3d = compute_3d_metric(r)

    g_rr_data, g_th_data, g_ph_data, g_d_data, floor_data, db_data = \
        compute_4d_metric(r, delta_b_data_interp)

    g_rr_eq, g_th_eq, g_ph_eq, g_d_eq, floor_eq, db_eq = \
        compute_4d_metric(r, delta_b_equation)

    # --- Entropy rates ---
    erate_3d = entropy_rate(r, g_rr_3d, g_th_3d, g_ph_3d)
    erate_data = entropy_rate(r, g_rr_data, g_th_data, g_ph_data, g_d_data)
    erate_eq = entropy_rate(r, g_rr_eq, g_th_eq, g_ph_eq, g_d_eq)

    # --- Coupling ---
    both_floor = floor_3d & floor_eq
    coupling_eq = np.full_like(r, np.nan)
    coupling_eq[both_floor] = 1.0 - (g_rr_eq[both_floor] / g_rr_3d[both_floor])

    both_floor_data = floor_3d & floor_data
    coupling_data = np.full_like(r, np.nan)
    coupling_data[both_floor_data] = 1.0 - (g_rr_data[both_floor_data] / g_rr_3d[both_floor_data])

    # --- Analysis ---
    r_floor_3d = FLOOR_VAL**0.25
    fi_eq = np.where(floor_eq)[0]
    r_floor_eq = r[fi_eq[0]] if len(fi_eq) > 0 else None
    fi_data = np.where(floor_data)[0]
    r_floor_data = r[fi_data[0]] if len(fi_data) > 0 else None

    print(f"\nBenford floor: {FLOOR_VAL:.6f}")
    print(f"3D floor at:   r = {r_floor_3d:.4f}")
    if r_floor_data: print(f"4D(data) at:   r ~ {r_floor_data:.4f}")
    if r_floor_eq:   print(f"4D(equation) at: r ~ {r_floor_eq:.4f}")

    vc_data = coupling_data[~np.isnan(coupling_data)]
    vc_eq = coupling_eq[~np.isnan(coupling_eq)]

    print(f"\nEmergent coupling (data):     {np.mean(vc_data)*100:.2f}% mean" if len(vc_data) > 0 else "")
    print(f"Emergent coupling (equation): {np.mean(vc_eq)*100:.2f}% mean" if len(vc_eq) > 0 else "")

    # How close is equation to data?
    diff_grr = np.abs(g_rr_eq - g_rr_data)
    diff_gd = np.abs(g_d_eq - g_d_data)
    diff_erate = np.abs(erate_eq - erate_data)

    print(f"\n--- Equation vs Data Agreement ---")
    print(f"g_rr mean |diff|:      {np.mean(diff_grr):.4f}")
    print(f"g_delta mean |diff|:   {np.mean(diff_gd):.4f}")
    print(f"entropy rate mean |diff|: {np.mean(diff_erate):.4f}")

    # Relative agreement where it matters (inside horizon)
    inside = r < 1.0
    if np.any(inside):
        rel_grr = np.mean(diff_grr[inside] / np.maximum(g_rr_data[inside], 1e-10)) * 100
        rel_gd = np.mean(diff_gd[inside] / np.maximum(g_d_data[inside], 1e-10)) * 100
        rel_erate = np.mean(diff_erate[inside] / np.maximum(erate_data[inside], 1e-10)) * 100
        print(f"\nInside horizon (r < 1):")
        print(f"  g_rr relative diff:      {rel_grr:.2f}%")
        print(f"  g_delta relative diff:   {rel_gd:.2f}%")
        print(f"  entropy rate rel diff:   {rel_erate:.2f}%")

    # Data points overlay
    r_pts = np.array([p[0] for p in CS_DATA])
    db_pts = np.array([p[1] for p in CS_DATA])
    gd_pts = np.log10(1.0 + 1.0 / db_pts)

    print("\n" + "=" * 70)

    # ===================================================================
    # VISUALIZATION — 4 panels
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
        for s in ax.spines.values(): s.set_color(gc)
        ax.grid(True, alpha=0.2, color=gc)
        ax.set_xlim(10, 0.01)
        ax.set_xscale('log')

    # --- Panel 1: All metric components (equation version) ---
    ax1 = fig.add_subplot(gs[0])
    # Radial: 3D vs 4D-equation
    ax1.plot(r, g_rr_3d, color='#ff4444', linewidth=1.5, alpha=0.4, label=r'$g_{rr}$ 3D')
    ax1.plot(r, g_rr_eq, color='#ff4444', linewidth=2.0, label=r'$g_{rr}$ 4D+CS (equation)')
    # Angular (same in both, just r²)
    ax1.plot(r, g_th_eq, color='#44ff44', linewidth=2.0, label=r'$g_{\theta\theta}$ = $g_{\phi\phi}$ = $r^2$')
    # CS dimension
    ax1.plot(r, g_d_eq, color='#ffff00', linewidth=2.0, label=r'$g_{\hat\delta}$ (CS)')
    ax1.axvline(x=1.0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.4, label='Horizon')
    ax1.axhline(y=FLOOR_VAL, color='#ff8800', linestyle=':', linewidth=1, alpha=0.5,
                label=f'Floor = {FLOOR_VAL:.4f}')
    ax1.set_yscale('log')
    ax1.set_ylim(1e-4, 1e9)
    style(ax1, 'Metric Component')
    ax1.set_title('Panel 1: All 4 Dimensions — The Full Metric', color=tc, fontsize=13, pad=10)
    ax1.legend(fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e', labelcolor=tc, loc='upper right')

    # Mark crossing points
    # Find where g_rr_eq crosses g_theta (angular)
    diff_rr_th = g_rr_eq - g_th_eq
    cross_indices = np.where(np.diff(np.sign(diff_rr_th)))[0]
    for ci in cross_indices:
        ax1.plot(r[ci], g_rr_eq[ci], 'o', color='#ffffff', markersize=8, zorder=10)
        ax1.annotate(f'r={r[ci]:.3f}', xy=(r[ci], g_rr_eq[ci]),
                     xytext=(10, 15), textcoords='offset points',
                     color='#ffffff', fontsize=9,
                     arrowprops=dict(arrowstyle='->', color='#ffffff', lw=0.8))

    # Find where g_delta crosses g_theta (angular)
    diff_gd_th = g_d_eq - g_th_eq
    cross_gd = np.where(np.diff(np.sign(diff_gd_th)))[0]
    for ci in cross_gd:
        ax1.plot(r[ci], g_d_eq[ci], 'o', color='#ffff00', markersize=8, zorder=10)
        ax1.annotate(f'r={r[ci]:.3f}', xy=(r[ci], g_d_eq[ci]),
                     xytext=(10, -20), textcoords='offset points',
                     color='#ffff00', fontsize=9,
                     arrowprops=dict(arrowstyle='->', color='#ffff00', lw=0.8))

    # --- Panel 2: g_delta comparison ---
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(r, g_d_data, color='#ffdd44', linewidth=1.8, label=r'$g_{\hat\delta}$ (data interpolation)')
    ax2.plot(r, g_d_eq, color='#00ffaa', linewidth=2.0, label=r'$g_{\hat\delta}$ (equation)', linestyle='--')
    ax2.scatter(r_pts, gd_pts, color='#ffdd44', s=25, zorder=5, edgecolors='#aa9900',
                linewidths=0.5, label='40 measured points')
    ax2.axvline(x=1.0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.4)
    style(ax2, r'$g_{\hat\delta}$')
    ax2.set_ylim(1.0, 3.0)
    ax2.set_title('Panel 2: CS Dimension — Data vs Equation', color=tc, fontsize=13, pad=10)
    ax2.legend(fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e', labelcolor=tc)

    # --- Panel 3: Entropy rate ---
    ax3 = fig.add_subplot(gs[2])
    ax3.plot(r, erate_3d, color='#888888', linewidth=1.5, label='3D entropy rate')
    ax3.plot(r, erate_data, color='#ff4444', linewidth=1.8, label='4D entropy rate (data)', alpha=0.7)
    ax3.plot(r, erate_eq, color='#00ffaa', linewidth=1.8, label='4D entropy rate (equation)', linestyle='--')
    ax3.axvline(x=1.0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.4)
    ax3.set_yscale('log')
    style(ax3, 'Entropy Rate')
    ax3.set_title('Panel 3: Entropy Rate — Data vs Equation', color=tc, fontsize=13, pad=10)
    ax3.legend(fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e', labelcolor=tc)

    # --- Panel 4: Emergent coupling ---
    ax4 = fig.add_subplot(gs[3])
    m_data = ~np.isnan(coupling_data)
    m_eq = ~np.isnan(coupling_eq)
    if np.any(m_data):
        ax4.plot(r[m_data], coupling_data[m_data]*100, color='#ff4444', linewidth=1.8,
                 label='Coupling (data)', alpha=0.7)
    if np.any(m_eq):
        ax4.plot(r[m_eq], coupling_eq[m_eq]*100, color='#00ffaa', linewidth=1.8,
                 label='Coupling (equation)', linestyle='--')
    ax4.axvline(x=1.0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.4)
    style(ax4, 'Coupling (%)')
    ax4.set_xlabel(r'$r / r_s$', color=tc, fontsize=12)
    ax4.set_title('Panel 4: Emergent Coupling — Data vs Equation', color=tc, fontsize=13, pad=10)
    ax4.legend(fontsize=9, facecolor='#2a2a4e', edgecolor='#3a3a5e', labelcolor=tc)

    fig.suptitle('4D Metric: Does the Equation Match the Data?',
                 color='#ffffff', fontsize=16, fontweight='bold')

    out = os.path.expanduser('~/Desktop/schwarzschild_benford_4d_equation.png')
    fig.savefig(out, dpi=200, facecolor=dark_bg)
    print(f"\nFigure saved to {out}")
    plt.close(fig)

if __name__ == '__main__':
    run()
