#!/usr/bin/env python3
"""
Test different Lp norms of the Benford probability vector as floor values.
Goal: determine if L2 (0.4068) is special, or if other norms also produce
reasonable black hole interior geometry.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Benford probability vector ──
digits = np.arange(1, 10)
P = np.log10(1 + 1/digits)

print("=" * 70)
print("Benford Probability Vector")
print("=" * 70)
for d, p in zip(digits, P):
    print(f"  P({d}) = {p:.6f}")

# ── Compute Lp norms for p = 1 through 10, plus infinity ──
print("\n" + "=" * 70)
print("Lp Norms of the Benford Probability Vector")
print("=" * 70)

norms = {}
for p_val in range(1, 11):
    lp = np.sum(P**p_val) ** (1/p_val)
    norms[f"L{p_val}"] = lp
    print(f"  L{p_val:>2d} = {lp:.6f}")

linf = np.max(P)
norms["L∞"] = linf
print(f"  L ∞ = {linf:.6f}")

# ── CS equation (same for all tests) ──
def delta_B(r):
    return 0.003389 * np.abs(np.log(r)) + 0.002508

def g_delta(r):
    db = delta_B(r)
    return np.log10(1 + 1/db)

# ── Run the 5D metric with each norm as floor ──
r = np.logspace(np.log10(0.01), np.log10(10), 2000)

# Select norms to test: L1, L2, L3, L4, L5, L∞
test_norms = ["L1", "L2", "L3", "L4", "L5", "L∞"]

print("\n" + "=" * 70)
print("5D Metric Results for Each Floor Value")
print("=" * 70)

results = {}
for name in test_norms:
    floor_val = norms[name]

    gd = g_delta(r)
    gtheta = r**2
    gphi = r**2
    gtime = -(1 - 1/r)

    det_spatial = r**4 * gd
    floor_active = det_spatial < floor_val

    grr = np.ones_like(r)
    grr[floor_active] = floor_val / (r[floor_active]**4 * gd[floor_active])

    # Find where floor activates
    activation_idx = np.where(floor_active)[0]
    if len(activation_idx) > 0:
        activation_r = r[activation_idx[-1]]
    else:
        activation_r = None

    # Entropy rate (all 5 components)
    dr = r * 0.001
    rlo = r - dr
    rhi = r + dr

    gd_lo = g_delta(rlo)
    gd_hi = g_delta(rhi)

    det_lo = rlo**4 * gd_lo
    det_hi = rhi**4 * gd_hi

    grr_lo = np.where(det_lo < floor_val, floor_val / (rlo**4 * gd_lo), 1.0)
    grr_hi = np.where(det_hi < floor_val, floor_val / (rhi**4 * gd_hi), 1.0)

    gtime_lo = -(1 - 1/rlo)
    gtime_hi = -(1 - 1/rhi)

    inv = 1 / (2 * dr)
    dgrr = (grr_hi - grr_lo) * inv
    dgtheta = (rhi**2 - rlo**2) * inv
    dgphi = dgtheta
    dgd = (gd_hi - gd_lo) * inv
    dgtime = (gtime_hi - gtime_lo) * inv

    entropy = np.sqrt(dgrr**2 + dgtheta**2 + dgphi**2 + dgd**2 + dgtime**2)

    # Max grr inside horizon (how much radial growth)
    inside = r < 1.0
    max_grr = np.max(grr[inside]) if np.any(inside) else 0
    max_entropy = np.max(entropy[inside]) if np.any(inside) else 0
    grr_at_001 = grr[np.argmin(np.abs(r - 0.01))]
    entropy_at_001 = entropy[np.argmin(np.abs(r - 0.01))]

    results[name] = {
        'floor': floor_val,
        'grr': grr,
        'entropy': entropy,
        'floor_active': floor_active,
        'activation_r': activation_r,
        'max_grr': max_grr,
        'max_entropy': max_entropy,
        'grr_at_001': grr_at_001,
        'entropy_at_001': entropy_at_001,
        'gtime': gtime,
        'gd': gd,
    }

    act_str = f"r = {activation_r:.4f} r_s" if activation_r else "never"
    print(f"\n  {name} (floor = {floor_val:.6f}):")
    print(f"    Floor activates at:  {act_str}")
    print(f"    g_rr at r=0.01:      {grr_at_001:.4f}")
    print(f"    Max g_rr inside:     {max_grr:.4f}")
    print(f"    Entropy at r=0.01:   {entropy_at_001:.4f}")
    print(f"    Max entropy inside:  {max_entropy:.4f}")

# ── Visualization: 6-panel comparison ──
fig, axes = plt.subplots(2, 3, figsize=(24, 15), facecolor='#1a1a2e')
fig.suptitle('Benford Floor: Which Lp Norm?',
             color='white', fontsize=22, fontweight='bold', y=0.98)

colors = {
    'L1': '#ff4444',
    'L2': '#00ff88',
    'L3': '#44aaff',
    'L4': '#ffaa00',
    'L5': '#cc44ff',
    'L∞': '#ff8888',
}

for idx, name in enumerate(test_norms):
    ax = axes.flat[idx]
    ax.set_facecolor('#1a1a2e')
    ax.tick_params(colors='#e0e0e0')
    for spine in ax.spines.values():
        spine.set_color('#333')

    res = results[name]

    # Plot g_rr
    ax.plot(r, res['grr'], color=colors[name], linewidth=2, label=f'g_rr')

    # Plot entropy rate (scaled to fit)
    entropy_scaled = res['entropy'] / max(res['max_entropy'], 1) * res['max_grr'] * 0.5
    ax.plot(r, entropy_scaled, color='#00ffff', linewidth=1, alpha=0.6, label='entropy (scaled)')

    # Horizon line
    ax.axvline(x=1.0, color='white', linestyle='--', linewidth=0.8, alpha=0.4)

    # Floor activation
    if res['activation_r']:
        ax.axvline(x=res['activation_r'], color=colors[name], linestyle=':',
                   linewidth=1, alpha=0.6)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.01, 10)
    ax.set_ylim(0.1, max(res['max_grr'] * 2, 10))
    ax.set_xlabel('r / r_s', color='#e0e0e0')
    ax.set_ylabel('g_rr', color='#e0e0e0')

    title = f'{name}  (floor = {res["floor"]:.4f})'
    if name == 'L2':
        title += '  ← CURRENT'
    ax.set_title(title, color=colors[name], fontsize=13, fontweight='bold')

    # Stats text — bottom right to avoid overlap with legend
    act_str = f'{res["activation_r"]:.3f}' if res['activation_r'] else 'never'
    stats = f'g_rr(0.01) = {res["grr_at_001"]:.1f}\nactivates: {act_str}'
    ax.text(0.97, 0.03, stats, transform=ax.transAxes, fontsize=9,
            color='#cccccc', ha='right', va='bottom', family='monospace',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#1a1a2e',
                     edgecolor='#444', alpha=0.9))

    ax.legend(fontsize=9, facecolor='#22223a', edgecolor='#444',
             labelcolor='#e0e0e0', loc='upper left')
    ax.grid(True, alpha=0.1, color='#444')

plt.tight_layout(rect=[0, 0.05, 1, 0.95])

# Add summary text at bottom
summary = "  |  ".join([f"{name}: floor={norms[name]:.4f}" for name in test_norms])
fig.text(0.5, 0.015, summary, ha='center', fontsize=10, color='#888888', family='monospace')

fig.savefig('/home/jackwayne/Desktop/benford_lp_norm_test.png', dpi=180,
            facecolor='#1a1a2e', bbox_inches='tight', pad_inches=0.3)
plt.close()

# ── Summary comparison ──
print("\n" + "=" * 70)
print("COMPARISON SUMMARY")
print("=" * 70)
print(f"{'Norm':<6} {'Floor':>8} {'Activates':>12} {'g_rr(0.01)':>12} {'Entropy(0.01)':>14} {'Interior'}")
print("-" * 70)
for name in test_norms:
    res = results[name]
    act_str = f"{res['activation_r']:.4f}" if res['activation_r'] else "never"

    # Classify interior behavior
    if res['grr_at_001'] < 2:
        interior = "FLAT (dead)"
    elif res['grr_at_001'] < 100:
        interior = "mild growth"
    elif res['grr_at_001'] < 10000:
        interior = "moderate growth"
    else:
        interior = "strong growth"

    print(f"{name:<6} {res['floor']:>8.4f} {act_str:>12} {res['grr_at_001']:>12.2f} {res['entropy_at_001']:>14.2f} {interior}")

print("\n" + "=" * 70)
print("Figure saved: ~/Desktop/benford_lp_norm_test.png")
print("=" * 70)
