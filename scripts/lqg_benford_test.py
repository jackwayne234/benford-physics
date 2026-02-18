#!/usr/bin/env python3
"""
Test: Do Loop Quantum Gravity eigenvalues follow Benford's Law?

LQG area eigenvalues: A_j = 8πγl_P² × √(j(j+1))
where j = 0, 1/2, 1, 3/2, 2, ... (spin quantum numbers)
γ = Barbero-Immirzi parameter ≈ 0.2375 (fixed by black hole entropy)
l_P = Planck length ≈ 1.616e-35 m

We also test:
- Single-spin eigenvalues: √(j(j+1)) for many j values
- Multi-spin (composite) areas: sums of √(j_i(j_i+1)) for random combinations
- Volume eigenvalues (more complex, approximate)

The key question: does the spectrum have Benford structure?
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter

# ── Benford reference ──
def benford_expected(d):
    return np.log10(1 + 1/d)

benford_probs = np.array([benford_expected(d) for d in range(1, 10)])

def first_digit(x):
    """Extract first significant digit of a positive number."""
    if x <= 0:
        return None
    s = f"{x:.15e}"
    return int(s[0])

def benford_analysis(values, label):
    """Run full Benford analysis on a set of values."""
    values = np.array([v for v in values if v > 0])
    n = len(values)

    # Extract first digits
    digits = [first_digit(v) for v in values]
    digits = [d for d in digits if d is not None and 1 <= d <= 9]

    # Count
    counts = Counter(digits)
    total = len(digits)

    # Observed distribution
    observed = np.array([counts.get(d, 0) / total for d in range(1, 10)])
    expected = benford_probs

    # Per-digit deviation
    epsilon = observed - expected

    # δ_B (L2 deviation)
    delta_b = np.sqrt(np.sum(epsilon**2))

    # MAD
    mad = np.mean(np.abs(epsilon))

    # Chi-squared
    expected_counts = np.array([total * benford_expected(d) for d in range(1, 10)])
    observed_counts = np.array([counts.get(d, 0) for d in range(1, 10)])
    chi2 = np.sum((observed_counts - expected_counts)**2 / expected_counts)

    # Classification
    if mad < 0.006:
        mad_class = "Close conformity"
    elif mad < 0.012:
        mad_class = "Acceptable conformity"
    elif mad < 0.015:
        mad_class = "Marginally acceptable"
    else:
        mad_class = "Nonconformity"

    if delta_b < 0.02:
        db_class = "Strong conformance"
    elif delta_b < 0.05:
        db_class = "Moderate conformance"
    elif delta_b < 0.10:
        db_class = "Weak conformance"
    else:
        db_class = "Significant deviation"

    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"  n = {n} values, {total} valid digits")
    print(f"{'='*60}")
    print(f"\n  Digit  Observed  Expected  Deviation")
    print(f"  {'─'*45}")
    for d in range(1, 10):
        print(f"    {d}     {observed[d-1]:.4f}    {expected[d-1]:.4f}   {epsilon[d-1]:+.4f}")

    print(f"\n  δ_B (L2 deviation) = {delta_b:.6f}  → {db_class}")
    print(f"  MAD               = {mad:.6f}  → {mad_class}")
    print(f"  χ² (df=8)         = {chi2:.2f}")

    return {
        'label': label,
        'observed': observed,
        'expected': expected,
        'epsilon': epsilon,
        'delta_b': delta_b,
        'mad': mad,
        'chi2': chi2,
        'mad_class': mad_class,
        'db_class': db_class,
        'n': n,
        'digits': digits,
    }


# ══════════════════════════════════════════════════════════════
# TEST 1: Single-spin area eigenvalues √(j(j+1))
# ══════════════════════════════════════════════════════════════

# Generate j values: 1/2, 1, 3/2, 2, ..., up to j_max
j_max = 5000
j_values = np.arange(0.5, j_max + 0.5, 0.5)  # half-integer steps
single_spin = np.sqrt(j_values * (j_values + 1))

result1 = benford_analysis(single_spin, f"TEST 1: Single-spin √(j(j+1)), j up to {j_max}")

# ══════════════════════════════════════════════════════════════
# TEST 2: Full area eigenvalues with physical constants
# A = 8πγl_P² × √(j(j+1))
# ══════════════════════════════════════════════════════════════

gamma = 0.2375  # Barbero-Immirzi parameter
lp = 1.616e-35  # Planck length in meters
prefactor = 8 * np.pi * gamma * lp**2

full_area = prefactor * single_spin

result2 = benford_analysis(full_area, f"TEST 2: Full area 8πγl_P²√(j(j+1)), j up to {j_max}")

# ══════════════════════════════════════════════════════════════
# TEST 3: Composite areas (sums of spins — realistic surfaces)
# A real surface is punctured by many spin network edges,
# each contributing √(j_i(j_i+1)) to the total area.
# ══════════════════════════════════════════════════════════════

np.random.seed(42)
n_surfaces = 10000
composite_areas = []

for _ in range(n_surfaces):
    # Random number of punctures (2 to 50)
    n_punctures = np.random.randint(2, 51)
    # Random spins (half-integers from 1/2 to 10)
    spins = np.random.choice(np.arange(0.5, 10.5, 0.5), size=n_punctures)
    # Total area
    area = prefactor * np.sum(np.sqrt(spins * (spins + 1)))
    composite_areas.append(area)

composite_areas = np.array(composite_areas)

result3 = benford_analysis(composite_areas, f"TEST 3: Composite areas (sums of {n_surfaces} random surfaces)")

# ══════════════════════════════════════════════════════════════
# TEST 4: Area ratios (how areas compare to each other)
# Scale-invariant quantities → strongest Benford signal
# ══════════════════════════════════════════════════════════════

# Ratios of consecutive eigenvalues
area_ratios = full_area[1:] / full_area[:-1]

result4 = benford_analysis(area_ratios, "TEST 4: Consecutive area eigenvalue ratios A(j+1/2)/A(j)")

# ══════════════════════════════════════════════════════════════
# TEST 5: Products of spins (multiplicative process)
# ══════════════════════════════════════════════════════════════

np.random.seed(123)
n_products = 10000
spin_products = []
for _ in range(n_products):
    n_factors = np.random.randint(2, 20)
    spins = np.random.choice(np.arange(0.5, 20.5, 0.5), size=n_factors)
    product = np.prod(np.sqrt(spins * (spins + 1)))
    spin_products.append(product)

spin_products = np.array(spin_products)

result5 = benford_analysis(spin_products, "TEST 5: Products of spin eigenvalues (multiplicative)")

# ══════════════════════════════════════════════════════════════
# TEST 6: The bare sequence j(j+1) — what the quantum numbers give you
# ══════════════════════════════════════════════════════════════

j_squared = j_values * (j_values + 1)

result6 = benford_analysis(j_squared, f"TEST 6: Bare j(j+1) values, j up to {j_max}")

# ══════════════════════════════════════════════════════════════
# VISUALIZATION
# ══════════════════════════════════════════════════════════════

all_results = [result1, result2, result3, result4, result5, result6]

fig, axes = plt.subplots(2, 3, figsize=(24, 15), facecolor='#1a1a2e')
fig.suptitle('Do LQG Eigenvalues Follow Benford\'s Law?',
             color='white', fontsize=22, fontweight='bold', y=0.98)

colors_bars = ['#ff6666', '#ff9944', '#ffcc44', '#88ff44', '#44ffaa',
               '#44ddff', '#4488ff', '#8844ff', '#cc44ff']

for idx, (ax, res) in enumerate(zip(axes.flat, all_results)):
    ax.set_facecolor('#1a1a2e')
    ax.tick_params(colors='#e0e0e0')
    for spine in ax.spines.values():
        spine.set_color('#333')

    x = np.arange(1, 10)
    width = 0.35

    # Observed bars
    bars = ax.bar(x - width/2, res['observed'], width, color=colors_bars,
                  alpha=0.8, label='Observed')

    # Expected bars (Benford)
    ax.bar(x + width/2, res['expected'], width, color='white',
           alpha=0.3, edgecolor='white', linewidth=1, label='Benford expected')

    ax.set_xlabel('First Digit', color='#e0e0e0')
    ax.set_ylabel('Frequency', color='#e0e0e0')
    ax.set_xticks(x)
    ax.set_ylim(0, 0.4)
    ax.grid(True, alpha=0.1, color='#444', axis='y')

    # Title with result
    db_color = '#00ff88' if res['delta_b'] < 0.02 else '#ffaa00' if res['delta_b'] < 0.05 else '#ff4444'
    ax.set_title(f"{res['label'].split(':')[0]}: δ_B = {res['delta_b']:.4f}",
                 color=db_color, fontsize=11, fontweight='bold')

    # Legend in upper left, stats in upper right (no overlap)
    ax.legend(fontsize=8, facecolor='#22223a', edgecolor='#444',
             labelcolor='#e0e0e0', loc='upper left')

    # Stats box — below the legend area
    stats = f"{res['db_class']}"
    ax.text(0.97, 0.85, stats, transform=ax.transAxes, fontsize=9,
            color='#cccccc', ha='right', va='top', family='monospace',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#1a1a2e',
                     edgecolor='#444', alpha=0.9))

plt.tight_layout(rect=[0, 0.05, 1, 0.95])

# Summary at bottom
summary_parts = []
for res in all_results:
    name = res['label'].split(':')[0]
    symbol = '✓' if res['delta_b'] < 0.05 else '✗'
    summary_parts.append(f"{name}: δ_B={res['delta_b']:.4f} {symbol}")
summary = "  |  ".join(summary_parts)
fig.text(0.5, 0.015, summary, ha='center', fontsize=9, color='#888888', family='monospace')

fig.savefig('/home/jackwayne/Desktop/lqg_benford_test.png', dpi=180,
            facecolor='#1a1a2e', bbox_inches='tight', pad_inches=0.3)
plt.close()

# ══════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("  FINAL SUMMARY: LQG Eigenvalues vs Benford's Law")
print("=" * 60)
print(f"\n  {'Test':<12} {'δ_B':>8} {'MAD':>8} {'Verdict'}")
print(f"  {'─'*50}")
for res in all_results:
    name = res['label'].split(':')[0]
    verdict = res['db_class']
    print(f"  {name:<12} {res['delta_b']:>8.4f} {res['mad']:>8.4f} {verdict}")

print(f"\n  Benford conformance thresholds:")
print(f"    δ_B < 0.02  → Strong conformance")
print(f"    δ_B < 0.05  → Moderate conformance")
print(f"    δ_B < 0.10  → Weak conformance")
print(f"    δ_B ≥ 0.10  → Significant deviation")

print(f"\n  Figure saved: ~/Desktop/lqg_benford_test.png")
print("=" * 60)
