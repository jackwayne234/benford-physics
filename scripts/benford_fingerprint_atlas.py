#!/usr/bin/env python3
"""
Benford Fingerprint Atlas
Compare per-digit deviation patterns across quantum statistics and LQG.
Goal: see if there's a pattern in HOW things conform or deviate.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Benford reference ──
benford = np.array([np.log10(1 + 1/d) for d in range(1, 10)])

# ══════════════════════════════════════════════════════════════
# DATA: Per-digit deviations ε(d) from existing results
# ══════════════════════════════════════════════════════════════

systems = {}

# Quantum statistics (from results/individual/ — 100k data points each)
systems['Bose-Einstein'] = {
    'epsilon': np.array([-0.00433, 0.003189, 0.001411, 0.00062, 0.000159,
                         -6.7e-05, -0.000232, -0.000363, -0.000387]),
    'delta_b': 0.005627,
    'type': 'Bosonic',
    'color': '#00ff88',
    'note': 'Completely monotonic → conforms',
}

systems['Fermi-Dirac'] = {
    'epsilon': np.array([-0.00659, 0.005209, 0.004961, 0.00515, -0.002391,
                         -0.001957, -0.001662, -0.001463, -0.001257]),
    'delta_b': 0.011736,
    'type': 'Fermionic',
    'color': '#ff4444',
    'note': 'Alternating coefficients → oscillatory deviation',
}

systems['Maxwell-Boltzmann'] = {
    'epsilon': np.array([-0.00917, 0.002319, 0.001641, 0.00128, 0.001039,
                         0.000873, 0.000768, 0.000667, 0.000583]),
    'delta_b': 0.009850,
    'type': 'Classical',
    'color': '#ffaa00',
    'note': 'Single exponential → monotonic deviation',
}

systems['Planck Spectrum'] = {
    'epsilon': np.array([0.02233, -0.014891, -0.006039, -0.00343, -0.001751,
                         -0.000457, 0.000468, 0.001467, 0.002303]),
    'delta_b': 0.027921,
    'type': 'Bosonic + prefactor',
    'color': '#ff88ff',
    'note': 'BE denominator masked by ν³ prefactor',
}

# LQG results (from today's test — regenerate the deviations)
def first_digit(x):
    if x <= 0: return None
    return int(f"{x:.15e}"[0])

def get_epsilon(values):
    digits = [first_digit(v) for v in values if v > 0]
    digits = [d for d in digits if d is not None and 1 <= d <= 9]
    total = len(digits)
    from collections import Counter
    counts = Counter(digits)
    observed = np.array([counts.get(d, 0) / total for d in range(1, 10)])
    return observed - benford

# Generate LQG data
j_values = np.arange(0.5, 5000.5, 0.5)
gamma = 0.2375
lp = 1.616e-35
prefactor = 8 * np.pi * gamma * lp**2

single_spin = np.sqrt(j_values * (j_values + 1))
full_area = prefactor * single_spin

np.random.seed(42)
composite = []
for _ in range(10000):
    n_p = np.random.randint(2, 51)
    spins = np.random.choice(np.arange(0.5, 10.5, 0.5), size=n_p)
    composite.append(prefactor * np.sum(np.sqrt(spins * (spins + 1))))
composite = np.array(composite)

np.random.seed(123)
products = []
for _ in range(10000):
    n_f = np.random.randint(2, 20)
    spins = np.random.choice(np.arange(0.5, 20.5, 0.5), size=n_f)
    products.append(np.prod(np.sqrt(spins * (spins + 1))))
products = np.array(products)

j_squared = j_values * (j_values + 1)

systems['LQG Single Spin'] = {
    'epsilon': get_epsilon(single_spin),
    'delta_b': 0.2035,
    'type': 'Quantum Gravity',
    'color': '#4488ff',
    'note': 'Additive sequence → deviates',
}

systems['LQG Full Area'] = {
    'epsilon': get_epsilon(full_area),
    'delta_b': 0.2101,
    'type': 'Quantum Gravity',
    'color': '#44aaff',
    'note': 'Physical constants × single spin',
}

systems['LQG Composite'] = {
    'epsilon': get_epsilon(composite),
    'delta_b': 0.1820,
    'type': 'Quantum Gravity',
    'color': '#44ccff',
    'note': 'Sums of random spins → still deviates',
}

systems['LQG Products'] = {
    'epsilon': get_epsilon(products),
    'delta_b': 0.0037,
    'type': 'Quantum Gravity (mult)',
    'color': '#00ffcc',
    'note': 'PRODUCTS of spins → strong Benford!',
}

systems['LQG j(j+1)'] = {
    'epsilon': get_epsilon(j_squared),
    'delta_b': 0.1027,
    'type': 'Quantum Gravity',
    'color': '#6688ff',
    'note': 'Bare quantum numbers',
}

# ══════════════════════════════════════════════════════════════
# FIGURE 1: Fingerprint Atlas — all deviation patterns
# ══════════════════════════════════════════════════════════════

names = list(systems.keys())
n_sys = len(names)

fig, axes = plt.subplots(3, 3, figsize=(24, 18), facecolor='#1a1a2e')
fig.suptitle('Benford Fingerprint Atlas\nPer-Digit Deviation Patterns ε(d) Across Quantum Systems',
             color='white', fontsize=20, fontweight='bold', y=0.98)

for idx in range(9):
    ax = axes.flat[idx]
    ax.set_facecolor('#1a1a2e')
    ax.tick_params(colors='#e0e0e0')
    for spine in ax.spines.values():
        spine.set_color('#333')

    if idx < n_sys:
        name = names[idx]
        sys = systems[name]
        eps = sys['epsilon']

        # Bar colors: green for positive, red for negative
        colors = ['#00ff88' if e >= 0 else '#ff4466' for e in eps]
        bars = ax.bar(range(1, 10), eps, color=colors, alpha=0.8, edgecolor='white', linewidth=0.5)

        ax.axhline(y=0, color='white', linewidth=0.5, alpha=0.5)
        ax.set_xticks(range(1, 10))
        ax.set_xlabel('Digit', color='#e0e0e0', fontsize=9)
        ax.set_ylabel('ε(d)', color='#e0e0e0', fontsize=9)

        # Set consistent y-limits based on max deviation
        max_dev = max(abs(eps.max()), abs(eps.min()))
        if max_dev < 0.01:
            ax.set_ylim(-0.012, 0.012)
        elif max_dev < 0.03:
            ax.set_ylim(-0.035, 0.035)
        else:
            ax.set_ylim(-max_dev * 1.3, max_dev * 1.3)

        # Color title by δ_B
        db = sys['delta_b']
        if db < 0.02:
            title_color = '#00ff88'
        elif db < 0.05:
            title_color = '#ffaa00'
        else:
            title_color = '#ff4444'

        ax.set_title(f"{name}\nδ_B = {db:.4f}  |  {sys['type']}",
                     color=title_color, fontsize=11, fontweight='bold')

        # Pattern description — positioned above x-axis label with padding
        signs = ''.join(['+' if e >= 0 else '−' for e in eps])
        ax.text(0.5, -0.18, f"Pattern: {signs}  |  {sys['note']}",
                transform=ax.transAxes, fontsize=8, color='#999999',
                ha='center', va='top', family='monospace')

        ax.grid(True, alpha=0.1, color='#444', axis='y')
    else:
        ax.set_visible(False)

plt.tight_layout(rect=[0, 0.02, 1, 0.93], h_pad=4.0, w_pad=2.0)

fig.savefig('/home/jackwayne/Desktop/benford_fingerprint_atlas.png', dpi=180,
            facecolor='#1a1a2e', bbox_inches='tight', pad_inches=0.3)
plt.close()

# ══════════════════════════════════════════════════════════════
# FIGURE 2: Overview — δ_B comparison + pattern classification
# ══════════════════════════════════════════════════════════════

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 10), facecolor='#1a1a2e')

# Left: δ_B bar chart sorted
sorted_names = sorted(names, key=lambda n: systems[n]['delta_b'])
dbs = [systems[n]['delta_b'] for n in sorted_names]
colors = [systems[n]['color'] for n in sorted_names]

bars = ax1.barh(range(len(sorted_names)), dbs, color=colors, alpha=0.8,
                edgecolor='white', linewidth=0.5)
ax1.set_yticks(range(len(sorted_names)))
ax1.set_yticklabels(sorted_names, color='#e0e0e0', fontsize=10)
ax1.set_xlabel('δ_B (Euclidean deviation from Benford)', color='#e0e0e0', fontsize=12)
ax1.set_title('Benford Conformance Ranking', color='white', fontsize=14, fontweight='bold')
ax1.set_facecolor('#1a1a2e')
ax1.tick_params(colors='#e0e0e0')
for spine in ax1.spines.values():
    spine.set_color('#333')

# Threshold lines
ax1.axvline(x=0.02, color='#00ff88', linestyle='--', linewidth=1, alpha=0.5)
ax1.text(0.02, len(sorted_names) - 0.5, ' Strong', color='#00ff88', fontsize=8, va='top')
ax1.axvline(x=0.05, color='#ffaa00', linestyle='--', linewidth=1, alpha=0.5)
ax1.text(0.05, len(sorted_names) - 0.5, ' Moderate', color='#ffaa00', fontsize=8, va='top')
ax1.axvline(x=0.10, color='#ff4444', linestyle='--', linewidth=1, alpha=0.5)
ax1.text(0.10, len(sorted_names) - 0.5, ' Weak', color='#ff4444', fontsize=8, va='top')

ax1.grid(True, alpha=0.1, color='#444', axis='x')

# Right: Overlaid deviation patterns for quantum stats only
ax2.set_facecolor('#1a1a2e')
ax2.tick_params(colors='#e0e0e0')
for spine in ax2.spines.values():
    spine.set_color('#333')

quantum_names = ['Bose-Einstein', 'Fermi-Dirac', 'Maxwell-Boltzmann', 'Planck Spectrum', 'LQG Products']
for name in quantum_names:
    sys = systems[name]
    ax2.plot(range(1, 10), sys['epsilon'], 'o-', color=sys['color'],
             linewidth=2, markersize=6, label=f"{name} (δ_B={sys['delta_b']:.4f})")

ax2.axhline(y=0, color='white', linewidth=0.5, alpha=0.5)
ax2.set_xticks(range(1, 10))
ax2.set_xlabel('Digit d', color='#e0e0e0', fontsize=12)
ax2.set_ylabel('ε(d) = P_obs(d) − P_benford(d)', color='#e0e0e0', fontsize=12)
ax2.set_title('Deviation Fingerprints Overlaid\n(Conforming + Near-Conforming Systems)',
              color='white', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10, facecolor='#22223a', edgecolor='#444', labelcolor='#e0e0e0',
           loc='upper right', bbox_to_anchor=(1.0, 1.0))
ax2.grid(True, alpha=0.15, color='#444')

plt.tight_layout()
fig2.savefig('/home/jackwayne/Desktop/benford_fingerprint_overview.png', dpi=180,
             facecolor='#1a1a2e', bbox_inches='tight', pad_inches=0.3)
plt.close()

# ══════════════════════════════════════════════════════════════
# ANALYSIS: Pattern classification
# ══════════════════════════════════════════════════════════════

print("=" * 70)
print("  BENFORD FINGERPRINT ATLAS — PATTERN ANALYSIS")
print("=" * 70)

print("\n  CONFORMING / NEAR-CONFORMING (δ_B < 0.05):")
print("  " + "─" * 66)
for name in sorted_names:
    sys = systems[name]
    if sys['delta_b'] < 0.05:
        eps = sys['epsilon']
        signs = ''.join(['+' if e >= 0 else '−' for e in eps])
        # Classify pattern shape
        sign_changes = sum(1 for i in range(len(eps)-1) if (eps[i] >= 0) != (eps[i+1] >= 0))

        if sign_changes <= 1:
            shape = "monotonic"
        elif sign_changes == 2:
            shape = "single oscillation"
        else:
            shape = f"oscillatory ({sign_changes} crossings)"

        print(f"\n  {name}")
        print(f"    δ_B = {sys['delta_b']:.6f}  |  Type: {sys['type']}")
        print(f"    Pattern: {signs}  →  {shape}")
        print(f"    ε(1) = {eps[0]:+.6f}  (digit 1 deviation)")
        print(f"    Note: {sys['note']}")

print("\n\n  DEVIATING (δ_B ≥ 0.05):")
print("  " + "─" * 66)
for name in sorted_names:
    sys = systems[name]
    if sys['delta_b'] >= 0.05:
        eps = sys['epsilon']
        signs = ''.join(['+' if e >= 0 else '−' for e in eps])
        sign_changes = sum(1 for i in range(len(eps)-1) if (eps[i] >= 0) != (eps[i+1] >= 0))
        shape = f"oscillatory ({sign_changes} crossings)" if sign_changes > 2 else "stepped" if sign_changes <= 1 else f"{sign_changes} crossings"

        print(f"\n  {name}")
        print(f"    δ_B = {sys['delta_b']:.6f}  |  Type: {sys['type']}")
        print(f"    Pattern: {signs}  →  {shape}")
        print(f"    Note: {sys['note']}")

# ── Key observations ──
print("\n\n" + "=" * 70)
print("  KEY OBSERVATIONS")
print("=" * 70)

# Compare ε(1) across conforming systems
print("\n  Digit 1 deviation across conforming systems:")
for name in ['LQG Products', 'Bose-Einstein', 'Maxwell-Boltzmann',
             'Fermi-Dirac', 'Planck Spectrum']:
    eps1 = systems[name]['epsilon'][0]
    db = systems[name]['delta_b']
    print(f"    {name:25s}  ε(1) = {eps1:+.6f}   δ_B = {db:.6f}")

# Dot products between deviation vectors (similarity)
print("\n  Deviation vector similarities (cosine):")
conf_names = ['Bose-Einstein', 'Fermi-Dirac', 'Maxwell-Boltzmann', 'Planck Spectrum', 'LQG Products']
for i in range(len(conf_names)):
    for j in range(i+1, len(conf_names)):
        v1 = systems[conf_names[i]]['epsilon']
        v2 = systems[conf_names[j]]['epsilon']
        cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        print(f"    {conf_names[i]:20s} vs {conf_names[j]:20s}  cosine = {cos:+.4f}")

# L2 norm of each deviation vector
print("\n  L2 norm of deviation vectors (= δ_B):")
for name in sorted_names:
    db = systems[name]['delta_b']
    db_calc = np.linalg.norm(systems[name]['epsilon'])
    print(f"    {name:25s}  δ_B = {db:.6f}  (computed: {db_calc:.6f})")

print(f"\n  Benford L2 norm (the floor): 0.4068")
print(f"  Any δ_B close to 0.4068? Let's see...")
for name in sorted_names:
    db = systems[name]['delta_b']
    ratio = db / 0.4068
    print(f"    {name:25s}  δ_B/0.4068 = {ratio:.4f}")

print("\n" + "=" * 70)
print("  Figures saved:")
print("    ~/Desktop/benford_fingerprint_atlas.png")
print("    ~/Desktop/benford_fingerprint_overview.png")
print("=" * 70)
