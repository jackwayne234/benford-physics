#!/usr/bin/env python3
"""
Matrix(CS) = f(r)?

Left side:  g_δ = log₁₀(1 + 1/δ_B)  using actual data at each r
Right side: g_δ = log₁₀(1 + 1/f(r))  using each candidate equation

Direct comparison — no matrix, no simulation.
Just: does the equation reproduce what we measured?
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ── The 40 Data Points ──
r_data = np.array([
    10.0, 7.0, 5.0, 3.0, 2.0, 1.5, 1.3, 1.2, 1.15, 1.1,
    1.08, 1.06, 1.04, 1.03, 1.02, 1.015, 1.01, 1.005, 1.002, 1.001,
    0.99, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,
    0.25, 0.2, 0.15, 0.12, 0.1, 0.08, 0.06, 0.04, 0.02, 0.01
])

db_data = np.array([
    0.027552, 0.010611, 0.004223, 0.002856, 0.002748, 0.005169,
    0.003184, 0.003569, 0.002757, 0.002136, 0.00255, 0.002722,
    0.003532, 0.003889, 0.003482, 0.003757, 0.004013, 0.0042,
    0.004049, 0.004166, 0.004267, 0.005472, 0.00286, 0.004293,
    0.003671, 0.006253, 0.006376, 0.00486, 0.007448, 0.01313,
    0.016065, 0.01205, 0.015154, 0.019625, 0.016711, 0.017486,
    0.01349, 0.017329, 0.018378, 0.014661
])

# ── Benford's equation ──
def B(x):
    """Benford's equation: log₁₀(1 + 1/x)"""
    return np.log10(1 + 1/np.maximum(x, 1e-10))

# ── LEFT SIDE: Matrix(CS) from data ──
g_delta_data = B(db_data)

# ── The 5 Candidate Equations ──

# 1. Damped oscillation (5 params)
def damped_osc(r, eq, A, omega, phi, gamma):
    x = np.log(r)
    return eq + A * np.exp(-gamma * np.abs(x)) * np.cos(omega * x + phi)

# 2. Horizon-centered (4 params)
def horizon_centered(r, eq, a, b, c):
    d = np.log(r)
    return eq * (1 + a * np.abs(d)**b) + c * np.exp(-5 * np.abs(d))

# 3. Benford attractor (3 params)
def benford_attractor(r, eq, a, b):
    return eq * (1 + a * np.abs(np.log(r))**b)

# 4. |ln(r)| linear (2 params)
def log_abs(r, a, b):
    return a * np.abs(np.log(r)) + b

# 5. Benford of r (2 params)
def benford_of_r(r, a, b):
    return a * np.log10(1 + 1/r) + b

candidates = [
    ("Damped Oscillation", damped_osc, [0.01, 0.01, 3.0, 0, 1.0], 5),
    ("Horizon-Centered",   horizon_centered, [0.01, 0.5, 1.0, 0.005], 4),
    ("Benford Attractor",  benford_attractor, [0.017, 0.5, 1.0], 3),
    ("|ln(r)| Linear",     log_abs, [0.005, 0.003], 2),
    ("Benford of r",       benford_of_r, [0.05, 0.001], 2),
]

# ── Fit each equation and compare ──
print("=" * 80)
print("Matrix(CS) = f(r) ?")
print("Comparing g_δ from data vs g_δ from each candidate equation")
print("=" * 80)

fitted = []
for name, func, p0, nparams in candidates:
    try:
        popt, _ = curve_fit(func, r_data, db_data, p0=p0, maxfev=20000)
        db_pred = func(r_data, *popt)
        # Clamp predicted δ_B to positive values
        db_pred = np.maximum(db_pred, 1e-10)
        g_delta_eq = B(db_pred)

        # How well does g_delta_eq match g_delta_data?
        rmse_g = np.sqrt(np.mean((g_delta_eq - g_delta_data)**2))
        max_err = np.max(np.abs(g_delta_eq - g_delta_data))
        mean_pct_err = np.mean(np.abs(g_delta_eq - g_delta_data) / g_delta_data) * 100

        fitted.append((name, popt, db_pred, g_delta_eq, rmse_g, max_err, mean_pct_err, nparams))

        print(f"\n{'─' * 60}")
        print(f"  {name}  ({nparams} parameters)")
        print(f"{'─' * 60}")
        print(f"  g_δ RMSE:        {rmse_g:.4f}")
        print(f"  g_δ Max Error:   {max_err:.4f}")
        print(f"  g_δ Mean % Err:  {mean_pct_err:.2f}%")

    except Exception as e:
        print(f"\n  {name}: FAILED ({e})")

# ── Sort by RMSE ──
fitted.sort(key=lambda x: x[4])

print("\n" + "=" * 80)
print("RANKINGS  (by g_δ RMSE — lower = better match to Matrix(CS))")
print("=" * 80)
for i, (name, popt, _, _, rmse_g, max_err, pct_err, nparams) in enumerate(fitted):
    print(f"\n  #{i+1}  {name}  ({nparams} params)")
    print(f"      g_δ RMSE = {rmse_g:.4f}    Max Err = {max_err:.4f}    Mean %Err = {pct_err:.2f}%")

# ── Detailed table for the top 3 ──
print("\n" + "=" * 80)
print("POINT-BY-POINT: Matrix(CS)  vs  Top 3 Equations")
print("=" * 80)

header = f"{'r':>8s}  {'δ_B data':>10s}  {'g_δ DATA':>10s}"
for i in range(min(3, len(fitted))):
    short = fitted[i][0][:12]
    header += f"  {'g_δ '+short:>16s}"
print(header)
print("─" * len(header))

for j in range(len(r_data)):
    row = f"{r_data[j]:8.3f}  {db_data[j]:10.6f}  {g_delta_data[j]:10.4f}"
    for i in range(min(3, len(fitted))):
        g_eq = fitted[i][3][j]
        diff = g_eq - g_delta_data[j]
        row += f"  {g_eq:10.4f} ({diff:+.2f})"
    print(row)

# ── Visualization ──
fig, axes = plt.subplots(2, 3, figsize=(20, 12), facecolor='#1a1a2e')
fig.suptitle('Matrix(CS) = f(r) ?    —    Does the equation match the data?',
             color='white', fontsize=16, fontweight='bold')

r_smooth = np.logspace(np.log10(0.01), np.log10(10), 500)

for idx in range(min(5, len(fitted))):
    ax = axes.flat[idx]
    name, popt, db_pred, g_delta_eq, rmse_g, max_err, pct_err, nparams = fitted[idx]

    ax.set_facecolor('#1a1a2e')
    ax.tick_params(colors='#e0e0e0')
    for spine in ax.spines.values():
        spine.set_color('#333')
    ax.grid(True, alpha=0.15, color='#444')
    ax.set_xscale('log')

    # Data: g_delta from actual measurements
    ax.scatter(r_data, g_delta_data, color='#ffdd44', s=40, zorder=5,
               label='Matrix(CS) from data', edgecolors='#aa9900', linewidths=0.5)

    # Equation: g_delta from the candidate
    ax.scatter(r_data, g_delta_eq, color='#00ffaa', s=25, zorder=4,
               marker='x', linewidths=1.5, label=f'{name}')

    # Smooth curve for the equation
    try:
        func = [c for c in candidates if c[0] == name][0][1]
        db_smooth = func(r_smooth, *popt)
        db_smooth = np.maximum(db_smooth, 1e-10)
        g_smooth = B(db_smooth)
        g_smooth = np.clip(g_smooth, 0, 3)
        ax.plot(r_smooth, g_smooth, color='#00ffaa', linewidth=1.5, alpha=0.6)
    except:
        pass

    # Horizon line
    ax.axvline(x=1.0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.3)

    ax.set_xlabel('r / r_s', color='#e0e0e0')
    ax.set_ylabel('g_δ̂', color='#e0e0e0')
    ax.set_title(f'#{idx+1}: {name} ({nparams}p)\nRMSE={rmse_g:.4f}  %Err={pct_err:.1f}%',
                 color='#e0e0e0', fontsize=10)
    ax.legend(fontsize=7, facecolor='#22223a', edgecolor='#444', labelcolor='#e0e0e0')
    ax.set_ylim(0.8, 2.8)

# 6th panel: overlay all on one chart
ax = axes.flat[5]
ax.set_facecolor('#1a1a2e')
ax.tick_params(colors='#e0e0e0')
for spine in ax.spines.values():
    spine.set_color('#333')
ax.grid(True, alpha=0.15, color='#444')
ax.set_xscale('log')

ax.scatter(r_data, g_delta_data, color='#ffdd44', s=50, zorder=5,
           label='Matrix(CS) DATA', edgecolors='#aa9900', linewidths=0.5)

colors = ['#00ffaa', '#ff6688', '#6688ff', '#ffaa00', '#aa66ff']
for idx in range(min(5, len(fitted))):
    name, popt, _, g_delta_eq, rmse_g, _, _, nparams = fitted[idx]
    try:
        func = [c for c in candidates if c[0] == name][0][1]
        db_smooth = func(r_smooth, *popt)
        db_smooth = np.maximum(db_smooth, 1e-10)
        g_smooth = B(db_smooth)
        g_smooth = np.clip(g_smooth, 0, 3)
        ax.plot(r_smooth, g_smooth, color=colors[idx], linewidth=1.5, alpha=0.7,
                label=f'{name} ({nparams}p)')
    except:
        pass

ax.axvline(x=1.0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.3)
ax.set_xlabel('r / r_s', color='#e0e0e0')
ax.set_ylabel('g_δ̂', color='#e0e0e0')
ax.set_title('All 5 equations vs Matrix(CS)', color='#e0e0e0', fontsize=11)
ax.legend(fontsize=7, facecolor='#22223a', edgecolor='#444', labelcolor='#e0e0e0')
ax.set_ylim(0.8, 2.8)

plt.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig('/home/jackwayne/Desktop/cs_matrix_vs_equations.png', dpi=200, facecolor='#1a1a2e')
print(f"\nFigure saved to ~/Desktop/cs_matrix_vs_equations.png")

# ── Final verdict ──
best = fitted[0]
print("\n" + "=" * 80)
print(f"VERDICT")
print("=" * 80)
print(f"\n  Best match to Matrix(CS):  {best[0]}")
print(f"  Parameters: {best[7]}")
print(f"  g_δ RMSE:      {best[4]:.4f}")
print(f"  g_δ Mean %Err: {best[6]:.2f}%")
print(f"\n  Can this equation REPLACE the 40 data points?")
if best[6] < 5:
    print(f"  → YES — {best[6]:.1f}% mean error is tight enough.")
elif best[6] < 15:
    print(f"  → MAYBE — {best[6]:.1f}% mean error. Captures the shape but not the noise.")
else:
    print(f"  → NOT YET — {best[6]:.1f}% mean error. The data has structure these equations miss.")
print("=" * 80)
