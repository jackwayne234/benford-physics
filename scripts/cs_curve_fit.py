#!/usr/bin/env python3
"""
Analyze the 40 CS delta_B data points to find the equation governing
how the Causal Set substrate responds to black hole geometry.

The goal: find δ_B(r) — the function that describes CS deviation
as a function of radial position through the black hole.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

# ── The Data ──
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

# Benford equilibrium (from the literature and Christopher's work)
DB_EQ = 0.017

print("=" * 70)
print("CS delta_B Curve Analysis")
print("=" * 70)

# ── Basic statistics by region ──
outside = r_data > 1.01
horizon = (r_data >= 0.99) & (r_data <= 1.01)
inside = r_data < 0.99

print(f"\nOutside horizon (r > 1.01): n={sum(outside)}")
print(f"  Mean δ_B = {db_data[outside].mean():.6f}")
print(f"  Min  δ_B = {db_data[outside].min():.6f} at r = {r_data[outside][db_data[outside].argmin()]:.4f}")

print(f"\nAt horizon (0.99-1.01): n={sum(horizon)}")
print(f"  Mean δ_B = {db_data[horizon].mean():.6f}")

print(f"\nInside horizon (r < 0.99): n={sum(inside)}")
print(f"  Mean δ_B = {db_data[inside].mean():.6f}")
print(f"  Approaches equilibrium: ~{DB_EQ}")

# ── Candidate Functions ──

# 1. Power law with offset: δ_B = a * |r - r0|^b + c
def power_offset(r, a, r0, b, c):
    return a * np.abs(r - r0)**b + c

# 2. Gravitational-inspired: δ_B = a / r^n + b * r^m + c
def grav_inspired(r, a, n, b, m, c):
    return a / r**n + b * r**m + c

# 3. Benford attractor: δ_B = eq * (1 + a * |ln(r)|^b * sign_term)
# The idea: equilibrium is 0.017, deviations are driven by |ln(r)|
def benford_attractor(r, eq, a, b):
    return eq * (1 + a * np.abs(np.log(r))**b)

# 4. Horizon-centered: deviation driven by distance from horizon in log space
def horizon_centered(r, eq, a, b, c):
    # distance from horizon in log space
    d = np.log(r)  # 0 at horizon, positive outside, negative inside
    return eq * (1 + a * np.abs(d)**b) + c * np.exp(-5 * np.abs(d))

# 5. Simple: two-regime piecewise — outside is 1/r decay, inside climbs to eq
def two_regime(r, a, b, eq, k):
    result = np.zeros_like(r)
    out = r >= 1
    inn = r < 1
    result[out] = a / r[out]**b
    result[inn] = eq * (1 - np.exp(-k * (1 - r[inn])))
    # smooth blend near horizon
    blend = np.exp(-20 * (r - 1)**2)
    return result * (1 - blend) + (a / 1.0**b) * blend  # not great, try anyway

# 6. Damped oscillation around equilibrium
def damped_osc(r, eq, A, omega, phi, gamma):
    x = np.log(r)  # log-radius as the "time" variable
    return eq + A * np.exp(-gamma * np.abs(x)) * np.cos(omega * x + phi)

# 7. Simplest possible: δ_B = a * |ln(r)| + b
def log_abs(r, a, b):
    return a * np.abs(np.log(r)) + b

# 8. Benford's equation applied to r itself: δ_B = a * log10(1 + 1/r) + b
def benford_of_r(r, a, b):
    return a * np.log10(1 + 1/r) + b

# 9. Benford of r with scaling: δ_B = a * log10(1 + c/r) + b
def benford_of_r_scaled(r, a, b, c):
    return a * np.log10(1 + c / r) + b

# ── Fit each candidate ──
results = []

# Fit 3: Benford attractor
try:
    popt, pcov = curve_fit(benford_attractor, r_data, db_data,
                           p0=[0.017, 0.5, 1.0], maxfev=10000)
    pred = benford_attractor(r_data, *popt)
    resid = np.sqrt(np.mean((pred - db_data)**2))
    results.append(('Benford attractor: eq*(1 + a*|ln(r)|^b)', popt, resid,
                     benford_attractor, ['eq', 'a', 'b']))
    print(f"\n3. Benford attractor: RMSE = {resid:.6f}")
    print(f"   eq={popt[0]:.6f}, a={popt[1]:.6f}, b={popt[2]:.6f}")
except Exception as e:
    print(f"\n3. Benford attractor: FAILED ({e})")

# Fit 6: Damped oscillation
try:
    popt, pcov = curve_fit(damped_osc, r_data, db_data,
                           p0=[0.01, 0.01, 3.0, 0, 1.0], maxfev=20000)
    pred = damped_osc(r_data, *popt)
    resid = np.sqrt(np.mean((pred - db_data)**2))
    results.append(('Damped oscillation', popt, resid, damped_osc,
                     ['eq', 'A', 'omega', 'phi', 'gamma']))
    print(f"\n6. Damped oscillation: RMSE = {resid:.6f}")
    print(f"   eq={popt[0]:.6f}, A={popt[1]:.6f}, omega={popt[2]:.6f}, phi={popt[3]:.6f}, gamma={popt[4]:.6f}")
except Exception as e:
    print(f"\n6. Damped oscillation: FAILED ({e})")

# Fit 7: Simple |ln(r)|
try:
    popt, pcov = curve_fit(log_abs, r_data, db_data, p0=[0.005, 0.003])
    pred = log_abs(r_data, *popt)
    resid = np.sqrt(np.mean((pred - db_data)**2))
    results.append(('|ln(r)| linear: a*|ln(r)| + b', popt, resid,
                     log_abs, ['a', 'b']))
    print(f"\n7. |ln(r)| linear: RMSE = {resid:.6f}")
    print(f"   a={popt[0]:.6f}, b={popt[1]:.6f}")
except Exception as e:
    print(f"\n7. |ln(r)| linear: FAILED ({e})")

# Fit 8: Benford's equation applied to r
try:
    popt, pcov = curve_fit(benford_of_r, r_data, db_data, p0=[0.05, 0.001])
    pred = benford_of_r(r_data, *popt)
    resid = np.sqrt(np.mean((pred - db_data)**2))
    results.append(('B(r): a*log10(1+1/r) + b', popt, resid,
                     benford_of_r, ['a', 'b']))
    print(f"\n8. Benford of r: RMSE = {resid:.6f}")
    print(f"   a={popt[0]:.6f}, b={popt[1]:.6f}")
except Exception as e:
    print(f"\n8. Benford of r: FAILED ({e})")

# Fit 9: Benford of r with scale
try:
    popt, pcov = curve_fit(benford_of_r_scaled, r_data, db_data,
                           p0=[0.05, 0.001, 1.0], maxfev=10000)
    pred = benford_of_r_scaled(r_data, *popt)
    resid = np.sqrt(np.mean((pred - db_data)**2))
    results.append(('B(r,c): a*log10(1+c/r) + b', popt, resid,
                     benford_of_r_scaled, ['a', 'b', 'c']))
    print(f"\n9. Benford of r (scaled): RMSE = {resid:.6f}")
    print(f"   a={popt[0]:.6f}, b={popt[1]:.6f}, c={popt[2]:.6f}")
except Exception as e:
    print(f"\n9. Benford of r (scaled): FAILED ({e})")

# Fit 4: Horizon-centered
try:
    popt, pcov = curve_fit(horizon_centered, r_data, db_data,
                           p0=[0.01, 0.5, 1.0, 0.005], maxfev=20000)
    pred = horizon_centered(r_data, *popt)
    resid = np.sqrt(np.mean((pred - db_data)**2))
    results.append(('Horizon-centered', popt, resid, horizon_centered,
                     ['eq', 'a', 'b', 'c']))
    print(f"\n4. Horizon-centered: RMSE = {resid:.6f}")
    print(f"   eq={popt[0]:.6f}, a={popt[1]:.6f}, b={popt[2]:.6f}, c={popt[3]:.6f}")
except Exception as e:
    print(f"\n4. Horizon-centered: FAILED ({e})")

# ── Sort by RMSE ──
results.sort(key=lambda x: x[2])

print("\n" + "=" * 70)
print("RANKINGS (by RMSE — lower is better)")
print("=" * 70)
for i, (name, popt, resid, func, pnames) in enumerate(results):
    print(f"\n{i+1}. {name}")
    print(f"   RMSE = {resid:.6f}")
    for pn, pv in zip(pnames, popt):
        print(f"   {pn} = {pv:.6f}")

# ── Visualization ──
fig, axes = plt.subplots(2, 2, figsize=(16, 12), facecolor='#1a1a2e')
fig.suptitle('CS δ_B Curve Fitting — Finding the CS Equation',
             color='white', fontsize=16, fontweight='bold')

r_smooth = np.logspace(np.log10(0.01), np.log10(10), 500)

for idx, ax in enumerate(axes.flat):
    ax.set_facecolor('#1a1a2e')
    ax.tick_params(colors='#e0e0e0')
    ax.set_xlabel('r / r_s', color='#e0e0e0')
    ax.set_ylabel('δ_B', color='#e0e0e0')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.15, color='#444')
    for spine in ax.spines.values():
        spine.set_color('#333')

    # Data points
    ax.scatter(r_data, db_data, color='#ffdd44', s=30, zorder=5, label='CS data')

    # Horizon line
    ax.axvline(x=1.0, color='#ffffff', linestyle='--', linewidth=0.8, alpha=0.4)

    # Benford equilibrium
    ax.axhline(y=DB_EQ, color='#ff8800', linestyle=':', linewidth=0.8, alpha=0.5,
               label=f'Benford eq = {DB_EQ}')

    if idx < len(results):
        name, popt, resid, func, pnames = results[idx]
        try:
            y_smooth = func(r_smooth, *popt)
            # Clip for display
            y_smooth = np.clip(y_smooth, 0, 0.05)
            ax.plot(r_smooth, y_smooth, color='#00ffaa', linewidth=2,
                    label=f'Fit (RMSE={resid:.5f})')
        except:
            pass
        ax.set_title(f'#{idx+1}: {name}', color='#e0e0e0', fontsize=11)
        ax.legend(fontsize=8, facecolor='#22223a', edgecolor='#444', labelcolor='#e0e0e0')
    else:
        ax.set_title('(no fit)', color='#666')

plt.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig('/home/jackwayne/Desktop/cs_curve_fit.png', dpi=200, facecolor='#1a1a2e')
print(f"\nFigure saved to ~/Desktop/cs_curve_fit.png")

# ── The winner's equation, printed cleanly ──
if results:
    best_name, best_popt, best_resid, best_func, best_pnames = results[0]
    print("\n" + "=" * 70)
    print(f"BEST FIT: {best_name}")
    print(f"RMSE: {best_resid:.6f}")
    print("Parameters:")
    for pn, pv in zip(best_pnames, best_popt):
        print(f"  {pn} = {pv:.6f}")
    print("=" * 70)
