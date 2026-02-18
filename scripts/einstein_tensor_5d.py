"""
Compute the Einstein tensor of the 5D metric with Euler prime substrate.

Metric: g_μν = diag(g_tt, g_rr, g_θθ, g_φφ, g_ξξ)
  g_tt(r)  = -(1 - r_s/r)           Schwarzschild time
  g_rr(r)  = 1 + r_s/r              Painlevé-Gullstrand radial
  g_θθ(r)  = r²                     Angular
  g_φφ(r)  = r²                     Angular (equatorial plane, sin²θ=1)
  g_ξξ(r)  = ln ζ(s(r))             Euler prime 5th dimension

s(r) = 1 + r/r_s  maps radius to zeta argument:
  singularity r→0: s→1, ζ diverges
  horizon r=r_s:   s=2, ζ=π²/6
  far away r→∞:    s→∞, ζ→1

All metric components depend only on r. Diagonal metric.
Compute at equatorial plane θ=π/2.
"""

import numpy as np
import json
import sys
from mpmath import zeta as mpzeta, log as mplog, mp

mp.dps = 30  # high precision

sys.path.insert(0, '/home/jackwayne/Desktop/Projects/Benford_Fun/scripts')
from benford_core import run_full_analysis, euclidean_deviation, observed_distribution, first_digit, per_digit_deviation

# ===========================================================================
# Metric definition
# ===========================================================================

R_S = 1.0  # Schwarzschild radius (geometric units)

def zeta_ln(s):
    """ln(ζ(s)) computed via mpmath."""
    if s <= 1.0:
        return float('inf')
    return float(mplog(mpzeta(s)))

# Precompute g_ξξ on a fine grid and interpolate for speed
print("Precomputing ζ values...")
_s_grid = np.linspace(1.005, 102.0, 20000)
_lnzeta_grid = np.array([zeta_ln(s) for s in _s_grid])
from scipy.interpolate import interp1d
_lnzeta_interp = interp1d(_s_grid, _lnzeta_grid, kind='cubic', fill_value='extrapolate')

def g_xixi(r):
    s = 1.0 + r / R_S
    return float(_lnzeta_interp(s))

def metric_diag(r):
    """Return diagonal metric components [g_tt, g_rr, g_θθ, g_φφ, g_ξξ]."""
    return np.array([
        -(1.0 - R_S / r),   # g_tt
        1.0 + R_S / r,       # g_rr
        r**2,                 # g_θθ
        r**2,                 # g_φφ (equatorial)
        g_xixi(r),            # g_ξξ
    ])

# ===========================================================================
# Numerical derivatives
# ===========================================================================

H = 1e-5  # step size for finite differences

def dg_dr(r):
    return (metric_diag(r + H) - metric_diag(r - H)) / (2 * H)

def d2g_dr2(r):
    return (metric_diag(r + H) - 2 * metric_diag(r) + metric_diag(r - H)) / (H * H)

# ===========================================================================
# Christoffel symbols for diagonal metric, all depending on x^1 = r
# ===========================================================================
# Non-zero:
#   Γ^1_{11} = g1'/(2*g1)
#   Γ^1_{aa} = -ga'/(2*g1)     for a≠1
#   Γ^a_{a1} = Γ^a_{1a} = ga'/(2*ga)  for a≠1
# ===========================================================================

N = 5  # dimensions

def compute_christoffel(r):
    g = metric_diag(r)
    dg = dg_dr(r)
    Gamma = np.zeros((N, N, N))

    # Γ^1_{11}
    Gamma[1, 1, 1] = dg[1] / (2 * g[1])

    for a in range(N):
        if a != 1:
            Gamma[1, a, a] = -dg[a] / (2 * g[1])
            Gamma[a, a, 1] = dg[a] / (2 * g[a])
            Gamma[a, 1, a] = dg[a] / (2 * g[a])

    return Gamma

# ===========================================================================
# Riemann → Ricci → Einstein tensor
# ===========================================================================

def compute_einstein_tensor(r):
    g = metric_diag(r)
    Gamma = compute_christoffel(r)

    # Derivative of Christoffel symbols w.r.t. r (the only non-trivial derivative)
    Gamma_p = compute_christoffel(r + H)
    Gamma_m = compute_christoffel(r - H)
    dGamma = (Gamma_p - Gamma_m) / (2 * H)

    # Riemann tensor: R^ρ_{σμν} = ∂_μ Γ^ρ_{νσ} - ∂_ν Γ^ρ_{μσ} + Γ^ρ_{μλ}Γ^λ_{νσ} - Γ^ρ_{νλ}Γ^λ_{μσ}
    # Only ∂_1 (∂_r) derivatives are non-zero
    Riemann = np.zeros((N, N, N, N))

    for rho in range(N):
        for sig in range(N):
            for mu in range(N):
                for nu in range(N):
                    val = 0.0
                    # ∂_μ Γ^ρ_{νσ}: non-zero only if μ=1
                    if mu == 1:
                        val += dGamma[rho, nu, sig]
                    # -∂_ν Γ^ρ_{μσ}: non-zero only if ν=1
                    if nu == 1:
                        val -= dGamma[rho, mu, sig]
                    # Γ^ρ_{μλ} Γ^λ_{νσ}
                    for lam in range(N):
                        val += Gamma[rho, mu, lam] * Gamma[lam, nu, sig]
                    # -Γ^ρ_{νλ} Γ^λ_{μσ}
                    for lam in range(N):
                        val -= Gamma[rho, nu, lam] * Gamma[lam, mu, sig]
                    Riemann[rho, sig, mu, nu] = val

    # Ricci tensor: R_{σν} = R^ρ_{σρν}
    Ricci = np.zeros((N, N))
    for sig in range(N):
        for nu in range(N):
            for rho in range(N):
                Ricci[sig, nu] += Riemann[rho, sig, rho, nu]

    # Ricci scalar: R = g^{μν} R_{μν} = Σ R_{μμ}/g_{μμ}
    R_scalar = sum(Ricci[mu, mu] / g[mu] for mu in range(N))

    # Einstein tensor: G_{μν} = R_{μν} - (1/2) R g_{μν}
    G = np.zeros((N, N))
    for mu in range(N):
        for nu in range(N):
            G[mu, nu] = Ricci[mu, nu] - 0.5 * R_scalar * (g[mu] if mu == nu else 0.0)

    return G, Ricci, R_scalar, Riemann

# ===========================================================================
# Main computation
# ===========================================================================

print("Computing Einstein tensor across radial range...")
print("=" * 70)

# Sample radial points: from just outside horizon to far away
# Avoid r < r_s + epsilon (coordinate issues) and r too close to 0
r_values = np.concatenate([
    np.linspace(1.01, 2.0, 200),     # near horizon
    np.linspace(2.0, 10.0, 300),      # intermediate
    np.linspace(10.0, 100.0, 500),    # far field
])

all_G_values = []       # all |G_{μν}| values
all_R_values = []       # all |R_{μν}| values
all_Riem_values = []    # all |R^ρ_{σμν}| values
diagonal_G = {i: [] for i in range(N)}  # G_{ii} by component

for idx, r in enumerate(r_values):
    if idx % 100 == 0:
        print(f"  r = {r:.3f} ({idx}/{len(r_values)})")

    G, Ricci, R_sc, Riemann = compute_einstein_tensor(r)

    # Collect all non-trivial Einstein tensor components
    for i in range(N):
        for j in range(N):
            val = abs(G[i, j])
            if val > 1e-15:
                all_G_values.append(val)
            if i == j:
                diagonal_G[i].append(abs(G[i, i]) if abs(G[i, i]) > 1e-15 else None)

    # Collect Ricci tensor values
    for i in range(N):
        for j in range(N):
            val = abs(Ricci[i, j])
            if val > 1e-15:
                all_R_values.append(val)

    # Collect Riemann tensor values
    for rho in range(N):
        for sig in range(N):
            for mu in range(N):
                for nu in range(N):
                    val = abs(Riemann[rho, sig, mu, nu])
                    if val > 1e-15:
                        all_Riem_values.append(val)

    # Also collect the scalar curvature
    if abs(R_sc) > 1e-15:
        all_G_values.append(abs(R_sc))

print(f"\nTotal Einstein tensor values collected: {len(all_G_values)}")
print(f"Total Ricci tensor values collected:    {len(all_R_values)}")
print(f"Total Riemann tensor values collected:  {len(all_Riem_values)}")

# ===========================================================================
# Benford analysis
# ===========================================================================

print("\n" + "=" * 70)
print("BENFORD ANALYSIS: Einstein Tensor G_μν")
print("=" * 70)

result_G = run_full_analysis(all_G_values)
print(f"  N values:      {result_G['n']}")
print(f"  VERDICT:       {result_G['verdict']}")
print(f"  δ_B:           {result_G['first_digit']['euclidean_deviation']}")
print(f"  MAD:           {result_G['first_digit']['mad']['mad']}")
print(f"  MAD class:     {result_G['first_digit']['mad']['classification']}")
print(f"  χ² p-value:    {result_G['first_digit']['chi_squared']['p_value']}")
print(f"  KS p-value:    {result_G['first_digit']['ks_test']['p_value']}")
print(f"\n  Observed distribution:")
for d in range(1, 10):
    obs = result_G['first_digit']['observed_distribution'][d]
    exp = result_G['first_digit']['expected_distribution'][d]
    dev = result_G['first_digit']['per_digit_deviation'][d]
    bar = "█" * int(obs * 200)
    print(f"    d={d}: obs={obs:.4f}  exp={exp:.4f}  ε={dev:+.4f}  {bar}")

print("\n" + "=" * 70)
print("BENFORD ANALYSIS: Ricci Tensor R_μν")
print("=" * 70)

result_R = run_full_analysis(all_R_values)
print(f"  N values:      {result_R['n']}")
print(f"  VERDICT:       {result_R['verdict']}")
print(f"  δ_B:           {result_R['first_digit']['euclidean_deviation']}")
print(f"  MAD:           {result_R['first_digit']['mad']['mad']}")
print(f"  χ² p-value:    {result_R['first_digit']['chi_squared']['p_value']}")

print("\n" + "=" * 70)
print("BENFORD ANALYSIS: Riemann Tensor R^ρ_σμν")
print("=" * 70)

result_Riem = run_full_analysis(all_Riem_values)
print(f"  N values:      {result_Riem['n']}")
print(f"  VERDICT:       {result_Riem['verdict']}")
print(f"  δ_B:           {result_Riem['first_digit']['euclidean_deviation']}")
print(f"  MAD:           {result_Riem['first_digit']['mad']['mad']}")
print(f"  χ² p-value:    {result_Riem['first_digit']['chi_squared']['p_value']}")

# ===========================================================================
# Per-component analysis
# ===========================================================================

print("\n" + "=" * 70)
print("PER-COMPONENT: Diagonal G_{ii}")
print("=" * 70)

labels = ['G_tt', 'G_rr', 'G_θθ', 'G_φφ', 'G_ξξ']
for i in range(N):
    vals = [v for v in diagonal_G[i] if v is not None]
    if len(vals) > 10:
        res = run_full_analysis(vals)
        print(f"  {labels[i]:6s}: n={res['n']:5d}  δ_B={res['first_digit']['euclidean_deviation']:.6f}  verdict={res['verdict']}")

# ===========================================================================
# Comparison: 4D Schwarzschild (no 5th dimension)
# ===========================================================================

print("\n" + "=" * 70)
print("CONTROL: 4D Schwarzschild Einstein Tensor (no 5th dimension)")
print("=" * 70)

def metric_4d(r):
    return np.array([
        -(1.0 - R_S / r),
        1.0 + R_S / r,
        r**2,
        r**2,
    ])

def compute_einstein_4d(r):
    N4 = 4
    h = H
    g = metric_4d(r)
    dg = (metric_4d(r + h) - metric_4d(r - h)) / (2 * h)

    def christoffel_4d(rv):
        gv = metric_4d(rv)
        dgv = (metric_4d(rv + h) - metric_4d(rv - h)) / (2 * h)
        G4 = np.zeros((N4, N4, N4))
        G4[1, 1, 1] = dgv[1] / (2 * gv[1])
        for a in range(N4):
            if a != 1:
                G4[1, a, a] = -dgv[a] / (2 * gv[1])
                G4[a, a, 1] = dgv[a] / (2 * gv[a])
                G4[a, 1, a] = dgv[a] / (2 * gv[a])
        return G4

    Gam = christoffel_4d(r)
    Gam_p = christoffel_4d(r + h)
    Gam_m = christoffel_4d(r - h)
    dGam = (Gam_p - Gam_m) / (2 * h)

    Riem = np.zeros((N4, N4, N4, N4))
    for rho in range(N4):
        for sig in range(N4):
            for mu in range(N4):
                for nu in range(N4):
                    val = 0.0
                    if mu == 1:
                        val += dGam[rho, nu, sig]
                    if nu == 1:
                        val -= dGam[rho, mu, sig]
                    for lam in range(N4):
                        val += Gam[rho, mu, lam] * Gam[lam, nu, sig]
                        val -= Gam[rho, nu, lam] * Gam[lam, mu, sig]
                    Riem[rho, sig, mu, nu] = val

    Ric = np.zeros((N4, N4))
    for sig in range(N4):
        for nu in range(N4):
            for rho in range(N4):
                Ric[sig, nu] += Riem[rho, sig, rho, nu]

    R_sc = sum(Ric[mu, mu] / g[mu] for mu in range(N4))

    Ein = np.zeros((N4, N4))
    for mu in range(N4):
        for nu in range(N4):
            Ein[mu, nu] = Ric[mu, nu] - 0.5 * R_sc * (g[mu] if mu == nu else 0.0)

    return Ein

all_G4_values = []
for r in r_values:
    G4 = compute_einstein_4d(r)
    for i in range(4):
        for j in range(4):
            val = abs(G4[i, j])
            if val > 1e-15:
                all_G4_values.append(val)

result_4d = run_full_analysis(all_G4_values)
print(f"  N values:      {result_4d['n']}")
print(f"  VERDICT:       {result_4d['verdict']}")
print(f"  δ_B:           {result_4d['first_digit']['euclidean_deviation']}")
print(f"  MAD:           {result_4d['first_digit']['mad']['mad']}")
print(f"  χ² p-value:    {result_4d['first_digit']['chi_squared']['p_value']}")

# ===========================================================================
# Save results
# ===========================================================================

output = {
    "description": "5D Einstein tensor with Euler prime substrate — Benford analysis",
    "metric": {
        "g_tt": "-(1 - r_s/r)",
        "g_rr": "1 + r_s/r",
        "g_thth": "r²",
        "g_phph": "r² (equatorial)",
        "g_xixi": "ln ζ(s), s = 1 + r/r_s",
    },
    "r_range": [float(r_values[0]), float(r_values[-1])],
    "n_radial_points": len(r_values),
    "einstein_tensor_5d": {
        "n_values": result_G['n'],
        "verdict": result_G['verdict'],
        "delta_b": result_G['first_digit']['euclidean_deviation'],
        "mad": result_G['first_digit']['mad']['mad'],
        "chi_sq_p": result_G['first_digit']['chi_squared']['p_value'],
        "ks_p": result_G['first_digit']['ks_test']['p_value'],
        "observed": result_G['first_digit']['observed_distribution'],
        "per_digit_deviation": result_G['first_digit']['per_digit_deviation'],
    },
    "ricci_tensor_5d": {
        "n_values": result_R['n'],
        "verdict": result_R['verdict'],
        "delta_b": result_R['first_digit']['euclidean_deviation'],
    },
    "riemann_tensor_5d": {
        "n_values": result_Riem['n'],
        "verdict": result_Riem['verdict'],
        "delta_b": result_Riem['first_digit']['euclidean_deviation'],
    },
    "control_4d_schwarzschild": {
        "n_values": result_4d['n'],
        "verdict": result_4d['verdict'],
        "delta_b": result_4d['first_digit']['euclidean_deviation'],
    },
}

with open('/home/jackwayne/Desktop/Projects/Benford_Fun/results/round_trip/einstein_5d_benford.json', 'w') as f:
    json.dump(output, f, indent=2)

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"  5D Einstein tensor G_μν:  δ_B = {result_G['first_digit']['euclidean_deviation']:.6f}  → {result_G['verdict']}")
print(f"  5D Ricci tensor R_μν:     δ_B = {result_R['first_digit']['euclidean_deviation']:.6f}  → {result_R['verdict']}")
print(f"  5D Riemann R^ρ_σμν:       δ_B = {result_Riem['first_digit']['euclidean_deviation']:.6f}  → {result_Riem['verdict']}")
print(f"  4D Schwarzschild (ctrl):  δ_B = {result_4d['first_digit']['euclidean_deviation']:.6f}  → {result_4d['verdict']}")
print(f"\nResults saved to results/round_trip/einstein_5d_benford.json")
