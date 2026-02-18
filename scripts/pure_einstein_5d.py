"""
Pure Einstein tensor from the 5D metric. No Benford analysis.
Just: metric in → tensor out → plot what you see.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpmath import zeta as mpzeta, log as mplog, mp
from scipy.interpolate import interp1d

mp.dps = 30

# ===========================================================================
# Metric (exactly as in the matrix)
# ===========================================================================

R_S = 1.0

# Precompute ln(ζ(s))
print("Precomputing ζ values...")
_s_grid = np.linspace(1.002, 150.0, 30000)
_lnzeta_grid = np.array([float(mplog(mpzeta(s))) for s in _s_grid])
_lnzeta_interp = interp1d(_s_grid, _lnzeta_grid, kind='cubic', fill_value='extrapolate')

# Also precompute ζ'(s)/ζ(s) = d/ds[ln ζ(s)] for reference
_dlnzeta_grid = np.gradient(_lnzeta_grid, _s_grid)
_dlnzeta_interp = interp1d(_s_grid, _dlnzeta_grid, kind='cubic', fill_value='extrapolate')

N = 5
H = 1e-5

def metric_diag(r):
    s = 1.0 + r / R_S
    return np.array([
        -(1.0 - R_S / r),
        1.0 + R_S / r,
        r**2,
        r**2,   # equatorial plane sin²θ = 1
        float(_lnzeta_interp(s)),
    ])

def compute_christoffel(r):
    g = metric_diag(r)
    gp = metric_diag(r + H)
    gm = metric_diag(r - H)
    dg = (gp - gm) / (2 * H)
    G = np.zeros((N, N, N))
    G[1, 1, 1] = dg[1] / (2 * g[1])
    for a in range(N):
        if a != 1:
            G[1, a, a] = -dg[a] / (2 * g[1])
            G[a, a, 1] = dg[a] / (2 * g[a])
            G[a, 1, a] = dg[a] / (2 * g[a])
    return G

def compute_tensors(r):
    g = metric_diag(r)
    Gamma = compute_christoffel(r)
    Gamma_p = compute_christoffel(r + H)
    Gamma_m = compute_christoffel(r - H)
    dGamma = (Gamma_p - Gamma_m) / (2 * H)

    Riemann = np.zeros((N, N, N, N))
    for rho in range(N):
        for sig in range(N):
            for mu in range(N):
                for nu in range(N):
                    val = 0.0
                    if mu == 1:
                        val += dGamma[rho, nu, sig]
                    if nu == 1:
                        val -= dGamma[rho, mu, sig]
                    for lam in range(N):
                        val += Gamma[rho, mu, lam] * Gamma[lam, nu, sig]
                        val -= Gamma[rho, nu, lam] * Gamma[lam, mu, sig]
                    Riemann[rho, sig, mu, nu] = val

    Ricci = np.zeros((N, N))
    for sig in range(N):
        for nu in range(N):
            for rho in range(N):
                Ricci[sig, nu] += Riemann[rho, sig, rho, nu]

    R_scalar = sum(Ricci[mu, mu] / g[mu] for mu in range(N))

    Einstein = np.zeros((N, N))
    for mu in range(N):
        for nu in range(N):
            Einstein[mu, nu] = Ricci[mu, nu] - 0.5 * R_scalar * (g[mu] if mu == nu else 0.0)

    return Einstein, Ricci, R_scalar, g

# ===========================================================================
# Compute across radial range
# ===========================================================================

print("Computing Einstein tensor...")
r_vals = np.concatenate([
    np.linspace(1.005, 2.0, 400),
    np.linspace(2.01, 10.0, 300),
    np.linspace(10.1, 50.0, 200),
    np.linspace(51.0, 100.0, 100),
])

G_tt = []
G_rr = []
G_thth = []
G_phph = []
G_xixi = []
R_scalars = []
g_xi_vals = []
det_spatial = []

for idx, r in enumerate(r_vals):
    if idx % 200 == 0:
        print(f"  r = {r:.3f}  ({idx}/{len(r_vals)})")
    E, Ric, Rsc, g = compute_tensors(r)
    G_tt.append(E[0, 0])
    G_rr.append(E[1, 1])
    G_thth.append(E[2, 2])
    G_phph.append(E[3, 3])
    G_xixi.append(E[4, 4])
    R_scalars.append(Rsc)
    g_xi_vals.append(g[4])
    det_spatial.append(g[1] * g[2] * g[3] * g[4])

G_tt = np.array(G_tt)
G_rr = np.array(G_rr)
G_thth = np.array(G_thth)
G_phph = np.array(G_phph)
G_xixi = np.array(G_xixi)
R_scalars = np.array(R_scalars)
g_xi_vals = np.array(g_xi_vals)
det_spatial = np.array(det_spatial)

# ===========================================================================
# Also compute the effective stress-energy from the 5th dimension
# T_μν^(eff) = G_μν^(5D) - G_μν^(4D)  (what the 5th dimension adds)
# ===========================================================================

print("Computing 4D control...")

def compute_einstein_4d(r):
    g = np.array([
        -(1.0 - R_S / r),
        1.0 + R_S / r,
        r**2,
        r**2,
    ])
    def christoffel_4d(rv):
        gv = np.array([-(1.0 - R_S/rv), 1.0 + R_S/rv, rv**2, rv**2])
        dgv = (np.array([-(1.0-R_S/(rv+H)), 1.0+R_S/(rv+H), (rv+H)**2, (rv+H)**2])
             - np.array([-(1.0-R_S/(rv-H)), 1.0+R_S/(rv-H), (rv-H)**2, (rv-H)**2])) / (2*H)
        G4 = np.zeros((4, 4, 4))
        G4[1,1,1] = dgv[1]/(2*gv[1])
        for a in range(4):
            if a != 1:
                G4[1,a,a] = -dgv[a]/(2*gv[1])
                G4[a,a,1] = dgv[a]/(2*gv[a])
                G4[a,1,a] = dgv[a]/(2*gv[a])
        return G4

    Gam = christoffel_4d(r)
    Gam_p = christoffel_4d(r + H)
    Gam_m = christoffel_4d(r - H)
    dGam = (Gam_p - Gam_m) / (2*H)

    Riem = np.zeros((4,4,4,4))
    for rho in range(4):
        for sig in range(4):
            for mu in range(4):
                for nu in range(4):
                    val = 0.0
                    if mu == 1: val += dGam[rho,nu,sig]
                    if nu == 1: val -= dGam[rho,mu,sig]
                    for lam in range(4):
                        val += Gam[rho,mu,lam]*Gam[lam,nu,sig]
                        val -= Gam[rho,nu,lam]*Gam[lam,mu,sig]
                    Riem[rho,sig,mu,nu] = val

    Ric = np.zeros((4,4))
    for sig in range(4):
        for nu in range(4):
            for rho in range(4):
                Ric[sig,nu] += Riem[rho,sig,rho,nu]

    Rsc = sum(Ric[mu,mu]/g[mu] for mu in range(4))
    Ein = np.zeros((4,4))
    for mu in range(4):
        for nu in range(4):
            Ein[mu,nu] = Ric[mu,nu] - 0.5*Rsc*(g[mu] if mu==nu else 0.0)
    return Ein

G4_tt = []
G4_rr = []
for r in r_vals:
    E4 = compute_einstein_4d(r)
    G4_tt.append(E4[0,0])
    G4_rr.append(E4[1,1])

G4_tt = np.array(G4_tt)
G4_rr = np.array(G4_rr)

# Effective stress-energy from the 5th dimension
T_eff_tt = G_tt - G4_tt
T_eff_rr = G_rr - G4_rr

# ===========================================================================
# Plot
# ===========================================================================

print("Generating figure...")

plt.style.use('dark_background')
fig = plt.figure(figsize=(16, 20))
gs = gridspec.GridSpec(4, 2, hspace=0.35, wspace=0.3,
                       left=0.08, right=0.95, top=0.93, bottom=0.04)

fig.suptitle(
    r"$G_{\mu\nu}^{(5)}$  —  Einstein Tensor from 5D Metric with Euler Prime Substrate"
    "\nPure output. No Benford analysis. Just what the metric produces.",
    fontsize=16, fontweight='bold', color='white', y=0.97
)

colors = {
    'tt': '#FF6B6B',    # red
    'rr': '#4ECDC4',    # teal
    'thth': '#45B7D1',  # blue
    'phph': '#96CEB4',  # green
    'xixi': '#FFEAA7',  # gold
    'scalar': '#DDA0DD', # plum
    'eff': '#FF9F43',   # orange
    'det': '#A29BFE',   # purple
}

# --- Panel 1: All diagonal Einstein tensor components ---
ax1 = fig.add_subplot(gs[0, :])
ax1.plot(r_vals, G_tt, color=colors['tt'], linewidth=1.5, label=r'$G_{tt}$')
ax1.plot(r_vals, G_rr, color=colors['rr'], linewidth=1.5, label=r'$G_{rr}$')
ax1.plot(r_vals, G_thth, color=colors['thth'], linewidth=1.2, label=r'$G_{\theta\theta}$')
ax1.plot(r_vals, G_phph, color=colors['phph'], linewidth=1.2, label=r'$G_{\phi\phi}$', linestyle='--')
ax1.plot(r_vals, G_xixi, color=colors['xixi'], linewidth=2.0, label=r'$G_{\hat{\zeta}\hat{\zeta}}$  (Euler prime)')
ax1.axhline(y=0, color='gray', linewidth=0.5, alpha=0.5)
ax1.axvline(x=R_S, color='red', linewidth=0.8, alpha=0.4, linestyle=':', label=r'$r = r_s$ (horizon)')
ax1.set_xlabel(r'$r / r_s$', fontsize=12)
ax1.set_ylabel(r'$G_{\mu\mu}$', fontsize=12)
ax1.set_title('Diagonal Einstein Tensor Components', fontsize=13, pad=10)
ax1.legend(fontsize=10, loc='upper right', ncol=3)
ax1.set_xlim(1.0, 20)
# Clip extreme values for readability
ymax = np.percentile(np.abs(np.concatenate([G_tt[:400], G_rr[:400], G_xixi[:400]])), 98)
ax1.set_ylim(-ymax*1.5, ymax*1.5)
ax1.grid(alpha=0.15)

# --- Panel 2: Near-horizon zoom ---
ax2 = fig.add_subplot(gs[1, 0])
mask_near = r_vals < 3.0
ax2.plot(r_vals[mask_near], G_tt[mask_near], color=colors['tt'], linewidth=1.5, label=r'$G_{tt}$')
ax2.plot(r_vals[mask_near], G_rr[mask_near], color=colors['rr'], linewidth=1.5, label=r'$G_{rr}$')
ax2.plot(r_vals[mask_near], G_xixi[mask_near], color=colors['xixi'], linewidth=2.0, label=r'$G_{\hat{\zeta}\hat{\zeta}}$')
ax2.axhline(y=0, color='gray', linewidth=0.5, alpha=0.5)
ax2.axvline(x=R_S, color='red', linewidth=0.8, alpha=0.4, linestyle=':')
ax2.set_xlabel(r'$r / r_s$', fontsize=11)
ax2.set_ylabel(r'$G_{\mu\mu}$', fontsize=11)
ax2.set_title('Near-Horizon Detail', fontsize=12, pad=8)
ax2.legend(fontsize=9)
ax2.grid(alpha=0.15)

# --- Panel 3: The 5th dimension metric component g_ξξ = ln ζ(s) ---
ax3 = fig.add_subplot(gs[1, 1])
ax3.plot(r_vals, g_xi_vals, color=colors['xixi'], linewidth=2.0)
# Mark key values
r_horizon = R_S
s_horizon = 1.0 + r_horizon / R_S  # s=2
g_at_horizon = float(_lnzeta_interp(s_horizon))
ax3.plot(r_horizon, g_at_horizon, 'o', color='red', markersize=8, zorder=5)
ax3.annotate(f'  horizon: ln ζ(2) = ln(π²/6) ≈ {g_at_horizon:.3f}',
             xy=(r_horizon, g_at_horizon), fontsize=9, color='red',
             xytext=(3, g_at_horizon+0.05))
ax3.set_xlabel(r'$r / r_s$', fontsize=11)
ax3.set_ylabel(r'$g_{\hat{\zeta}} = \ln\,\zeta(s)$', fontsize=11)
ax3.set_title(r'5th Dimension: $\ln\,\zeta(s)$,   $s = 1 + r/r_s$', fontsize=12, pad=8)
ax3.set_xlim(1.0, 50)
ax3.grid(alpha=0.15)

# --- Panel 4: Ricci scalar ---
ax4 = fig.add_subplot(gs[2, 0])
ax4.plot(r_vals, R_scalars, color=colors['scalar'], linewidth=1.5)
ax4.axhline(y=0, color='gray', linewidth=0.5, alpha=0.5)
ax4.axvline(x=R_S, color='red', linewidth=0.8, alpha=0.4, linestyle=':')
ax4.set_xlabel(r'$r / r_s$', fontsize=11)
ax4.set_ylabel(r'$R$', fontsize=11)
ax4.set_title(r'Ricci Scalar $R = g^{\mu\nu} R_{\mu\nu}$', fontsize=12, pad=8)
ax4.set_xlim(1.0, 20)
ax4.grid(alpha=0.15)

# --- Panel 5: Effective stress-energy from the 5th dimension ---
ax5 = fig.add_subplot(gs[2, 1])
ax5.plot(r_vals, T_eff_tt, color=colors['eff'], linewidth=1.5, label=r'$\Delta G_{tt} = G_{tt}^{(5D)} - G_{tt}^{(4D)}$')
ax5.plot(r_vals, T_eff_rr, color=colors['rr'], linewidth=1.5, label=r'$\Delta G_{rr} = G_{rr}^{(5D)} - G_{rr}^{(4D)}$')
ax5.axhline(y=0, color='gray', linewidth=0.5, alpha=0.5)
ax5.axvline(x=R_S, color='red', linewidth=0.8, alpha=0.4, linestyle=':')
ax5.set_xlabel(r'$r / r_s$', fontsize=11)
ax5.set_ylabel(r'$\Delta G_{\mu\mu}$', fontsize=11)
ax5.set_title('What the 5th Dimension Adds (5D minus 4D)', fontsize=12, pad=8)
ax5.legend(fontsize=9)
ax5.set_xlim(1.0, 20)
ax5.grid(alpha=0.15)

# --- Panel 6: Spatial determinant ---
ax6 = fig.add_subplot(gs[3, 0])
ax6.plot(r_vals, det_spatial, color=colors['det'], linewidth=1.5)
ax6.axhline(y=0.4068, color=colors['xixi'], linewidth=1.0, linestyle='--', alpha=0.7,
            label=r'$\|P\|_2^2 = 0.4068$ (Benford floor — not enforced)')
ax6.axvline(x=R_S, color='red', linewidth=0.8, alpha=0.4, linestyle=':')
ax6.set_xlabel(r'$r / r_s$', fontsize=11)
ax6.set_ylabel(r'$\det(g_{\mathrm{spatial}})$', fontsize=11)
ax6.set_title(r'Spatial Determinant: $g_{rr} \cdot r^2 \cdot r^2 \cdot g_{\hat{\zeta}}$', fontsize=12, pad=8)
ax6.legend(fontsize=9)
ax6.set_xlim(1.0, 10)
ax6.set_yscale('log')
ax6.grid(alpha=0.15)

# --- Panel 7: First-digit distribution of |G_μν| values (raw histogram, no Benford overlay) ---
ax7 = fig.add_subplot(gs[3, 1])

# Collect all |G_μν| values
all_G = []
for arr in [G_tt, G_rr, G_thth, G_phph, G_xixi]:
    for v in arr:
        av = abs(v)
        if av > 1e-15:
            all_G.append(av)

# Extract first digits
def first_digit(x):
    if x <= 0 or not np.isfinite(x):
        return None
    return int(f"{x:.15e}"[0])

digits = [first_digit(v) for v in all_G]
digits = [d for d in digits if d is not None]

counts = {d: 0 for d in range(1, 10)}
for d in digits:
    counts[d] += 1
total = len(digits)
freqs = {d: counts[d]/total for d in range(1, 10)}

bar_colors = ['#FF6B6B', '#FF9F43', '#FFEAA7', '#96CEB4', '#4ECDC4',
              '#45B7D1', '#A29BFE', '#DDA0DD', '#FF85A1']
bars = ax7.bar(range(1, 10), [freqs[d] for d in range(1, 10)],
               color=bar_colors, alpha=0.85, edgecolor='white', linewidth=0.5)

# Add percentage labels
for d in range(1, 10):
    ax7.text(d, freqs[d] + 0.005, f'{freqs[d]:.1%}', ha='center', va='bottom',
             fontsize=8, color='white')

ax7.set_xlabel('First Significant Digit', fontsize=11)
ax7.set_ylabel('Frequency', fontsize=11)
ax7.set_title(f'First-Digit Distribution of |G_μν| (n={total})', fontsize=12, pad=8)
ax7.set_xticks(range(1, 10))
ax7.grid(alpha=0.15, axis='y')

# --- Info box ---
fig.text(0.5, 0.005,
         r"Metric: $g_{\mu\nu}^{(5)} = \mathrm{diag}\!\left(-(1-r_s/r),\; 1+r_s/r,\; r^2,\; r^2,\; \ln\zeta(s)\right)$"
         r"     |     $s(r) = 1 + r/r_s$"
         r"     |     $\zeta(s) = \prod_{p\;\mathrm{prime}} \frac{1}{1-p^{-s}} = \sum_{n=1}^{\infty} \frac{1}{n^s}$",
         ha='center', fontsize=10, color='#888888',
         fontstyle='italic')

outpath = '/home/jackwayne/Desktop/Projects/Benford_Fun/results/einstein_tensor_5d_pure.png'
plt.savefig(outpath, dpi=180, facecolor='#0D1117', edgecolor='none')
print(f"Saved: {outpath}")
plt.close()
