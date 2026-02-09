#!/usr/bin/env python3
"""
Schwarzschild-Benford Simulation in Painlevé-Gullstrand Coordinates

Uses Painlevé-Gullstrand (PG) coordinates instead of Schwarzschild.
On a constant Painlevé-time slice the spatial metric is:

    ds² = dr² + r² dΩ²          (flat Euclidean 3-space)

    g_rr   = 1       (everywhere — no horizon singularity)
    g_θθ   = r²
    g_φφ   = r²      (equatorial)

The gravitational physics lives entirely in the time-space cross terms
(g_Tr = √(r_s/r)), which we discard: time is emergent in our framework.

Advantage over Schwarzschild: the event horizon is unremarkable in the
spatial metric.  There is no coordinate spike.  The chart shows clean
physics from far away through the horizon to the singularity.

Reference: C. Riner (2026)
"""

import numpy as np
import matplotlib.pyplot as plt
import json
import os

# ── Benford Distribution ────────────────────────────────────────────

def benford_probs():
    """P(d) = log10(1 + 1/d) for d = 1,...,9"""
    return np.array([np.log10(1.0 + 1.0 / d) for d in range(1, 10)])


def benford_floor_value():
    """L2 norm of the Benford probability vector."""
    P = benford_probs()
    return float(np.sqrt(np.sum(P ** 2)))


# ── Painlevé-Gullstrand Spatial Metric (r_s = 1, equatorial) ────────

def pg_metric(r):
    """
    Spatial metric on a T_PG = const slice.  Flat Euclidean.
    Returns g_rr, g_theta, g_phi.
    """
    g_rr = np.ones_like(r)      # identically 1
    g_theta = r ** 2
    g_phi = r ** 2
    return g_rr, g_theta, g_phi


def pg_derivatives(r):
    """Exact dg/dr for the PG spatial metric."""
    dg_rr = np.zeros_like(r)    # d/dr(1) = 0
    dg_theta = 2.0 * r
    dg_phi = 2.0 * r
    return dg_rr, dg_theta, dg_phi


def composite_rate(dg_rr, dg_theta, dg_phi):
    """L2-norm of the vector of rates — emergent time."""
    return np.sqrt(dg_rr ** 2 + dg_theta ** 2 + dg_phi ** 2)


# ── Benford Floor Modification ──────────────────────────────────────

def apply_floor(g_rr, g_theta, g_phi, floor_val):
    """
    Where |det(g_spatial)| < floor, uniformly scale component magnitudes
    so the determinant equals the floor.  Signs preserved.
    """
    det = g_rr * g_theta * g_phi
    abs_det = np.abs(det)

    m_rr = g_rr.copy()
    m_theta = g_theta.copy()
    m_phi = g_phi.copy()

    below = abs_det < floor_val
    if np.any(below):
        s = (floor_val / abs_det[below]) ** (1.0 / 3.0)
        m_rr[below] *= s
        m_theta[below] *= s
        m_phi[below] *= s

    return m_rr, m_theta, m_phi


# ── Load Causal-Set Delta_B Data ────────────────────────────────────

def load_cs_data(path):
    """Return (r_array, delta_b_array) for positive-r CS infalling points."""
    with open(path) as f:
        data = json.load(f)
    entries = data["infalling_observer"]["Causal Set"]
    r_vals, db_vals = [], []
    for e in entries:
        if e["r_ratio"] > 0:
            r_vals.append(e["r_ratio"])
            db_vals.append(e["delta_b"])
    return np.array(r_vals), np.array(db_vals)


# ── Main ────────────────────────────────────────────────────────────

def main():
    floor_val = benford_floor_value()
    P = benford_probs()

    print("Benford probabilities P(d):")
    for d, p in enumerate(P, 1):
        print(f"  d={d}: {p:.6f}")
    print(f"\nBenford floor (L2 norm of P): {floor_val:.6f}")

    # ── radial grid ─────────────────────────────────────────────
    # No need to skip r=1 — PG coordinates are smooth there.
    r = np.logspace(np.log10(10.0), np.log10(0.01), 1200)

    # ── standard PG metric ──────────────────────────────────────
    grr, gtt, gpp = pg_metric(r)
    dgrr, dgtt, dgpp = pg_derivatives(r)
    rate_std = composite_rate(dgrr, dgtt, dgpp)
    # Analytic:  rate_std = sqrt(0 + 4r² + 4r²) = 2√2 · r
    det_std = grr * gtt * gpp  # = r^4

    # ── Benford-modified metric ─────────────────────────────────
    mgrr, mgtt, mgpp = apply_floor(grr, gtt, gpp, floor_val)

    # numerical derivatives for modified components
    dm_rr = np.gradient(mgrr, r)
    dm_tt = np.gradient(mgtt, r)
    dm_pp = np.gradient(mgpp, r)
    rate_mod = composite_rate(dm_rr, dm_tt, dm_pp)

    # ── CS data ─────────────────────────────────────────────────
    cs_r, cs_db = load_cs_data("results/round_trip/black_hole_wall.json")

    # ── where does the floor activate? ──────────────────────────
    # det = r^4 < floor_val  =>  r < floor_val^(1/4)
    floor_onset = floor_val ** 0.25
    print(f"\nPG spatial determinant: det = r⁴")
    print(f"Benford floor activates at r/r_s = {floor_onset:.4f}")
    print(f"  (Schwarzschild version: r ≈ 0.665)")

    # ── interior masks ──────────────────────────────────────────
    far = r > 5
    near_h = (r > 0.95) & (r < 1.05)
    near_s = r < 0.05
    interior = r < 0.95

    # ── analysis ────────────────────────────────────────────────
    print(f"\n{'='*65}")
    print("ANALYSIS — Painlevé-Gullstrand Coordinates")
    print(f"{'='*65}")
    print(f"CS δ_B range : {cs_db.min():.5f} – {cs_db.max():.5f}")
    print(f"CS δ_B mean  : {np.mean(cs_db):.5f}")

    print(f"\nEmergent time — STANDARD PG metric (= 2√2·r):")
    print(f"  Far (r>5)            : mean rate = {rate_std[far].mean():.4f}")
    print(f"  At horizon (r≈1)     : mean rate = {rate_std[near_h].mean():.4f}")
    print(f"  Near singularity     : mean rate = {rate_std[near_s].mean():.4f}")
    print(f"  Horizon is smooth?   : "
          f"{'YES' if rate_std[near_h].std() / rate_std[near_h].mean() < 0.1 else 'NO'}")

    print(f"\nEmergent time — BENFORD-MODIFIED PG:")
    print(f"  Far (r>5)            : mean rate = {rate_mod[far].mean():.4f}")
    print(f"  At horizon (r≈1)     : mean rate = {rate_mod[near_h].mean():.4f}")
    print(f"  Near singularity     : mean rate = {rate_mod[near_s].mean():.4f}")

    # Does time still drop in standard? Rise in modified?
    r_050 = np.argmin(np.abs(r - 0.50))
    r_005 = np.argmin(np.abs(r - 0.05))
    print(f"\nStandard rate at r=0.50: {rate_std[r_050]:.4f}")
    print(f"Standard rate at r=0.05: {rate_std[r_005]:.4f}")
    print(f"  → {'Drops' if rate_std[r_005] < rate_std[r_050] else 'Rises'} "
          f"toward singularity (time {'freezes' if rate_std[r_005] < rate_std[r_050] else 'accelerates'})")

    print(f"\nModified rate at r=0.50: {rate_mod[r_050]:.4f}")
    print(f"Modified rate at r=0.05: {rate_mod[r_005]:.4f}")
    print(f"  → {'Drops' if rate_mod[r_005] < rate_mod[r_050] else 'Rises'} "
          f"toward singularity (time {'freezes' if rate_mod[r_005] < rate_mod[r_050] else 'accelerates'})")

    # Check if there's a dip at the CS delta_B scale
    cs_db_max = cs_db.max()
    std_crosses_db = np.any(rate_std[interior] < cs_db_max)
    mod_crosses_db = np.any(rate_mod[interior] < cs_db_max)
    print(f"\nCS δ_B max = {cs_db_max:.5f}")
    print(f"Standard rate drops below CS δ_B max? "
          f"{'YES' if std_crosses_db else 'NO'}")
    print(f"Modified rate drops below CS δ_B max? "
          f"{'YES' if mod_crosses_db else 'NO'}")

    if std_crosses_db:
        cross_r = r[interior][rate_std[interior] < cs_db_max]
        print(f"  Standard crosses below δ_B at r ≈ {cross_r[0]:.4f}")

    print(f"{'='*65}")

    # ═══════════════════════════════════════════════════════════════
    #  VISUALIZATION — 4 panels, dark background
    # ═══════════════════════════════════════════════════════════════

    BG = "#1a1a2e"
    fig, (ax1, ax2, ax3a, ax3b) = plt.subplots(
        4, 1, figsize=(14, 26),
        gridspec_kw={"hspace": 0.18, "top": 0.94, "bottom": 0.05,
                     "height_ratios": [1, 1, 0.85, 0.85]}
    )
    fig.patch.set_facecolor(BG)

    def style_ax(ax):
        ax.set_facecolor(BG)
        ax.tick_params(colors="white", which="both", labelsize=10)
        for spine in ax.spines.values():
            spine.set_color("#555")
        ax.yaxis.label.set_color("white")
        ax.xaxis.label.set_color("white")
        ax.title.set_color("white")
        ax.grid(True, which="both", color="#333", linewidth=0.4, alpha=0.5)

    for a in (ax1, ax2, ax3a, ax3b):
        style_ax(a)

    lw = 1.6
    al = 0.92

    def set_full_xaxis(ax, label=False):
        ax.set_xscale("log")
        ax.set_xlim(10, 0.01)
        if label:
            ax.set_xlabel(
                "r / r_s     Far away  ───────────────────►  "
                "Event Horizon  ───────────────────►  Singularity",
                fontsize=11, color="white")

    # ── Panel 1 — Standard PG Spatial Metric ─────────────────────
    ax1.set_title("Panel 1 — Painlevé-Gullstrand Spatial Metric (flat space)",
                  fontsize=13, fontweight="bold", pad=10)

    ax1.semilogy(r, grr, color="red", lw=lw, alpha=al,
                 label="g_rr = 1  (flat)")
    ax1.semilogy(r, gtt, color="#00ff00", lw=lw, alpha=al,
                 label="g_θθ = r²")
    ax1.semilogy(r, gpp, color="#4488ff", lw=lw, alpha=al,
                 ls="--", label="g_φφ = r²")
    ax1.axvline(1.0, color="white", ls="--", lw=1, alpha=0.55,
                label="Event horizon r = r_s")
    ax1.set_ylabel("Component magnitude", fontsize=11)
    ax1.set_ylim(1e-5, 1e3)
    ax1.legend(loc="upper right", fontsize=9, facecolor="#2a2a3e",
               edgecolor="#555", labelcolor="white")
    set_full_xaxis(ax1)

    # ── Panel 2 — Benford-modified PG + CS overlay ───────────────
    ax2.set_title(
        f"Panel 2 — Benford-Modified PG Metric  "
        f"(floor = L² norm = {floor_val:.4f})",
        fontsize=13, fontweight="bold", pad=10)

    ax2.semilogy(r, mgrr, color="red", lw=lw, alpha=al,
                 label="g_rr modified")
    ax2.semilogy(r, mgtt, color="#00ff00", lw=lw, alpha=al,
                 label="g_θθ modified")
    ax2.semilogy(r, mgpp, color="#4488ff", lw=lw, alpha=al,
                 ls="--", label="g_φφ modified")

    ax2.semilogy(cs_r, cs_db, color="#ffff00", lw=2.2, marker="o", ms=4,
                 alpha=0.95, zorder=5, label="CS δ_B (measured)")
    ax2.axhline(floor_val, color="orange", ls="--", lw=1.8, alpha=0.85,
                label=f"Benford floor = {floor_val:.4f}")
    ax2.axvline(1.0, color="white", ls="--", lw=1, alpha=0.55)

    ax2.set_ylabel("Magnitude / δ_B", fontsize=11)
    ax2.set_ylim(1e-5, 1e3)
    ax2.legend(loc="upper right", fontsize=9, facecolor="#2a2a3e",
               edgecolor="#555", labelcolor="white")
    set_full_xaxis(ax2)

    # ── Panel 3a — Emergent time (full range) ────────────────────
    ax3a.set_title("Panel 3a — Emergent Time — Full Range  (no horizon spike!)",
                   fontsize=13, fontweight="bold", pad=10)

    ax3a.semilogy(r, rate_std, color="#888888", lw=lw, alpha=0.8,
                  label="Standard PG  (= 2√2·r)")
    ax3a.semilogy(r, rate_mod, color="white", lw=2.2, alpha=0.95,
                  label="Benford-modified")
    ax3a.axvline(1.0, color="white", ls="--", lw=1, alpha=0.55,
                 label="Event horizon")
    ax3a.axhline(floor_val, color="orange", ls="--", lw=1.5, alpha=0.7,
                 label=f"Benford floor = {floor_val:.4f}")

    ax3a.set_ylabel("Composite rate  √(Σ(dg/dr)²)", fontsize=11)
    ax3a.legend(loc="upper right", fontsize=9, facecolor="#2a2a3e",
                edgecolor="#555", labelcolor="white")
    set_full_xaxis(ax3a)

    # ── Panel 3b — Emergent time (interior zoom) ─────────────────
    ax3b.set_title("Panel 3b — Emergent Time — Interior Only  "
                   "(Painlevé-Gullstrand, no artifacts)",
                   fontsize=13, fontweight="bold", pad=10)

    int_mask = r <= 0.95
    ax3b.semilogy(r[int_mask], rate_std[int_mask], color="#888888", lw=lw,
                  alpha=0.8, label="Standard  (time freezes → 0)")
    ax3b.semilogy(r[int_mask], rate_mod[int_mask], color="white", lw=2.2,
                  alpha=0.95, label="Benford-modified  (time accelerates)")

    # Benford floor
    ax3b.axhline(floor_val, color="orange", ls="--", lw=1.8, alpha=0.85,
                 label=f"Benford floor = {floor_val:.4f}")

    # CS delta_B overlay (interior points)
    cs_int = cs_r <= 0.95
    ax3b.semilogy(cs_r[cs_int], cs_db[cs_int], color="#ffff00", lw=2.2,
                  marker="o", ms=4, alpha=0.95, zorder=5,
                  label="CS δ_B (measured)")

    ax3b.set_xscale("log")
    ax3b.set_xlim(0.95, 0.01)
    ax3b.set_ylabel("Composite rate  √(Σ(dg/dr)²)", fontsize=11)
    ax3b.set_xlabel(
        "r / r_s     Just inside horizon  ──────────────────────────►"
        "  Singularity",
        fontsize=11, color="white")
    ax3b.legend(loc="upper right", fontsize=9, facecolor="#2a2a3e",
                edgecolor="#555", labelcolor="white")

    fig.suptitle(
        "Schwarzschild-Benford in Painlevé-Gullstrand Coordinates\n"
        "Flat spatial metric — no horizon artifact — time emergent from dg/dr",
        fontsize=15, fontweight="bold", color="white"
    )

    outpath = "results/figures/schwarzschild_benford_ef.png"
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=200, bbox_inches="tight", facecolor=BG)
    plt.close()
    print(f"\nFigure saved → {outpath}")


if __name__ == "__main__":
    main()
