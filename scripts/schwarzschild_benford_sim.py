#!/usr/bin/env python3
"""
Schwarzschild Metric with Benford Floor Simulation

Visualizes the 3 spatial Schwarzschild metric components from r/r_s = 10
down to r/r_s = 0.01 with two modifications:
  1. The determinant of the spatial metric is floored at the L2 norm of
     the Benford distribution — no hard-coded constant.
  2. Time is NOT a basis vector.  It emerges as the composite rate of
     change of the spatial components (dg/dr).

Overlays the Causal Set model's measured delta_B values from the
black_hole_wall experiment for comparison.

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
    """
    L2 norm of the Benford probability vector.
    This is the natural geometric scale derived from Benford's law —
    NOT a hard-coded constant.
    """
    P = benford_probs()
    return float(np.sqrt(np.sum(P ** 2)))


# ── Schwarzschild Spatial Metric (r_s = 1, equatorial) ──────────────

def metric_components(r):
    """Return g_rr, g_theta, g_phi at radial coordinate r (arrays ok)."""
    g_rr = 1.0 / (1.0 - 1.0 / r)       # radial — diverges at r=1
    g_theta = r ** 2                      # angular
    g_phi = r ** 2                        # angular
    return g_rr, g_theta, g_phi


def analytical_derivatives(r):
    """Exact dg/dr for the standard Schwarzschild spatial metric."""
    dg_rr = -1.0 / (r - 1.0) ** 2
    dg_theta = 2.0 * r
    dg_phi = 2.0 * r
    return dg_rr, dg_theta, dg_phi


def composite_rate(dg_rr, dg_theta, dg_phi):
    """L2-norm of the vector of rates — this IS emergent time."""
    return np.sqrt(dg_rr ** 2 + dg_theta ** 2 + dg_phi ** 2)


# ── Benford Floor Modification ──────────────────────────────────────

def apply_floor(g_rr, g_theta, g_phi, floor_val):
    """
    Where |det(g_spatial)| drops below the Benford floor, uniformly
    scale component magnitudes so the determinant equals the floor.
    Signs are preserved.
    """
    det = g_rr * g_theta * g_phi
    abs_det = np.abs(det)

    m_rr = g_rr.copy()
    m_theta = g_theta.copy()
    m_phi = g_phi.copy()

    below = abs_det < floor_val
    if np.any(below):
        # uniform cube-root scaling:  |det| * s^3 = floor  =>  s = (floor/|det|)^{1/3}
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


# ── Radial Grid (avoids the coordinate singularity at r=1) ──────────

def build_grid(n_per_segment=600):
    """Log-spaced grid from r=10 → 1.0005 (outside) + 0.9995 → 0.01 (inside)."""
    r_out = np.logspace(np.log10(10.0), np.log10(1.0005), n_per_segment)
    r_in = np.logspace(np.log10(0.9995), np.log10(0.01), n_per_segment)
    return r_out, r_in


# ── Main ────────────────────────────────────────────────────────────

def main():
    floor_val = benford_floor_value()
    P = benford_probs()

    print("Benford probabilities P(d):")
    for d, p in enumerate(P, 1):
        print(f"  d={d}: {p:.6f}")
    print(f"\nBenford floor (L2 norm of P): {floor_val:.6f}")

    # ── grids ────────────────────────────────────────────────────
    r_out, r_in = build_grid()
    r_full = np.concatenate([r_out, r_in])

    # ── standard metric ──────────────────────────────────────────
    grr_o, gtt_o, gpp_o = metric_components(r_out)
    grr_i, gtt_i, gpp_i = metric_components(r_in)

    # analytical derivatives (separate segments to avoid r=1)
    dgrr_o, dgtt_o, dgpp_o = analytical_derivatives(r_out)
    dgrr_i, dgtt_i, dgpp_i = analytical_derivatives(r_in)
    rate_std_o = composite_rate(dgrr_o, dgtt_o, dgpp_o)
    rate_std_i = composite_rate(dgrr_i, dgtt_i, dgpp_i)

    # concatenate for plotting
    grr_std = np.concatenate([grr_o, grr_i])
    gtt_std = np.concatenate([gtt_o, gtt_i])
    gpp_std = np.concatenate([gpp_o, gpp_i])
    rate_std = np.concatenate([rate_std_o, rate_std_i])

    # ── Benford-modified metric ──────────────────────────────────
    mgrr_o, mgtt_o, mgpp_o = apply_floor(grr_o, gtt_o, gpp_o, floor_val)
    mgrr_i, mgtt_i, mgpp_i = apply_floor(grr_i, gtt_i, gpp_i, floor_val)

    # numerical derivatives for modified (no closed form)
    def num_rate(r_seg, a, b, c):
        da = np.gradient(a, r_seg)
        db = np.gradient(b, r_seg)
        dc = np.gradient(c, r_seg)
        return composite_rate(da, db, dc)

    rate_mod_o = num_rate(r_out, mgrr_o, mgtt_o, mgpp_o)
    rate_mod_i = num_rate(r_in, mgrr_i, mgtt_i, mgpp_i)

    mgrr = np.concatenate([mgrr_o, mgrr_i])
    mgtt = np.concatenate([mgtt_o, mgtt_i])
    mgpp = np.concatenate([mgpp_o, mgpp_i])
    rate_mod = np.concatenate([rate_mod_o, rate_mod_i])

    # ── CS data ──────────────────────────────────────────────────
    cs_r, cs_db = load_cs_data("results/round_trip/black_hole_wall.json")

    # ── where does the floor kick in? ────────────────────────────
    det_std_full = grr_std * gtt_std * gpp_std
    floor_radius_candidates = r_full[np.abs(det_std_full) < floor_val]
    if len(floor_radius_candidates) > 0:
        floor_onset_r = floor_radius_candidates[0]  # largest r where floor activates
    else:
        floor_onset_r = None

    # ── print analysis ───────────────────────────────────────────
    print(f"\n{'='*65}")
    print("ANALYSIS")
    print(f"{'='*65}")
    if floor_onset_r is not None:
        print(f"Benford floor activates at r/r_s ≈ {floor_onset_r:.4f}")
    print(f"CS δ_B range : {cs_db.min():.5f} – {cs_db.max():.5f}")
    print(f"CS δ_B mean  : {np.mean(cs_db):.5f}")
    print(f"CS δ_B near singularity (r<0.1): mean = "
          f"{cs_db[cs_r < 0.1].mean():.5f}")

    far = r_out > 5
    near_h = (r_out > 1.001) & (r_out < 1.05)
    near_s = r_in < 0.05

    print(f"\nEmergent time — STANDARD metric:")
    print(f"  Far (r>5)            : mean rate = {rate_std_o[far].mean():.2f}")
    print(f"  Near horizon (≈1.02) : mean rate = {rate_std_o[near_h].mean():.2f}")
    print(f"  Near singularity     : mean rate = {rate_std_i[near_s].mean():.2f}")

    print(f"\nEmergent time — BENFORD-MODIFIED:")
    print(f"  Far (r>5)            : mean rate = {rate_mod_o[far].mean():.2f}")
    print(f"  Near horizon (≈1.02) : mean rate = {rate_mod_o[near_h].mean():.2f}")
    print(f"  Near singularity     : mean rate = {rate_mod_i[near_s].mean():.2f}")

    ratio_horizon = rate_std_o[near_h].mean() / rate_std_o[far].mean()
    ratio_singularity_std = rate_std_i[near_s].mean() / rate_std_o[far].mean()
    ratio_singularity_mod = rate_mod_i[near_s].mean() / rate_mod_o[far].mean()
    print(f"\n  Horizon / far ratio  (standard)  : {ratio_horizon:.1f}×")
    print(f"  Singularity / far    (standard)  : {ratio_singularity_std:.1f}×")
    print(f"  Singularity / far    (modified)  : {ratio_singularity_mod:.1f}×")

    # ── does emergent time show expected behaviour? ───────────────
    # The rate spikes at the horizon (dg_rr/dr diverges) and near
    # the singularity.  Check if the Benford floor tames the
    # singularity spike.
    print(f"\n  Standard  singularity peak rate  : {rate_std_i[near_s].max():.2f}")
    print(f"  Modified  singularity peak rate  : {rate_mod_i[near_s].max():.2f}")
    if rate_mod_i[near_s].max() < rate_std_i[near_s].max():
        print("  → Benford floor REDUCES the singularity rate spike.")
    print(f"{'='*65}")

    # ═══════════════════════════════════════════════════════════════
    #  VISUALIZATION — 4 stacked panels, dark background
    #  v2: x-axis flipped (far left → singularity right),
    #      Panel 3 split into full-range + zoomed interior
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
    alpha = 0.92

    # x-axis: far away (10) on LEFT → singularity (0.01) on RIGHT
    def set_full_xaxis(ax, label=False):
        ax.set_xscale("log")
        ax.set_xlim(10, 0.01)
        if label:
            ax.set_xlabel(
                "r / r_s     Far away  ───────────────────►  "
                "Event Horizon  ───────────────────►  Singularity",
                fontsize=11, color="white")

    # ── Panel 1 — Standard Schwarzschild ─────────────────────────
    ax1.set_title("Panel 1 — Standard Schwarzschild Spatial Metric",
                  fontsize=13, fontweight="bold", pad=10)

    ax1.semilogy(r_full, np.abs(grr_std), color="red", lw=lw, alpha=alpha,
                 label="|g_rr| radial")
    ax1.semilogy(r_full, gtt_std, color="#00ff00", lw=lw, alpha=alpha,
                 label="g_θθ angular")
    ax1.semilogy(r_full, gpp_std, color="#4488ff", lw=lw, alpha=alpha,
                 ls="--", label="g_φφ angular")
    ax1.axvline(1.0, color="white", ls="--", lw=1, alpha=0.55,
                label="Event horizon r = r_s")
    ax1.set_ylabel("Component magnitude", fontsize=11)
    ax1.set_ylim(1e-5, 1e5)
    ax1.legend(loc="upper right", fontsize=9, facecolor="#2a2a3e",
               edgecolor="#555", labelcolor="white")
    set_full_xaxis(ax1)

    # ── Panel 2 — Benford-modified + CS overlay ──────────────────
    ax2.set_title(
        f"Panel 2 — Benford-Modified Metric  "
        f"(floor = L² norm = {floor_val:.4f})",
        fontsize=13, fontweight="bold", pad=10)

    ax2.semilogy(r_full, np.abs(mgrr), color="red", lw=lw, alpha=alpha,
                 label="|g_rr| modified")
    ax2.semilogy(r_full, mgtt, color="#00ff00", lw=lw, alpha=alpha,
                 label="g_θθ modified")
    ax2.semilogy(r_full, mgpp, color="#4488ff", lw=lw, alpha=alpha,
                 ls="--", label="g_φφ modified")

    # CS delta_B overlay
    ax2.semilogy(cs_r, cs_db, color="#ffff00", lw=2.2, marker="o", ms=4,
                 alpha=0.95, zorder=5, label="CS δ_B (measured)")

    # Benford floor line
    ax2.axhline(floor_val, color="orange", ls="--", lw=1.8, alpha=0.85,
                label=f"Benford floor = {floor_val:.4f}")
    ax2.axvline(1.0, color="white", ls="--", lw=1, alpha=0.55)

    ax2.set_ylabel("Magnitude / δ_B", fontsize=11)
    ax2.set_ylim(1e-5, 1e5)
    ax2.legend(loc="upper right", fontsize=9, facecolor="#2a2a3e",
               edgecolor="#555", labelcolor="white")
    set_full_xaxis(ax2)

    # ── Panel 3a — Emergent time (full range) ────────────────────
    ax3a.set_title("Panel 3a — Emergent Time — Full Range",
                   fontsize=13, fontweight="bold", pad=10)

    ax3a.semilogy(r_full, rate_std, color="#888888", lw=lw, alpha=0.8,
                  label="Standard metric")
    ax3a.semilogy(r_full, rate_mod, color="white", lw=2.2, alpha=0.95,
                  label="Benford-modified")
    ax3a.axvline(1.0, color="white", ls="--", lw=1, alpha=0.55,
                 label="Event horizon")

    ax3a.set_ylabel("Composite rate  √(Σ(dg/dr)²)", fontsize=11)
    ax3a.legend(loc="upper right", fontsize=9, facecolor="#2a2a3e",
                edgecolor="#555", labelcolor="white")
    set_full_xaxis(ax3a)

    # ── Panel 3b — Emergent time (interior ZOOM) ─────────────────
    ax3b.set_title("Panel 3b — Emergent Time — Interior Only "
                   "(horizon spike removed)",
                   fontsize=13, fontweight="bold", pad=10)

    # plot only interior data
    ax3b.semilogy(r_in, rate_std_i, color="#888888", lw=lw, alpha=0.8,
                  label="Standard  (time freezes)")
    ax3b.semilogy(r_in, rate_mod_i, color="white", lw=2.2, alpha=0.95,
                  label="Benford-modified  (time accelerates)")

    # Benford floor line
    ax3b.axhline(floor_val, color="orange", ls="--", lw=1.8, alpha=0.85,
                 label=f"Benford floor = {floor_val:.4f}")

    # CS delta_B overlay (only interior points: r < 0.95)
    cs_interior = cs_r <= 0.95
    ax3b.semilogy(cs_r[cs_interior], cs_db[cs_interior],
                  color="#ffff00", lw=2.2, marker="o", ms=4,
                  alpha=0.95, zorder=5, label="CS δ_B (measured)")

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
        "Schwarzschild Metric with Benford Floor\n"
        "Time as Emergent Rate of Change — no g_tt",
        fontsize=16, fontweight="bold", color="white"
    )

    # ── save v2 ──────────────────────────────────────────────────
    outpath = "results/figures/schwarzschild_benford_floor_v3.png"
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=200, bbox_inches="tight", facecolor=BG)
    plt.close()
    print(f"\nFigure saved → {outpath}")


if __name__ == "__main__":
    main()
