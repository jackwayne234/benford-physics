#!/usr/bin/env python3
"""Visualizations for the Round-Trip Experiments.

Generates figures showing:
  1. ε(d) fingerprints — the reference atlas (BE, FD, MB, Planck)
  2. Dimension sweep calibration curve with inversion
  3. Eta recovery curve with inversion
  4. The Planck Wall — all 5 QG models across temperature
  5. The Whiteboard — all exotic candidates ranked by δ_B
  6. Anyon interpolation — smooth path from BE to FD
  7. Mass dial — δ_B vs m² including tachyonic territory

Saves all figures to results/round_trip/figures/
"""

import json
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

FIG_DIR = "results/round_trip/figures"

# Style
plt.rcParams.update({
    "figure.facecolor": "#0d1117",
    "axes.facecolor": "#161b22",
    "axes.edgecolor": "#30363d",
    "axes.labelcolor": "#c9d1d9",
    "text.color": "#c9d1d9",
    "xtick.color": "#8b949e",
    "ytick.color": "#8b949e",
    "grid.color": "#21262d",
    "grid.alpha": 0.8,
    "font.family": "monospace",
    "font.size": 11,
})

COLORS = {
    "BE": "#58a6ff",
    "FD": "#f78166",
    "MB": "#7ee787",
    "Planck": "#d2a8ff",
    "Standard": "#58a6ff",
    "LQG": "#7ee787",
    "GUP": "#f78166",
    "DSR": "#d2a8ff",
    "Hagedorn": "#ffd700",
}


def load_json(path):
    with open(path) as f:
        return json.load(f)


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1: ε(d) Reference Fingerprints
# ═══════════════════════════════════════════════════════════════════════════

def fig_fingerprints():
    """Plot ε(d) for the four fundamental distributions."""
    refs = {
        "Bose-Einstein": load_json("results/individual/bose_einstein_numerical.json"),
        "Fermi-Dirac": load_json("results/individual/fermi_dirac_numerical.json"),
        "Maxwell-Boltzmann": load_json("results/individual/maxwell_boltzmann_numerical.json"),
        "Planck": load_json("results/individual/planck_radiation_spectrum.json"),
    }
    colors = ["#58a6ff", "#f78166", "#7ee787", "#d2a8ff"]
    digits = list(range(1, 10))

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
    fig.suptitle("ε(d) Fingerprints — The Four Fundamental Distributions",
                 fontsize=16, fontweight="bold", y=0.97)

    for ax, (name, data), color in zip(axes.flat, refs.items(), colors):
        eps = [data["per_digit_deviation"][str(d)] for d in digits]
        ax.bar(digits, eps, color=color, alpha=0.85, edgecolor=color, linewidth=0.5)
        ax.axhline(0, color="#8b949e", linewidth=0.8, linestyle="--")
        ax.set_title(f"{name}  (δ_B = {data['delta_b']:.4f})", fontsize=13, color=color)
        ax.set_ylabel("ε(d)")
        ax.set_ylim(-0.03, 0.03)
        ax.grid(True, axis="y")

        # Sign pattern
        signs = "".join("+" if e > 0 else "−" for e in eps)
        ax.text(0.98, 0.95, signs, transform=ax.transAxes, ha="right", va="top",
                fontsize=14, color=color, fontfamily="monospace")

    for ax in axes[1]:
        ax.set_xlabel("First Digit d")
        ax.set_xticks(digits)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    path = os.path.join(FIG_DIR, "01_fingerprints.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2: Dimension Sweep Calibration Curve
# ═══════════════════════════════════════════════════════════════════════════

def fig_dimension_sweep():
    """Plot δ_B vs exponent n with inversion arrow."""
    data = load_json("results/round_trip/dimension_sweep.json")
    sweep = data["sweep"]
    ns = [r["exponent"] for r in sweep]
    dbs = [r["delta_b"] for r in sweep]

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.plot(ns, dbs, "o-", color="#d2a8ff", markersize=10, linewidth=2.5, label="δ_B(n)")
    ax.fill_between(ns, dbs, alpha=0.15, color="#d2a8ff")

    # Mark n=0 (BE) and n=3 (Planck)
    ax.plot(0, dbs[0], "s", color="#58a6ff", markersize=14, zorder=5, label=f"n=0: BE (δ_B={dbs[0]:.4f})")
    planck_db = [r["delta_b"] for r in sweep if r["exponent"] == 3][0]
    ax.plot(3, planck_db, "D", color="#ffd700", markersize=14, zorder=5, label=f"n=3: Planck (δ_B={planck_db:.4f})")

    # Inversion arrow
    ax.annotate("", xy=(3, 0.002), xytext=(3, planck_db - 0.002),
                arrowprops=dict(arrowstyle="->", color="#ffd700", lw=2))
    ax.annotate("", xy=(0.05, planck_db), xytext=(2.95, planck_db),
                arrowprops=dict(arrowstyle="->", color="#ffd700", lw=2))
    ax.axhline(planck_db, color="#ffd700", linewidth=1, linestyle=":", alpha=0.5)
    ax.axvline(3, color="#ffd700", linewidth=1, linestyle=":", alpha=0.5)
    ax.text(1.5, planck_db + 0.001, f"δ_B = {planck_db:.4f} → n = 3.000",
            color="#ffd700", fontsize=12, ha="center")

    ax.set_xlabel("Exponent n  (spatial dimensions)", fontsize=13)
    ax.set_ylabel("δ_B  (Euclidean deviation)", fontsize=13)
    ax.set_title("Dimension Sweep: B(x) = x^n / (e^x − 1)\nδ_B = 0.028 inverts to n = 3 (exact)",
                 fontsize=15, fontweight="bold")
    ax.legend(fontsize=11, loc="upper left")
    ax.grid(True, alpha=0.5)
    ax.set_xticks(ns)

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "02_dimension_sweep.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 3: Eta Recovery Curve
# ═══════════════════════════════════════════════════════════════════════════

def fig_eta_recovery():
    """Plot δ_B vs alpha with inversion."""
    data = load_json("results/round_trip/eta_recovery.json")
    sweep = data["sweep"]
    alphas = [r["alpha"] for r in sweep]
    dbs = [r["delta_b"] for r in sweep]

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.plot(alphas, dbs, "o-", color="#f78166", markersize=10, linewidth=2.5, label="δ_B(α)")
    ax.fill_between(alphas, dbs, alpha=0.15, color="#f78166")

    # Mark α=0 (BE) and α=1 (FD)
    ax.plot(0, dbs[0], "s", color="#58a6ff", markersize=14, zorder=5, label=f"α=0: BE (δ_B={dbs[0]:.4f})")
    fd_db = [r["delta_b"] for r in sweep if r["alpha"] == 1.0][0]
    ax.plot(1.0, fd_db, "D", color="#f78166", markersize=14, zorder=5, label=f"α=1: FD (δ_B={fd_db:.4f})")

    # Inversion annotation
    ax.axhline(fd_db, color="#f78166", linewidth=1, linestyle=":", alpha=0.5)
    ax.axvline(1.0, color="#f78166", linewidth=1, linestyle=":", alpha=0.5)
    ax.text(0.5, fd_db + 0.0005, f"δ_B = {fd_db:.4f} → α = 1.000\nη(1) = ln(2) = 0.6931",
            color="#f78166", fontsize=11, ha="center")

    ax.set_xlabel("α  (Dirichlet eta modulation strength)", fontsize=13)
    ax.set_ylabel("δ_B  (Euclidean deviation)", fontsize=13)
    ax.set_title("Eta Recovery: n(x) = 1/(e^x−1) − α·2/(e^{2x}−1)\nδ_B = 0.0117 inverts to α = 1 (exact), confirming η(1) = ln(2)",
                 fontsize=15, fontweight="bold")
    ax.legend(fontsize=11, loc="upper left")
    ax.grid(True, alpha=0.5)

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "03_eta_recovery.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 4: The Planck Wall
# ═══════════════════════════════════════════════════════════════════════════

def fig_planck_wall():
    """Plot δ_B vs T/T_P for all five QG models."""
    # Use high-res data if available
    hires_path = "results/round_trip/planck_wall_hires.json"
    lores_path = "results/round_trip/planck_wall.json"
    data = load_json(hires_path if os.path.exists(hires_path) else lores_path)

    fig, ax = plt.subplots(figsize=(14, 8))

    model_styles = {
        "Standard": {"color": "#58a6ff", "marker": "o", "ls": "-"},
        "LQG": {"color": "#7ee787", "marker": "s", "ls": "--"},
        "GUP": {"color": "#f78166", "marker": "^", "ls": "-."},
        "DSR": {"color": "#d2a8ff", "marker": "D", "ls": ":"},
        "Hagedorn": {"color": "#ffd700", "marker": "*", "ls": "-"},
    }

    for model_name, style in model_styles.items():
        entries = data["results"][model_name]
        temps = [e["T_planck_units"] for e in entries if e["computable"]]
        dbs = [e["delta_b"] for e in entries if e["computable"]]
        is_hires = len(temps) > 30
        ms = (6 if model_name == "Hagedorn" else 3) if is_hires else (12 if model_name == "Hagedorn" else 8)
        lw = 2.5 if model_name == "Hagedorn" else 1.5
        ax.plot(temps, dbs, marker=style["marker"], linestyle=style["ls"],
                color=style["color"], markersize=ms, linewidth=lw,
                label=model_name, alpha=0.9)

    # The wall / Big Bang
    ax.axvline(1.0, color="#ff4444", linewidth=3, linestyle="-", alpha=0.6, label="◉ BIG BANG")
    ax.axvspan(0.8, 1.2, alpha=0.08, color="#ff4444")
    ax.text(1.15, 0.90, "◉ BIG BANG", color="#ff4444",
            fontsize=14, fontweight="bold", alpha=0.85, transform=ax.get_xaxis_transform(), va="top")

    # Timeline context
    ax.text(0.003, 0.90, "Our universe\n(after Big Bang)", color="#8b949e",
            fontsize=10, ha="center", va="top", transform=ax.get_xaxis_transform())
    ax.text(30, 0.90, "Before the\nBig Bang?", color="#c9d1d9",
            fontsize=11, ha="center", fontweight="bold", va="top", transform=ax.get_xaxis_transform())

    # Conformance zones
    ax.axhspan(0, 0.025, alpha=0.06, color="#7ee787")
    ax.text(0.002, 0.012, "CONFORMS", color="#7ee787", fontsize=9, alpha=0.7)
    ax.axhspan(0.025, 0.1, alpha=0.04, color="#ffd700")
    ax.text(0.002, 0.06, "MARGINAL", color="#ffd700", fontsize=9, alpha=0.7)

    ax.set_xscale("log")
    ax.set_xlabel("T / T_Planck  (higher T = further back → through Big Bang → before)", fontsize=12)
    ax.set_ylabel("δ_B  (Euclidean deviation)", fontsize=14)
    ax.set_title("The Planck Wall: Which Quantum Gravity Models Survive the Big Bang?\nδ_B across temperature for five proposals",
                 fontsize=16, fontweight="bold")
    ax.legend(fontsize=12, loc="upper right", framealpha=0.3)
    ax.grid(True, alpha=0.4)
    ax.set_ylim(-0.02, 0.95)
    ax.set_xlim(0.0008, 150)

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "04_planck_wall.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 5: The Whiteboard — All candidates ranked
# ═══════════════════════════════════════════════════════════════════════════

def fig_whiteboard():
    """Horizontal bar chart of all candidates sorted by δ_B."""
    data = load_json("results/round_trip/whiteboard.json")
    candidates = data["candidates"]

    # Separate EXISTS and UNDEFINED
    exists = [c for c in candidates if c["status"] == "EXISTS"]
    undefined = [c for c in candidates if c["status"] == "UNDEFINED"]

    # Sort exists by δ_B
    exists.sort(key=lambda x: x["delta_b"])

    names = [c["name"] for c in exists]
    dbs = [c["delta_b"] for c in exists]

    # Color by category
    def get_color(c):
        cat = c.get("category", "")
        if "Fractional" in cat:
            return "#d2a8ff"
        if "Exotic" in cat:
            return "#f78166"
        if "Black Hole" in cat:
            return "#ffd700"
        if "Dark" in cat:
            return "#7ee787"
        if "BSM" in cat or "QFT" in cat or "Gravity" in cat:
            return "#58a6ff"
        return "#8b949e"

    colors = [get_color(c) for c in exists]

    fig, ax = plt.subplots(figsize=(14, 12))

    bars = ax.barh(range(len(names)), dbs, color=colors, alpha=0.85, edgecolor="#30363d")

    # Add δ_B labels
    for i, (db, name) in enumerate(zip(dbs, names)):
        if db < 0.3:
            ax.text(db + 0.005, i, f"{db:.4f}", va="center", fontsize=9, color="#c9d1d9")
        else:
            ax.text(db - 0.01, i, f"{db:.3f}", va="center", ha="right", fontsize=9, color="#0d1117")

    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=10)
    ax.invert_yaxis()

    # Zone markers
    ax.axvline(0.025, color="#7ee787", linewidth=2, linestyle="--", alpha=0.5, label="CONFORMS threshold")
    ax.axvline(0.1, color="#ffd700", linewidth=2, linestyle="--", alpha=0.5, label="DEVIATES threshold")

    # Add UNDEFINED entries as text
    if undefined:
        undef_text = "UNDEFINED (do not exist):\n" + "\n".join(f"  ✗ {c['name']}" for c in undefined)
        ax.text(0.55, 0.02, undef_text, transform=ax.transAxes, fontsize=10,
                color="#f78166", va="bottom", fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.5", facecolor="#1a1a2e", edgecolor="#f78166", alpha=0.8))

    ax.set_xlabel("δ_B  (Euclidean deviation from Benford)", fontsize=13)
    ax.set_title("The Whiteboard: Exotic Physics Existence Filter\nAll candidates ranked by δ_B",
                 fontsize=16, fontweight="bold")
    ax.legend(fontsize=11, loc="lower right")
    ax.grid(True, axis="x", alpha=0.4)
    ax.set_xlim(0, max(dbs) * 1.15)

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "05_whiteboard.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 6: Anyon Interpolation — ε(d) evolution
# ═══════════════════════════════════════════════════════════════════════════

def fig_anyons():
    """Show ε(d) fingerprints evolving from BE to FD through anyonic statistics."""
    data = load_json("results/round_trip/whiteboard.json")
    candidates = data["candidates"]

    anyons = [c for c in candidates if "Anyon" in c["name"] and c["status"] == "EXISTS"]
    anyons.sort(key=lambda c: float(c["name"].split("g=")[1].split(" ")[0]))

    digits = list(range(1, 10))
    colors_grad = ["#58a6ff", "#7eb8ff", "#a4c8ff", "#d4a0ff", "#f78166"]

    fig, axes = plt.subplots(1, 5, figsize=(18, 5), sharey=True)
    fig.suptitle("Anyon Interpolation: ε(d) evolves smoothly from BE → FD",
                 fontsize=16, fontweight="bold", y=1.02)

    for ax, anyon, color in zip(axes, anyons, colors_grad):
        eps = [anyon["epsilon_d"][str(d)] for d in digits]
        ax.bar(digits, eps, color=color, alpha=0.85, edgecolor=color, linewidth=0.5)
        ax.axhline(0, color="#8b949e", linewidth=0.8, linestyle="--")
        g_val = anyon["name"].split("g=")[1].split(" ")[0]
        label = anyon["name"].split("(")[1].rstrip(")")
        ax.set_title(f"g={g_val}\n({label})\nδ_B={anyon['delta_b']:.4f}", fontsize=10, color=color)
        ax.set_xlabel("Digit d")
        ax.set_xticks(digits)
        ax.grid(True, axis="y", alpha=0.4)

    axes[0].set_ylabel("ε(d)")
    axes[0].set_ylim(-0.015, 0.015)

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "06_anyons.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 7: Planck Wall — Hagedorn spotlight
# ═══════════════════════════════════════════════════════════════════════════

def fig_hagedorn_spotlight():
    """Spotlight on Hagedorn: the phase transition through the wall."""
    hires_path = "results/round_trip/planck_wall_hires.json"
    lores_path = "results/round_trip/planck_wall.json"
    data = load_json(hires_path if os.path.exists(hires_path) else lores_path)
    entries = data["results"]["Hagedorn"]

    temps = [e["T_planck_units"] for e in entries if e["computable"]]
    dbs = [e["delta_b"] for e in entries if e["computable"]]

    fig, ax = plt.subplots(figsize=(14, 7))

    # Fill regions
    ax.axvspan(0.0005, 1.0, alpha=0.06, color="#f78166", label="Pre-wall")
    ax.axvspan(1.0, 150, alpha=0.06, color="#7ee787", label="Post-wall")

    # Plot — smaller markers for high-res data
    ms = 5 if len(temps) > 30 else 12
    ax.plot(temps, dbs, "o-", color="#ffd700", markersize=ms, linewidth=2.5, zorder=5)

    # Annotate key points
    peak_idx = np.argmax(dbs)
    ax.annotate(f"Peak chaos\nδ_B = {dbs[peak_idx]:.3f}", xy=(temps[peak_idx], dbs[peak_idx]),
                xytext=(temps[peak_idx] * 3, dbs[peak_idx]),
                fontsize=11, color="#ffd700",
                arrowprops=dict(arrowstyle="->", color="#ffd700", lw=1.5))

    # Post-wall conformance
    post = [(t, d) for t, d in zip(temps, dbs) if t > 2]
    if post:
        mean_post = np.mean([d for _, d in post])
        ax.annotate(f"Settles to δ_B ≈ {mean_post:.3f}\n(near-perfect Benford)", xy=(10, mean_post),
                    xytext=(20, 0.1), fontsize=11, color="#7ee787",
                    arrowprops=dict(arrowstyle="->", color="#7ee787", lw=1.5))

    # The wall / Big Bang
    ax.axvline(1.0, color="#ff4444", linewidth=3, alpha=0.6)
    ax.text(1.05, 0.50, "◉ BIG BANG\n   (T = T_Planck)", color="#ff4444",
            fontsize=13, fontweight="bold", alpha=0.85)

    # Timeline annotations
    ax.annotate("Today's universe\n(T << T_P)", xy=(0.001, 0.01), fontsize=10,
                color="#8b949e", ha="center")
    ax.annotate("Before the\nBig Bang?", xy=(50, 0.01), fontsize=11,
                color="#7ee787", ha="center", fontweight="bold")

    # Arrow showing direction of time (right to left = going back in time)
    ax.annotate("", xy=(0.003, -0.015), xytext=(0.5, -0.015),
                arrowprops=dict(arrowstyle="<-", color="#8b949e", lw=1.5))
    ax.text(0.03, -0.013, "← direction of time", color="#8b949e", fontsize=9)
    ax.annotate("", xy=(2, -0.015), xytext=(100, -0.015),
                arrowprops=dict(arrowstyle="<-", color="#7ee787", lw=1.5))
    ax.text(10, -0.013, "← pre-Big-Bang era", color="#7ee787", fontsize=9)

    # Conformance threshold
    ax.axhline(0.025, color="#7ee787", linewidth=1.5, linestyle="--", alpha=0.5)
    ax.text(0.002, 0.027, "CONFORMS threshold", color="#7ee787", fontsize=9, alpha=0.7)

    ax.set_xscale("log")
    ax.set_xlabel("T / T_Planck  (higher T = further back in time → through the Big Bang → before)",
                  fontsize=12)
    ax.set_ylabel("δ_B  (Euclidean deviation)", fontsize=14)
    ax.set_title("Hagedorn (String Theory): Through the Big Bang and Out the Other Side\nThe only model that emerges cleaner after the singularity",
                 fontsize=15, fontweight="bold")
    ax.grid(True, alpha=0.4)
    ax.set_ylim(-0.02, max(dbs) * 1.15)
    ax.set_xlim(0.0008, 150)

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "07_hagedorn_spotlight.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 8: Black Hole Wall — Infalling observer
# ═══════════════════════════════════════════════════════════════════════════

def fig_black_hole_wall():
    """Plot δ_B vs r/r_s — full journey through a black hole including post-singularity."""
    bh_path = "results/round_trip/black_hole_wall.json"
    if not os.path.exists(bh_path):
        print(f"  Skipping black hole wall — {bh_path} not found")
        return
    data = load_json(bh_path)

    model_styles = {
        "Standard":     {"color": "#58a6ff", "ls": "-",  "lw": 2.0},
        "LQG":          {"color": "#7ee787", "ls": "--", "lw": 2.0},
        "GUP":          {"color": "#f78166", "ls": "-.", "lw": 2.0},
        "DSR":          {"color": "#d2a8ff", "ls": ":",  "lw": 2.0},
        "Hagedorn":     {"color": "#ffd700", "ls": "-",  "lw": 3.0},
        "Causal Set":   {"color": "#00d4ff", "ls": "-",  "lw": 3.0},
        "Asym. Safety": {"color": "#ff69b4", "ls": "--", "lw": 2.0},
        "Horava-Lif.":  {"color": "#ff6347", "ls": "-.", "lw": 2.0},
        "Noncommut.":   {"color": "#98fb98", "ls": ":",  "lw": 2.0},
        "CDT":          {"color": "#dda0dd", "ls": "--", "lw": 2.0},
    }

    # Single wide figure: full journey
    fig, ax = plt.subplots(figsize=(22, 10))
    fig.suptitle("Black Hole: Through the Horizon, Through the Singularity, Out the Other Side",
                 fontsize=20, fontweight="bold", y=0.97)

    for model_name, style in model_styles.items():
        entries = data["infalling_observer"][model_name]
        pts = [(e["r_ratio"], e["delta_b"]) for e in entries if e.get("computable")]
        if pts:
            rs, dbs = zip(*sorted(pts, reverse=True))
            ax.plot(rs, dbs, color=style["color"], linestyle=style["ls"],
                    linewidth=style["lw"], marker="o", markersize=5,
                    label=model_name, alpha=0.9)

    # Event horizon
    ax.axvline(1.0, color="#ff4444", linewidth=4, alpha=0.6)
    ax.axvspan(0.95, 1.05, alpha=0.08, color="#ff4444")
    ax.text(1.0, 0.92, "EVENT\nHORIZON", color="#ff4444",
            fontsize=16, fontweight="bold", ha="center",
            transform=ax.get_xaxis_transform())

    # Singularity
    ax.axvline(0.0, color="#ffffff", linewidth=3, alpha=0.5)
    ax.axvspan(-0.03, 0.03, alpha=0.08, color="#ffffff")
    ax.text(0.0, 0.92, "SINGULARITY", color="#ffffff",
            fontsize=16, fontweight="bold", ha="center", alpha=0.8,
            transform=ax.get_xaxis_transform())

    # Zone labels
    ax.text(5.0, 0.82, "Outside\n(approaching BH)", color="#8b949e",
            fontsize=14, ha="center", transform=ax.get_xaxis_transform())
    ax.text(0.5, 0.82, "Inside\n(falling to singularity)", color="#f78166",
            fontsize=14, ha="center", transform=ax.get_xaxis_transform())
    ax.text(-0.5, 0.82, "Other Side?\n(post-singularity bounce)", color="#7ee787",
            fontsize=14, ha="center", fontweight="bold",
            transform=ax.get_xaxis_transform())

    # Conformance zone
    ax.axhspan(0, 0.025, alpha=0.05, color="#7ee787")
    ax.axhline(0.025, color="#7ee787", linewidth=1.5, linestyle="--", alpha=0.4)
    ax.text(10, 0.030, "CONFORMS", color="#7ee787", fontsize=13, alpha=0.6)

    ax.set_xlabel("r / r_s  (far away → horizon → singularity → other side)",
                  fontsize=15)
    ax.set_ylabel("δ_B  (Euclidean deviation)", fontsize=15)
    ax.set_xlim(10.5, -1.1)  # Far left to post-singularity right
    ax.set_ylim(-0.02, 0.95)
    ax.grid(True, alpha=0.4)
    ax.tick_params(labelsize=13)

    # Legend outside the plot area
    ax.legend(fontsize=13, loc="upper left", bbox_to_anchor=(0.01, 0.78),
              framealpha=0.4, ncol=2)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    path = os.path.join(FIG_DIR, "08_black_hole_wall.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 10: Hawking Radiation vs Causal Set — focused comparison
# ═══════════════════════════════════════════════════════════════════════════

def fig_hawking_vs_causal_set():
    """Focused comparison: Hawking radiation and Causal Set Theory through a black hole."""
    bh_path = "results/round_trip/black_hole_wall.json"
    wb_path = "results/round_trip/whiteboard.json"
    if not os.path.exists(bh_path) or not os.path.exists(wb_path):
        print(f"  Skipping Hawking vs CS — missing data files")
        return

    bh_data = load_json(bh_path)
    wb_data = load_json(wb_path)

    # ── Pull Hawking radiation ε(d) signatures from the whiteboard ──
    hawking_variants = {}
    for c in wb_data["candidates"]:
        if "Hawking" in c.get("name", "") and c["status"] == "EXISTS":
            hawking_variants[c["name"]] = c

    # ── Pull Causal Set BH journey ──
    cs_entries = bh_data["infalling_observer"]["Causal Set"]
    cs_pts = [(e["r_ratio"], e["delta_b"]) for e in cs_entries if e.get("computable")]
    cs_pts.sort(key=lambda p: -p[0])  # far → close

    # Hawking color map
    hawk_colors = {
        "0.5": "#ffd700", "1.0": "#ffaa00",
        "2.0": "#ff7700", "5.0": "#ff4400",
    }

    # ── Build figure: 2 rows × 2 cols ──
    fig, ((ax_top_l, ax_top_r), (ax_bot_l, ax_bot_r)) = plt.subplots(
        2, 2, figsize=(24, 14),
        gridspec_kw={"width_ratios": [2, 1], "hspace": 0.32, "wspace": 0.25}
    )

    fig.suptitle("Hawking Radiation  vs  Causal Set Theory",
                 fontsize=24, fontweight="bold", y=0.98, color="#00d4ff")

    # ════════════════════════════════════════════════════════════
    # TOP LEFT — Full BH journey, zoomed to comparison zone
    # ════════════════════════════════════════════════════════════
    cs_r, cs_db = zip(*cs_pts)
    ax_top_l.plot(cs_r, cs_db, color="#00d4ff", linewidth=3.5, marker="o",
                  markersize=7, label="Causal Set (infalling)", alpha=0.95, zorder=5)

    # Hawking δ_B reference lines — stagger labels to avoid overlap
    hawk_sorted = sorted(hawking_variants.items(),
                         key=lambda x: x[1]["delta_b"])
    label_x_positions = [8.0, 6.0, 4.0, 2.0]  # stagger along x
    for i, (name, cdata) in enumerate(hawk_sorted):
        wc = name.split("=")[1]
        color = hawk_colors.get(wc, "#ffd700")
        db = cdata["delta_b"]
        ax_top_l.axhline(db, color=color, linewidth=2.5, linestyle="--", alpha=0.6)
        lx = label_x_positions[i] if i < len(label_x_positions) else 5.0
        ax_top_l.annotate(f"Hawking ω_c={wc}", xy=(lx, db),
                          xytext=(lx, db + 0.008),
                          color=color, fontsize=13, fontweight="bold",
                          ha="center", va="bottom",
                          arrowprops=dict(arrowstyle="-", color=color, alpha=0.5),
                          bbox=dict(boxstyle="round,pad=0.15", facecolor="#161b22",
                                    edgecolor=color, alpha=0.85))

    # Event horizon + singularity
    ax_top_l.axvline(1.0, color="#ff4444", linewidth=4, alpha=0.5)
    ax_top_l.axvspan(0.95, 1.05, alpha=0.06, color="#ff4444")
    ax_top_l.text(1.0, 0.93, "EVENT\nHORIZON", color="#ff4444",
                  fontsize=16, fontweight="bold", ha="center",
                  transform=ax_top_l.get_xaxis_transform())

    ax_top_l.axvline(0.0, color="#ffffff", linewidth=3, alpha=0.4)
    ax_top_l.text(0.0, 0.93, "SINGULARITY", color="#ffffff",
                  fontsize=14, fontweight="bold", ha="center", alpha=0.7,
                  transform=ax_top_l.get_xaxis_transform())

    # Conformance zone
    ax_top_l.axhspan(0, 0.025, alpha=0.05, color="#7ee787")
    ax_top_l.axhline(0.025, color="#7ee787", linewidth=1.5, linestyle="--", alpha=0.3)
    ax_top_l.text(8.5, 0.027, "Benford conformance zone", color="#7ee787",
                  fontsize=13, alpha=0.6)

    ax_top_l.set_xlabel("r / r_s  (far  →  horizon  →  singularity  →  other side)",
                        fontsize=16)
    ax_top_l.set_ylabel("δ_B", fontsize=16)
    ax_top_l.set_title("Causal Set Journey Through the Black Hole\n"
                        "with Hawking radiation δ_B levels overlaid",
                        fontsize=16, fontweight="bold")
    ax_top_l.set_xlim(10.5, -1.1)
    ax_top_l.set_ylim(-0.01, 0.15)  # Zoomed to comparison zone
    ax_top_l.grid(True, alpha=0.4)
    ax_top_l.tick_params(labelsize=14)
    ax_top_l.legend(fontsize=14, loc="upper left", bbox_to_anchor=(0.55, 0.99),
                    framealpha=0.5)

    # ════════════════════════════════════════════════════════════
    # TOP RIGHT — Causal Set ε(d) at key radii
    # ════════════════════════════════════════════════════════════
    digits = list(range(1, 10))
    key_radii = [
        ("Far out (r = 10 r_s)", 10.0, "#00d4ff"),
        ("At horizon (r ≈ r_s)", 1.001, "#00ffff"),
        ("Deep inside (r = 0.1 r_s)", 0.1, "#00aa88"),
    ]
    for label, target_r, color in key_radii:
        entry = min(cs_entries, key=lambda e: abs(e["r_ratio"] - target_r))
        if entry.get("computable") and entry.get("epsilon_d"):
            eps = [entry["epsilon_d"][str(d)] for d in digits]
            ax_top_r.plot(digits, eps, "o-", color=color, linewidth=2.5,
                          markersize=8, label=f"{label}  δ_B={entry['delta_b']:.3f}",
                          alpha=0.9)
    ax_top_r.axhline(0, color="#8b949e", linewidth=1, linestyle="--")
    ax_top_r.set_xlabel("First Digit d", fontsize=15)
    ax_top_r.set_ylabel("ε(d)", fontsize=15)
    ax_top_r.set_title("Causal Set Fingerprints\nat 3 positions", fontsize=16,
                        fontweight="bold", color="#00d4ff")
    ax_top_r.set_xticks(digits)
    ax_top_r.tick_params(labelsize=13)
    ax_top_r.grid(True, alpha=0.4)
    ax_top_r.legend(fontsize=11, framealpha=0.4, loc="upper right")

    # ════════════════════════════════════════════════════════════
    # BOTTOM LEFT — Full BH journey, FULL scale (shows the spike)
    # ════════════════════════════════════════════════════════════
    ax_bot_l.plot(cs_r, cs_db, color="#00d4ff", linewidth=3, marker="o",
                  markersize=6, label="Causal Set", alpha=0.95, zorder=5)

    # Same Hawking lines
    for name, cdata in sorted(hawking_variants.items()):
        wc = name.split("=")[1]
        color = hawk_colors.get(wc, "#ffd700")
        db = cdata["delta_b"]
        ax_bot_l.axhline(db, color=color, linewidth=2, linestyle="--", alpha=0.5)

    ax_bot_l.axvline(1.0, color="#ff4444", linewidth=4, alpha=0.5)
    ax_bot_l.axvline(0.0, color="#ffffff", linewidth=3, alpha=0.4)

    # Zoom box showing where top-left panel is
    from matplotlib.patches import Rectangle
    rect = Rectangle((10.5, -0.01), -(10.5 + 1.1), 0.16, linewidth=2,
                      edgecolor="#00d4ff", facecolor="#00d4ff", alpha=0.08)
    ax_bot_l.add_patch(rect)
    ax_bot_l.text(5, 0.14, "← zoomed view above", color="#00d4ff",
                  fontsize=13, ha="center", alpha=0.7)

    ax_bot_l.set_xlabel("r / r_s", fontsize=16)
    ax_bot_l.set_ylabel("δ_B", fontsize=16)
    ax_bot_l.set_title("Full Scale — Singularity Spike",
                        fontsize=16, fontweight="bold")
    ax_bot_l.set_xlim(10.5, -1.1)
    max_db = max(cs_db)
    ax_bot_l.set_ylim(-0.02, max_db * 1.15)
    ax_bot_l.grid(True, alpha=0.4)
    ax_bot_l.tick_params(labelsize=14)
    ax_bot_l.legend(fontsize=14, loc="upper left", bbox_to_anchor=(0.55, 0.99),
                    framealpha=0.5)

    # ════════════════════════════════════════════════════════════
    # BOTTOM RIGHT — Hawking ε(d) at 4 greybody cutoffs
    # ════════════════════════════════════════════════════════════
    for name, cdata in sorted(hawking_variants.items()):
        wc = name.split("=")[1]
        color = hawk_colors.get(wc, "#ffd700")
        eps = [cdata["epsilon_d"][str(d)] for d in digits]
        ax_bot_r.plot(digits, eps, "s-", color=color, linewidth=2.5,
                      markersize=8, label=f"ω_c = {wc}  δ_B={cdata['delta_b']:.3f}",
                      alpha=0.9)
    ax_bot_r.axhline(0, color="#8b949e", linewidth=1, linestyle="--")
    ax_bot_r.set_xlabel("First Digit d", fontsize=15)
    ax_bot_r.set_ylabel("ε(d)", fontsize=15)
    ax_bot_r.set_title("Hawking Radiation Fingerprints\nat 4 greybody cutoffs",
                        fontsize=16, fontweight="bold", color="#ffd700")
    ax_bot_r.set_xticks(digits)
    ax_bot_r.tick_params(labelsize=13)
    ax_bot_r.grid(True, alpha=0.4)
    ax_bot_r.legend(fontsize=11, framealpha=0.4, loc="upper right")

    fig.subplots_adjust(left=0.06, right=0.97, bottom=0.06, top=0.92)
    path = os.path.join(FIG_DIR, "10_hawking_vs_causal_set.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 9: Big Bang vs Black Hole comparison
# ═══════════════════════════════════════════════════════════════════════════

def fig_bb_vs_bh():
    """Side-by-side: Big Bang wall vs Black Hole wall for the top models."""
    bh_path = "results/round_trip/black_hole_wall.json"
    bb_path = "results/round_trip/planck_wall_extended.json"
    if not os.path.exists(bh_path) or not os.path.exists(bb_path):
        print(f"  Skipping BB vs BH — missing data files")
        return

    bh_data = load_json(bh_path)
    bb_data = load_json(bb_path)

    # Focus on the top 4 most interesting models
    spotlight = ["Hagedorn", "Causal Set", "Standard", "LQG"]
    colors = {"Hagedorn": "#ffd700", "Causal Set": "#00d4ff",
              "Standard": "#58a6ff", "LQG": "#7ee787"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    fig.suptitle("Same Physics, Two Walls: Big Bang vs Black Hole\n"
                 "Do QG models behave the same at both singularities?",
                 fontsize=16, fontweight="bold", y=0.98)

    # ── Left: Big Bang (T/T_P) ──
    ax1.set_title("Big Bang Wall\n(T sweeps through T_Planck)", fontsize=13, color="#ff4444")
    ax1.axvline(1.0, color="#ff4444", linewidth=3, alpha=0.5)
    ax1.text(1.1, 0.42, "BIG\nBANG", color="#ff4444", fontsize=12, fontweight="bold")

    for name in spotlight:
        entries = bb_data["results"][name]
        pts = [(e["T_planck_units"], e["delta_b"]) for e in entries if e.get("computable")]
        if pts:
            ts, dbs = zip(*sorted(pts))
            ax1.plot(ts, dbs, "o-", color=colors[name], markersize=4,
                     linewidth=2, label=name, alpha=0.9)

    ax1.set_xscale("log")
    ax1.set_xlabel("T / T_Planck", fontsize=12)
    ax1.set_ylabel("δ_B  (Euclidean deviation)", fontsize=13)
    ax1.legend(fontsize=11, framealpha=0.3)
    ax1.grid(True, alpha=0.4)
    ax1.set_ylim(-0.02, 0.5)
    ax1.set_xlim(0.0008, 150)
    ax1.axhspan(0, 0.025, alpha=0.04, color="#7ee787")
    ax1.axhline(0.025, color="#7ee787", linewidth=1, linestyle="--", alpha=0.4)

    # ── Right: Black Hole (r/r_s, infalling) ──
    ax2.set_title("Black Hole Wall\n(r sweeps through event horizon to singularity)",
                  fontsize=13, color="#ff4444")
    ax2.axvline(1.0, color="#ff4444", linewidth=3, alpha=0.5)
    ax2.text(0.85, 0.42, "EVENT\nHORIZON", color="#ff4444", fontsize=12,
             fontweight="bold", ha="right")
    ax2.axvline(0.0, color="#ffffff", linewidth=2, alpha=0.3)
    ax2.text(0.02, 0.42, "SINGULARITY", color="#ffffff", fontsize=10, alpha=0.6)

    for name in spotlight:
        entries = bh_data["infalling_observer"][name]
        pts = [(e["r_ratio"], e["delta_b"]) for e in entries if e.get("computable")]
        if pts:
            rs, dbs = zip(*sorted(pts, reverse=True))
            ax2.plot(rs, dbs, "o-", color=colors[name], markersize=5,
                     linewidth=2, label=name, alpha=0.9)

    ax2.set_xlabel("r / r_s  (→ horizon → singularity)", fontsize=12)
    ax2.set_xlim(10.5, -0.02)  # Reversed: far on left, singularity on right
    ax2.legend(fontsize=11, framealpha=0.3)
    ax2.grid(True, alpha=0.4)
    ax2.set_ylim(-0.02, 0.5)
    ax2.axhspan(0, 0.025, alpha=0.04, color="#7ee787")
    ax2.axhline(0.025, color="#7ee787", linewidth=1, linestyle="--", alpha=0.4)

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    path = os.path.join(FIG_DIR, "09_bb_vs_bh.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("Generating visualizations...")
    print("=" * 60)
    print()

    fig_fingerprints()
    fig_dimension_sweep()
    fig_eta_recovery()
    fig_planck_wall()
    fig_whiteboard()
    fig_anyons()
    fig_hagedorn_spotlight()
    fig_black_hole_wall()
    fig_hawking_vs_causal_set()
    fig_bb_vs_bh()

    print()
    print(f"All figures saved to {FIG_DIR}/")
    print()

    # List generated files
    files = sorted(os.listdir(FIG_DIR))
    for f in files:
        size = os.path.getsize(os.path.join(FIG_DIR, f))
        print(f"  {f}  ({size // 1024} KB)")


if __name__ == "__main__":
    main()
