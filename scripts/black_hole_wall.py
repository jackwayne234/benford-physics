#!/usr/bin/env python3
"""Experiment 7: Black Hole Wall — 10 QG Models Through an Event Horizon

Same idea as the Planck Wall, but the "wall" is a black hole event horizon
instead of the Big Bang singularity.

Two perspectives modeled:

  A) STATIC OBSERVER (hovering outside):
     Local temperature = T_H / √(1 - r_s/r)  (Tolman redshift)
     Diverges at the horizon — this IS the wall.
     Only exists for r > r_s (can't hover inside).

  B) INFALLING OBSERVER (free-falling through):
     Effective tidal temperature ~ T_H × (r_s/r)^{3/2}  (curvature scale)
     Smooth at the horizon (equivalence principle — nothing special happens).
     Diverges at the singularity r → 0 — that's the REAL wall inside.

The experiment:
  - Sweep from far outside (r = 10 r_s) through the horizon down to
    near the singularity (r = 0.01 r_s)
  - At each radius, compute the effective temperature under both models
  - Run all 10 QG spectrum generators at that temperature
  - Track δ_B as a function of r/r_s
  - The horizon (r/r_s = 1) and singularity (r → 0) are the two walls

All quantities in Planck units.
"""

import json
import math
import os
import sys
import numpy as np
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benford_core import run_full_analysis


# ── Constants ────────────────────────────────────────────────────────────
T_PLANCK_KELVIN = 1.416784e32
E_PLANCK_GEV = 1.220910e19


# ── QG Model Generators (same as planck_wall_extended.py) ────────────────

def standard_spectrum(T, k):
    E = k
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    s = np.zeros_like(k)
    s[valid] = k[valid]**2 / denom[valid]
    return s

def lqg_spectrum(T, k):
    vk = k[k < np.pi]
    if len(vk) == 0:
        return np.array([])
    E = 2.0 * np.abs(np.sin(vk / 2.0))
    E = np.maximum(E, 1e-10)
    exp = np.minimum(E / T, 500)
    d = np.exp(exp) - 1.0
    v = d > 0
    s = np.zeros_like(vk)
    s[v] = vk[v]**2 / d[v]
    return s

def gup_spectrum(T, k):
    E = k * np.sqrt(1.0 + k**2)
    g = k**2 / (1.0 + k**2)
    exp = np.minimum(E / T, 500)
    d = np.exp(exp) - 1.0
    v = d > 0
    s = np.zeros_like(k)
    s[v] = g[v] / d[v]
    return s

def dsr_spectrum(T, k):
    E = 1.0 - np.exp(-k)
    E = np.maximum(E, 1e-10)
    exp = np.minimum(E / T, 500)
    d = np.exp(exp) - 1.0
    v = d > 0
    s = np.zeros_like(k)
    s[v] = k[v]**2 / d[v]
    return s

def hagedorn_spectrum(T, k):
    E = k
    growth = np.exp(np.minimum(k, 500))
    exp = np.minimum(E / T, 500)
    d = np.exp(exp) - 1.0
    v = d > 0
    s = np.zeros_like(k)
    s[v] = k[v]**2 * growth[v] / d[v]
    return s

def causal_set_spectrum(T, k):
    E = k
    suppression = np.exp(-k**2)
    g = k**2 * suppression
    exp = np.minimum(E / T, 500)
    d = np.exp(exp) - 1.0
    v = d > 0
    s = np.zeros_like(k)
    s[v] = g[v] / d[v]
    return s

def asymptotic_safety_spectrum(T, k):
    E = k
    d_s = 2.0 + 2.0 / (1.0 + k**2)
    g = k**(d_s - 1.0)
    exp = np.minimum(E / T, 500)
    d = np.exp(exp) - 1.0
    v = d > 0
    s = np.zeros_like(k)
    s[v] = g[v] / d[v]
    return s

def horava_lifshitz_spectrum(T, k):
    E = np.sqrt(k**2 + k**4 + k**6)
    g = k**2
    exp = np.minimum(E / T, 500)
    d = np.exp(exp) - 1.0
    v = d > 0
    s = np.zeros_like(k)
    s[v] = g[v] / d[v]
    return s

def noncommutative_spectrum(T, k):
    E = np.sqrt(k**2 + k**4)
    g = k**2 * (1.0 + k**2)
    exp = np.minimum(E / T, 500)
    d = np.exp(exp) - 1.0
    v = d > 0
    s = np.zeros_like(k)
    s[v] = g[v] / d[v]
    return s

def cdt_spectrum(T, k):
    E = k
    d_s = 2.0 + 2.0 / (1.0 + k**4)
    g = k**(d_s - 1.0)
    exp = np.minimum(E / T, 500)
    d = np.exp(exp) - 1.0
    v = d > 0
    s = np.zeros_like(k)
    s[v] = g[v] / d[v]
    return s


# ── Analysis ─────────────────────────────────────────────────────────────

def analyze_spectrum(spectrum_values):
    values = [float(v) for v in spectrum_values if v > 0 and np.isfinite(v)]
    if len(values) < 50:
        return None
    result = run_full_analysis(values)
    fd = result["first_digit"]
    return {
        "delta_b": fd["euclidean_deviation"],
        "epsilon_d": {str(k): v for k, v in fd["per_digit_deviation"].items()},
        "mad": fd["mad"]["mad"],
        "n": result["n"],
        "verdict": result["verdict"],
    }


# ── Temperature models ───────────────────────────────────────────────────

def static_temperature(r_ratio, T_H):
    """Tolman-redshifted temperature for a static (hovering) observer.
    Only valid outside the horizon (r_ratio > 1).
    Diverges as r_ratio → 1."""
    if r_ratio <= 1.0:
        return None  # Can't hover inside
    return T_H / math.sqrt(1.0 - 1.0 / r_ratio)


def infalling_temperature(r_ratio, T_H):
    """Effective tidal temperature for an infalling observer.
    Based on the curvature scale: T ~ (r_s/r)^{3/2}.
    Smooth at the horizon, diverges at the singularity (r → 0).

    Normalized so that at r = r_s: T_eff = T_H (matches Hawking temp).
    At r → 0: T_eff → ∞ (singularity).
    At r → ∞: T_eff → 0 (flat space).

    For r_ratio < 0 (post-singularity / bounce):
    The "bounce" models (LQG, string theory) predict the singularity
    resolves and spacetime continues. Temperature mirrors the approach:
    T_eff = T_H × (r_s/|r|)^{3/2}, decreasing as you move away from
    the singularity on the other side.
    """
    if r_ratio == 0:
        return None
    # Use |r| for both sides of the singularity (bounce symmetry)
    return T_H * (1.0 / abs(r_ratio)) ** 1.5


# ── Main ─────────────────────────────────────────────────────────────────

def run_black_hole_wall():

    # Black hole parameters
    # T_H = 0.05 T_P gives a nice range:
    #   - Static observer hits T_P near the horizon
    #   - Infalling observer hits T_P at r ≈ 0.14 r_s
    T_H = 0.05  # Hawking temperature in Planck units

    # Radial grid: outside → horizon → singularity
    # Outside (r > r_s): from far to very close to horizon
    outside = [10.0, 7.0, 5.0, 3.0, 2.0, 1.5, 1.3, 1.2, 1.15, 1.10,
               1.08, 1.06, 1.04, 1.03, 1.02, 1.015, 1.01, 1.005, 1.002, 1.001]
    # Inside (r < r_s): from just inside to near singularity
    inside = [0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40,
              0.30, 0.25, 0.20, 0.15, 0.12, 0.10, 0.08, 0.06, 0.04, 0.02, 0.01]
    # Post-singularity (bounce): mirror of inside, representing "the other side"
    # Negative r values = past the singularity
    post_singularity = [-0.01, -0.02, -0.04, -0.06, -0.08, -0.10, -0.12, -0.15,
                        -0.20, -0.25, -0.30, -0.40, -0.50, -0.70, -1.0]

    all_r = sorted(outside + inside + post_singularity, reverse=True)  # Far → singularity → past

    k = np.linspace(0.001, 50, 100000)

    models = {
        "Standard":      lambda T: standard_spectrum(T, k),
        "LQG":           lambda T: lqg_spectrum(T, k),
        "GUP":           lambda T: gup_spectrum(T, k),
        "DSR":           lambda T: dsr_spectrum(T, k),
        "Hagedorn":      lambda T: hagedorn_spectrum(T, k),
        "Causal Set":    lambda T: causal_set_spectrum(T, k),
        "Asym. Safety":  lambda T: asymptotic_safety_spectrum(T, k),
        "Horava-Lif.":   lambda T: horava_lifshitz_spectrum(T, k),
        "Noncommut.":    lambda T: noncommutative_spectrum(T, k),
        "CDT":           lambda T: cdt_spectrum(T, k),
    }

    print("=" * 84)
    print("EXPERIMENT 7: Black Hole Wall — 10 QG Models Through an Event Horizon")
    print(f"  Hawking temperature: T_H = {T_H} T_P")
    print(f"  {len(all_r)} radial points from r = {max(all_r)} r_s to r = {min(all_r)} r_s")
    print("  Two perspectives: static observer (outside only) + infalling observer (all r)")
    print("=" * 84)
    print()

    # ── Compute temperature profiles ─────────────────────────────────────
    print("── TEMPERATURE PROFILE ──")
    print(f"  {'r/r_s':>8s}  {'Zone':>8s}  {'T_static':>12s}  {'T_infall':>12s}")
    print(f"  {'─'*8}  {'─'*8}  {'─'*12}  {'─'*12}")

    for r in all_r:
        zone = "outside" if r > 1.0 else "inside"
        T_s = static_temperature(r, T_H)
        T_i = infalling_temperature(r, T_H)
        T_s_str = f"{T_s:.6f}" if T_s is not None else "N/A"
        T_i_str = f"{T_i:.6f}" if T_i is not None else "N/A"
        marker = ""
        if T_s is not None and abs(T_s - 1.0) < 0.05:
            marker = " ← T = T_Planck (static)"
        if T_i is not None and abs(T_i - 1.0) < 0.05:
            marker = " ← T = T_Planck (infall)"
        if abs(r - 1.0) < 0.005:
            marker = " ← HORIZON"
        print(f"  {r:8.3f}  {zone:>8s}  {T_s_str:>12s}  {T_i_str:>12s}{marker}")

    print()

    # ── Run infalling observer perspective (all radii) ───────────────────
    print("=" * 84)
    print("INFALLING OBSERVER — Through the horizon to the singularity")
    print("  T_eff(r) = T_H × (r_s/r)^{3/2}")
    print("=" * 84)
    print()

    infall_results = {name: [] for name in models}

    for model_name, func in models.items():
        print(f"  {model_name:<16s}", end="", flush=True)

        for r in all_r:
            T_eff = infalling_temperature(r, T_H)
            if T_eff is None or T_eff < 1e-6:
                entry = {"r_ratio": r, "T_eff": None, "computable": False,
                         "zone": "outside" if r > 1.0 else "inside"}
                infall_results[model_name].append(entry)
                continue

            spectrum = func(T_eff)
            result = analyze_spectrum(spectrum)

            if r > 1.0:
                zone = "outside"
            elif r > 0:
                zone = "inside"
            else:
                zone = "post-singularity"

            entry = {
                "r_ratio": r,
                "T_eff": round(float(T_eff), 8),
                "zone": zone,
            }

            if result:
                entry.update(result)
                entry["computable"] = True
            else:
                entry["delta_b"] = None
                entry["computable"] = False
                n_pos = len([v for v in spectrum if v > 0 and np.isfinite(v)])
                entry["n_positive"] = n_pos

            infall_results[model_name].append(entry)

        n_comp = sum(1 for e in infall_results[model_name] if e["computable"])
        n_total = len(all_r)
        print(f"{n_comp}/{n_total} computable")

    print()

    # ── Run static observer perspective (outside only) ───────────────────
    print("=" * 84)
    print("STATIC OBSERVER — Hovering outside, approaching the horizon")
    print("  T_local(r) = T_H / √(1 - r_s/r)")
    print("=" * 84)
    print()

    static_results = {name: [] for name in models}
    outside_r = [r for r in all_r if r > 1.0]

    for model_name, func in models.items():
        print(f"  {model_name:<16s}", end="", flush=True)

        for r in outside_r:
            T_local = static_temperature(r, T_H)

            spectrum = func(T_local)
            result = analyze_spectrum(spectrum)

            entry = {
                "r_ratio": r,
                "T_local": round(float(T_local), 8),
            }

            if result:
                entry.update(result)
                entry["computable"] = True
            else:
                entry["delta_b"] = None
                entry["computable"] = False
                n_pos = len([v for v in spectrum if v > 0 and np.isfinite(v)])
                entry["n_positive"] = n_pos

            static_results[model_name].append(entry)

        n_comp = sum(1 for e in static_results[model_name] if e["computable"])
        print(f"{n_comp}/{len(outside_r)} computable")

    print()

    # ══════════════════════════════════════════════════════════════════════
    # INFALLING OBSERVER: δ_B vs r/r_s
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 84)
    print("INFALLING OBSERVER: δ_B approaching the singularity")
    print("=" * 84)
    print()

    # Print in two blocks of 5
    for block_name, block_models in [("ORIGINAL 5", list(models.keys())[:5]),
                                      ("NEW 5", list(models.keys())[5:])]:
        print(f"── {block_name} (infalling) ──")
        print(f"  {'r/r_s':>8s}  {'zone':>8s}  {'T_eff':>8s}", end="")
        for name in block_models:
            print(f"  {name:>12s}", end="")
        print()
        print(f"  {'─'*8}  {'─'*8}  {'─'*8}", end="")
        for _ in block_models:
            print(f"  {'─'*12}", end="")
        print()

        display_r = [10.0, 5.0, 2.0, 1.5, 1.1, 1.01, 0.99, 0.90, 0.70,
                     0.50, 0.30, 0.15, 0.10, 0.04, 0.01,
                     -0.01, -0.04, -0.10, -0.15, -0.30, -0.50, -1.0]

        for r in display_r:
            T_eff = infalling_temperature(r, T_H)
            zone = "outside" if r > 1.0 else "inside"
            line = f"  {r:8.3f}  {zone:>8s}  {T_eff:8.4f}"

            for name in block_models:
                entry = next((e for e in infall_results[name]
                              if abs(e["r_ratio"] - r) < 0.0005), None)
                if entry and entry["computable"]:
                    line += f"  {entry['delta_b']:12.6f}"
                else:
                    line += f"  {'UNDEF':>12s}"

            if abs(r - 1.01) < 0.005 or abs(r - 0.99) < 0.005:
                line += "  ← HORIZON"
            if abs(r - 0.01) < 0.005:
                line += "  ← SINGULARITY"
            if abs(r - (-0.01)) < 0.005:
                line += "  ← POST-SINGULARITY"
            print(line)
        print()

    # ══════════════════════════════════════════════════════════════════════
    # STATIC OBSERVER: δ_B approaching the horizon
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 84)
    print("STATIC OBSERVER: δ_B hovering near the event horizon")
    print("=" * 84)
    print()

    for block_name, block_models in [("ORIGINAL 5", list(models.keys())[:5]),
                                      ("NEW 5", list(models.keys())[5:])]:
        print(f"── {block_name} (static) ──")
        print(f"  {'r/r_s':>8s}  {'T_local':>10s}", end="")
        for name in block_models:
            print(f"  {name:>12s}", end="")
        print()
        print(f"  {'─'*8}  {'─'*10}", end="")
        for _ in block_models:
            print(f"  {'─'*12}", end="")
        print()

        static_display = [10.0, 5.0, 2.0, 1.5, 1.2, 1.1, 1.06, 1.04,
                          1.02, 1.01, 1.005, 1.002, 1.001]

        for r in static_display:
            T_local = static_temperature(r, T_H)
            line = f"  {r:8.3f}  {T_local:10.4f}"

            for name in block_models:
                entry = next((e for e in static_results[name]
                              if abs(e["r_ratio"] - r) < 0.0005), None)
                if entry and entry["computable"]:
                    line += f"  {entry['delta_b']:12.6f}"
                else:
                    line += f"  {'UNDEF':>12s}"

            # Mark where T_local crosses T_P
            if abs(T_local - 1.0) < 0.1:
                line += f"  ← T ≈ T_P"
            print(line)
        print()

    # ══════════════════════════════════════════════════════════════════════
    # RANKINGS
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 84)
    print("RANKING: Post-horizon quality (infalling observer, r < r_s)")
    print("=" * 84)
    print()

    rankings = []
    for name in models:
        inside_entries = [e for e in infall_results[name]
                         if 0 < e["r_ratio"] < 1.0 and e["computable"]]
        outside_entries = [e for e in infall_results[name]
                          if e["r_ratio"] > 1.0 and e["computable"]]
        post_entries = [e for e in infall_results[name]
                        if e["r_ratio"] < 0 and e["computable"]]

        if inside_entries:
            mean_in = np.mean([e["delta_b"] for e in inside_entries])
            min_in = min(e["delta_b"] for e in inside_entries)
            max_in = max(e["delta_b"] for e in inside_entries)
        else:
            mean_in = float("inf")
            min_in = float("inf")
            max_in = float("inf")

        if outside_entries:
            mean_out = np.mean([e["delta_b"] for e in outside_entries])
        else:
            mean_out = float("inf")

        n_undef = sum(1 for e in infall_results[name]
                      if e["r_ratio"] < 1.0 and not e["computable"])

        rankings.append({
            "name": name,
            "mean_inside": mean_in, "min_inside": min_in, "max_inside": max_in,
            "mean_outside": mean_out, "n_undefined_inside": n_undef,
        })

    rankings.sort(key=lambda x: x["mean_inside"])

    print(f"  {'Rank':>4s}  {'Model':>16s}  {'Inside δ_B':>12s}  "
          f"{'Outside δ_B':>12s}  {'Undef':>5s}  {'Character':>20s}")
    print(f"  {'─'*4}  {'─'*16}  {'─'*12}  {'─'*12}  {'─'*5}  {'─'*20}")

    for i, r in enumerate(rankings):
        if r["mean_inside"] == float("inf"):
            char = "BREAKS inside"
        elif r["mean_inside"] < r["mean_outside"] * 0.5:
            char = "IMPROVES inside"
        elif r["mean_inside"] > r["mean_outside"] * 2.0:
            char = "DEGRADES inside"
        elif abs(r["mean_inside"] - r["mean_outside"]) / max(r["mean_outside"], 0.001) < 0.3:
            char = "FLAT"
        else:
            char = "MIXED"

        inside_str = f"{r['mean_inside']:12.6f}" if r["mean_inside"] != float("inf") else f"{'NO DATA':>12s}"
        outside_str = f"{r['mean_outside']:12.6f}" if r["mean_outside"] != float("inf") else f"{'NO DATA':>12s}"

        print(f"  {i+1:4d}  {r['name']:>16s}  {inside_str}  "
              f"{outside_str}  {r['n_undefined_inside']:5d}  {char:>20s}")

    print()

    # ══════════════════════════════════════════════════════════════════════
    # COMPARISON: Big Bang Wall vs Black Hole Wall
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 84)
    print("BIG BANG vs BLACK HOLE: Do models behave the same at both walls?")
    print("=" * 84)
    print()

    # Load Big Bang results for comparison
    bb_path = "results/round_trip/planck_wall_extended.json"
    if os.path.exists(bb_path):
        with open(bb_path) as f:
            bb_data = json.load(f)

        print(f"  {'Model':>16s}  {'BB post-wall':>14s}  {'BH inside':>14s}  "
              f"{'Ratio':>8s}  {'Same behavior?':>16s}")
        print(f"  {'─'*16}  {'─'*14}  {'─'*14}  {'─'*8}  {'─'*16}")

        for r in rankings:
            name = r["name"]
            bb_entries = bb_data["results"].get(name, [])
            bb_post = [e for e in bb_entries
                       if e.get("T_planck_units", 0) > 1.0 and e.get("computable", False)]

            if bb_post and r["mean_inside"] != float("inf"):
                bb_mean = np.mean([e["delta_b"] for e in bb_post])
                ratio = r["mean_inside"] / bb_mean if bb_mean > 0 else float("inf")
                same = "YES" if 0.5 < ratio < 2.0 else "NO"
                print(f"  {name:>16s}  {bb_mean:14.6f}  {r['mean_inside']:14.6f}  "
                      f"{ratio:8.2f}  {same:>16s}")
            else:
                print(f"  {name:>16s}  {'N/A':>14s}  {'N/A':>14s}  "
                      f"{'---':>8s}  {'N/A':>16s}")

        print()
        print("  Ratio ≈ 1.0 means the model behaves identically at both walls.")
        print("  Ratio >> 1 means worse inside a black hole than after the Big Bang.")
        print("  Ratio << 1 means better inside a black hole than after the Big Bang.")
    else:
        print("  (Big Bang results not found — run planck_wall_extended.py first)")

    print()

    # ══════════════════════════════════════════════════════════════════════
    # SAVE
    # ══════════════════════════════════════════════════════════════════════

    output = {
        "experiment": "black_hole_wall",
        "description": "10 QG models through a black hole event horizon to the singularity",
        "black_hole": {
            "T_hawking_planck": T_H,
            "T_hawking_kelvin": float(T_H * T_PLANCK_KELVIN),
        },
        "radial_points": all_r,
        "models_tested": list(models.keys()),
        "infalling_observer": {name: entries for name, entries in infall_results.items()},
        "static_observer": {name: entries for name, entries in static_results.items()},
        "rankings": [
            {"model": r["name"],
             "mean_inside_delta_b": round(r["mean_inside"], 6) if r["mean_inside"] != float("inf") else None,
             "mean_outside_delta_b": round(r["mean_outside"], 6) if r["mean_outside"] != float("inf") else None,
             "n_undefined_inside": r["n_undefined_inside"]}
            for r in rankings
        ],
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = "results/round_trip/black_hole_wall.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Results saved to {out_path}")

    return output


if __name__ == "__main__":
    run_black_hole_wall()
