#!/usr/bin/env python3
"""Experiment 6b: Planck Wall — HIGH RESOLUTION

Same five QG models as planck_wall.py, but with ~80 temperature points
concentrated near the Planck wall (T/T_P = 1.0) to pinpoint exactly where
phase transitions, inflection points, and gradient changes happen.

Temperature grid:
  - 0.001 to 0.5:   coarse (background)
  - 0.5 to 2.0:     fine (0.02 T_P steps — the critical zone)
  - 2.0 to 100:     medium (post-wall behavior)
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


# ── Model Generators (same as planck_wall.py) ────────────────────────────

def standard_spectrum(T, k):
    E = k
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = k[valid]**2 / denom[valid]
    return spectrum


def lqg_spectrum(T, k):
    valid_k = k[k < np.pi]
    if len(valid_k) == 0:
        return np.array([])
    E = 2.0 * np.abs(np.sin(valid_k / 2.0))
    E = np.maximum(E, 1e-10)
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(valid_k)
    spectrum[valid] = valid_k[valid]**2 / denom[valid]
    return spectrum


def gup_spectrum(T, k, beta=1.0):
    E = k * np.sqrt(1.0 + beta * k**2)
    g = k**2 / (1.0 + beta * k**2)
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = g[valid] / denom[valid]
    return spectrum


def dsr_spectrum(T, k):
    E = 1.0 - np.exp(-k)
    E = np.maximum(E, 1e-10)
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = k[valid]**2 / denom[valid]
    return spectrum


def hagedorn_spectrum(T, k, T_H=1.0):
    E = k
    growth = np.exp(np.minimum(k / T_H, 500))
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = k[valid]**2 * growth[valid] / denom[valid]
    return spectrum


# ── Analysis helper ──────────────────────────────────────────────────────

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


# ── Gradient analysis ────────────────────────────────────────────────────

def find_transitions(temps, deltas):
    """Find inflection points, steepest gradients, and local extrema."""
    if len(temps) < 3:
        return {}

    t = np.array(temps)
    d = np.array(deltas)

    # Numerical gradient: dδ/dT
    grad = np.gradient(d, t)

    # Second derivative: d²δ/dT²
    grad2 = np.gradient(grad, t)

    # Steepest ascent/descent
    max_grad_idx = int(np.argmax(np.abs(grad)))

    # Local extrema (peaks and valleys)
    peaks = []
    valleys = []
    for i in range(1, len(d) - 1):
        if d[i] > d[i-1] and d[i] > d[i+1]:
            peaks.append({"T": float(t[i]), "delta_b": float(d[i])})
        if d[i] < d[i-1] and d[i] < d[i+1]:
            valleys.append({"T": float(t[i]), "delta_b": float(d[i])})

    # Sign changes in second derivative (inflection points)
    inflections = []
    for i in range(1, len(grad2)):
        if grad2[i-1] * grad2[i] < 0:
            # Linear interpolation of the zero crossing
            frac = abs(grad2[i-1]) / (abs(grad2[i-1]) + abs(grad2[i]))
            T_inflect = float(t[i-1] + frac * (t[i] - t[i-1]))
            d_inflect = float(np.interp(T_inflect, t, d))
            inflections.append({"T": round(T_inflect, 4), "delta_b": round(d_inflect, 6)})

    return {
        "steepest_gradient": {
            "T": float(t[max_grad_idx]),
            "delta_b": float(d[max_grad_idx]),
            "gradient": float(grad[max_grad_idx]),
        },
        "peaks": peaks,
        "valleys": valleys,
        "inflection_points": inflections,
        "gradient_at_wall": float(np.interp(1.0, t, grad)),
    }


# ── Main ─────────────────────────────────────────────────────────────────

def run_planck_wall_hires():
    """High-resolution sweep across the Planck wall."""

    # Build temperature grid — dense near T_P = 1.0
    coarse_low = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4]
    fine = list(np.arange(0.50, 2.02, 0.02))  # 0.50, 0.52, ..., 2.00
    coarse_high = [2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 50.0, 100.0]

    temperatures = sorted(set([round(t, 4) for t in coarse_low + fine + coarse_high]))

    k = np.linspace(0.001, 50, 100000)

    models = {
        "Standard":  lambda T: standard_spectrum(T, k),
        "LQG":       lambda T: lqg_spectrum(T, k),
        "GUP":       lambda T: gup_spectrum(T, k),
        "DSR":       lambda T: dsr_spectrum(T, k),
        "Hagedorn":  lambda T: hagedorn_spectrum(T, k),
    }

    results = {name: [] for name in models}

    print("=" * 80)
    print("EXPERIMENT 6b: Planck Wall — HIGH RESOLUTION")
    print(f"  {len(temperatures)} temperature points ({len(fine)} in the critical zone 0.50–2.00 T_P)")
    print("  All quantities in Planck units (E_P = T_P = 1)")
    print("=" * 80)
    print()

    for model_name, func in models.items():
        print(f"── {model_name.upper()} ", end="")
        sys.stdout.flush()

        for T in temperatures:
            spectrum = func(T)
            result = analyze_spectrum(spectrum)

            entry = {
                "T_planck_units": T,
                "T_kelvin": float(T * T_PLANCK_KELVIN),
                "T_GeV": float(T * E_PLANCK_GEV),
            }

            if result:
                entry.update(result)
                entry["computable"] = True
            else:
                entry["delta_b"] = None
                entry["computable"] = False
                n_pos = len([v for v in spectrum if v > 0 and np.isfinite(v)])
                entry["n_positive"] = n_pos

            results[model_name].append(entry)

        n_comp = sum(1 for e in results[model_name] if e["computable"])
        print(f"({n_comp}/{len(temperatures)} computable)")

    print()

    # ── Transition analysis ──────────────────────────────────────────────
    print("=" * 80)
    print("TRANSITION ANALYSIS")
    print("=" * 80)
    print()

    transitions = {}

    for model_name in models:
        entries = results[model_name]
        computable = [(e["T_planck_units"], e["delta_b"])
                      for e in entries if e["computable"]]

        if len(computable) < 3:
            transitions[model_name] = {"error": "too few computable points"}
            continue

        temps, deltas = zip(*computable)
        trans = find_transitions(list(temps), list(deltas))
        transitions[model_name] = trans

        print(f"── {model_name.upper()} ──")
        print(f"  Steepest gradient at T/T_P = {trans['steepest_gradient']['T']:.4f}")
        print(f"    δ_B = {trans['steepest_gradient']['delta_b']:.6f}, "
              f"dδ/dT = {trans['steepest_gradient']['gradient']:.6f}")
        print(f"  Gradient at the wall (T=1): dδ/dT = {trans['gradient_at_wall']:.6f}")

        if trans["peaks"]:
            print(f"  Peaks ({len(trans['peaks'])}):")
            for p in trans["peaks"]:
                print(f"    T/T_P = {p['T']:.4f}, δ_B = {p['delta_b']:.6f}")

        if trans["valleys"]:
            print(f"  Valleys ({len(trans['valleys'])}):")
            for v in trans["valleys"]:
                print(f"    T/T_P = {v['T']:.4f}, δ_B = {v['delta_b']:.6f}")

        if trans["inflection_points"]:
            print(f"  Inflection points ({len(trans['inflection_points'])}):")
            for ip in trans["inflection_points"][:5]:  # Show first 5
                print(f"    T/T_P = {ip['T']:.4f}, δ_B = {ip['delta_b']:.6f}")

        print()

    # ── High-res comparison table near the wall ──────────────────────────
    print("=" * 80)
    print("δ_B NEAR THE WALL (T/T_P = 0.80 to 1.20, step 0.04)")
    print("=" * 80)
    print()

    wall_temps = [t for t in temperatures if 0.78 <= t <= 1.22]
    # Thin to every-other for display readability
    display_temps = wall_temps[::2]

    print(f"  {'T/T_P':>7s}", end="")
    for name in models:
        print(f"  {name:>10s}", end="")
    print()
    print(f"  {'─'*7}", end="")
    for _ in models:
        print(f"  {'─'*10}", end="")
    print()

    for T in display_temps:
        line = f"  {T:7.2f}"
        is_wall = abs(T - 1.0) < 0.005
        for name in models:
            entry = next((e for e in results[name]
                          if abs(e["T_planck_units"] - T) < 0.005), None)
            if entry and entry["computable"]:
                line += f"  {entry['delta_b']:10.6f}"
            else:
                line += f"  {'UNDEF':>10s}"
        if is_wall:
            line += "  ◄ WALL"
        print(line)

    print()

    # ── Hagedorn spotlight ───────────────────────────────────────────────
    print("=" * 80)
    print("HAGEDORN SPOTLIGHT: The Phase Transition")
    print("=" * 80)
    print()

    hag_entries = [e for e in results["Hagedorn"] if e["computable"]]
    if hag_entries:
        # Find the peak (maximum chaos)
        peak_entry = max(hag_entries, key=lambda e: e["delta_b"])
        # Find pre-wall minimum
        pre_wall = [e for e in hag_entries if e["T_planck_units"] <= 0.5]
        post_wall = [e for e in hag_entries if e["T_planck_units"] >= 2.0]

        print(f"  Peak chaos: T/T_P = {peak_entry['T_planck_units']:.4f}, "
              f"δ_B = {peak_entry['delta_b']:.6f}")

        if pre_wall:
            pre_min = min(pre_wall, key=lambda e: e["delta_b"])
            print(f"  Pre-wall minimum: T/T_P = {pre_min['T_planck_units']:.4f}, "
                  f"δ_B = {pre_min['delta_b']:.6f}")

        if post_wall:
            post_min = min(post_wall, key=lambda e: e["delta_b"])
            print(f"  Post-wall minimum: T/T_P = {post_min['T_planck_units']:.4f}, "
                  f"δ_B = {post_min['delta_b']:.6f}")

        print()
        print("  Hagedorn δ_B through the wall:")
        hag_wall_zone = [e for e in hag_entries
                         if 0.5 <= e["T_planck_units"] <= 3.0]
        for e in hag_wall_zone:
            bar_len = int(e["delta_b"] * 200)
            bar = "█" * bar_len
            marker = ""
            if abs(e["T_planck_units"] - 1.0) < 0.005:
                marker = " ◄ WALL"
            elif abs(e["T_planck_units"] - peak_entry["T_planck_units"]) < 0.005:
                marker = " ◄ PEAK CHAOS"
            print(f"    T={e['T_planck_units']:5.2f}  δ_B={e['delta_b']:.6f}  {bar}{marker}")

    print()

    # ── Model quality ranking ────────────────────────────────────────────
    print("=" * 80)
    print("MODEL RANKING: Post-wall statistical quality")
    print("=" * 80)
    print()

    rankings = []
    for name in models:
        post = [e for e in results[name]
                if e["T_planck_units"] > 1.0 and e["computable"]]
        if post:
            mean_db = np.mean([e["delta_b"] for e in post])
            min_db = min(e["delta_b"] for e in post)
            max_db = max(e["delta_b"] for e in post)
            rankings.append((name, mean_db, min_db, max_db, len(post)))

    rankings.sort(key=lambda x: x[1])

    print(f"  {'Rank':>4s}  {'Model':>10s}  {'Mean δ_B':>10s}  "
          f"{'Min δ_B':>10s}  {'Max δ_B':>10s}  {'Points':>6s}")
    print(f"  {'─'*4}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*6}")

    for i, (name, mean_db, min_db, max_db, n) in enumerate(rankings):
        print(f"  {i+1:4d}  {name:>10s}  {mean_db:10.6f}  "
              f"{min_db:10.6f}  {max_db:10.6f}  {n:6d}")

    print()
    if rankings:
        winner = rankings[0]
        print(f"  WINNER: {winner[0]} (mean post-wall δ_B = {winner[1]:.6f})")
        print(f"  → Produces the most natural statistical structure beyond the Planck wall")

    print()

    # ── Save ─────────────────────────────────────────────────────────────
    output = {
        "experiment": "planck_wall_hires",
        "description": "High-resolution Planck wall sweep — pinpointing phase transitions",
        "n_temperature_points": len(temperatures),
        "critical_zone": "0.50 to 2.00 T_P in steps of 0.02",
        "models_tested": list(models.keys()),
        "temperatures_planck_units": temperatures,
        "results": {name: entries for name, entries in results.items()},
        "transitions": transitions,
        "rankings": [
            {"model": name, "mean_post_wall_delta_b": round(mean_db, 6),
             "min_post_wall_delta_b": round(min_db, 6),
             "max_post_wall_delta_b": round(max_db, 6)}
            for name, mean_db, min_db, max_db, _ in rankings
        ],
        "physical_scales": {
            "T_Planck_kelvin": T_PLANCK_KELVIN,
            "E_Planck_GeV": E_PLANCK_GEV,
        },
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = "results/round_trip/planck_wall_hires.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Results saved to {out_path}")

    return output


if __name__ == "__main__":
    run_planck_wall_hires()
