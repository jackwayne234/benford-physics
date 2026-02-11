#!/usr/bin/env python3
"""Robustness Testing: How stable are δ_B values under different computational parameters?

Varies:
  1. n_modes:   10k, 50k, 100k (baseline), 200k, 500k
  2. k_max:     20, 50 (baseline), 100, 200
  3. k_min:     0.0001, 0.001 (baseline), 0.01
  4. Grid type: linspace (baseline), logspace

Tests all 10 models at 8 key radial positions through the black hole.
Reports whether rankings and key values (CS equilibrium, horizon invisibility) are stable.
"""

import json
import math
import os
import sys
import time
import numpy as np
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benford_core import run_full_analysis


# ── Spectrum generators (from black_hole_wall.py) ────────────────────────

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
        "n": result["n"],
        "verdict": result["verdict"],
    }


def infalling_temperature(r_ratio, T_H):
    if r_ratio == 0:
        return None
    return T_H * (1.0 / abs(r_ratio)) ** 1.5


# ── Models dict builder ──────────────────────────────────────────────────

ALL_MODELS = {
    "Standard":      standard_spectrum,
    "LQG":           lqg_spectrum,
    "GUP":           gup_spectrum,
    "DSR":           dsr_spectrum,
    "Hagedorn":      hagedorn_spectrum,
    "Causal Set":    causal_set_spectrum,
    "Asym. Safety":  asymptotic_safety_spectrum,
    "Horava-Lif.":   horava_lifshitz_spectrum,
    "Noncommut.":    noncommutative_spectrum,
    "CDT":           cdt_spectrum,
}

# Key radial positions (representative subset)
KEY_POSITIONS = [10.0, 1.1, 1.001, 0.99, 0.5, 0.1, 0.04, 0.01]

T_H = 0.05  # Hawking temperature (Planck units)


def run_single_config(k_grid, models, positions, label):
    """Run all models at all positions for a given k grid. Returns dict of results."""
    results = {}
    for model_name, func in models.items():
        model_results = []
        for r in positions:
            T_eff = infalling_temperature(r, T_H)
            if T_eff is None or T_eff < 1e-6:
                model_results.append({"r_ratio": r, "delta_b": None})
                continue

            spectrum = func(T_eff, k_grid)
            result = analyze_spectrum(spectrum)

            if result:
                model_results.append({
                    "r_ratio": r,
                    "delta_b": result["delta_b"],
                    "n": result["n"],
                    "verdict": result["verdict"],
                })
            else:
                model_results.append({"r_ratio": r, "delta_b": None})

        results[model_name] = model_results
    return results


def compute_rankings(results, positions):
    """Compute mean interior δ_B and rank models."""
    rankings = []
    for model_name, entries in results.items():
        inside = [e["delta_b"] for e in entries
                  if e["r_ratio"] < 1.0 and e["delta_b"] is not None]
        outside = [e["delta_b"] for e in entries
                   if e["r_ratio"] > 1.0 and e["delta_b"] is not None]
        mean_in = np.mean(inside) if inside else float("inf")
        mean_out = np.mean(outside) if outside else float("inf")
        rankings.append({
            "model": model_name,
            "mean_inside": round(float(mean_in), 6),
            "mean_outside": round(float(mean_out), 6),
        })
    rankings.sort(key=lambda x: x["mean_inside"])
    return rankings


def get_cs_at_position(results, r_target):
    """Get Causal Set δ_B at a specific radial position."""
    for entry in results.get("Causal Set", []):
        if abs(entry["r_ratio"] - r_target) < 0.0005:
            return entry.get("delta_b")
    return None


def horizon_jump(results):
    """Compute |δ_B(1.001) - δ_B(0.99)| for Causal Set."""
    outside = get_cs_at_position(results, 1.001)
    inside = get_cs_at_position(results, 0.99)
    if outside is not None and inside is not None:
        return abs(outside - inside)
    return None


# ══════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════

def main():
    start_time = time.time()
    all_sweep_results = {}

    # ── SWEEP 1: Number of modes ─────────────────────────────────────────
    print("=" * 80)
    print("SWEEP 1: Varying number of momentum modes (k_min=0.001, k_max=50, linspace)")
    print("=" * 80)

    n_modes_values = [10000, 50000, 100000, 200000, 500000]
    sweep1 = {}

    for n_modes in n_modes_values:
        label = f"n_modes={n_modes}"
        print(f"\n  Running {label}...", end=" ", flush=True)
        k = np.linspace(0.001, 50, n_modes)
        t0 = time.time()
        results = run_single_config(k, ALL_MODELS, KEY_POSITIONS, label)
        rankings = compute_rankings(results, KEY_POSITIONS)
        dt = time.time() - t0
        print(f"done ({dt:.1f}s)")

        cs_eq = get_cs_at_position(results, 0.04)
        cs_horizon_outside = get_cs_at_position(results, 1.001)
        cs_horizon_inside = get_cs_at_position(results, 0.99)
        h_jump = horizon_jump(results)

        sweep1[n_modes] = {
            "rankings": rankings,
            "cs_equilibrium_0.04": cs_eq,
            "cs_horizon_outside": cs_horizon_outside,
            "cs_horizon_inside": cs_horizon_inside,
            "horizon_jump": h_jump,
            "first_place": rankings[0]["model"] if rankings else None,
            "all_positions": {
                model: {str(e["r_ratio"]): e.get("delta_b") for e in entries}
                for model, entries in results.items()
            },
        }

        # Print summary
        print(f"    Rank 1: {rankings[0]['model']} (mean inside δ_B = {rankings[0]['mean_inside']:.4f})")
        print(f"    Rank 2: {rankings[1]['model']} (mean inside δ_B = {rankings[1]['mean_inside']:.4f})")
        print(f"    CS equilibrium (r=0.04): {cs_eq}")
        print(f"    CS horizon jump: {h_jump}")

    all_sweep_results["n_modes"] = sweep1

    # ── SWEEP 2: k_max ───────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("SWEEP 2: Varying k_max (n_modes=100000, k_min=0.001, linspace)")
    print("=" * 80)

    k_max_values = [20, 50, 100, 200]
    sweep2 = {}

    for k_max in k_max_values:
        label = f"k_max={k_max}"
        print(f"\n  Running {label}...", end=" ", flush=True)
        k = np.linspace(0.001, k_max, 100000)
        t0 = time.time()
        results = run_single_config(k, ALL_MODELS, KEY_POSITIONS, label)
        rankings = compute_rankings(results, KEY_POSITIONS)
        dt = time.time() - t0
        print(f"done ({dt:.1f}s)")

        cs_eq = get_cs_at_position(results, 0.04)
        h_jump = horizon_jump(results)

        sweep2[k_max] = {
            "rankings": rankings,
            "cs_equilibrium_0.04": cs_eq,
            "horizon_jump": h_jump,
            "first_place": rankings[0]["model"] if rankings else None,
            "all_positions": {
                model: {str(e["r_ratio"]): e.get("delta_b") for e in entries}
                for model, entries in results.items()
            },
        }

        print(f"    Rank 1: {rankings[0]['model']} (mean inside δ_B = {rankings[0]['mean_inside']:.4f})")
        print(f"    Rank 2: {rankings[1]['model']} (mean inside δ_B = {rankings[1]['mean_inside']:.4f})")
        print(f"    CS equilibrium (r=0.04): {cs_eq}")
        print(f"    CS horizon jump: {h_jump}")

    all_sweep_results["k_max"] = sweep2

    # ── SWEEP 3: k_min ───────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("SWEEP 3: Varying k_min (n_modes=100000, k_max=50, linspace)")
    print("=" * 80)

    k_min_values = [0.0001, 0.001, 0.01, 0.1]
    sweep3 = {}

    for k_min in k_min_values:
        label = f"k_min={k_min}"
        print(f"\n  Running {label}...", end=" ", flush=True)
        k = np.linspace(k_min, 50, 100000)
        t0 = time.time()
        results = run_single_config(k, ALL_MODELS, KEY_POSITIONS, label)
        rankings = compute_rankings(results, KEY_POSITIONS)
        dt = time.time() - t0
        print(f"done ({dt:.1f}s)")

        cs_eq = get_cs_at_position(results, 0.04)
        h_jump = horizon_jump(results)

        sweep3[k_min] = {
            "rankings": rankings,
            "cs_equilibrium_0.04": cs_eq,
            "horizon_jump": h_jump,
            "first_place": rankings[0]["model"] if rankings else None,
            "all_positions": {
                model: {str(e["r_ratio"]): e.get("delta_b") for e in entries}
                for model, entries in results.items()
            },
        }

        print(f"    Rank 1: {rankings[0]['model']} (mean inside δ_B = {rankings[0]['mean_inside']:.4f})")
        print(f"    Rank 2: {rankings[1]['model']} (mean inside δ_B = {rankings[1]['mean_inside']:.4f})")
        print(f"    CS equilibrium (r=0.04): {cs_eq}")
        print(f"    CS horizon jump: {h_jump}")

    all_sweep_results["k_min"] = sweep3

    # ── SWEEP 4: Grid type ───────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("SWEEP 4: linspace vs logspace (n_modes=100000, k_range=[0.001, 50])")
    print("=" * 80)

    grids = {
        "linspace": np.linspace(0.001, 50, 100000),
        "logspace": np.logspace(np.log10(0.001), np.log10(50), 100000),
    }
    sweep4 = {}

    for grid_name, k in grids.items():
        label = f"grid={grid_name}"
        print(f"\n  Running {label}...", end=" ", flush=True)
        t0 = time.time()
        results = run_single_config(k, ALL_MODELS, KEY_POSITIONS, label)
        rankings = compute_rankings(results, KEY_POSITIONS)
        dt = time.time() - t0
        print(f"done ({dt:.1f}s)")

        cs_eq = get_cs_at_position(results, 0.04)
        h_jump = horizon_jump(results)

        sweep4[grid_name] = {
            "rankings": rankings,
            "cs_equilibrium_0.04": cs_eq,
            "horizon_jump": h_jump,
            "first_place": rankings[0]["model"] if rankings else None,
            "all_positions": {
                model: {str(e["r_ratio"]): e.get("delta_b") for e in entries}
                for model, entries in results.items()
            },
        }

        print(f"    Rank 1: {rankings[0]['model']} (mean inside δ_B = {rankings[0]['mean_inside']:.4f})")
        print(f"    Rank 2: {rankings[1]['model']} (mean inside δ_B = {rankings[1]['mean_inside']:.4f})")
        print(f"    CS equilibrium (r=0.04): {cs_eq}")
        print(f"    CS horizon jump: {h_jump}")

    all_sweep_results["grid_type"] = sweep4

    # ══════════════════════════════════════════════════════════════════════
    # SUMMARY
    # ══════════════════════════════════════════════════════════════════════

    total_time = time.time() - start_time

    print("\n" + "=" * 80)
    print("ROBUSTNESS SUMMARY")
    print("=" * 80)

    # Sweep 1 summary
    print("\n── Sweep 1: Number of Modes ──")
    print(f"  {'n_modes':>10s}  {'1st place':>14s}  {'CS inside':>10s}  {'CS equil':>10s}  {'H jump':>10s}")
    print(f"  {'─'*10}  {'─'*14}  {'─'*10}  {'─'*10}  {'─'*10}")
    for n, data in sweep1.items():
        cs_in = data["rankings"][0]["mean_inside"] if data["first_place"] == "Causal Set" else "—"
        cs_in_str = f"{cs_in:.4f}" if isinstance(cs_in, float) else cs_in
        # Find CS specifically
        cs_rank = next((r for r in data["rankings"] if r["model"] == "Causal Set"), None)
        cs_inside = f"{cs_rank['mean_inside']:.4f}" if cs_rank else "—"
        cs_eq = f"{data['cs_equilibrium_0.04']:.4f}" if data["cs_equilibrium_0.04"] else "—"
        hj = f"{data['horizon_jump']:.4f}" if data["horizon_jump"] else "—"
        print(f"  {n:>10d}  {data['first_place']:>14s}  {cs_inside:>10s}  {cs_eq:>10s}  {hj:>10s}")

    # Sweep 2 summary
    print("\n── Sweep 2: Momentum Range (k_max) ──")
    print(f"  {'k_max':>10s}  {'1st place':>14s}  {'CS inside':>10s}  {'CS equil':>10s}  {'H jump':>10s}")
    print(f"  {'─'*10}  {'─'*14}  {'─'*10}  {'─'*10}  {'─'*10}")
    for kmax, data in sweep2.items():
        cs_rank = next((r for r in data["rankings"] if r["model"] == "Causal Set"), None)
        cs_inside = f"{cs_rank['mean_inside']:.4f}" if cs_rank else "—"
        cs_eq = f"{data['cs_equilibrium_0.04']:.4f}" if data["cs_equilibrium_0.04"] else "—"
        hj = f"{data['horizon_jump']:.4f}" if data["horizon_jump"] else "—"
        print(f"  {kmax:>10d}  {data['first_place']:>14s}  {cs_inside:>10s}  {cs_eq:>10s}  {hj:>10s}")

    # Sweep 3 summary
    print("\n── Sweep 3: Momentum Range (k_min) ──")
    print(f"  {'k_min':>10s}  {'1st place':>14s}  {'CS inside':>10s}  {'CS equil':>10s}  {'H jump':>10s}")
    print(f"  {'─'*10}  {'─'*14}  {'─'*10}  {'─'*10}  {'─'*10}")
    for kmin, data in sweep3.items():
        cs_rank = next((r for r in data["rankings"] if r["model"] == "Causal Set"), None)
        cs_inside = f"{cs_rank['mean_inside']:.4f}" if cs_rank else "—"
        cs_eq = f"{data['cs_equilibrium_0.04']:.4f}" if data["cs_equilibrium_0.04"] else "—"
        hj = f"{data['horizon_jump']:.4f}" if data["horizon_jump"] else "—"
        print(f"  {kmin:>10.4f}  {data['first_place']:>14s}  {cs_inside:>10s}  {cs_eq:>10s}  {hj:>10s}")

    # Sweep 4 summary
    print("\n── Sweep 4: Grid Type ──")
    print(f"  {'grid':>10s}  {'1st place':>14s}  {'CS inside':>10s}  {'CS equil':>10s}  {'H jump':>10s}")
    print(f"  {'─'*10}  {'─'*14}  {'─'*10}  {'─'*10}  {'─'*10}")
    for gtype, data in sweep4.items():
        cs_rank = next((r for r in data["rankings"] if r["model"] == "Causal Set"), None)
        cs_inside = f"{cs_rank['mean_inside']:.4f}" if cs_rank else "—"
        cs_eq = f"{data['cs_equilibrium_0.04']:.4f}" if data["cs_equilibrium_0.04"] else "—"
        hj = f"{data['horizon_jump']:.4f}" if data["horizon_jump"] else "—"
        print(f"  {gtype:>10s}  {data['first_place']:>14s}  {cs_inside:>10s}  {cs_eq:>10s}  {hj:>10s}")

    # ── Overall verdict ──────────────────────────────────────────────────
    print("\n── Key Claims Stability ──")

    # Claim: CS ranks first
    all_firsts = []
    for sweep_name, sweep_data in all_sweep_results.items():
        for param, data in sweep_data.items():
            all_firsts.append(data["first_place"])
    cs_first_count = sum(1 for f in all_firsts if f == "Causal Set")
    print(f"  CS ranks 1st: {cs_first_count}/{len(all_firsts)} configurations ({100*cs_first_count/len(all_firsts):.0f}%)")

    # Claim: CS equilibrium ~0.015-0.017
    all_eq = []
    for sweep_name, sweep_data in all_sweep_results.items():
        for param, data in sweep_data.items():
            eq = data.get("cs_equilibrium_0.04")
            if eq is not None:
                all_eq.append(eq)
    if all_eq:
        print(f"  CS equilibrium (r=0.04): min={min(all_eq):.4f}, max={max(all_eq):.4f}, "
              f"mean={np.mean(all_eq):.4f}, std={np.std(all_eq):.4f}")

    # Claim: horizon invisible (jump < 0.002)
    all_jumps = []
    for sweep_name, sweep_data in all_sweep_results.items():
        for param, data in sweep_data.items():
            hj = data.get("horizon_jump")
            if hj is not None:
                all_jumps.append(hj)
    if all_jumps:
        print(f"  Horizon jump |δ_B(1.001)−δ_B(0.99)|: min={min(all_jumps):.4f}, "
              f"max={max(all_jumps):.4f}, mean={np.mean(all_jumps):.4f}")
        invisible_count = sum(1 for j in all_jumps if j < 0.005)
        print(f"  Horizon invisible (jump < 0.005): {invisible_count}/{len(all_jumps)} configurations")

    print(f"\n  Total runtime: {total_time:.1f}s")

    # ── Save full results ────────────────────────────────────────────────
    # Convert numpy types for JSON
    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, dict):
            return {str(k): convert(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [convert(v) for v in obj]
        return obj

    output = {
        "experiment": "robustness_test",
        "description": "Parameter sensitivity analysis for black hole wall experiment",
        "baseline_params": {
            "n_modes": 100000,
            "k_min": 0.001,
            "k_max": 50,
            "grid_type": "linspace",
            "T_H": 0.05,
        },
        "key_positions": KEY_POSITIONS,
        "sweeps": convert(all_sweep_results),
        "summary": {
            "cs_first_place_fraction": cs_first_count / len(all_firsts),
            "cs_equilibrium_range": [min(all_eq), max(all_eq)] if all_eq else None,
            "cs_equilibrium_mean": float(np.mean(all_eq)) if all_eq else None,
            "horizon_jump_range": [min(all_jumps), max(all_jumps)] if all_jumps else None,
            "horizon_invisible_fraction": invisible_count / len(all_jumps) if all_jumps else None,
        },
        "total_runtime_seconds": round(total_time, 1),
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = "results/round_trip/robustness_test.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n  Results saved to {out_path}")


if __name__ == "__main__":
    main()
