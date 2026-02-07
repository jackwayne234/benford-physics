#!/usr/bin/env python3
"""Experiment 6: The Planck Wall — Which quantum gravity proposals survive?

General relativity breaks down at the Planck scale (T ~ T_P ≈ 1.4 × 10³² K).
Different quantum gravity proposals handle this differently:

  1. Standard (GR + QFT): no modification to dispersion, no UV cutoff
  2. Loop Quantum Gravity: polymer dispersion E = (2/l_P)|sin(l_P k/2)|,
     maximum energy, discrete spectrum
  3. GUP (Generalized Uncertainty Principle): modified dispersion
     E = k√(1 + β k²), minimum length, steeper UV behavior
  4. DSR (Doubly Special Relativity): energy saturates at E_P,
     E = E_P(1 - e^{-k/E_P})
  5. String/Hagedorn: exponential growth of states, maximum temperature T_H

For each model, we generate the thermal distribution at temperatures from
T << T_P to T >> T_P and run the Benford analysis.  Models that maintain
computable δ_B through the Planck wall describe physics on the other side.
Models that go UNDEFINED hit a genuine singularity.

All quantities in Planck units: E_P = T_P = k_P = l_P = 1.
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

# Planck units: E_P = T_P = l_P = 1
# Physical values (for display only):
T_PLANCK_KELVIN = 1.416784e32   # K
E_PLANCK_GEV = 1.220910e19     # GeV


# ── Model Generators ─────────────────────────────────────────────────────

def standard_spectrum(T, k):
    """Standard GR+QFT: E = k, density of states g(k) = k².
    Spectrum: S(k) = k² / (e^{k/T} - 1).
    No modification — should always be computable.
    """
    E = k
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    # Avoid division by zero
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = k[valid]**2 / denom[valid]
    return spectrum


def lqg_spectrum(T, k):
    """Loop Quantum Gravity: polymer dispersion.
    E(k) = 2|sin(k/2)| (in Planck units with polymer scale = 1).
    Maximum energy E_max = 2 at k = π.
    Modes only exist for 0 < k < π (first Brillouin zone).
    """
    # Truncate to first Brillouin zone
    valid_k = k[k < np.pi]
    if len(valid_k) == 0:
        return np.array([])
    E = 2.0 * np.abs(np.sin(valid_k / 2.0))
    # Avoid E = 0 at k = 0
    E = np.maximum(E, 1e-10)
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(valid_k)
    # Use same density of states k² within the zone
    spectrum[valid] = valid_k[valid]**2 / denom[valid]
    return spectrum


def gup_spectrum(T, k, beta=1.0):
    """GUP (Generalized Uncertainty Principle): modified dispersion.
    E = k · √(1 + β k²).  At low k: E ≈ k.  At high k: E ≈ √β k².
    Modified density of states: g(k) = k² / (1 + β k²).
    """
    E = k * np.sqrt(1.0 + beta * k**2)
    g = k**2 / (1.0 + beta * k**2)
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = g[valid] / denom[valid]
    return spectrum


def dsr_spectrum(T, k):
    """DSR (Doubly Special Relativity): energy saturates at E_P = 1.
    E = 1 - e^{-k} (approaches 1 asymptotically).
    At low k: E ≈ k.  As k → ∞: E → 1.
    """
    E = 1.0 - np.exp(-k)
    # Avoid E = 0 exactly
    E = np.maximum(E, 1e-10)
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = k[valid]**2 / denom[valid]
    return spectrum


def hagedorn_spectrum(T, k, T_H=1.0):
    """String theory / Hagedorn: exponential growth of states.
    Density of states: g(E) ~ E^{-5/2} e^{E/T_H} (for bosonic string, D=26).
    In k-space with E = k: g(k) = k² · e^{k/T_H}.
    Spectrum: S(k) = k² · e^{k/T_H} / (e^{k/T} - 1).

    At T < T_H: converges (Boltzmann suppression wins).
    At T = T_H: diverges (exponentials cancel).
    At T > T_H: violently diverges → UNDEFINED.
    """
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
    """Run Benford analysis on spectrum values. Returns dict or None."""
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


# ── Main experiment ──────────────────────────────────────────────────────

def run_planck_wall():
    """Sweep temperature from sub-Planck to super-Planck for all models."""

    # Temperature grid (in Planck units)
    temperatures = [
        0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5,
        0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1,
        1.3, 1.5, 2.0, 3.0, 5.0, 10.0, 50.0, 100.0
    ]

    # Momentum grid (in Planck units)
    k = np.linspace(0.001, 50, 100000)

    models = {
        "Standard":  {"func": lambda T: standard_spectrum(T, k), "color": "blue"},
        "LQG":       {"func": lambda T: lqg_spectrum(T, k), "color": "green"},
        "GUP":       {"func": lambda T: gup_spectrum(T, k), "color": "orange"},
        "DSR":       {"func": lambda T: dsr_spectrum(T, k), "color": "red"},
        "Hagedorn":  {"func": lambda T: hagedorn_spectrum(T, k), "color": "purple"},
    }

    results = {name: [] for name in models}

    print("=" * 78)
    print("EXPERIMENT 6: The Planck Wall")
    print("  Pushing thermal distributions through T_Planck under five QG proposals")
    print("  All quantities in Planck units (E_P = T_P = 1)")
    print("=" * 78)
    print()

    # Run sweep for each model
    for model_name, model_info in models.items():
        print(f"─── {model_name.upper()} ───")
        func = model_info["func"]

        for T in temperatures:
            spectrum = func(T)
            result = analyze_spectrum(spectrum)

            T_physical = T * T_PLANCK_KELVIN
            T_gev = T * E_PLANCK_GEV

            entry = {
                "T_planck_units": T,
                "T_kelvin": float(T_physical),
                "T_GeV": float(T_gev),
            }

            if result:
                entry.update(result)
                entry["computable"] = True
                status = f"δ_B={result['delta_b']:.6f}  n={result['n']:6d}  {result['verdict']}"
            else:
                entry["delta_b"] = None
                entry["computable"] = False
                n_pos = len([v for v in spectrum if v > 0 and np.isfinite(v)])
                entry["n_positive"] = n_pos
                status = f"UNDEFINED ({n_pos} valid modes)"

            results[model_name].append(entry)

            # Print key temperatures
            if T in [0.001, 0.1, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0, 5.0, 10.0, 100.0]:
                label = ""
                if abs(T - 1.0) < 0.001:
                    label = " ← PLANCK WALL"
                print(f"  T/T_P = {T:7.3f}  {status}{label}")

        print()

    # ══════════════════════════════════════════════════════════════════════
    # THE WALL
    # ══════════════════════════════════════════════════════════════════════

    print()
    print("█" * 78)
    print("█" + " " * 76 + "█")
    print("█" + "         T H E   P L A N C K   W A L L".center(76) + "█")
    print("█" + "    Which models survive beyond T_Planck?".center(76) + "█")
    print("█" + " " * 76 + "█")
    print("█" * 78)
    print()

    # Comparison table at key temperatures
    key_temps = [0.001, 0.1, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0, 5.0, 10.0, 100.0]

    print(f"  {'T/T_P':>7s}", end="")
    for name in models:
        print(f"  {name:>12s}", end="")
    print()
    print(f"  {'─'*7}", end="")
    for _ in models:
        print(f"  {'─'*12}", end="")
    print()

    for T in key_temps:
        line = f"  {T:7.3f}"
        for name in models:
            entry = next((e for e in results[name] if abs(e["T_planck_units"] - T) < 0.0001), None)
            if entry and entry["computable"]:
                line += f"  {entry['delta_b']:12.6f}"
            else:
                line += f"  {'UNDEFINED':>12s}"
        if abs(T - 1.0) < 0.001:
            line += "  ← PLANCK WALL"
        print(line)

    print()

    # ── Survival analysis ────────────────────────────────────────────────
    print("─" * 78)
    print("SURVIVAL ANALYSIS")
    print("─" * 78)
    print()

    for name in models:
        entries = results[name]
        computable = [e for e in entries if e["computable"]]
        undefined = [e for e in entries if not e["computable"]]

        # Find max computable temperature
        max_T = max(e["T_planck_units"] for e in computable) if computable else 0

        # Find where it breaks (first UNDEFINED)
        break_T = None
        for e in entries:
            if not e["computable"]:
                break_T = e["T_planck_units"]
                break

        print(f"  {name}:")
        print(f"    Computable points: {len(computable)}/{len(entries)}")
        print(f"    Max computable T/T_P: {max_T:.3f}", end="")
        if max_T >= 100:
            print("  (survives to 100 T_P and beyond)")
        elif max_T >= 1.0:
            print(f"  (survives through the wall)")
        else:
            print(f"  (breaks BEFORE the wall)")

        if break_T is not None:
            print(f"    First breakdown at T/T_P = {break_T:.3f}", end="")
            physical_T = break_T * T_PLANCK_KELVIN
            physical_E = break_T * E_PLANCK_GEV
            print(f"  ({physical_E:.2e} GeV)")
        else:
            print(f"    No breakdown detected in tested range")

        # δ_B trajectory
        if computable:
            pre_wall = [e for e in computable if e["T_planck_units"] <= 1.0]
            post_wall = [e for e in computable if e["T_planck_units"] > 1.0]
            if pre_wall:
                pre_dbs = [e["delta_b"] for e in pre_wall]
                print(f"    Pre-wall δ_B range:  [{min(pre_dbs):.6f}, {max(pre_dbs):.6f}]")
            if post_wall:
                post_dbs = [e["delta_b"] for e in post_wall]
                print(f"    Post-wall δ_B range: [{min(post_dbs):.6f}, {max(post_dbs):.6f}]")

        print()

    # ── Phase diagram ────────────────────────────────────────────────────
    print("─" * 78)
    print("PHASE DIAGRAM: Computability across the Planck Wall")
    print("─" * 78)
    print()

    for name in models:
        entries = results[name]
        bar = ""
        for e in entries:
            if e["computable"]:
                if e["T_planck_units"] < 0.9:
                    bar += "░"  # pre-wall, computable
                elif e["T_planck_units"] <= 1.1:
                    bar += "▓"  # at the wall
                else:
                    bar += "█"  # post-wall, computable
            else:
                bar += "✗"
        print(f"  {name:<12s} [{bar}]")

    print()
    print("  Legend: ░ = pre-wall  ▓ = at wall  █ = post-wall  ✗ = UNDEFINED")
    print(f"  Temperature: {temperatures[0]} → {temperatures[-1]} T_P ({len(temperatures)} points)")
    print()

    # ── Verdict ──────────────────────────────────────────────────────────
    print("─" * 78)
    print("VERDICT")
    print("─" * 78)
    print()

    survivors = []
    breakers = []
    for name in models:
        entries = results[name]
        post_wall = [e for e in entries if e["T_planck_units"] > 1.0 and e["computable"]]
        if post_wall:
            survivors.append(name)
        else:
            breakers.append(name)

    if survivors:
        print(f"  SURVIVE the Planck Wall ({len(survivors)} models):")
        for name in survivors:
            entries = results[name]
            max_T = max(e["T_planck_units"] for e in entries if e["computable"])
            post_wall = [e for e in entries if e["T_planck_units"] > 1.0 and e["computable"]]
            mean_db = np.mean([e["delta_b"] for e in post_wall])
            print(f"    ✓ {name}: computable to T = {max_T:.0f} T_P, "
                  f"mean post-wall δ_B = {mean_db:.6f}")
        print()

    if breakers:
        print(f"  BREAK at the Planck Wall ({len(breakers)} models):")
        for name in breakers:
            entries = results[name]
            break_T = None
            for e in entries:
                if not e["computable"]:
                    break_T = e["T_planck_units"]
                    break
            if break_T:
                print(f"    ✗ {name}: breaks at T = {break_T:.3f} T_P")
            else:
                max_T = max(e["T_planck_units"] for e in entries if e["computable"])
                print(f"    ✗ {name}: last computable at T = {max_T:.3f} T_P")
        print()

    print("  INTERPRETATION:")
    if survivors and breakers:
        print(f"  Models that maintain computable Benford fingerprints beyond T_Planck")
        print(f"  describe a universe where the singularity is resolved — physics")
        print(f"  continues through the wall. Models that go UNDEFINED hit a genuine")
        print(f"  breakdown — the distribution ceases to exist at the singularity.")
    print()

    # ── Save ─────────────────────────────────────────────────────────────
    output = {
        "experiment": "planck_wall",
        "description": "Benford existence filter applied to quantum gravity proposals at the Planck scale",
        "models_tested": list(models.keys()),
        "temperatures_planck_units": temperatures,
        "results": {name: entries for name, entries in results.items()},
        "survival": {
            "survivors": survivors,
            "breakers": breakers,
        },
        "physical_scales": {
            "T_Planck_kelvin": T_PLANCK_KELVIN,
            "E_Planck_GeV": E_PLANCK_GEV,
        },
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = "results/round_trip/planck_wall.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Results saved to {out_path}")

    return output


if __name__ == "__main__":
    run_planck_wall()
