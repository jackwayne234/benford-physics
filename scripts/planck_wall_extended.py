#!/usr/bin/env python3
"""Experiment 6c: Planck Wall — EXTENDED (10 Quantum Gravity Models)

Original 5 models from planck_wall.py:
  1. Standard (GR + QFT)
  2. Loop Quantum Gravity (LQG)
  3. GUP (Generalized Uncertainty Principle)
  4. DSR (Doubly Special Relativity)
  5. Hagedorn (String Theory)

New 5 models:
  6. Causal Set Theory — random discrete spacetime, Gaussian UV suppression
  7. Asymptotic Safety — running spectral dimension (4 → 2)
  8. Horava-Lifshitz — anisotropic scaling, E² = k² + k⁴ + k⁶
  9. Non-commutative Geometry — minimum area, E² = k² + θk⁴
  10. Causal Dynamical Triangulations (CDT) — sharp dimensional reduction (4 → 2)

Uses the same high-res temperature grid as planck_wall_hires.py.
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


# ══════════════════════════════════════════════════════════════════════════
# ORIGINAL 5 MODELS
# ══════════════════════════════════════════════════════════════════════════

def standard_spectrum(T, k):
    """Standard GR+QFT: E = k, g(k) = k²."""
    E = k
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = k[valid]**2 / denom[valid]
    return spectrum


def lqg_spectrum(T, k):
    """Loop Quantum Gravity: E = 2|sin(k/2)|, first Brillouin zone only."""
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
    """GUP: E = k√(1 + βk²), g(k) = k²/(1 + βk²)."""
    E = k * np.sqrt(1.0 + beta * k**2)
    g = k**2 / (1.0 + beta * k**2)
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = g[valid] / denom[valid]
    return spectrum


def dsr_spectrum(T, k):
    """DSR: E = 1 - e^{-k}, energy saturates at E_P = 1."""
    E = 1.0 - np.exp(-k)
    E = np.maximum(E, 1e-10)
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = k[valid]**2 / denom[valid]
    return spectrum


def hagedorn_spectrum(T, k, T_H=1.0):
    """String/Hagedorn: g(k) = k² × e^{k/T_H}, exponential state growth."""
    E = k
    growth = np.exp(np.minimum(k / T_H, 500))
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = k[valid]**2 * growth[valid] / denom[valid]
    return spectrum


# ══════════════════════════════════════════════════════════════════════════
# NEW 5 MODELS
# ══════════════════════════════════════════════════════════════════════════

def causal_set_spectrum(T, k):
    """Causal Set Theory: discrete random spacetime.

    The fundamental discreteness introduces a nonlocal d'Alembertian which
    produces a Gaussian UV suppression: modes above the Planck scale are
    exponentially damped.  The sprinkling is random (Poisson process), so
    mode density fluctuates — but in the mean-field limit the effect is:

        S(k) = k² × exp(-k²) / (e^{k/T} - 1)

    At low k: standard (exp(-k²) ≈ 1)
    At k ~ 1 (Planck): modes are suppressed by 1/e
    At k >> 1: modes vanish exponentially — spacetime runs out of points
    """
    E = k
    # Gaussian UV suppression from discreteness
    suppression = np.exp(-k**2)
    g = k**2 * suppression
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = g[valid] / denom[valid]
    return spectrum


def asymptotic_safety_spectrum(T, k):
    """Asymptotic Safety: running spectral dimension 4 → 2.

    The gravitational coupling runs to a UV fixed point.  The key observable
    prediction: the spectral dimension d_s(k) runs smoothly from 4 at low
    energy to 2 at the Planck scale and above.

        d_s(k) = 2 + 2 / (1 + (k/k_P)²)

    The density of states follows g(k) = k^{d_s - 1}:
        Low k:  d_s ≈ 4 → g(k) ≈ k³  (standard 3+1D)
        High k: d_s ≈ 2 → g(k) ≈ k¹  (effective 1+1D)

    Space itself loses two dimensions near the singularity.
    """
    E = k
    # Running spectral dimension
    d_s = 2.0 + 2.0 / (1.0 + k**2)  # k in Planck units, so k_P = 1
    g = k**(d_s - 1.0)
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = g[valid] / denom[valid]
    return spectrum


def horava_lifshitz_spectrum(T, k):
    """Horava-Lifshitz gravity: anisotropic scaling.

    Time and space scale differently at high energy.  The dispersion relation
    picks up higher spatial derivatives:

        E² = k² + k⁴ + k⁶   (in Planck units)

    At low k: E ≈ k (standard GR)
    At high k: E ≈ k³ (much steeper — modes cost much more energy)

    The k⁶ term is what makes gravity power-counting renormalizable.
    Density of states: standard k² (the modification is in dispersion, not DoS).
    """
    E_squared = k**2 + k**4 + k**6
    E = np.sqrt(E_squared)
    g = k**2
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = g[valid] / denom[valid]
    return spectrum


def noncommutative_spectrum(T, k, theta=1.0):
    """Non-commutative Geometry: spacetime coordinates don't commute.

    [x_μ, x_ν] = iθ_μν introduces a minimum area ~ θ.
    Modified dispersion: E² = k² + θk⁴
    UV/IR mixing modifies the density of states: g(k) = k²(1 + θk²)

    At low k: E ≈ k, g ≈ k²  (standard)
    At high k: E ≈ √θ k², g ≈ θk⁴  (modes get denser AND more expensive)

    The competing effects (more modes but higher energy per mode) make this
    model's behavior near the wall nontrivial.
    """
    E_squared = k**2 + theta * k**4
    E = np.sqrt(E_squared)
    g = k**2 * (1.0 + theta * k**2)
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = g[valid] / denom[valid]
    return spectrum


def cdt_spectrum(T, k):
    """Causal Dynamical Triangulations (CDT): sharp dimensional reduction.

    Like asymptotic safety, predicts spectral dimension drops from 4 to 2.
    But the mechanism is different: causal gluing of simplices with a
    preferred foliation creates a SHARPER transition than the smooth
    running of asymptotic safety.

        d_s(k) = 2 + 2 / (1 + (k/k_P)⁴)

    The k⁴ exponent (vs k² in AS) makes the transition more abrupt —
    a narrow crossover region rather than a gradual fade.

    At low k: d_s ≈ 4 → g(k) ≈ k³
    At high k: d_s ≈ 2 → g(k) ≈ k¹
    Transition width: narrower than asymptotic safety
    """
    E = k
    # Sharp dimensional reduction (k⁴ vs k² gives narrower transition)
    d_s = 2.0 + 2.0 / (1.0 + k**4)
    g = k**(d_s - 1.0)
    exponent = np.minimum(E / T, 500)
    denom = np.exp(exponent) - 1.0
    valid = denom > 0
    spectrum = np.zeros_like(k)
    spectrum[valid] = g[valid] / denom[valid]
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
    if len(temps) < 3:
        return {}
    t = np.array(temps)
    d = np.array(deltas)
    grad = np.gradient(d, t)
    grad2 = np.gradient(grad, t)
    max_grad_idx = int(np.argmax(np.abs(grad)))
    peaks = []
    valleys = []
    for i in range(1, len(d) - 1):
        if d[i] > d[i-1] and d[i] > d[i+1]:
            peaks.append({"T": float(t[i]), "delta_b": float(d[i])})
        if d[i] < d[i-1] and d[i] < d[i+1]:
            valleys.append({"T": float(t[i]), "delta_b": float(d[i])})
    return {
        "steepest_gradient": {
            "T": float(t[max_grad_idx]),
            "delta_b": float(d[max_grad_idx]),
            "gradient": float(grad[max_grad_idx]),
        },
        "peaks": peaks,
        "valleys": valleys,
        "gradient_at_wall": float(np.interp(1.0, t, grad)),
    }


# ── Main ─────────────────────────────────────────────────────────────────

def run_extended():
    """Run all 10 QG models through the Planck wall at high resolution."""

    # Same high-res temperature grid
    coarse_low = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4]
    fine = list(np.arange(0.50, 2.02, 0.02))
    coarse_high = [2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 50.0, 100.0]
    temperatures = sorted(set([round(t, 4) for t in coarse_low + fine + coarse_high]))

    k = np.linspace(0.001, 50, 100000)

    # All 10 models
    models = {
        # Original 5
        "Standard":      lambda T: standard_spectrum(T, k),
        "LQG":           lambda T: lqg_spectrum(T, k),
        "GUP":           lambda T: gup_spectrum(T, k),
        "DSR":           lambda T: dsr_spectrum(T, k),
        "Hagedorn":      lambda T: hagedorn_spectrum(T, k),
        # New 5
        "Causal Set":    lambda T: causal_set_spectrum(T, k),
        "Asym. Safety":  lambda T: asymptotic_safety_spectrum(T, k),
        "Horava-Lif.":   lambda T: horava_lifshitz_spectrum(T, k),
        "Noncommut.":    lambda T: noncommutative_spectrum(T, k),
        "CDT":           lambda T: cdt_spectrum(T, k),
    }

    results = {name: [] for name in models}

    print("=" * 82)
    print("EXPERIMENT 6c: Planck Wall — EXTENDED (10 Quantum Gravity Models)")
    print(f"  {len(temperatures)} temperature points, {len(models)} models")
    print("  All quantities in Planck units (E_P = T_P = 1)")
    print("=" * 82)
    print()

    for model_name, func in models.items():
        print(f"  {model_name:<16s}", end="", flush=True)

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
        n_undef = len(temperatures) - n_comp
        status = f"{n_comp}/{len(temperatures)} computable"
        if n_undef > 0:
            status += f"  ({n_undef} UNDEFINED)"
        print(status)

    print()

    # ══════════════════════════════════════════════════════════════════════
    # COMPARISON TABLE AT THE WALL
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 82)
    print("δ_B AT KEY TEMPERATURES")
    print("=" * 82)
    print()

    key_temps = [0.001, 0.1, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0, 5.0, 10.0, 100.0]

    # Print in two blocks of 5 models each for readability
    for block_name, block_models in [("ORIGINAL 5", list(models.keys())[:5]),
                                      ("NEW 5", list(models.keys())[5:])]:
        print(f"── {block_name} ──")
        print(f"  {'T/T_P':>7s}", end="")
        for name in block_models:
            print(f"  {name:>14s}", end="")
        print()
        print(f"  {'─'*7}", end="")
        for _ in block_models:
            print(f"  {'─'*14}", end="")
        print()

        for T in key_temps:
            line = f"  {T:7.3f}"
            for name in block_models:
                entry = next((e for e in results[name]
                              if abs(e["T_planck_units"] - T) < 0.0001), None)
                if entry and entry["computable"]:
                    line += f"  {entry['delta_b']:14.6f}"
                else:
                    line += f"  {'UNDEFINED':>14s}"
            if abs(T - 1.0) < 0.001:
                line += "  ← WALL"
            print(line)
        print()

    # ══════════════════════════════════════════════════════════════════════
    # TRANSITION ANALYSIS (new 5 only — original 5 already analyzed)
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 82)
    print("TRANSITION ANALYSIS — NEW MODELS")
    print("=" * 82)
    print()

    transitions = {}
    new_models = list(models.keys())[5:]

    for model_name in new_models:
        entries = results[model_name]
        computable = [(e["T_planck_units"], e["delta_b"])
                      for e in entries if e["computable"]]

        if len(computable) < 3:
            transitions[model_name] = {"error": "too few computable points"}
            print(f"── {model_name}: TOO FEW COMPUTABLE POINTS ──")
            print()
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
            top_peaks = sorted(trans["peaks"], key=lambda p: p["delta_b"], reverse=True)[:3]
            print(f"  Top peaks:")
            for p in top_peaks:
                print(f"    T/T_P = {p['T']:.4f}, δ_B = {p['delta_b']:.6f}")

        if trans["valleys"]:
            bot_valleys = sorted(trans["valleys"], key=lambda v: v["delta_b"])[:3]
            print(f"  Deepest valleys:")
            for v in bot_valleys:
                print(f"    T/T_P = {v['T']:.4f}, δ_B = {v['delta_b']:.6f}")

        print()

    # ══════════════════════════════════════════════════════════════════════
    # FULL RANKING — ALL 10 MODELS
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 82)
    print("FULL RANKING: All 10 Models — Post-Wall Statistical Quality")
    print("=" * 82)
    print()

    rankings = []
    for name in models:
        entries = results[name]
        post = [e for e in entries if e["T_planck_units"] > 1.0 and e["computable"]]
        pre = [e for e in entries if e["T_planck_units"] <= 1.0 and e["computable"]]
        all_comp = [e for e in entries if e["computable"]]

        if post:
            mean_post = np.mean([e["delta_b"] for e in post])
            min_post = min(e["delta_b"] for e in post)
            max_post = max(e["delta_b"] for e in post)
        else:
            mean_post = float("inf")
            min_post = float("inf")
            max_post = float("inf")

        if pre:
            mean_pre = np.mean([e["delta_b"] for e in pre])
        else:
            mean_pre = float("inf")

        n_comp = len(all_comp)
        n_undef = len(entries) - n_comp

        rankings.append({
            "name": name,
            "mean_post": mean_post,
            "min_post": min_post,
            "max_post": max_post,
            "mean_pre": mean_pre,
            "n_computable": n_comp,
            "n_undefined": n_undef,
            "n_post": len(post),
        })

    rankings.sort(key=lambda x: x["mean_post"])

    print(f"  {'Rank':>4s}  {'Model':>16s}  {'Mean δ_B':>10s}  {'Min δ_B':>10s}  "
          f"{'Max δ_B':>10s}  {'Comp':>5s}  {'Undef':>5s}")
    print(f"  {'─'*4}  {'─'*16}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*5}  {'─'*5}")

    for i, r in enumerate(rankings):
        if r["mean_post"] == float("inf"):
            print(f"  {i+1:4d}  {r['name']:>16s}  {'NO DATA':>10s}  {'':>10s}  "
                  f"{'':>10s}  {r['n_computable']:5d}  {r['n_undefined']:5d}")
        else:
            print(f"  {i+1:4d}  {r['name']:>16s}  {r['mean_post']:10.6f}  {r['min_post']:10.6f}  "
                  f"{r['max_post']:10.6f}  {r['n_computable']:5d}  {r['n_undefined']:5d}")

    print()

    # ══════════════════════════════════════════════════════════════════════
    # CHARACTER ANALYSIS — What kind of physics does each model describe?
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 82)
    print("CHARACTER ANALYSIS: How each model behaves at the wall")
    print("=" * 82)
    print()

    for name in models:
        entries = results[name]
        computable = [e for e in entries if e["computable"]]
        pre = [e for e in computable if e["T_planck_units"] <= 1.0]
        post = [e for e in computable if e["T_planck_units"] > 1.0]
        at_wall = next((e for e in computable
                        if abs(e["T_planck_units"] - 1.0) < 0.005), None)

        if not computable:
            print(f"  {name}: ALL UNDEFINED — does not survive at any temperature")
            print()
            continue

        pre_mean = np.mean([e["delta_b"] for e in pre]) if pre else None
        post_mean = np.mean([e["delta_b"] for e in post]) if post else None
        wall_db = at_wall["delta_b"] if at_wall else None

        # Classify behavior
        if post_mean is not None and pre_mean is not None:
            if post_mean < pre_mean * 0.5:
                character = "IMPROVES after the wall"
            elif post_mean > pre_mean * 2.0:
                character = "DEGRADES after the wall"
            elif abs(post_mean - pre_mean) / max(pre_mean, 0.001) < 0.3:
                character = "FLAT — wall doesn't matter"
            else:
                character = "MIXED behavior"
        elif post_mean is None:
            character = "BREAKS at or before the wall"
        else:
            character = "UNKNOWN"

        print(f"  {name}:")
        if pre_mean is not None:
            print(f"    Pre-wall mean δ_B:  {pre_mean:.6f}")
        if wall_db is not None:
            print(f"    At the wall:         {wall_db:.6f}")
        if post_mean is not None:
            print(f"    Post-wall mean δ_B: {post_mean:.6f}")
        print(f"    Character: {character}")

        n_undef = len(entries) - len(computable)
        if n_undef > 0:
            undef_temps = [e["T_planck_units"] for e in entries if not e["computable"]]
            print(f"    UNDEFINED at {n_undef} temperatures: "
                  f"T/T_P = {min(undef_temps):.3f}–{max(undef_temps):.3f}")

        print()

    # ══════════════════════════════════════════════════════════════════════
    # HEAD-TO-HEAD: Asymptotic Safety vs CDT (both predict d_s = 4 → 2)
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 82)
    print("HEAD-TO-HEAD: Asymptotic Safety vs CDT")
    print("  Both predict dimensional reduction 4 → 2, but differ in sharpness")
    print("=" * 82)
    print()

    as_entries = {e["T_planck_units"]: e for e in results["Asym. Safety"] if e["computable"]}
    cdt_entries = {e["T_planck_units"]: e for e in results["CDT"] if e["computable"]}

    compare_temps = [0.1, 0.5, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 5.0, 10.0]
    print(f"  {'T/T_P':>7s}  {'Asym.Safety':>12s}  {'CDT':>12s}  {'Difference':>12s}  {'Closer to Benford':>18s}")
    print(f"  {'─'*7}  {'─'*12}  {'─'*12}  {'─'*12}  {'─'*18}")

    for T in compare_temps:
        as_e = next((as_entries[t] for t in as_entries if abs(t - T) < 0.005), None)
        cdt_e = next((cdt_entries[t] for t in cdt_entries if abs(t - T) < 0.005), None)

        if as_e and cdt_e:
            diff = as_e["delta_b"] - cdt_e["delta_b"]
            winner = "Asym. Safety" if as_e["delta_b"] < cdt_e["delta_b"] else "CDT"
            marker = " ← WALL" if abs(T - 1.0) < 0.005 else ""
            print(f"  {T:7.1f}  {as_e['delta_b']:12.6f}  {cdt_e['delta_b']:12.6f}  "
                  f"{diff:+12.6f}  {winner:>18s}{marker}")

    print()

    # ══════════════════════════════════════════════════════════════════════
    # SAVE
    # ══════════════════════════════════════════════════════════════════════

    # Build rankings for JSON
    rankings_json = []
    for r in rankings:
        entry = {"model": r["name"], "n_computable": r["n_computable"],
                 "n_undefined": r["n_undefined"]}
        if r["mean_post"] != float("inf"):
            entry["mean_post_wall_delta_b"] = round(r["mean_post"], 6)
            entry["min_post_wall_delta_b"] = round(r["min_post"], 6)
            entry["max_post_wall_delta_b"] = round(r["max_post"], 6)
        rankings_json.append(entry)

    output = {
        "experiment": "planck_wall_extended",
        "description": "10 quantum gravity models through the Planck wall at high resolution",
        "n_temperature_points": len(temperatures),
        "models_tested": list(models.keys()),
        "model_descriptions": {
            "Standard": "GR+QFT, no modification. E=k, g(k)=k².",
            "LQG": "Polymer dispersion E=2|sin(k/2)|, Brillouin zone k<π.",
            "GUP": "Minimum length. E=k√(1+k²), g(k)=k²/(1+k²).",
            "DSR": "Energy saturation. E=1-e^{-k}, max energy = E_P.",
            "Hagedorn": "String theory. g(k)=k²·e^{k/T_H}, exponential state growth.",
            "Causal Set": "Discrete random spacetime. Gaussian UV suppression exp(-k²).",
            "Asym. Safety": "Running spectral dim 4→2. d_s=2+2/(1+k²).",
            "Horava-Lif.": "Anisotropic scaling. E²=k²+k⁴+k⁶.",
            "Noncommut.": "Minimum area. E²=k²+k⁴, g(k)=k²(1+k²).",
            "CDT": "Sharp dimensional reduction 4→2. d_s=2+2/(1+k⁴).",
        },
        "temperatures_planck_units": temperatures,
        "results": {name: entries for name, entries in results.items()},
        "transitions_new_models": transitions,
        "rankings": rankings_json,
        "physical_scales": {
            "T_Planck_kelvin": T_PLANCK_KELVIN,
            "E_Planck_GeV": E_PLANCK_GEV,
        },
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = "results/round_trip/planck_wall_extended.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Results saved to {out_path}")

    return output


if __name__ == "__main__":
    run_extended()
