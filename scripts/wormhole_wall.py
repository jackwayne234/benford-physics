#!/usr/bin/env python3
"""Experiment 8: Wormhole Wall — 10 QG Models Through a Morris-Thorne/Ellis Wormhole

The "wall" here is the wormhole throat at l = 0, where tidal forces peak.
Unlike the black hole and Big Bang, there is no singularity and no horizon —
the throat is traversable and symmetric.

Geometry (Morris-Thorne / Ellis):
  r(l) = sqrt(l² + b₀²)
  - l = proper distance from the throat
  - b₀ = throat radius (shape parameter)
  - Throat at l = 0: r = b₀ (minimum radius)
  - No singularity, no horizon, fully symmetric: l → -l

Temperature model (tidal):
  T(l) = (1/2π) × b₀² / (l² + b₀²)^{3/2}
  - Peaks at the throat: T(0) = 1/(2π b₀)
  - Falls off as |l|^{-3} far from throat
  - Same on both sides (symmetric)

Two observer types:
  A) TRAVERSING — standard spectra at tidal temperature
  B) CASIMIR — restricted modes due to throat geometry:
     casimir_factor = 1 - exp(-(k r_local / π)²)
     Suppresses long-wavelength modes that don't fit through the throat.
     Plus Pure Casimir: standard BE + mode restriction only (no QG).

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


# ── QG Model Generators (same as black_hole_wall.py) ────────────────

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


# ── Wormhole geometry ───────────────────────────────────────────────────

def tidal_temperature(l, b_0):
    """Tidal temperature at proper distance l from the throat.
    T(l) = (1/2π) × b₀² / (l² + b₀²)^{3/2}
    Peaks at throat (l=0): T(0) = 1/(2π b₀)
    """
    return (1.0 / (2.0 * math.pi)) * b_0**2 / (l**2 + b_0**2)**1.5


def casimir_modified_spectrum(base_func, T, k, r_local):
    """Apply Casimir mode restriction to a base spectrum.
    Suppresses long-wavelength modes that don't fit through the throat.
    """
    base = base_func(T, k)
    # Use matching k subset (handles LQG's Brillouin zone cutoff)
    k_eff = k[:len(base)]
    casimir_factor = 1.0 - np.exp(-(k_eff * r_local / np.pi)**2)
    return base * casimir_factor


# ── Main ─────────────────────────────────────────────────────────────────

def run_wormhole_wall():

    b_0 = 0.1  # Throat radius in Planck units

    # Sweep grid: 43 points, symmetric about the throat
    right_side = [0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07,
                  0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0]
    left_side = [-x for x in right_side]
    all_l = sorted(left_side + [0.0] + right_side)

    # Log-spaced k grid to handle wide temperature range
    # (throat is hot ~1.6 T_P, far regions are cold ~10^-6 T_P)
    k = np.logspace(-5, np.log10(50), 100000)

    models = {
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

    print("=" * 84)
    print("EXPERIMENT 8: Wormhole Wall — 10 QG Models Through a Morris-Thorne/Ellis Wormhole")
    print(f"  Throat radius: b_0 = {b_0} l_P")
    print(f"  Throat temperature: T(0) = {tidal_temperature(0, b_0):.6f} T_P")
    print(f"  {len(all_l)} proper-distance points from l = {min(all_l)} to l = {max(all_l)}")
    print("  Two observer types: traversing (standard) + Casimir (mode-restricted)")
    print("=" * 84)
    print()

    # ── Temperature profile ──────────────────────────────────────────────
    print("── TEMPERATURE PROFILE ──")
    print(f"  {'l/b_0':>8s}  {'l':>10s}  {'r_local':>10s}  {'T_tidal':>12s}")
    print(f"  {'─'*8}  {'─'*10}  {'─'*10}  {'─'*12}")

    for l in all_l:
        r_local = math.sqrt(l**2 + b_0**2)
        T = tidal_temperature(l, b_0)
        marker = ""
        if abs(l) < 0.0005:
            marker = " <- THROAT"
        print(f"  {l/b_0:8.2f}  {l:10.4f}  {r_local:10.4f}  {T:12.8f}{marker}")

    print()

    # ── Run traversing observer (all l) ──────────────────────────────────
    print("=" * 84)
    print("TRAVERSING OBSERVER — Standard spectra at tidal temperature")
    print("=" * 84)
    print()

    trav_results = {name: [] for name in models}

    for model_name, func in models.items():
        print(f"  {model_name:<16s}", end="", flush=True)

        for l in all_l:
            T = tidal_temperature(l, b_0)
            r_local = math.sqrt(l**2 + b_0**2)

            spectrum = func(T, k)
            result = analyze_spectrum(spectrum)

            entry = {
                "l": l,
                "l_over_b0": round(l / b_0, 6),
                "r_local": round(r_local, 8),
                "T_tidal": round(T, 8),
            }

            if result:
                entry.update(result)
                entry["computable"] = True
            else:
                entry["delta_b"] = None
                entry["computable"] = False
                n_pos = len([v for v in spectrum if v > 0 and np.isfinite(v)])
                entry["n_positive"] = n_pos

            trav_results[model_name].append(entry)

        n_comp = sum(1 for e in trav_results[model_name] if e["computable"])
        print(f"{n_comp}/{len(all_l)} computable")

    print()

    # ── Run Casimir observer (all l, mode-restricted) ────────────────────
    print("=" * 84)
    print("CASIMIR OBSERVER — Mode-restricted spectra")
    print("=" * 84)
    print()

    cas_results = {name: [] for name in models}

    for model_name, func in models.items():
        print(f"  {model_name:<16s}", end="", flush=True)

        for l in all_l:
            T = tidal_temperature(l, b_0)
            r_local = math.sqrt(l**2 + b_0**2)

            spectrum = casimir_modified_spectrum(func, T, k, r_local)
            result = analyze_spectrum(spectrum)

            entry = {
                "l": l,
                "l_over_b0": round(l / b_0, 6),
                "r_local": round(r_local, 8),
                "T_tidal": round(T, 8),
            }

            if result:
                entry.update(result)
                entry["computable"] = True
            else:
                entry["delta_b"] = None
                entry["computable"] = False
                n_pos = len([v for v in spectrum if v > 0 and np.isfinite(v)])
                entry["n_positive"] = n_pos

            cas_results[model_name].append(entry)

        n_comp = sum(1 for e in cas_results[model_name] if e["computable"])
        print(f"{n_comp}/{len(all_l)} computable")

    print()

    # ── Pure Casimir (standard BE + mode restriction, no QG) ─────────────
    print("── PURE CASIMIR (standard + mode restriction) ──")

    pure_casimir_results = []
    for l in all_l:
        T = tidal_temperature(l, b_0)
        r_local = math.sqrt(l**2 + b_0**2)

        spectrum = casimir_modified_spectrum(standard_spectrum, T, k, r_local)
        result = analyze_spectrum(spectrum)

        entry = {
            "l": l,
            "l_over_b0": round(l / b_0, 6),
            "r_local": round(r_local, 8),
            "T_tidal": round(T, 8),
        }

        if result:
            entry.update(result)
            entry["computable"] = True
        else:
            entry["delta_b"] = None
            entry["computable"] = False
            n_pos = len([v for v in spectrum if v > 0 and np.isfinite(v)])
            entry["n_positive"] = n_pos

        pure_casimir_results.append(entry)

    n_comp = sum(1 for e in pure_casimir_results if e["computable"])
    print(f"  Pure Casimir: {n_comp}/{len(all_l)} computable")
    print()

    # ══════════════════════════════════════════════════════════════════════
    # TRAVERSING OBSERVER: δ_B vs l/b_0
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 84)
    print("TRAVERSING OBSERVER: delta_B across the wormhole")
    print("=" * 84)
    print()

    for block_name, block_models in [("ORIGINAL 5", list(models.keys())[:5]),
                                      ("NEW 5", list(models.keys())[5:])]:
        print(f"── {block_name} (traversing) ──")
        print(f"  {'l/b_0':>8s}  {'T_tidal':>10s}", end="")
        for name in block_models:
            print(f"  {name:>12s}", end="")
        print()
        print(f"  {'─'*8}  {'─'*10}", end="")
        for _ in block_models:
            print(f"  {'─'*12}", end="")
        print()

        display_l = [-10.0, -5.0, -1.0, -0.5, -0.1, -0.01,
                     0.0,
                     0.01, 0.1, 0.5, 1.0, 5.0, 10.0]

        for l in display_l:
            T = tidal_temperature(l, b_0)
            line = f"  {l/b_0:8.1f}  {T:10.6f}"

            for name in block_models:
                entry = next((e for e in trav_results[name]
                              if abs(e["l"] - l) < 0.0005), None)
                if entry and entry["computable"]:
                    line += f"  {entry['delta_b']:12.6f}"
                else:
                    line += f"  {'UNDEF':>12s}"

            if abs(l) < 0.005:
                line += "  <- THROAT"
            print(line)
        print()

    # ══════════════════════════════════════════════════════════════════════
    # ANALYSIS 1: SYMMETRY TEST
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 84)
    print("SYMMETRY TEST: delta_B(l) = delta_B(-l) ?")
    print("=" * 84)
    print()

    symmetry_test = {}

    for name in models:
        max_asym = 0.0

        for l in right_side:
            entry_pos = next((e for e in trav_results[name]
                              if abs(e["l"] - l) < 1e-8), None)
            entry_neg = next((e for e in trav_results[name]
                              if abs(e["l"] - (-l)) < 1e-8), None)

            if (entry_pos and entry_pos["computable"] and
                entry_neg and entry_neg["computable"]):
                asym = abs(entry_pos["delta_b"] - entry_neg["delta_b"])
                max_asym = max(max_asym, asym)

        is_symmetric = max_asym < 1e-10
        status = "SYMMETRIC" if is_symmetric else "ASYMMETRIC"
        print(f"  {name:<16s}  max_asym = {max_asym:.2e}  -> {status}")

        symmetry_test[name] = {
            "max_asymmetry": float(max_asym),
            "is_symmetric": is_symmetric,
        }

    print()

    # ══════════════════════════════════════════════════════════════════════
    # ANALYSIS 2: CASIMIR COMPARISON (at throat)
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 84)
    print("CASIMIR EFFECT: delta_B with vs without mode restriction (at throat)")
    print("=" * 84)
    print()

    casimir_effect = {}

    print(f"  {'Model':>16s}  {'Standard':>12s}  {'Casimir':>12s}  {'Delta':>10s}  {'%':>8s}")
    print(f"  {'─'*16}  {'─'*12}  {'─'*12}  {'─'*10}  {'─'*8}")

    for name in models:
        trav_throat = next((e for e in trav_results[name]
                           if abs(e["l"]) < 1e-8), None)
        cas_throat = next((e for e in cas_results[name]
                          if abs(e["l"]) < 1e-8), None)

        if (trav_throat and trav_throat["computable"] and
            cas_throat and cas_throat["computable"]):
            db_std = trav_throat["delta_b"]
            db_cas = cas_throat["delta_b"]
            diff = db_cas - db_std
            pct = 100.0 * diff / db_std if db_std > 0 else 0.0
            print(f"  {name:>16s}  {db_std:12.6f}  {db_cas:12.6f}  {diff:+10.6f}  {pct:+8.1f}%")

            casimir_effect[name] = {
                "delta_b_standard": round(float(db_std), 6),
                "delta_b_casimir": round(float(db_cas), 6),
                "difference": round(float(diff), 6),
                "percent_change": round(float(pct), 2),
            }

    # Pure Casimir at throat
    pure_throat = next((e for e in pure_casimir_results
                       if abs(e["l"]) < 1e-8), None)
    if pure_throat and pure_throat["computable"]:
        print(f"  {'Pure Casimir':>16s}  {'---':>12s}  {pure_throat['delta_b']:12.6f}  {'---':>10s}  {'---':>8s}")
        casimir_effect["Pure Casimir"] = {
            "delta_b_casimir": round(float(pure_throat["delta_b"]), 6),
        }

    print()

    # ══════════════════════════════════════════════════════════════════════
    # ANALYSIS 3: THROAT SIZE SWEEP
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 84)
    print("THROAT SIZE SWEEP: delta_B at throat (l=0) for various b_0")
    print("=" * 84)
    print()

    b0_values = [0.01, 0.05, 0.1, 0.5, 1.0]
    throat_sweep = {}

    print(f"  {'b_0':>6s}  {'T(0)':>10s}", end="")
    for name in list(models.keys())[:5]:
        print(f"  {name:>12s}", end="")
    print()
    print(f"  {'─'*6}  {'─'*10}", end="")
    for _ in list(models.keys())[:5]:
        print(f"  {'─'*12}", end="")
    print()

    for b in b0_values:
        T_throat = tidal_temperature(0, b)
        line = f"  {b:6.2f}  {T_throat:10.6f}"
        throat_sweep[str(b)] = {"T_throat": round(float(T_throat), 8), "models": {}}

        for name, func in models.items():
            spectrum = func(T_throat, k)
            result = analyze_spectrum(spectrum)
            if result:
                db = result["delta_b"]
                throat_sweep[str(b)]["models"][name] = round(float(db), 6)
                if name in list(models.keys())[:5]:
                    line += f"  {db:12.6f}"
            else:
                throat_sweep[str(b)]["models"][name] = None
                if name in list(models.keys())[:5]:
                    line += f"  {'UNDEF':>12s}"

        print(line)

    print()

    # Print second block (models 6-10)
    print(f"  {'b_0':>6s}  {'T(0)':>10s}", end="")
    for name in list(models.keys())[5:]:
        print(f"  {name:>12s}", end="")
    print()
    print(f"  {'─'*6}  {'─'*10}", end="")
    for _ in list(models.keys())[5:]:
        print(f"  {'─'*12}", end="")
    print()

    for b in b0_values:
        T_throat = tidal_temperature(0, b)
        line = f"  {b:6.2f}  {T_throat:10.6f}"
        for name in list(models.keys())[5:]:
            db = throat_sweep[str(b)]["models"].get(name)
            if db is not None:
                line += f"  {db:12.6f}"
            else:
                line += f"  {'UNDEF':>12s}"
        print(line)

    print()

    # ══════════════════════════════════════════════════════════════════════
    # ANALYSIS 4: THREE-WALL COMPARISON
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 84)
    print("THREE-WALL COMPARISON: Wormhole vs Black Hole vs Big Bang")
    print("=" * 84)
    print()

    three_wall = {}

    # Load Big Bang and Black Hole data
    bb_path = "results/round_trip/planck_wall_extended.json"
    bh_path = "results/round_trip/black_hole_wall.json"
    bb_data = None
    bh_data = None

    if os.path.exists(bb_path):
        with open(bb_path) as f:
            bb_data = json.load(f)
    else:
        print("  (Big Bang results not found — run planck_wall_extended.py first)")

    if os.path.exists(bh_path):
        with open(bh_path) as f:
            bh_data = json.load(f)
    else:
        print("  (Black Hole results not found — run black_hole_wall.py first)")

    if bb_data and bh_data:
        print(f"  {'Model':>16s}  {'Wormhole':>12s}  {'Black Hole':>12s}  {'Big Bang':>12s}  {'WH<BH<BB?':>10s}")
        print(f"  {'─'*16}  {'─'*12}  {'─'*12}  {'─'*12}  {'─'*10}")

        for name in models:
            # Wormhole: mean delta_B across all computable traversing entries
            wh_entries = [e for e in trav_results[name] if e["computable"]]
            wh_mean = float(np.mean([e["delta_b"] for e in wh_entries])) if wh_entries else None

            # Black Hole: mean delta_B inside (r < r_s)
            bh_entries = bh_data["infalling_observer"].get(name, [])
            bh_inside = [e for e in bh_entries
                         if 0 < e.get("r_ratio", 0) < 1.0 and e.get("computable", False)]
            bh_mean = float(np.mean([e["delta_b"] for e in bh_inside])) if bh_inside else None

            # Big Bang: mean delta_B post-wall (T > T_P)
            bb_entries = bb_data["results"].get(name, [])
            bb_post = [e for e in bb_entries
                       if e.get("T_planck_units", 0) > 1.0 and e.get("computable", False)]
            bb_mean = float(np.mean([e["delta_b"] for e in bb_post])) if bb_post else None

            wh_str = f"{wh_mean:12.6f}" if wh_mean is not None else f"{'N/A':>12s}"
            bh_str = f"{bh_mean:12.6f}" if bh_mean is not None else f"{'N/A':>12s}"
            bb_str = f"{bb_mean:12.6f}" if bb_mean is not None else f"{'N/A':>12s}"

            ordering = "---"
            if wh_mean is not None and bh_mean is not None and bb_mean is not None:
                if wh_mean < bh_mean < bb_mean:
                    ordering = "YES"
                else:
                    ordering = "no"

            print(f"  {name:>16s}  {wh_str}  {bh_str}  {bb_str}  {ordering:>10s}")

            three_wall[name] = {
                "wormhole_mean": round(float(wh_mean), 6) if wh_mean is not None else None,
                "black_hole_mean": round(float(bh_mean), 6) if bh_mean is not None else None,
                "big_bang_mean": round(float(bb_mean), 6) if bb_mean is not None else None,
                "wh_lt_bh_lt_bb": ordering == "YES",
            }

        print()
        print("  WH<BH<BB = wormhole delta_B < black hole delta_B < big bang delta_B")
        print("  (expected: wormhole is the gentlest wall — no singularity, no horizon)")

    print()

    # ══════════════════════════════════════════════════════════════════════
    # RANKINGS
    # ══════════════════════════════════════════════════════════════════════

    print("=" * 84)
    print("RANKING: Mean delta_B through the wormhole (traversing observer)")
    print("=" * 84)
    print()

    rankings = []
    for name in models:
        entries = [e for e in trav_results[name] if e["computable"]]
        if entries:
            mean_db = float(np.mean([e["delta_b"] for e in entries]))
            min_db = float(min(e["delta_b"] for e in entries))
            max_db = float(max(e["delta_b"] for e in entries))
        else:
            mean_db = float("inf")
            min_db = float("inf")
            max_db = float("inf")

        throat_entry = next((e for e in trav_results[name]
                            if abs(e["l"]) < 1e-8 and e["computable"]), None)
        throat_db = float(throat_entry["delta_b"]) if throat_entry else None

        rankings.append({
            "name": name,
            "mean_delta_b": mean_db,
            "min_delta_b": min_db,
            "max_delta_b": max_db,
            "throat_delta_b": throat_db,
            "n_computable": len(entries),
            "n_total": len(all_l),
        })

    rankings.sort(key=lambda x: x["mean_delta_b"])

    print(f"  {'Rank':>4s}  {'Model':>16s}  {'Mean delta_B':>12s}  {'Throat delta_B':>14s}  {'Min':>10s}  {'Max':>10s}")
    print(f"  {'─'*4}  {'─'*16}  {'─'*12}  {'─'*14}  {'─'*10}  {'─'*10}")

    for i, r in enumerate(rankings):
        mean_str = f"{r['mean_delta_b']:12.6f}" if r["mean_delta_b"] != float("inf") else f"{'NO DATA':>12s}"
        throat_str = f"{r['throat_delta_b']:14.6f}" if r["throat_delta_b"] is not None else f"{'N/A':>14s}"
        min_str = f"{r['min_delta_b']:10.6f}" if r["min_delta_b"] != float("inf") else f"{'---':>10s}"
        max_str = f"{r['max_delta_b']:10.6f}" if r["max_delta_b"] != float("inf") else f"{'---':>10s}"
        print(f"  {i+1:4d}  {r['name']:>16s}  {mean_str}  {throat_str}  {min_str}  {max_str}")

    print()

    # Total computable count
    total_comp = sum(sum(1 for e in trav_results[name] if e["computable"])
                     for name in models)
    total_entries = len(models) * len(all_l)
    print(f"  Total: {total_comp}/{total_entries} entries computable")
    print()

    # ══════════════════════════════════════════════════════════════════════
    # SAVE
    # ══════════════════════════════════════════════════════════════════════

    output = {
        "experiment": "wormhole_wall",
        "description": "10 QG models through a Morris-Thorne/Ellis traversable wormhole",
        "wormhole": {
            "type": "Morris-Thorne / Ellis",
            "b_0": b_0,
            "b_0_description": "throat radius in Planck units",
            "T_throat": round(float(tidal_temperature(0, b_0)), 8),
            "geometry": "r(l) = sqrt(l^2 + b_0^2)",
            "properties": "no singularity, no horizon, symmetric",
        },
        "proper_distance_points": all_l,
        "models_tested": list(models.keys()),
        "traversing_observer": {name: entries for name, entries in trav_results.items()},
        "casimir_observer": {name: entries for name, entries in cas_results.items()},
        "pure_casimir": pure_casimir_results,
        "rankings": [
            {"model": r["name"],
             "mean_delta_b": round(r["mean_delta_b"], 6) if r["mean_delta_b"] != float("inf") else None,
             "throat_delta_b": round(r["throat_delta_b"], 6) if r["throat_delta_b"] is not None else None,
             "n_computable": r["n_computable"]}
            for r in rankings
        ],
        "symmetry_test": symmetry_test,
        "casimir_effect": casimir_effect,
        "throat_size_sweep": throat_sweep,
        "three_wall_comparison": three_wall,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    os.makedirs("results/round_trip", exist_ok=True)
    out_path = "results/round_trip/wormhole_wall.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Results saved to {out_path}")

    return output


if __name__ == "__main__":
    run_wormhole_wall()
