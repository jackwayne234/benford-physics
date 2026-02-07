#!/usr/bin/env python3
"""Experiment 5: The Whiteboard — Existence Filter for Exotic Physics.

Run every candidate through the Benford filter and report:
  - EXISTS:    computable δ_B with characteristic fingerprint
  - MARGINAL:  computable but degrading (losing modes)
  - UNDEFINED: no valid statistical sample — does not exist as physical distribution

Candidates:
  1. Anyons          — fractional exclusion statistics (between BE and FD)
  2. Negative-mass bosons   — negative energy occupation
  3. Negative-mass fermions — negative energy, Pauli exclusion
  4. Phantom energy   — wrong-sign kinetic term (w < -1 dark energy)
  5. Hawking radiation — greybody-modified black hole spectrum
  6. Unruh radiation  — accelerating observer thermal bath
  7. Gravitons        — massless spin-2 bosons
  8. Majorana fermions — self-conjugate (their own antiparticle)
  9. Axion-like boson — ultralight dark matter candidate
  10. Sterile neutrino — non-thermally produced fermion
"""

import json
import math
import os
import sys
import numpy as np
from scipy.optimize import brentq
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benford_core import run_full_analysis, euclidean_deviation, per_digit_deviation, observed_distribution, first_digit


# ── Helpers ──────────────────────────────────────────────────────────────

def analyze_values(values):
    """Run Benford analysis on a list of values, return summary dict or None."""
    positive = [v for v in values if isinstance(v, (int, float)) and v > 0 and math.isfinite(v)]
    if len(positive) < 50:
        return None
    result = run_full_analysis(positive)
    fd = result["first_digit"]
    eps = fd["per_digit_deviation"]
    return {
        "delta_b": fd["euclidean_deviation"],
        "epsilon_d": {str(k): v for k, v in eps.items()},
        "mad": fd["mad"]["mad"],
        "chi_sq_p": fd["chi_squared"]["p_value"],
        "n": result["n"],
        "verdict": result["verdict"],
    }


def classify_shape(eps_dict):
    """Classify ε(d) vector shape."""
    eps = [eps_dict[str(d)] for d in range(1, 10)]
    if all(abs(e) < 0.002 for e in eps):
        return "flat"
    tail = eps[1:]
    sign_changes = sum(
        1 for i in range(len(tail) - 1)
        if tail[i] * tail[i + 1] < 0 and abs(tail[i]) > 0.0005
    )
    if sign_changes >= 2:
        return "oscillatory"
    if abs(eps[0]) > 0.002 and abs(eps[-1]) > 0.002 and eps[0] * eps[-1] < 0:
        return "sweep"
    diffs = [tail[i + 1] - tail[i] for i in range(len(tail) - 1)]
    if all(d >= -0.0005 for d in diffs):
        return "monotone"
    if all(d <= 0.0005 for d in diffs):
        return "monotone"
    return "mixed"


def eps_distance(eps1, eps2):
    """L2 distance between two ε(d) dicts."""
    return math.sqrt(sum((eps1[str(d)] - eps2[str(d)])**2 for d in range(1, 10)))


# ── Reference fingerprints (from completed experiments) ─────────────────

REFERENCE = {
    "BE":    {"delta_b": 0.005627, "epsilon_d": {"1": -0.00433, "2": 0.003189, "3": 0.001411, "4": 0.00062, "5": 0.000159, "6": -6.7e-05, "7": -0.000232, "8": -0.000363, "9": -0.000387}},
    "FD":    {"delta_b": 0.011736, "epsilon_d": {"1": -0.00659, "2": 0.00521, "3": 0.00496, "4": 0.00515, "5": -0.00239, "6": -0.00196, "7": -0.00166, "8": -0.00146, "9": -0.00126}},
    "MB":    {"delta_b": 0.00985,  "epsilon_d": {"1": -0.00917, "2": 0.00232, "3": 0.00164, "4": 0.00128, "5": 0.00104, "6": 0.00087, "7": 0.00077, "8": 0.00067, "9": 0.00058}},
    "Planck": {"delta_b": 0.027921, "epsilon_d": {"1": 0.02233, "2": -0.014891, "3": -0.006039, "4": -0.00343, "5": -0.001751, "6": -0.000457, "7": 0.000468, "8": 0.001467, "9": 0.002303}},
}


def find_nearest_reference(eps_dict):
    """Find the closest known distribution by ε(d) distance."""
    best_name, best_dist = None, float("inf")
    for name, ref in REFERENCE.items():
        d = eps_distance(eps_dict, ref["epsilon_d"])
        if d < best_dist:
            best_dist = d
            best_name = name
    return best_name, best_dist


# ══════════════════════════════════════════════════════════════════════════
# CANDIDATE GENERATORS
# ══════════════════════════════════════════════════════════════════════════

def generate_anyon(g, n_points=100000):
    """Haldane fractional exclusion statistics.

    Solve: w^g · (1+w)^{1-g} = e^x  for w, then n = 1/(w + g).
    g=0: bosons (BE), g=0.5: semions, g=1: fermions (FD).
    """
    x = np.linspace(0.001, 50, n_points)
    z = np.exp(x)
    values = []

    for i in range(n_points):
        xi = x[i]
        zi = z[i]
        # f(w) = g*ln(w) + (1-g)*ln(1+w) - xi = 0, w > 0
        try:
            if g == 0:
                w = zi - 1.0
            elif g == 1:
                w = zi
            else:
                def eq(w):
                    return g * math.log(w) + (1 - g) * math.log(1 + w) - xi
                # Bracket: w is between something small and e^x
                w_lo = max(1e-15, zi * 0.01)
                w_hi = zi * 10
                # Ensure bracket
                if eq(w_lo) * eq(w_hi) > 0:
                    w_lo = 1e-15
                    w_hi = zi * 100
                w = brentq(eq, w_lo, w_hi, xtol=1e-12)
            n_val = 1.0 / (w + g)
            if n_val > 0 and math.isfinite(n_val):
                values.append(n_val)
        except (ValueError, ZeroDivisionError):
            continue

    return values


def generate_negative_mass_boson(m_abs=5.0, n_points=100000):
    """Negative-mass boson: energy E = -√(k² + m²), occupation = 1/(e^{E/T} - 1).

    With E < 0, the exponential e^{E/T} < 1, so denominator < 0, making n < 0.
    Expect: all negative → UNDEFINED.
    """
    k = np.linspace(0.001, 50, n_points)
    E = -np.sqrt(k**2 + m_abs**2)
    # n = 1/(e^E - 1), with E < 0
    occupation = 1.0 / (np.exp(np.maximum(E, -500)) - 1.0)
    return [float(v) for v in occupation if v > 0 and math.isfinite(v)]


def generate_negative_mass_fermion(m_abs=5.0, n_points=100000):
    """Negative-mass fermion: energy E = -√(k² + m²), occupation = 1/(e^{E/T} + 1).

    With E < 0, e^{E/T} < 1, but denominator = e^{E/T} + 1 > 1, so n > 0.
    This gives n = 1 - n_FD(|E|) — an inverted Fermi-Dirac.
    """
    k = np.linspace(0.001, 50, n_points)
    E = -np.sqrt(k**2 + m_abs**2)
    occupation = 1.0 / (np.exp(np.maximum(E, -500)) + 1.0)
    return [float(v) for v in occupation if v > 0 and math.isfinite(v)]


def generate_phantom_energy(n_points=100000):
    """Phantom energy field (w < -1): wrong-sign kinetic term.

    The field has negative kinetic energy, so mode energies are E_k = -ω_k.
    Bosonic occupation: n = 1/(e^{-ω/T} - 1).
    Since e^{-ω} < 1 for ω > 0, denominator < 0 → n < 0.
    Expect: UNDEFINED.
    """
    omega = np.linspace(0.001, 50, n_points)
    occupation = 1.0 / (np.exp(-omega) - 1.0)
    return [float(v) for v in occupation if v > 0 and math.isfinite(v)]


def generate_hawking_greybody(omega_c=2.0, n_points=100000):
    """Hawking radiation with greybody factor.

    Spectrum: N(ω) ∝ ω² · Γ(ω) / (e^{ω/T} - 1)
    Greybody: Γ(ω) = 1 - exp(-(ω/ω_c)²)

    This suppresses low-frequency modes (unlike pure Planck).
    At high ω, Γ → 1 and spectrum ≈ Planck.
    """
    omega = np.linspace(0.001, 50, n_points)
    greybody = 1.0 - np.exp(-(omega / omega_c)**2)
    spectrum = omega**2 * greybody / (np.exp(np.minimum(omega, 500)) - 1.0)
    return [float(v) for v in spectrum if v > 0 and math.isfinite(v)]


def generate_unruh(n_points=100000):
    """Unruh radiation: pure Planck spectrum at Unruh temperature.

    Identical to thermal radiation — an accelerating observer sees the
    Minkowski vacuum as a thermal bath.  n(ω) = ω³ / (e^ω - 1) (Planck).
    Should match Planck fingerprint exactly.
    """
    omega = np.linspace(0.001, 50, n_points)
    spectrum = omega**3 / (np.exp(np.minimum(omega, 500)) - 1.0)
    return [float(v) for v in spectrum if v > 0 and math.isfinite(v)]


def generate_graviton_thermal(n_points=100000):
    """Thermal gravitons: massless spin-2 bosons.

    Same Planck spectrum as photons in 3+1D (2 polarization states,
    density of states ∝ ω²). The occupation per mode is pure BE.
    Spectrum: N(ω) ∝ ω² · 1/(e^ω - 1) ∝ ω³/(e^ω - 1) (with density of states).
    Overall multiplicative constant (2 polarizations) doesn't affect first digits.
    Should match Planck fingerprint.
    """
    omega = np.linspace(0.001, 50, n_points)
    # ω² density of states × ω energy per mode × 1/(e^ω - 1) = ω³/(e^ω - 1)
    spectrum = omega**3 / (np.exp(np.minimum(omega, 500)) - 1.0)
    return [float(v) for v in spectrum if v > 0 and math.isfinite(v)]


def generate_majorana(n_points=100000):
    """Majorana fermion: self-conjugate (own antiparticle).

    Occupation per mode is still Fermi-Dirac: n = 1/(e^ε + 1).
    The Majorana condition halves degrees of freedom vs Dirac,
    but an overall factor doesn't change first digits.
    Should match FD fingerprint exactly.
    """
    x = np.linspace(0.001, 50, n_points)
    occupation = 1.0 / (np.exp(x) + 1.0)
    return [float(v) for v in occupation if v > 0 and math.isfinite(v)]


def generate_axion(m_over_T=0.001, n_points=100000):
    """Axion-like ultralight boson.

    Relativistic BE with very small mass: E = √(k² + m²), m/T << 1.
    Should be very close to massless BE.
    """
    k = np.linspace(0.001, 50, n_points)
    E = np.sqrt(k**2 + m_over_T**2)
    occupation = 1.0 / (np.exp(np.minimum(E, 500)) - 1.0)
    return [float(v) for v in occupation if v > 0 and math.isfinite(v)]


def generate_sterile_neutrino(n_points=100000):
    """Sterile neutrino (Dodelson-Widrow-like non-thermal production).

    The DW spectrum is approximately: f(p) ∝ p · n_FD(p)
    (one extra power of momentum from production rate scaling).
    This is FD with a momentum-dependent modification.
    """
    p = np.linspace(0.001, 50, n_points)
    spectrum = p / (np.exp(np.minimum(p, 500)) + 1.0)
    return [float(v) for v in spectrum if v > 0 and math.isfinite(v)]


# ══════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════

def run_whiteboard():
    """Run all candidates through the existence filter."""

    print("=" * 78)
    print("THE WHITEBOARD: Exotic Physics Existence Filter")
    print("  Running all candidates through δ_B / ε(d) analysis")
    print("=" * 78)
    print()

    candidates = []

    # ── 1. Anyons (fractional exclusion) ─────────────────────────────────
    print("─── ANYONS (Haldane Fractional Exclusion Statistics) ───")
    print("  g=0: boson,  g=0.5: semion,  g=1: fermion")
    print()

    for g in [0, 0.25, 0.5, 0.75, 1.0]:
        label = {0: "boson", 0.25: "1/4-on", 0.5: "semion", 0.75: "3/4-on", 1.0: "fermion"}.get(g, f"g={g}")
        print(f"  Computing g={g:.2f} ({label})...", end=" ", flush=True)
        values = generate_anyon(g)
        result = analyze_values(values)
        if result:
            shape = classify_shape(result["epsilon_d"])
            nearest, dist = find_nearest_reference(result["epsilon_d"])
            print(f"δ_B={result['delta_b']:.6f}  shape={shape}  nearest={nearest} (d={dist:.4f})")
            candidates.append({
                "name": f"Anyon g={g} ({label})",
                "category": "Fractional Statistics",
                "physics": f"Haldane exclusion g={g}",
                "status": "EXISTS",
                **result,
                "shape": shape,
                "nearest_known": nearest,
                "nearest_distance": round(dist, 6),
            })
        else:
            print("UNDEFINED")
            candidates.append({
                "name": f"Anyon g={g} ({label})",
                "category": "Fractional Statistics",
                "physics": f"Haldane exclusion g={g}",
                "status": "UNDEFINED",
                "delta_b": None, "epsilon_d": None, "n": len(values),
            })

    print()

    # ── 2. Negative-mass bosons ──────────────────────────────────────────
    print("─── NEGATIVE-MASS BOSONS ───")
    print("  E = -√(k²+m²),  n = 1/(e^{E/T} - 1)")
    print()

    for m in [1.0, 5.0, 10.0]:
        print(f"  |m|/T = {m:.1f}...", end=" ", flush=True)
        values = generate_negative_mass_boson(m)
        result = analyze_values(values)
        if result:
            shape = classify_shape(result["epsilon_d"])
            nearest, dist = find_nearest_reference(result["epsilon_d"])
            print(f"δ_B={result['delta_b']:.6f}  shape={shape}  n={result['n']}")
            candidates.append({
                "name": f"Negative-mass boson |m|/T={m}",
                "category": "Exotic Mass",
                "physics": "Negative energy bosonic occupation",
                "status": "EXISTS",
                **result,
                "shape": shape,
                "nearest_known": nearest,
                "nearest_distance": round(dist, 6),
            })
        else:
            n_positive = len(values)
            print(f"UNDEFINED ({n_positive} positive values out of 100000)")
            candidates.append({
                "name": f"Negative-mass boson |m|/T={m}",
                "category": "Exotic Mass",
                "physics": "Negative energy bosonic occupation",
                "status": "UNDEFINED",
                "delta_b": None, "epsilon_d": None, "n": n_positive,
                "reason": "All occupation numbers negative (e^{E/T} < 1 → denominator < 0)",
            })

    print()

    # ── 3. Negative-mass fermions ────────────────────────────────────────
    print("─── NEGATIVE-MASS FERMIONS ───")
    print("  E = -√(k²+m²),  n = 1/(e^{E/T} + 1)")
    print()

    for m in [1.0, 5.0, 10.0]:
        print(f"  |m|/T = {m:.1f}...", end=" ", flush=True)
        values = generate_negative_mass_fermion(m)
        result = analyze_values(values)
        if result:
            shape = classify_shape(result["epsilon_d"])
            nearest, dist = find_nearest_reference(result["epsilon_d"])
            print(f"δ_B={result['delta_b']:.6f}  shape={shape}  n={result['n']}  nearest={nearest}")
            candidates.append({
                "name": f"Negative-mass fermion |m|/T={m}",
                "category": "Exotic Mass",
                "physics": "Negative energy fermionic occupation (inverted FD)",
                "status": "EXISTS",
                **result,
                "shape": shape,
                "nearest_known": nearest,
                "nearest_distance": round(dist, 6),
            })
        else:
            print(f"UNDEFINED ({len(values)} positive values)")
            candidates.append({
                "name": f"Negative-mass fermion |m|/T={m}",
                "category": "Exotic Mass",
                "physics": "Negative energy fermionic occupation",
                "status": "UNDEFINED",
                "delta_b": None, "epsilon_d": None, "n": len(values),
            })

    print()

    # ── 4. Phantom energy ────────────────────────────────────────────────
    print("─── PHANTOM ENERGY (w < -1) ───")
    print("  Wrong-sign kinetic term: E_k = -ω_k,  n = 1/(e^{-ω/T} - 1)")
    print()

    print("  Computing...", end=" ", flush=True)
    values = generate_phantom_energy()
    result = analyze_values(values)
    if result:
        shape = classify_shape(result["epsilon_d"])
        print(f"δ_B={result['delta_b']:.6f}  shape={shape}")
        candidates.append({
            "name": "Phantom energy",
            "category": "Dark Energy",
            "physics": "Wrong-sign kinetic term (w < -1)",
            "status": "EXISTS",
            **result,
            "shape": shape,
        })
    else:
        print(f"UNDEFINED ({len(values)} positive values out of 100000)")
        candidates.append({
            "name": "Phantom energy",
            "category": "Dark Energy",
            "physics": "Wrong-sign kinetic term (w < -1)",
            "status": "UNDEFINED",
            "delta_b": None, "epsilon_d": None, "n": len(values),
            "reason": "e^{-ω/T} < 1 for ω > 0 → denominator < 0 → all occupation numbers negative",
        })

    print()

    # ── 5. Hawking radiation ─────────────────────────────────────────────
    print("─── HAWKING RADIATION (greybody-modified) ───")
    print("  N(ω) ∝ ω² · Γ(ω) / (e^ω - 1),  Γ = 1 - e^{-(ω/ω_c)²}")
    print()

    for omega_c in [0.5, 1.0, 2.0, 5.0]:
        print(f"  ω_c = {omega_c}...", end=" ", flush=True)
        values = generate_hawking_greybody(omega_c)
        result = analyze_values(values)
        if result:
            shape = classify_shape(result["epsilon_d"])
            nearest, dist = find_nearest_reference(result["epsilon_d"])
            print(f"δ_B={result['delta_b']:.6f}  shape={shape}  nearest={nearest} (d={dist:.4f})")
            candidates.append({
                "name": f"Hawking radiation ω_c={omega_c}",
                "category": "Black Holes",
                "physics": f"Greybody-modified Planck (ω_c={omega_c})",
                "status": "EXISTS",
                **result,
                "shape": shape,
                "nearest_known": nearest,
                "nearest_distance": round(dist, 6),
            })
        else:
            print("UNDEFINED")
            candidates.append({
                "name": f"Hawking radiation ω_c={omega_c}",
                "category": "Black Holes",
                "physics": f"Greybody-modified Planck",
                "status": "UNDEFINED",
                "delta_b": None, "epsilon_d": None, "n": len(values),
            })

    print()

    # ── 6. Unruh radiation ───────────────────────────────────────────────
    print("─── UNRUH RADIATION ───")
    print("  Pure Planck spectrum (accelerating observer thermal bath)")
    print()

    print("  Computing...", end=" ", flush=True)
    values = generate_unruh()
    result = analyze_values(values)
    if result:
        shape = classify_shape(result["epsilon_d"])
        nearest, dist = find_nearest_reference(result["epsilon_d"])
        print(f"δ_B={result['delta_b']:.6f}  shape={shape}  nearest={nearest} (d={dist:.6f})")
        candidates.append({
            "name": "Unruh radiation",
            "category": "Relativistic QFT",
            "physics": "Accelerating-observer thermal bath (= Planck)",
            "status": "EXISTS",
            **result,
            "shape": shape,
            "nearest_known": nearest,
            "nearest_distance": round(dist, 6),
        })

    print()

    # ── 7. Gravitons ─────────────────────────────────────────────────────
    print("─── THERMAL GRAVITONS ───")
    print("  Massless spin-2: same spectrum as photons (Planck)")
    print()

    print("  Computing...", end=" ", flush=True)
    values = generate_graviton_thermal()
    result = analyze_values(values)
    if result:
        shape = classify_shape(result["epsilon_d"])
        nearest, dist = find_nearest_reference(result["epsilon_d"])
        print(f"δ_B={result['delta_b']:.6f}  shape={shape}  nearest={nearest} (d={dist:.6f})")
        candidates.append({
            "name": "Gravitons (thermal)",
            "category": "Quantum Gravity",
            "physics": "Massless spin-2 thermal bosons",
            "status": "EXISTS",
            **result,
            "shape": shape,
            "nearest_known": nearest,
            "nearest_distance": round(dist, 6),
        })

    print()

    # ── 8. Majorana fermions ─────────────────────────────────────────────
    print("─── MAJORANA FERMIONS ───")
    print("  Self-conjugate: own antiparticle, but same FD occupation per mode")
    print()

    print("  Computing...", end=" ", flush=True)
    values = generate_majorana()
    result = analyze_values(values)
    if result:
        shape = classify_shape(result["epsilon_d"])
        nearest, dist = find_nearest_reference(result["epsilon_d"])
        print(f"δ_B={result['delta_b']:.6f}  shape={shape}  nearest={nearest} (d={dist:.6f})")
        candidates.append({
            "name": "Majorana fermion",
            "category": "BSM Particles",
            "physics": "Self-conjugate fermion (= FD per mode)",
            "status": "EXISTS",
            **result,
            "shape": shape,
            "nearest_known": nearest,
            "nearest_distance": round(dist, 6),
        })

    print()

    # ── 9. Axion-like boson ──────────────────────────────────────────────
    print("─── AXION-LIKE BOSON ───")
    print("  Ultralight boson: m/T << 1, nearly massless BE")
    print()

    for m_over_T in [0.001, 0.01, 0.1]:
        print(f"  m/T = {m_over_T}...", end=" ", flush=True)
        values = generate_axion(m_over_T)
        result = analyze_values(values)
        if result:
            shape = classify_shape(result["epsilon_d"])
            nearest, dist = find_nearest_reference(result["epsilon_d"])
            print(f"δ_B={result['delta_b']:.6f}  shape={shape}  nearest={nearest} (d={dist:.4f})")
            candidates.append({
                "name": f"Axion m/T={m_over_T}",
                "category": "Dark Matter",
                "physics": f"Ultralight boson m/T={m_over_T}",
                "status": "EXISTS",
                **result,
                "shape": shape,
                "nearest_known": nearest,
                "nearest_distance": round(dist, 6),
            })

    print()

    # ── 10. Sterile neutrino ─────────────────────────────────────────────
    print("─── STERILE NEUTRINO (Dodelson-Widrow) ───")
    print("  Non-thermal production: f(p) ∝ p / (e^p + 1)")
    print()

    print("  Computing...", end=" ", flush=True)
    values = generate_sterile_neutrino()
    result = analyze_values(values)
    if result:
        shape = classify_shape(result["epsilon_d"])
        nearest, dist = find_nearest_reference(result["epsilon_d"])
        print(f"δ_B={result['delta_b']:.6f}  shape={shape}  nearest={nearest} (d={dist:.4f})")
        candidates.append({
            "name": "Sterile neutrino (DW)",
            "category": "Dark Matter",
            "physics": "Non-thermal FD with momentum weighting",
            "status": "EXISTS",
            **result,
            "shape": shape,
            "nearest_known": nearest,
            "nearest_distance": round(dist, 6),
        })

    print()

    # ══════════════════════════════════════════════════════════════════════
    # THE WHITEBOARD
    # ══════════════════════════════════════════════════════════════════════

    print()
    print("█" * 78)
    print("█" + " " * 76 + "█")
    print("█" + "           T H E   W H I T E B O A R D".center(76) + "█")
    print("█" + "       Exotic Physics Existence Filter Results".center(76) + "█")
    print("█" + " " * 76 + "█")
    print("█" * 78)
    print()

    # Group by status
    exists = [c for c in candidates if c["status"] == "EXISTS"]
    undefined = [c for c in candidates if c["status"] == "UNDEFINED"]

    # ── EXISTS ───────────────────────────────────────────────────────────
    print(f"  EXIST ({len(exists)} candidates)")
    print(f"  {'─'*74}")
    print(f"  {'Candidate':<35s}  {'δ_B':>8s}  {'Shape':<12s}  {'Nearest':>8s}  {'Dist':>8s}")
    print(f"  {'─'*35}  {'─'*8}  {'─'*12}  {'─'*8}  {'─'*8}")

    for c in sorted(exists, key=lambda x: x["delta_b"]):
        nearest = c.get("nearest_known", "—")
        dist = c.get("nearest_distance", 0)
        shape = c.get("shape", "—")
        match_marker = ""
        if dist < 0.001:
            match_marker = " ≡"  # exact match to known
        elif dist < 0.005:
            match_marker = " ≈"  # close match
        print(f"  {c['name']:<35s}  {c['delta_b']:8.4f}  {shape:<12s}  "
              f"{nearest:>8s}  {dist:8.4f}{match_marker}")

    print()

    # ── UNDEFINED ────────────────────────────────────────────────────────
    print(f"  DO NOT EXIST ({len(undefined)} candidates)")
    print(f"  {'─'*74}")

    for c in undefined:
        reason = c.get("reason", "Insufficient positive-valued modes")
        print(f"  ✗ {c['name']}")
        print(f"    {reason}")

    print()

    # ── Key findings ─────────────────────────────────────────────────────
    print("─" * 78)
    print("KEY FINDINGS")
    print("─" * 78)

    # Find the anyons
    anyon_entries = [c for c in exists if "Anyon" in c["name"]]
    if anyon_entries:
        print()
        print("  ANYONS: Fractional exclusion statistics produce a smooth")
        print("  interpolation between BE and FD fingerprints:")
        for a in anyon_entries:
            nearest = a.get("nearest_known", "—")
            print(f"    {a['name']:<30s}  δ_B={a['delta_b']:.4f}  nearest={nearest}")

    # Negative mass
    neg_boson = [c for c in candidates if "Negative-mass boson" in c["name"]]
    neg_fermion = [c for c in exists if "Negative-mass fermion" in c["name"]]
    if neg_boson and all(c["status"] == "UNDEFINED" for c in neg_boson):
        print()
        print("  NEGATIVE-MASS BOSONS: All UNDEFINED.")
        print("  Negative energy → e^{E/T} < 1 → denominator (e^{E/T} - 1) < 0")
        print("  → all occupation numbers negative. Cannot exist thermally.")

    if neg_fermion:
        print()
        print("  NEGATIVE-MASS FERMIONS: Computable! (Pauli exclusion saves them)")
        print("  n = 1/(e^{-|E|/T} + 1) = 1 - n_FD(|E|): inverted Fermi-Dirac.")
        for f in neg_fermion:
            print(f"    {f['name']:<35s}  δ_B={f['delta_b']:.4f}")

    phantom = [c for c in candidates if "Phantom" in c["name"]]
    if phantom and phantom[0]["status"] == "UNDEFINED":
        print()
        print("  PHANTOM ENERGY: UNDEFINED.")
        print("  Same mechanism as negative-mass bosons — wrong-sign kinetic term")
        print("  makes all occupation numbers negative. Thermodynamically unviable.")

    # Hawking / Graviton / Unruh verification
    verification = [c for c in exists if any(x in c["name"] for x in ["Unruh", "Graviton", "Majorana"])]
    if verification:
        print()
        print("  VERIFICATION (should match known distributions):")
        for v in verification:
            nearest = v.get("nearest_known", "—")
            dist = v.get("nearest_distance", 0)
            marker = "EXACT MATCH" if dist < 0.001 else f"close (d={dist:.4f})"
            print(f"    {v['name']:<30s}  → {nearest}: {marker}")

    hawking = [c for c in exists if "Hawking" in c["name"]]
    if hawking:
        print()
        print("  HAWKING RADIATION (greybody-modified):")
        print("  Departs from pure Planck as greybody cutoff ω_c decreases:")
        for h in hawking:
            nearest = h.get("nearest_known", "—")
            dist = h.get("nearest_distance", 0)
            print(f"    {h['name']:<35s}  δ_B={h['delta_b']:.4f}  nearest={nearest} (d={dist:.4f})")

    print()

    # ── Save ─────────────────────────────────────────────────────────────
    output = {
        "experiment": "whiteboard",
        "description": "Exotic physics existence filter — all candidates tested",
        "candidates": candidates,
        "summary": {
            "total_tested": len(candidates),
            "exist": len(exists),
            "undefined": len(undefined),
            "exist_names": [c["name"] for c in exists],
            "undefined_names": [c["name"] for c in undefined],
        },
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = "results/round_trip/whiteboard.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Results saved to {out_path}")

    return output


if __name__ == "__main__":
    run_whiteboard()
