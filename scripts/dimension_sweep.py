#!/usr/bin/env python3
"""Experiment 1: Dimensionality Sweep — δ_B as a calibration curve.

The Planck formula uses ν³ because space has 3 dimensions (density of states
∝ ν^{D-1}, energy per mode ∝ ν, so the spectral density prefactor ∝ ν^D = ν³).

This script sweeps the exponent n in B(ν) = ν^n / (e^ν - 1) and measures δ_B
at each n.  If δ_B forms a smooth calibration curve, we can invert it: given
the measured Planck δ_B ≈ 0.028, can we recover n = 3?

Success criterion: inverting δ_B = 0.028 recovers n = 3.0 ± 0.1
"""

import json
import os
import sys
import numpy as np
from datetime import datetime, timezone

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benford_core import run_full_analysis
from scripts.data_fetchers import generate_be_with_prefactor


def run_dimension_sweep():
    """Sweep n = 0 to 5 in steps of 0.5, measuring δ_B and ε(d) at each."""
    exponents = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
    sweep_results = []

    print("=" * 70)
    print("EXPERIMENT 1: Dimensionality Sweep")
    print("  Formula: B(x) = x^n / (e^x - 1)")
    print("  Sweeping n = 0, 0.5, 1, ..., 5")
    print("=" * 70)
    print()

    for n in exponents:
        values = generate_be_with_prefactor(n)
        analysis = run_full_analysis(values)

        fd = analysis["first_digit"]
        delta_b = fd["euclidean_deviation"]
        epsilon_d = fd["per_digit_deviation"]
        mad = fd["mad"]["mad"]
        chi_sq_p = fd["chi_squared"]["p_value"]

        entry = {
            "exponent": n,
            "delta_b": delta_b,
            "epsilon_d": {str(k): v for k, v in epsilon_d.items()},
            "mad": mad,
            "chi_squared_p": chi_sq_p,
            "n_values": analysis["n"],
            "verdict": analysis["verdict"],
        }
        sweep_results.append(entry)

        # Format ε(d) compactly
        eps_str = " ".join(f"{epsilon_d[d]:+.5f}" for d in range(1, 10))
        print(f"  n={n:4.1f}  δ_B={delta_b:.6f}  MAD={mad:.6f}  "
              f"χ²p={chi_sq_p:.4f}  {analysis['verdict']}")

    print()

    # --- Calibration curve summary ---
    print("-" * 70)
    print("CALIBRATION CURVE: n → δ_B")
    print("-" * 70)
    print(f"  {'n':>5s}  {'δ_B':>10s}  {'MAD':>10s}  {'Verdict':>10s}  ε(d) signature")
    print(f"  {'─'*5}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*40}")
    for r in sweep_results:
        eps = r["epsilon_d"]
        # Classify shape
        eps_vals = [eps[str(d)] for d in range(1, 10)]
        shape = classify_epsilon_shape(eps_vals)
        print(f"  {r['exponent']:5.1f}  {r['delta_b']:10.6f}  {r['mad']:10.6f}  "
              f"{r['verdict']:>10s}  {shape}")

    print()

    # --- Inverse: find n at δ_B = 0.028 ---
    target_delta_b = 0.028
    ns = np.array([r["exponent"] for r in sweep_results])
    dbs = np.array([r["delta_b"] for r in sweep_results])

    # Linear interpolation to find n where δ_B = target
    recovered_n = np.interp(target_delta_b, dbs, ns)

    # Also find n for the exact measured Planck δ_B
    planck_measured = 0.027921
    recovered_n_exact = np.interp(planck_measured, dbs, ns)

    print("-" * 70)
    print("ROUND-TRIP INVERSION")
    print("-" * 70)
    print(f"  Target δ_B = {target_delta_b:.6f}  →  recovered n = {recovered_n:.4f}")
    print(f"  Planck δ_B = {planck_measured:.6f}  →  recovered n = {recovered_n_exact:.4f}")
    print(f"  True value: n = 3 (3 spatial dimensions)")
    print(f"  Error: |n_recovered - 3| = {abs(recovered_n - 3):.4f}")
    print()

    success = abs(recovered_n - 3.0) < 0.1
    if success:
        print("  ✓ SUCCESS: δ_B = 0.028 inverts to n = 3.0 ± 0.1")
    else:
        print(f"  ✗ FAILED: recovered n = {recovered_n:.4f}, outside ±0.1 of 3.0")

    print()

    # --- Sanity checks ---
    print("-" * 70)
    print("SANITY CHECKS")
    print("-" * 70)
    be_delta = sweep_results[0]["delta_b"]  # n=0
    planck_delta = [r for r in sweep_results if r["exponent"] == 3][0]["delta_b"]

    be_check = abs(be_delta - 0.006) < 0.002
    planck_check = abs(planck_delta - 0.028) < 0.002

    print(f"  n=0 (pure BE):  δ_B = {be_delta:.6f}  "
          f"(expected ~0.006)  {'✓' if be_check else '✗'}")
    print(f"  n=3 (Planck):   δ_B = {planck_delta:.6f}  "
          f"(expected ~0.028)  {'✓' if planck_check else '✗'}")
    print()

    # --- Save results ---
    output = {
        "experiment": "dimension_sweep",
        "description": "δ_B as a function of density-of-states exponent n in x^n/(e^x-1)",
        "formula": "B(x) = x^n / (exp(x) - 1)",
        "sweep": sweep_results,
        "inversion": {
            "target_delta_b": target_delta_b,
            "recovered_exponent": round(float(recovered_n), 6),
            "planck_measured_delta_b": planck_measured,
            "recovered_exponent_exact": round(float(recovered_n_exact), 6),
            "true_exponent": 3,
            "error": round(float(abs(recovered_n - 3)), 6),
            "success": bool(success),
        },
        "sanity_checks": {
            "n0_matches_BE": bool(be_check),
            "n0_delta_b": be_delta,
            "n3_matches_Planck": bool(planck_check),
            "n3_delta_b": planck_delta,
        },
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = "results/round_trip/dimension_sweep.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Results saved to {out_path}")

    return output


def classify_epsilon_shape(eps_vals):
    """Classify an ε(d) vector into a shape type.

    Categories:
      - "flat":       all |ε| < 0.002 (noise-level)
      - "monotone":   monotonically increasing or decreasing after digit 1
      - "oscillatory": sign changes in the tail (digits 4-9)
      - "sweep":      smooth transition from positive to negative (or vice versa)
    """
    if all(abs(e) < 0.002 for e in eps_vals):
        return "flat (noise)"

    # Count sign changes in digits 2-9 (index 1-8)
    tail = eps_vals[1:]  # digits 2-9
    sign_changes = sum(
        1 for i in range(len(tail) - 1)
        if tail[i] * tail[i + 1] < 0 and abs(tail[i]) > 0.0005
    )

    if sign_changes >= 2:
        return "oscillatory"

    # Check for sweep: eps_vals goes from one sign to another smoothly
    if eps_vals[0] > 0.002 and eps_vals[-1] > 0.002:
        return "monotone (positive)"
    if eps_vals[0] < -0.002 and eps_vals[-1] < -0.002:
        return "monotone (negative)"
    if eps_vals[0] * eps_vals[-1] < 0:
        return "sweep"

    # Check monotonicity of tail
    diffs = [tail[i + 1] - tail[i] for i in range(len(tail) - 1)]
    if all(d >= -0.0005 for d in diffs):
        return "monotone (rising)"
    if all(d <= 0.0005 for d in diffs):
        return "monotone (falling)"

    return "mixed"


if __name__ == "__main__":
    run_dimension_sweep()
