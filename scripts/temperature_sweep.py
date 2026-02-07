#!/usr/bin/env python3
"""Temperature Sweep Analysis for Quantum Statistical Distributions.

Evaluates BE, FD, and MB at many temperatures, computing delta_B at each.
Validates three quantitative predictions from Riner (2026):
  Prediction 1: FD deviation period = 1 in log10(T)
  Prediction 2: FD amplitude = |eta(s)| = 1.054x MB baseline
  Prediction 3: delta_B(T) ~ delta_max * |cos(2*pi*log10(T) + phi)|

Outputs results/temperature_sweep.json for the dashboard.
"""

import json
import os
import sys
import numpy as np
from datetime import datetime, timezone

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benford_core import (
    first_digit, observed_distribution, euclidean_deviation,
    per_digit_deviation, BENFORD_EXPECTED
)


def compute_delta_b_at_temperature(dist_type, T, n_points=50000):
    """Compute delta_B for a quantum distribution at temperature T.

    dist_type: 'BE', 'FD', or 'MB'
    T: temperature in natural units (kT)
    n_points: number of energy grid points
    """
    # Fixed uniform x-grid (dimensionless energy).
    # At temperature T, the PHYSICAL quantity is T * n(x) — the thermal
    # energy contribution. Multiplying by T shifts all significands by
    # log10(T), producing the period-1 oscillation in log10(T) that the
    # paper predicts from the Fourier harmonic e^(2*pi*i*log10(T)).
    x = np.linspace(0.001, 50, n_points)

    if dist_type == "BE":
        n = 1.0 / (np.exp(np.minimum(x, 500)) - 1.0)
    elif dist_type == "FD":
        n = 1.0 / (np.exp(np.minimum(x, 500)) + 1.0)
    elif dist_type == "MB":
        n = np.exp(-np.minimum(x, 500))
    else:
        raise ValueError(f"Unknown dist_type: {dist_type}")

    # Scale by T: physical observable = kT * n(eps/kT)
    vals_raw = T * n

    # Filter to positive finite values
    mask = (vals_raw > 0) & np.isfinite(vals_raw)
    vals = vals_raw[mask]

    if len(vals) < 50:
        return None, None, None

    # Extract first digits
    digits = []
    for v in vals:
        d = first_digit(v)
        if d is not None:
            digits.append(d)

    if len(digits) < 50:
        return None, None, None

    obs = observed_distribution(digits)
    delta_b = euclidean_deviation(obs)
    eps_d = per_digit_deviation(obs)

    return delta_b, eps_d, len(digits)


def run_sweep():
    """Run the full temperature sweep and save results."""
    # Temperature range: 8 decades
    log_T_values = np.linspace(-3, 5, 200)
    T_values = 10.0 ** log_T_values

    results = {
        "analysis": "temperature_sweep",
        "description": "delta_B as a function of temperature for BE, FD, MB distributions",
        "paper_predictions": {
            "prediction_1": "FD deviation oscillates with period 1 in log10(T)",
            "prediction_2": "FD amplitude |eta(s)| = 1.054x MB baseline",
            "prediction_3": "delta_B(T) ~ delta_max * |cos(2*pi*log10(T) + phi)|"
        },
        "dirichlet_eta_values": {
            "zeta_s_magnitude": 0.4214,
            "exclusion_factor": 2.5021,
            "eta_s_magnitude": 1.0545,
            "s": "2*pi*i / ln(10)"
        },
        "temperatures": [],
        "log10_T": [],
        "BE": {"delta_b": [], "label": "Bose-Einstein (all + coefficients)"},
        "FD": {"delta_b": [], "label": "Fermi-Dirac (alternating coefficients)"},
        "MB": {"delta_b": [], "label": "Maxwell-Boltzmann (single exponential)"},
        "FD_over_MB_ratio": [],
        "timestamp": datetime.now(timezone.utc).isoformat()
    }

    print(f"Running temperature sweep across {len(T_values)} temperatures...")
    print(f"  log10(T) range: [{log_T_values[0]:.1f}, {log_T_values[-1]:.1f}]")
    print()

    for i, T in enumerate(T_values):
        log_t = log_T_values[i]
        if i % 20 == 0:
            print(f"  [{i+1}/{len(T_values)}] T = {T:.3e} (log10(T) = {log_t:.2f})")

        results["temperatures"].append(float(T))
        results["log10_T"].append(float(log_t))

        for dist in ["BE", "FD", "MB"]:
            db, _, _ = compute_delta_b_at_temperature(dist, T, n_points=30000)
            if db is None:
                db = float("nan")
            results[dist]["delta_b"].append(float(db))

        # FD/MB ratio
        fd_db = results["FD"]["delta_b"][-1]
        mb_db = results["MB"]["delta_b"][-1]
        if mb_db > 0 and not (np.isnan(fd_db) or np.isnan(mb_db)):
            results["FD_over_MB_ratio"].append(float(fd_db / mb_db))
        else:
            results["FD_over_MB_ratio"].append(float("nan"))

    # Compute summary statistics
    fd_vals = [v for v in results["FD"]["delta_b"] if not np.isnan(v)]
    mb_vals = [v for v in results["MB"]["delta_b"] if not np.isnan(v)]
    be_vals = [v for v in results["BE"]["delta_b"] if not np.isnan(v)]
    ratio_vals = [v for v in results["FD_over_MB_ratio"] if not np.isnan(v)]

    results["summary"] = {
        "BE_mean_delta_b": float(np.mean(be_vals)) if be_vals else None,
        "BE_max_delta_b": float(np.max(be_vals)) if be_vals else None,
        "FD_mean_delta_b": float(np.mean(fd_vals)) if fd_vals else None,
        "FD_max_delta_b": float(np.max(fd_vals)) if fd_vals else None,
        "MB_mean_delta_b": float(np.mean(mb_vals)) if mb_vals else None,
        "MB_max_delta_b": float(np.max(mb_vals)) if mb_vals else None,
        "mean_FD_MB_ratio": float(np.mean(ratio_vals)) if ratio_vals else None,
        "predicted_FD_MB_ratio": 1.054,
        "total_temperatures": len(T_values)
    }

    # Save
    out_path = "results/temperature_sweep.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")
    print(f"\nSummary:")
    print(f"  BE  mean delta_B: {results['summary']['BE_mean_delta_b']:.6f}")
    print(f"  FD  mean delta_B: {results['summary']['FD_mean_delta_b']:.6f}")
    print(f"  MB  mean delta_B: {results['summary']['MB_mean_delta_b']:.6f}")
    print(f"  FD/MB ratio:      {results['summary']['mean_FD_MB_ratio']:.4f}  (predicted: 1.054)")


if __name__ == "__main__":
    run_sweep()
