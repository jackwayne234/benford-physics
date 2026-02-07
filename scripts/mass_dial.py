#!/usr/bin/env python3
"""Experiment 4: Mass Dial — δ_B across the mass spectrum, including tachyons.

Uses the relativistic dispersion relation E = √(k² + m²) in natural units
(T=1) to generate bosonic occupation numbers:

    n(k) = 1 / (e^{√(k² + m²)} - 1)

Three regimes:
    m² > 0  :  Normal massive particle (real mass)
    m² = 0  :  Massless boson (photon-like, pure BE)
    m² < 0  :  Tachyonic (imaginary mass, E only real for k > |m|)

For tachyons, low-momentum modes have imaginary energy — the occupation
number is undefined there.  As |m²| grows, fewer valid modes survive.
If δ_B becomes uncomputable (too few valid points), the distribution
is physically unviable.

The experiment:
  1. Sweep m²/T² from -25 (deep tachyonic) to +25 (heavy)
  2. Measure δ_B, ε(d), and the fraction of surviving valid modes at each point
  3. Invert: pick random δ_B values and recover masses
  4. Report what happens in tachyonic territory
"""

import json
import math
import os
import sys
import random
import numpy as np
from datetime import datetime, timezone

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benford_core import run_full_analysis


def generate_relativistic_be(m_squared, n_points=100000):
    """Generate bosonic occupation numbers with relativistic dispersion.

    Parameters
    ----------
    m_squared : float
        m²/T² in natural units.  Positive = normal mass, negative = tachyonic.
    n_points : int
        Number of momentum grid points.

    Returns
    -------
    values : list of float
        Positive finite occupation numbers.
    n_valid : int
        Number of momentum modes with real energy.
    n_total : int
        Total momentum grid points.
    """
    k = np.linspace(0.001, 50, n_points)

    # E² = k² + m²
    E_squared = k**2 + m_squared

    # Only modes with E² > 0 have real energy
    valid_mask = E_squared > 0
    n_valid = int(np.sum(valid_mask))
    n_total = n_points

    if n_valid < 50:
        return [], n_valid, n_total

    k_valid = k[valid_mask]
    E = np.sqrt(E_squared[valid_mask])

    # Bosonic occupation: n(k) = 1/(e^E - 1)
    # Clip E to avoid overflow
    E_clipped = np.minimum(E, 500)
    occupation = 1.0 / (np.exp(E_clipped) - 1.0)

    # Filter to positive finite values
    out = []
    for v in occupation:
        fv = float(v)
        if fv > 0 and math.isfinite(fv):
            out.append(fv)

    return out, n_valid, n_total


def run_mass_dial():
    """Sweep m² from tachyonic to massive, measuring δ_B at each point."""

    # --- Define sweep points ---
    # Tachyonic region: m² < 0
    tachyonic_m2 = [-25, -20, -16, -12, -9, -6, -4, -2, -1, -0.5, -0.1]
    # Massless
    zero = [0]
    # Normal mass region: m² > 0
    normal_m2 = [0.1, 0.5, 1, 2, 4, 6, 9, 12, 16, 20, 25]

    all_m2 = tachyonic_m2 + zero + normal_m2
    sweep_results = []

    print("=" * 74)
    print("EXPERIMENT 4: Mass Dial")
    print("  Dispersion: E = sqrt(k² + m²),  Occupation: n = 1/(e^E - 1)")
    print("  Sweeping m²/T² from -25 (tachyonic) through 0 (massless) to +25 (massive)")
    print("=" * 74)
    print()
    print(f"  {'m²/T²':>8s}  {'regime':>10s}  {'valid%':>7s}  {'n_pts':>6s}  "
          f"{'δ_B':>10s}  {'MAD':>10s}  {'verdict':>10s}")
    print(f"  {'─'*8}  {'─'*10}  {'─'*7}  {'─'*6}  {'─'*10}  {'─'*10}  {'─'*10}")

    for m2 in all_m2:
        values, n_valid, n_total = generate_relativistic_be(m2)
        valid_frac = n_valid / n_total

        if m2 < 0:
            regime = "tachyonic"
        elif m2 == 0:
            regime = "massless"
        else:
            regime = "massive"

        entry = {
            "m_squared": m2,
            "regime": regime,
            "valid_fraction": round(valid_frac, 4),
            "n_valid_modes": n_valid,
            "n_total_modes": n_total,
        }

        if len(values) >= 50:
            analysis = run_full_analysis(values)
            fd = analysis["first_digit"]
            delta_b = fd["euclidean_deviation"]
            epsilon_d = fd["per_digit_deviation"]
            mad = fd["mad"]["mad"]

            entry["delta_b"] = delta_b
            entry["epsilon_d"] = {str(k): v for k, v in epsilon_d.items()}
            entry["mad"] = mad
            entry["n_analyzed"] = analysis["n"]
            entry["verdict"] = analysis["verdict"]
            entry["computable"] = True

            print(f"  {m2:8.1f}  {regime:>10s}  {valid_frac:6.1%}  "
                  f"{analysis['n']:6d}  {delta_b:10.6f}  {mad:10.6f}  "
                  f"{analysis['verdict']:>10s}")
        else:
            entry["delta_b"] = None
            entry["epsilon_d"] = None
            entry["mad"] = None
            entry["n_analyzed"] = len(values)
            entry["verdict"] = "UNDEFINED"
            entry["computable"] = False

            print(f"  {m2:8.1f}  {regime:>10s}  {valid_frac:6.1%}  "
                  f"{len(values):6d}  {'UNDEFINED':>10s}  {'---':>10s}  "
                  f"{'UNDEFINED':>10s}")

        sweep_results.append(entry)

    print()

    # --- Tachyonic analysis ---
    print("-" * 74)
    print("TACHYONIC TERRITORY ANALYSIS")
    print("-" * 74)

    tachyonic_entries = [r for r in sweep_results if r["regime"] == "tachyonic"]
    computable_tachyons = [r for r in tachyonic_entries if r["computable"]]
    undefined_tachyons = [r for r in tachyonic_entries if not r["computable"]]

    print(f"  Total tachyonic points tested: {len(tachyonic_entries)}")
    print(f"  Computable (enough valid modes): {len(computable_tachyons)}")
    print(f"  Undefined (too few valid modes): {len(undefined_tachyons)}")
    print()

    if undefined_tachyons:
        print("  Undefined at m² = " +
              ", ".join(str(r["m_squared"]) for r in undefined_tachyons))
        print("  → These distributions cannot produce enough real-valued modes")
        print("    to compute meaningful first-digit statistics.")
        print("  → In the Benford framework: they don't exist as physical distributions.")
        print()

    if computable_tachyons:
        print("  Computable tachyonic points:")
        for r in computable_tachyons:
            print(f"    m²={r['m_squared']:6.1f}  valid={r['valid_fraction']:5.1%}  "
                  f"δ_B={r['delta_b']:.6f}  verdict={r['verdict']}")
        print()
        print("  Note: Even when computable, tachyonic distributions are sampling")
        print("  only high-momentum modes (k > |m|). The missing low-k modes")
        print("  represent the tachyonic instability — physics that can't sustain itself.")

    print()

    # --- Phase transition boundary ---
    print("-" * 74)
    print("PHASE BOUNDARY: where does computability break down?")
    print("-" * 74)

    for r in tachyonic_entries:
        status = "computable" if r["computable"] else "UNDEFINED"
        bar = "█" * int(r["valid_fraction"] * 50) if r["valid_fraction"] > 0 else ""
        print(f"  m²={r['m_squared']:6.1f}  valid={r['valid_fraction']:5.1%}  "
              f"{bar}  {status}")

    print()

    # --- Normal mass calibration and random inversions ---
    print("-" * 74)
    print("NORMAL MASS CALIBRATION CURVE")
    print("-" * 74)

    # Collect computable entries with m² >= 0
    normal_entries = [r for r in sweep_results
                      if r["m_squared"] >= 0 and r["computable"]]

    m2_arr = np.array([r["m_squared"] for r in normal_entries])
    db_arr = np.array([r["delta_b"] for r in normal_entries])

    print(f"  {'m²/T²':>8s}  {'m/T':>8s}  {'δ_B':>10s}")
    print(f"  {'─'*8}  {'─'*8}  {'─'*10}")
    for r in normal_entries:
        m_over_t = math.sqrt(abs(r["m_squared"])) if r["m_squared"] >= 0 else 0
        print(f"  {r['m_squared']:8.1f}  {m_over_t:8.3f}  {r['delta_b']:10.6f}")

    print()

    # --- Random δ_B inversions ---
    print("-" * 74)
    print("RANDOM δ_B INVERSIONS → MASS")
    print("-" * 74)

    db_min = float(db_arr.min())
    db_max = float(db_arr.max())
    print(f"  Valid δ_B range for inversion: [{db_min:.6f}, {db_max:.6f}]")
    print()

    random.seed(42)
    random_dbs = sorted([round(random.uniform(db_min, db_max), 6) for _ in range(6)])

    inversions = []
    print(f"  {'δ_B (random)':>14s}  {'m²/T²':>10s}  {'m/T':>10s}  "
          f"{'At T=2.725K':>14s}  {'At T=150MeV':>14s}")
    print(f"  {'─'*14}  {'─'*10}  {'─'*10}  {'─'*14}  {'─'*14}")

    # Physical constants for mass conversion
    k_B = 8.617333e-5  # eV/K
    T_cmb = 2.725      # K (CMB temperature)
    T_qcd = 150e6      # eV (QCD transition ~150 MeV)

    for db in random_dbs:
        # Interpolate to find m²
        recovered_m2 = float(np.interp(db, db_arr, m2_arr))
        recovered_m = math.sqrt(abs(recovered_m2)) if recovered_m2 >= 0 else 0

        # Convert to physical masses
        # m_physical = (m/T) * k_B * T
        mass_cmb_eV = recovered_m * k_B * T_cmb
        mass_qcd_eV = recovered_m * T_qcd

        # Format mass strings
        if mass_cmb_eV < 1e-3:
            cmb_str = f"{mass_cmb_eV*1e6:.2f} μeV"
        elif mass_cmb_eV < 1:
            cmb_str = f"{mass_cmb_eV*1e3:.2f} meV"
        else:
            cmb_str = f"{mass_cmb_eV:.2f} eV"

        if mass_qcd_eV < 1e6:
            qcd_str = f"{mass_qcd_eV/1e3:.1f} keV"
        elif mass_qcd_eV < 1e9:
            qcd_str = f"{mass_qcd_eV/1e6:.1f} MeV"
        else:
            qcd_str = f"{mass_qcd_eV/1e9:.2f} GeV"

        print(f"  {db:14.6f}  {recovered_m2:10.4f}  {recovered_m:10.4f}  "
              f"{cmb_str:>14s}  {qcd_str:>14s}")

        inversions.append({
            "delta_b_input": db,
            "recovered_m_squared": round(recovered_m2, 6),
            "recovered_m_over_T": round(recovered_m, 6),
            "mass_at_CMB_temp_eV": float(mass_cmb_eV),
            "mass_at_QCD_temp_eV": float(mass_qcd_eV),
        })

    print()
    print("  (T_CMB = 2.725 K ≈ 0.235 meV;  T_QCD ≈ 150 MeV)")
    print()

    # --- Summary ---
    print("=" * 74)
    print("MASS DIAL SUMMARY")
    print("=" * 74)

    n_computable = sum(1 for r in sweep_results if r["computable"])
    n_undefined = sum(1 for r in sweep_results if not r["computable"])

    print(f"  Points swept: {len(sweep_results)}")
    print(f"  Computable: {n_computable}")
    print(f"  Undefined (tachyonic breakdown): {n_undefined}")
    print()
    print(f"  Massless (m²=0): δ_B = {[r['delta_b'] for r in sweep_results if r['m_squared'] == 0][0]:.6f}")
    print(f"  Heaviest tested (m²=25): δ_B = {[r['delta_b'] for r in sweep_results if r['m_squared'] == 25][0]:.6f}")
    print()

    if n_undefined > 0:
        print("  TACHYON VERDICT: Deep tachyonic distributions (m² << 0) produce")
        print("  UNDEFINED δ_B — not zero, not negative, simply not computable.")
        print("  The distribution has too few real-valued modes to form a")
        print("  meaningful statistical sample. In the Benford framework,")
        print("  these distributions do not exist as physical objects.")
    else:
        print("  All points were computable (no tachyonic breakdown in tested range).")

    print()

    # --- Save ---
    output = {
        "experiment": "mass_dial",
        "description": "δ_B across the mass spectrum from tachyonic (m²<0) to massive (m²>0)",
        "formula": "n(k) = 1 / (e^{sqrt(k² + m²)} - 1)",
        "sweep": sweep_results,
        "tachyonic_analysis": {
            "total_tested": len(tachyonic_entries),
            "computable": len(computable_tachyons),
            "undefined": len(undefined_tachyons),
            "undefined_at_m2": [r["m_squared"] for r in undefined_tachyons],
        },
        "inversions": inversions,
        "physical_scales": {
            "T_CMB_kelvin": T_cmb,
            "T_CMB_eV": k_B * T_cmb,
            "T_QCD_eV": T_qcd,
        },
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = "results/round_trip/mass_dial.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Results saved to {out_path}")

    return output


if __name__ == "__main__":
    run_mass_dial()
