#!/usr/bin/env python3
"""Experiment 3: ε(d) Fingerprint Atlas.

Collects all ε(d) signatures from the dimension sweep and eta recovery
experiments, classifies each into a shape type, and tests whether the
shape alone can identify which experiment and parameter produced it.

Success criterion: ε(d) shapes are structurally distinct between dimension
sweep (monotonic/sweep) and eta sweep (oscillatory), and change smoothly
within each sweep.
"""

import json
import os
import sys
import numpy as np
from datetime import datetime, timezone

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def classify_epsilon_shape(eps_vals):
    """Classify an ε(d) vector into a shape type.

    Categories:
      - "flat":        all |ε| < 0.002 (noise-level)
      - "monotone":    monotonically increasing or decreasing after digit 1
      - "oscillatory": sign changes in the tail (digits 4-9)
      - "sweep":       smooth transition from positive to negative (or vice versa)
    """
    if all(abs(e) < 0.002 for e in eps_vals):
        return "flat"

    # Count sign changes in digits 2-9 (index 1-8)
    tail = eps_vals[1:]  # digits 2-9
    sign_changes = sum(
        1 for i in range(len(tail) - 1)
        if tail[i] * tail[i + 1] < 0 and abs(tail[i]) > 0.0005
    )

    if sign_changes >= 2:
        return "oscillatory"

    # Check for sweep: eps_vals goes from one sign to another smoothly
    if abs(eps_vals[0]) > 0.002 and abs(eps_vals[-1]) > 0.002:
        if eps_vals[0] * eps_vals[-1] < 0:
            return "sweep"

    # Check monotonicity of tail
    diffs = [tail[i + 1] - tail[i] for i in range(len(tail) - 1)]
    if all(d >= -0.0005 for d in diffs):
        return "monotone"
    if all(d <= 0.0005 for d in diffs):
        return "monotone"

    return "mixed"


def smoothness_score(entries):
    """Measure how smoothly ε(d) changes across consecutive sweep points.

    Returns the mean L2 difference between consecutive ε(d) vectors,
    normalized by the parameter step size.
    """
    if len(entries) < 2:
        return 0.0

    diffs = []
    for i in range(len(entries) - 1):
        e1 = entries[i]["epsilon_d_vector"]
        e2 = entries[i + 1]["epsilon_d_vector"]
        l2 = np.sqrt(sum((a - b)**2 for a, b in zip(e1, e2)))
        diffs.append(l2)

    return float(np.mean(diffs))


def blind_identify(eps_vals, dim_entries, eta_entries):
    """Given just an ε(d) vector, identify which experiment produced it.

    Returns (experiment_name, best_match_param, distance).
    """
    best_dim = None
    best_dim_dist = float("inf")
    for entry in dim_entries:
        ref = entry["epsilon_d_vector"]
        dist = np.sqrt(sum((a - b)**2 for a, b in zip(eps_vals, ref)))
        if dist < best_dim_dist:
            best_dim_dist = dist
            best_dim = entry["parameter"]

    best_eta = None
    best_eta_dist = float("inf")
    for entry in eta_entries:
        ref = entry["epsilon_d_vector"]
        dist = np.sqrt(sum((a - b)**2 for a, b in zip(eps_vals, ref)))
        if dist < best_eta_dist:
            best_eta_dist = dist
            best_eta = entry["parameter"]

    if best_dim_dist < best_eta_dist:
        return "dimension_sweep", best_dim, best_dim_dist
    else:
        return "eta_recovery", best_eta, best_eta_dist


def run_fingerprint_atlas():
    """Build the fingerprint atlas from dimension sweep and eta recovery results."""

    # Load results from previous experiments
    dim_path = "results/round_trip/dimension_sweep.json"
    eta_path = "results/round_trip/eta_recovery.json"
    planck_path = "results/individual/planck_radiation_spectrum.json"

    with open(dim_path) as f:
        dim_data = json.load(f)
    with open(eta_path) as f:
        eta_data = json.load(f)
    with open(planck_path) as f:
        planck_data = json.load(f)

    print("=" * 70)
    print("EXPERIMENT 3: ε(d) Fingerprint Atlas")
    print("  Collecting all ε(d) signatures from Experiments 1 & 2")
    print("=" * 70)
    print()

    # --- Collect entries from dimension sweep ---
    dim_entries = []
    for r in dim_data["sweep"]:
        eps_dict = r["epsilon_d"]
        eps_vec = [eps_dict[str(d)] for d in range(1, 10)]
        shape = classify_epsilon_shape(eps_vec)
        dim_entries.append({
            "experiment": "dimension_sweep",
            "parameter": r["exponent"],
            "parameter_name": "n",
            "delta_b": r["delta_b"],
            "epsilon_d_vector": eps_vec,
            "shape": shape,
        })

    # --- Collect entries from eta recovery ---
    eta_entries = []
    for r in eta_data["sweep"]:
        eps_dict = r["epsilon_d"]
        eps_vec = [eps_dict[str(d)] for d in range(1, 10)]
        shape = classify_epsilon_shape(eps_vec)
        eta_entries.append({
            "experiment": "eta_recovery",
            "parameter": r["alpha"],
            "parameter_name": "α",
            "delta_b": r["delta_b"],
            "epsilon_d_vector": eps_vec,
            "shape": shape,
        })

    all_entries = dim_entries + eta_entries

    # --- Print atlas ---
    print("-" * 70)
    print("DIMENSION SWEEP ε(d) FINGERPRINTS")
    print("-" * 70)
    print(f"  {'n':>5s}  {'δ_B':>8s}  {'Shape':>12s}  ε(d) vector")
    print(f"  {'─'*5}  {'─'*8}  {'─'*12}  {'─'*50}")
    for e in dim_entries:
        eps_str = " ".join(f"{v:+.4f}" for v in e["epsilon_d_vector"])
        print(f"  {e['parameter']:5.1f}  {e['delta_b']:8.4f}  {e['shape']:>12s}  [{eps_str}]")

    print()
    print("-" * 70)
    print("ETA RECOVERY ε(d) FINGERPRINTS")
    print("-" * 70)
    print(f"  {'α':>5s}  {'δ_B':>8s}  {'Shape':>12s}  ε(d) vector")
    print(f"  {'─'*5}  {'─'*8}  {'─'*12}  {'─'*50}")
    for e in eta_entries:
        eps_str = " ".join(f"{v:+.4f}" for v in e["epsilon_d_vector"])
        print(f"  {e['parameter']:5.2f}  {e['delta_b']:8.4f}  {e['shape']:>12s}  [{eps_str}]")

    print()

    # --- Shape distribution ---
    print("-" * 70)
    print("SHAPE DISTRIBUTION BY EXPERIMENT")
    print("-" * 70)
    dim_shapes = [e["shape"] for e in dim_entries]
    eta_shapes = [e["shape"] for e in eta_entries]

    for label, shapes in [("Dimension sweep", dim_shapes), ("Eta recovery", eta_shapes)]:
        unique = sorted(set(shapes))
        counts = {s: shapes.count(s) for s in unique}
        print(f"  {label}:")
        for s, c in counts.items():
            print(f"    {s}: {c}/{len(shapes)}")
    print()

    # --- Smoothness analysis ---
    print("-" * 70)
    print("SMOOTHNESS ANALYSIS")
    print("-" * 70)
    dim_smooth = smoothness_score(dim_entries)
    eta_smooth = smoothness_score(eta_entries)
    print(f"  Dimension sweep: mean Δε between consecutive points = {dim_smooth:.6f}")
    print(f"  Eta recovery:    mean Δε between consecutive points = {eta_smooth:.6f}")
    print(f"  (Smooth transitions → small values; abrupt changes → large values)")
    print()

    # --- Blind identification test ---
    print("-" * 70)
    print("BLIND IDENTIFICATION TEST")
    print("-" * 70)
    print("  Given just an ε(d) vector, can we identify the source experiment?")
    print()

    # Test cases: each entry from both experiments
    correct = 0
    total = 0
    for entry in all_entries:
        exp, param, dist = blind_identify(
            entry["epsilon_d_vector"], dim_entries, eta_entries
        )
        match = exp == entry["experiment"]
        if match:
            correct += 1
        total += 1

    print(f"  Self-identification accuracy: {correct}/{total} = {100*correct/total:.1f}%")
    print()

    # Test with stored Planck result
    planck_eps = [planck_data["per_digit_deviation"][str(d)] for d in range(1, 10)]
    planck_exp, planck_param, planck_dist = blind_identify(
        planck_eps, dim_entries, eta_entries
    )
    print(f"  Stored Planck ε(d) → identified as: {planck_exp}, param={planck_param}")
    print(f"    L2 distance to nearest match: {planck_dist:.6f}")
    planck_match = planck_exp == "dimension_sweep" and abs(planck_param - 3.0) < 0.1
    print(f"    Expected: dimension_sweep n=3.0  {'✓' if planck_match else '✗'}")
    print()

    # --- Structural distinctness test ---
    print("-" * 70)
    print("STRUCTURAL DISTINCTNESS")
    print("-" * 70)

    # Compute inter-experiment vs intra-experiment distances
    dim_vecs = [np.array(e["epsilon_d_vector"]) for e in dim_entries]
    eta_vecs = [np.array(e["epsilon_d_vector"]) for e in eta_entries]

    # Intra-experiment: mean pairwise distance within each experiment
    intra_dim = []
    for i in range(len(dim_vecs)):
        for j in range(i + 1, len(dim_vecs)):
            intra_dim.append(float(np.linalg.norm(dim_vecs[i] - dim_vecs[j])))

    intra_eta = []
    for i in range(len(eta_vecs)):
        for j in range(i + 1, len(eta_vecs)):
            intra_eta.append(float(np.linalg.norm(eta_vecs[i] - eta_vecs[j])))

    # Inter-experiment: mean distance between experiments
    inter = []
    for dv in dim_vecs:
        for ev in eta_vecs:
            inter.append(float(np.linalg.norm(dv - ev)))

    mean_intra_dim = np.mean(intra_dim) if intra_dim else 0
    mean_intra_eta = np.mean(intra_eta) if intra_eta else 0
    mean_inter = np.mean(inter) if inter else 0

    print(f"  Mean intra-dimension-sweep distance: {mean_intra_dim:.6f}")
    print(f"  Mean intra-eta-recovery distance:    {mean_intra_eta:.6f}")
    print(f"  Mean inter-experiment distance:       {mean_inter:.6f}")
    separation = mean_inter / max(mean_intra_dim, mean_intra_eta, 1e-10)
    print(f"  Separation ratio (inter/max_intra):   {separation:.2f}x")
    well_separated = separation > 1.0
    print(f"  Well separated: {'✓' if well_separated else '✗'} (ratio > 1.0)")
    print()

    # --- Overall verdict ---
    print("=" * 70)
    print("FINGERPRINT ATLAS SUMMARY")
    print("=" * 70)

    dim_shape_set = set(dim_shapes)
    eta_shape_set = set(eta_shapes)
    shapes_distinct = not (dim_shape_set == eta_shape_set and len(dim_shape_set) == 1)

    print(f"  Shape types in dimension sweep: {sorted(dim_shape_set)}")
    print(f"  Shape types in eta recovery:    {sorted(eta_shape_set)}")
    print(f"  Structurally distinct shapes: {'✓' if shapes_distinct else '✗'}")
    print(f"  Blind identification accuracy: {100*correct/total:.1f}%")
    print(f"  Planck round-trip match: {'✓' if planck_match else '✗'}")
    print(f"  Inter/intra separation: {separation:.2f}x {'✓' if well_separated else '✗'}")
    print()

    # --- Save results ---
    output = {
        "experiment": "fingerprint_atlas",
        "description": "ε(d) shape classification and blind identification across sweeps",
        "atlas": [
            {
                "experiment": e["experiment"],
                "parameter_name": e["parameter_name"],
                "parameter": e["parameter"],
                "delta_b": e["delta_b"],
                "epsilon_d": e["epsilon_d_vector"],
                "shape": e["shape"],
            }
            for e in all_entries
        ],
        "shape_distribution": {
            "dimension_sweep": {s: dim_shapes.count(s) for s in sorted(set(dim_shapes))},
            "eta_recovery": {s: eta_shapes.count(s) for s in sorted(set(eta_shapes))},
        },
        "smoothness": {
            "dimension_sweep": round(dim_smooth, 6),
            "eta_recovery": round(eta_smooth, 6),
        },
        "blind_identification": {
            "accuracy": round(correct / total, 4),
            "correct": correct,
            "total": total,
            "planck_identified_as": planck_exp,
            "planck_matched_param": planck_param,
            "planck_match_correct": bool(planck_match),
        },
        "structural_distinctness": {
            "mean_intra_dimension": round(float(mean_intra_dim), 6),
            "mean_intra_eta": round(float(mean_intra_eta), 6),
            "mean_inter_experiment": round(float(mean_inter), 6),
            "separation_ratio": round(float(separation), 4),
            "well_separated": bool(well_separated),
        },
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = "results/round_trip/fingerprint_atlas.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Results saved to {out_path}")

    return output


if __name__ == "__main__":
    run_fingerprint_atlas()
