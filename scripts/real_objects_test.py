#!/usr/bin/env python3
"""
Test the 4D Benford metric against real objects:
Earth, Jupiter, Neutron Star, Quasar (SMBH)

For each: compute metric components, floor status, CS deviation,
entropy rate at the object's surface and a few points approaching it.
"""

import math
import numpy as np

# Constants
G = 6.674e-11       # m^3 kg^-1 s^-2
c = 2.998e8          # m/s
M_sun = 1.989e30     # kg

# CS Equation
CS_A = 0.003389
CS_B = 0.002508
def delta_b(r_ratio):
    return CS_A * abs(math.log(r_ratio)) + CS_B

# Benford floor
probs = [math.log10(1 + 1/d) for d in range(1, 10)]
FLOOR_VAL = math.sqrt(sum(p**2 for p in probs))

def g_delta(r_ratio):
    db = delta_b(r_ratio)
    return math.log10(1 + 1/max(db, 1e-10))

def metric_4d(r_ratio):
    gtt = r_ratio**2
    gpp = r_ratio**2
    gd = g_delta(r_ratio)
    det = r_ratio**4 * gd
    floor_active = det < FLOOR_VAL
    if floor_active:
        grr = FLOOR_VAL / (gtt * gpp * gd)
    else:
        grr = 1.0
    return {
        'grr': grr, 'gtt': gtt, 'gpp': gpp, 'gd': gd,
        'det': grr * gtt * gpp * gd,
        'delta_b': delta_b(r_ratio),
        'floor_active': floor_active
    }

def r_schwarzschild(mass_kg):
    return 2 * G * mass_kg / (c**2)

# ── Objects ──
objects = [
    {
        'name': 'Earth',
        'mass_kg': 5.972e24,
        'radius_m': 6.371e6,
        'description': 'Rocky planet, weak gravity'
    },
    {
        'name': 'Jupiter',
        'mass_kg': 1.898e27,
        'radius_m': 6.9911e7,
        'description': 'Gas giant, moderate gravity'
    },
    {
        'name': 'Neutron Star (typical)',
        'mass_kg': 1.4 * M_sun,
        'radius_m': 10e3,  # 10 km
        'description': '1.4 solar masses, 10 km radius'
    },
    {
        'name': 'Neutron Star (massive)',
        'mass_kg': 2.1 * M_sun,
        'radius_m': 10e3,
        'description': '2.1 solar masses, near collapse'
    },
    {
        'name': 'Sagittarius A* (Milky Way SMBH)',
        'mass_kg': 4.0e6 * M_sun,
        'radius_m': None,  # it's a black hole
        'description': '4 million solar masses, our galactic center'
    },
    {
        'name': 'TON 618 (Quasar)',
        'mass_kg': 66e9 * M_sun,
        'radius_m': None,  # it's a black hole
        'description': '66 billion solar masses, one of the largest known'
    },
]

print("=" * 70)
print("4D Benford Metric — Real Objects Test")
print(f"CS Equation: δ_B(r) = {CS_A} × |ln(r)| + {CS_B}")
print(f"Benford Floor: {FLOOR_VAL:.6f}")
print("=" * 70)

for obj in objects:
    r_s = r_schwarzschild(obj['mass_kg'])
    print(f"\n{'─' * 70}")
    print(f"  {obj['name']}")
    print(f"  {obj['description']}")
    print(f"{'─' * 70}")
    print(f"  Mass: {obj['mass_kg']:.3e} kg")
    print(f"  Schwarzschild radius: {r_s:.4e} m", end='')

    if r_s < 1:
        print(f" ({r_s*1000:.4f} mm)")
    elif r_s < 1000:
        print(f" ({r_s:.2f} m)")
    elif r_s < 1e6:
        print(f" ({r_s/1000:.2f} km)")
    else:
        print(f" ({r_s/1e9:.2f} billion m)")

    is_bh = obj['radius_m'] is None

    if not is_bh:
        surface_r_ratio = obj['radius_m'] / r_s
        print(f"  Physical radius: {obj['radius_m']:.4e} m")
        print(f"  Surface at r/r_s = {surface_r_ratio:.2f}")

        # Test at surface and a few distances
        test_points = [
            ('10× radius', 10 * surface_r_ratio),
            ('2× radius', 2 * surface_r_ratio),
            ('Surface', surface_r_ratio),
        ]

        print(f"\n  {'Location':<20s}  {'r/r_s':>14s}  {'g_rr':>10s}  {'g_θθ':>12s}  {'g_δ̂':>8s}  {'δ_B':>10s}  {'Floor?':>8s}")
        print(f"  {'─'*20}  {'─'*14}  {'─'*10}  {'─'*12}  {'─'*8}  {'─'*10}  {'─'*8}")

        for label, r_ratio in test_points:
            m = metric_4d(r_ratio)
            print(f"  {label:<20s}  {r_ratio:>14.2f}  {m['grr']:>10.6f}  {m['gtt']:>12.2f}  {m['gd']:>8.4f}  {m['delta_b']:>10.6f}  {'YES' if m['floor_active'] else 'no':>8s}")

        # Key question: does the floor ever activate?
        m_surface = metric_4d(surface_r_ratio)
        det_surface = surface_r_ratio**4 * g_delta(surface_r_ratio)
        print(f"\n  Determinant at surface: {det_surface:.4e}")
        print(f"  Floor value:           {FLOOR_VAL:.4e}")
        print(f"  Ratio (det/floor):     {det_surface/FLOOR_VAL:.2e}")
        if det_surface > FLOOR_VAL:
            print(f"  → Floor NEVER activates. Standard GR is sufficient.")
            print(f"  → Metric is effectively flat: g_rr = 1, g_θθ = r², no Benford effects.")
        else:
            print(f"  → Floor IS active at the surface!")

    else:
        # Black hole — test at various r/r_s including inside
        print(f"  This is a black hole — no physical surface")

        test_points = [
            ('Far away', 100),
            ('10 r_s', 10),
            ('3 r_s', 3),
            ('Horizon', 1.0),
            ('Inside (0.5)', 0.5),
            ('Inside (0.1)', 0.1),
            ('Deep (0.01)', 0.01),
        ]

        print(f"\n  {'Location':<20s}  {'r/r_s':>10s}  {'g_rr':>12s}  {'g_θθ':>10s}  {'g_δ̂':>8s}  {'δ_B':>10s}  {'Floor?':>8s}")
        print(f"  {'─'*20}  {'─'*10}  {'─'*12}  {'─'*10}  {'─'*8}  {'─'*10}  {'─'*8}")

        for label, r_ratio in test_points:
            m = metric_4d(r_ratio)
            grr_str = f"{m['grr']:.2f}" if m['grr'] < 1e6 else f"{m['grr']:.2e}"
            print(f"  {label:<20s}  {r_ratio:>10.2f}  {grr_str:>12s}  {m['gtt']:>10.4f}  {m['gd']:>8.4f}  {m['delta_b']:>10.6f}  {'YES' if m['floor_active'] else 'no':>8s}")

        # Floor activation point
        for r_test in np.logspace(np.log10(10), np.log10(0.001), 10000):
            m = metric_4d(r_test)
            if m['floor_active']:
                print(f"\n  Floor activates at r/r_s ≈ {r_test:.4f}")
                print(f"  → Same as every other black hole (scale-free)")
                break

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
The 4D Benford metric produces:
  - Earth:        Flat. g_rr = 1.000000. No floor. No Benford effects. ✓
  - Jupiter:      Flat. g_rr = 1.000000. No floor. No Benford effects. ✓
  - Neutron star:  Flat (barely). Surface is at r/r_s ≈ 2-3. No floor.
                   But CS deviation is measurable (δ_B ~ 0.005).
                   Model correctly predicts strong-but-not-extreme gravity. ✓
  - Black holes:  Full geometry redistribution. Floor activates at r/r_s ≈ 0.64.
                   Scale-free — same r/r_s behavior for any mass. ✓

Key result: the model is BORING for boring objects and DRAMATIC for
extreme objects. It doesn't overpredict. It scales correctly.
""")
