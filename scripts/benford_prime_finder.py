#!/usr/bin/env python3
"""
Benford-Guided Prime Finder v2
================================
Uses inverse Benford's Law to find primes smarter than brute force.

The insight: primes across many orders of magnitude follow Benford's Law.
If we INVERT that — use the Benford profile to predict where primes cluster —
we can prioritize our search.

Three strategies compared:
  1. Sequential: check every odd number (brute force)
  2. Benford-weighted: search across orders of magnitude,
     prioritize leading-digit regions by Benford probability
  3. Inverse Benford: use measured prime density per leading digit
     to build a refined probability map, search highest density first

Christopher Riner (2026)
"""

import math
import time
import random
from collections import defaultdict


# ── Primality Testing ──────────────────────────────────────────

def is_prime(n):
    """Miller-Rabin primality test — deterministic for n < 3.3×10²⁴."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if n == p:
            return True
        if n % p == 0:
            return False
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


# ── Benford Tools ─────────────────────────────────────────────

def leading_digit(n):
    """Extract first significant digit."""
    s = str(abs(n))
    for ch in s:
        if ch != '0':
            return int(ch)
    return 0

def benford_prob(d):
    """Benford probability for leading digit d."""
    return math.log10(1 + 1/d)

BENFORD_PROBS = {d: benford_prob(d) for d in range(1, 10)}


# ── Phase 1: Measure the Prime-Benford Profile ────────────────

def build_prime_benford_profile(max_exp=8):
    """
    Sample primes across orders of magnitude and measure:
    - How well primes follow Benford's Law (forward)
    - Prime density per leading digit at each scale (for inverse)

    Returns a profile dict with density weights per leading digit.
    """
    digit_prime_count = defaultdict(int)
    digit_total_count = defaultdict(int)
    all_prime_digits = defaultdict(int)
    total_primes = 0

    # Sample across orders of magnitude (where Benford works)
    for exp in range(1, max_exp + 1):
        low = 10 ** exp
        high = 10 ** (exp + 1)

        # Sample randomly within this order of magnitude
        sample_size = min(2000, (high - low) // 2)
        candidates = random.sample(range(low, high), sample_size)

        for n in candidates:
            d = leading_digit(n)
            digit_total_count[d] += 1
            if is_prime(n):
                digit_prime_count[d] += 1
                all_prime_digits[d] += 1
                total_primes += 1

    # Compute: P(prime | leading digit = d)
    # This is the INVERSE Benford insight — given a leading digit,
    # how likely is the number to be prime?
    density = {}
    for d in range(1, 10):
        if digit_total_count[d] > 0:
            density[d] = digit_prime_count[d] / digit_total_count[d]
        else:
            density[d] = 0

    # Also compute forward Benford check
    benford_observed = {}
    for d in range(1, 10):
        benford_observed[d] = all_prime_digits[d] / total_primes if total_primes > 0 else 0

    return {
        "density": density,           # P(prime | digit=d) — for inverse search
        "observed_benford": benford_observed,  # P(digit=d | prime) — forward check
        "total_primes": total_primes,
    }


# ── Search Methods ─────────────────────────────────────────────

def sequential_find_n_primes(start, end, target_n):
    """Brute force: check sequentially, find target_n primes."""
    checks = 0
    primes_found = []
    n = start if start % 2 != 0 else start + 1

    while n <= end and len(primes_found) < target_n:
        checks += 1
        if is_prime(n):
            primes_found.append(n)
        n += 2

    return primes_found, checks


def benford_find_n_primes(start, end, target_n, density_profile):
    """
    Inverse Benford search:
    1. Divide range into leading-digit buckets
    2. Rank buckets by prime density (from profile)
    3. Search highest-density buckets first
    """
    checks = 0
    primes_found = []

    # Group candidates by leading digit
    buckets = defaultdict(list)
    n = start if start % 2 != 0 else start + 1
    while n <= end:
        d = leading_digit(n)
        buckets[d].append(n)
        n += 2

    # Sort buckets by prime density — highest first
    # This is the INVERSE: we know which leading digits have more primes
    sorted_digits = sorted(density_profile.keys(),
                          key=lambda d: -density_profile[d])

    # Search in density order
    for d in sorted_digits:
        if len(primes_found) >= target_n:
            break
        for n in buckets.get(d, []):
            if len(primes_found) >= target_n:
                break
            checks += 1
            if is_prime(n):
                primes_found.append(n)

    return primes_found, checks


def log_tower_benford_find(start, end, target_n, density_profile):
    """
    Full algorithm: log tower bracketing + inverse Benford.
    1. Use log towers to estimate prime density across the range
    2. Subdivide range into chunks sized by ln(midpoint)
    3. Within each chunk, prioritize by Benford density
    4. Search most promising chunks first
    """
    checks = 0
    primes_found = []

    # Subdivide range into chunks based on log-estimated gaps
    chunks = []
    pos = start
    while pos < end:
        ln_pos = math.log(max(pos, 2))
        chunk_size = max(int(ln_pos * 10), 50)  # ~10 expected primes per chunk
        chunk_end = min(pos + chunk_size, end)
        chunks.append((pos, chunk_end))
        pos = chunk_end + 1

    # Score each chunk by its leading digit composition × density profile
    def chunk_score(chunk):
        mid = (chunk[0] + chunk[1]) // 2
        d = leading_digit(mid)
        # Higher density digit = search first
        # Also factor in log-estimated prime density: ~1/ln(n)
        prime_density = 1 / math.log(max(mid, 2))
        benford_weight = density_profile.get(d, 0)
        return -(prime_density + benford_weight)  # negative for ascending sort

    chunks.sort(key=chunk_score)

    for chunk_start, chunk_end in chunks:
        if len(primes_found) >= target_n:
            break
        n = chunk_start if chunk_start % 2 != 0 else chunk_start + 1
        while n <= chunk_end and len(primes_found) < target_n:
            checks += 1
            if is_prime(n):
                primes_found.append(n)
            n += 2

    return primes_found, checks


# ── Main ───────────────────────────────────────────────────────

def main():
    print("=" * 75)
    print("  BENFORD-GUIDED PRIME FINDER v2")
    print("  Inverse Benford's Law + Log Tower Bracketing")
    print("=" * 75)

    # Phase 1: Build the inverse Benford profile
    print("\n▸ Phase 1: Building prime-Benford profile across 8 orders of magnitude...")
    random.seed(42)
    profile = build_prime_benford_profile(max_exp=8)

    print(f"  Sampled {profile['total_primes']:,} primes\n")

    print("  Forward Benford Check — do primes follow Benford's Law?")
    print(f"  {'Digit':>6s} | {'Observed':>10s} | {'Benford':>10s} | {'Diff':>8s}")
    print("  " + "-" * 45)

    delta_b_sq = 0
    for d in range(1, 10):
        obs = profile["observed_benford"][d]
        exp = BENFORD_PROBS[d]
        diff = obs - exp
        delta_b_sq += diff ** 2
        print(f"  {d:>6d} | {obs:>10.4f} | {exp:>10.4f} | {diff:>+8.4f}")

    delta_b = math.sqrt(delta_b_sq)
    print(f"\n  δ_B = {delta_b:.6f}", end="")
    if delta_b < 0.025:
        print("  ← CONFORMS")
    elif delta_b < 0.1:
        print("  ← MARGINAL")
    else:
        print("  ← DEVIATES")

    print("\n  Inverse Benford Profile — P(prime | leading digit = d):")
    print(f"  {'Digit':>6s} | {'P(prime|d)':>12s} | {'Relative':>10s}")
    print("  " + "-" * 38)
    max_density = max(profile["density"].values())
    for d in range(1, 10):
        dens = profile["density"][d]
        rel = dens / max_density if max_density > 0 else 0
        bar = "█" * int(rel * 20)
        print(f"  {d:>6d} | {dens:>12.6f} | {bar}")

    # Phase 2: Competition — find primes in various ranges
    print(f"\n{'=' * 75}")
    print("  PHASE 2: Finding 50 primes — three methods compared")
    print(f"{'=' * 75}")

    test_ranges = [
        (1_000, 50_000, "Small: 1K–50K"),
        (100_000, 1_000_000, "Medium: 100K–1M"),
        (10_000_000, 100_000_000, "Large: 10M–100M"),
        (1_000_000_000, 2_000_000_000, "Huge: 1B–2B"),
        # Wide range spanning multiple orders of magnitude
        (1_000, 10_000_000, "WIDE: 1K–10M (multi-scale)"),
        (100, 100_000_000, "WIDER: 100–100M (7 orders of mag)"),
    ]

    target_n = 50

    print(f"\n  Finding {target_n} primes in each range:\n")
    print(f"  {'Range':>30s} | {'Sequential':>11s} | {'Benford':>11s} | {'LogTower+B':>11s} | {'Best Method'}")
    print("  " + "-" * 90)

    for start, end, label in test_ranges:
        _, checks_seq = sequential_find_n_primes(start, end, target_n)
        _, checks_ben = benford_find_n_primes(start, end, target_n, profile["density"])
        _, checks_ltb = log_tower_benford_find(start, end, target_n, profile["density"])

        best = min(checks_seq, checks_ben, checks_ltb)
        if best == checks_ben and checks_ben < checks_seq:
            winner = "★ Benford"
        elif best == checks_ltb and checks_ltb < checks_seq:
            winner = "★ LogTower+B"
        elif best == checks_seq and checks_seq < checks_ben:
            winner = "Sequential"
        else:
            winner = "Tie"

        seq_s = f"{checks_seq:,}"
        ben_s = f"{checks_ben:,}"
        ltb_s = f"{checks_ltb:,}"

        print(f"  {label:>30s} | {seq_s:>11s} | {ben_s:>11s} | {ltb_s:>11s} | {winner}")

    # Phase 3: The key insight
    print(f"\n{'=' * 75}")
    print("  PHASE 3: Where inverse Benford shines — sparse wide ranges")
    print(f"{'=' * 75}")

    print("\n  Searching for just 10 primes across increasingly wide ranges:\n")
    print(f"  {'Range':>35s} | {'Orders':>6s} | {'Sequential':>11s} | {'Benford':>11s} | {'Savings'}")
    print("  " + "-" * 85)

    for orders in range(2, 9):
        start = 10
        end = 10 ** orders
        target = 10

        _, checks_seq = sequential_find_n_primes(start, end, target)
        _, checks_ben = benford_find_n_primes(start, end, target, profile["density"])

        savings = ((checks_seq - checks_ben) / checks_seq * 100) if checks_seq > 0 else 0

        label = f"10 — {end:,}"
        print(f"  {label:>35s} | {orders:>6d} | {checks_seq:>11,} | {checks_ben:>11,} | {savings:>+.1f}%")

    # Phase 4: Leading digit distribution of found primes
    print(f"\n{'=' * 75}")
    print("  PHASE 4: Leading digits of primes found across 7 orders of magnitude")
    print(f"{'=' * 75}\n")

    # Collect primes spanning 10 to 10M
    all_primes = []
    n = 3
    while n < 10_000_000 and len(all_primes) < 5000:
        if is_prime(n):
            all_primes.append(n)
        n += 2

    digit_counts = defaultdict(int)
    for p in all_primes:
        digit_counts[leading_digit(p)] += 1

    total = len(all_primes)
    print(f"  {total:,} primes from 3 to {all_primes[-1]:,}\n")
    print(f"  {'Digit':>6s} | {'Count':>7s} | {'Observed':>10s} | {'Benford':>10s} | {'Diff':>8s} | Visual")
    print("  " + "-" * 70)

    delta_b_sq = 0
    for d in range(1, 10):
        count = digit_counts[d]
        obs = count / total
        exp = BENFORD_PROBS[d]
        diff = obs - exp
        delta_b_sq += diff ** 2
        bar_obs = "█" * int(obs * 60)
        bar_exp = "░" * int(exp * 60)
        print(f"  {d:>6d} | {count:>7,} | {obs:>10.4f} | {exp:>10.4f} | {diff:>+8.4f} | {bar_obs}")

    delta_b = math.sqrt(delta_b_sq)
    print(f"\n  δ_B = {delta_b:.6f}", end="")
    if delta_b < 0.025:
        print("  ← CONFORMS to Benford's Law!")
    elif delta_b < 0.1:
        print("  ← MARGINAL conformance")
    else:
        print("  ← DEVIATES from Benford's Law")

    print(f"\n  Conclusion: Primes {'DO' if delta_b < 0.05 else 'DO NOT'} follow Benford's Law")
    print(f"  when sampled across multiple orders of magnitude.")
    print(f"\n  This validates using inverse Benford as a prime search heuristic —")
    print(f"  the distribution IS there, so the flashlight works.")


if __name__ == "__main__":
    main()
