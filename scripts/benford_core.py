"""
benford_core.py -- Core analysis engine for Benford's Law testing.

Pure computation module: no file I/O, no plotting, no side effects.
All public functions operate on in-memory data structures.
"""

import math
from typing import Dict, List, Optional, Union

import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Reference distributions
# ---------------------------------------------------------------------------

BENFORD_EXPECTED: Dict[int, float] = {
    d: math.log10(1 + 1 / d) for d in range(1, 10)
}

BENFORD_FIRST_TWO: Dict[int, float] = {
    dd: math.log10(1 + 1 / dd) for dd in range(10, 100)
}

# ---------------------------------------------------------------------------
# Digit extraction helpers
# ---------------------------------------------------------------------------


def first_digit(x: Union[int, float]) -> Optional[int]:
    """Return the first significant digit (1-9) of a positive number.

    Returns None for values <= 0, NaN, or Inf.
    """
    if x is None:
        return None
    try:
        val = float(x)
    except (TypeError, ValueError):
        return None
    if val <= 0 or math.isnan(val) or math.isinf(val):
        return None
    return int(f"{val:.15e}"[0])


def first_two_digits(x: Union[int, float]) -> Optional[int]:
    """Return the first two significant digits (10-99) of a positive number.

    Returns None for values <= 0, NaN, or Inf.
    """
    if x is None:
        return None
    try:
        val = float(x)
    except (TypeError, ValueError):
        return None
    if val <= 0 or math.isnan(val) or math.isinf(val):
        return None
    return int(f"{val:.15e}"[0:4].replace(".", "")[:2])


# ---------------------------------------------------------------------------
# Distribution computation
# ---------------------------------------------------------------------------


def observed_distribution(digits: List[int]) -> Dict[int, float]:
    """Compute the observed proportion for each digit 1-9.

    Parameters
    ----------
    digits : list[int]
        First-digit values, each expected to be in 1..9.

    Returns
    -------
    dict mapping digit (1-9) to observed proportion (float).
    All nine keys are always present; missing digits get 0.0.
    """
    if not digits:
        return {d: 0.0 for d in range(1, 10)}
    n = len(digits)
    counts: Dict[int, int] = {d: 0 for d in range(1, 10)}
    for d in digits:
        if d in counts:
            counts[d] += 1
    return {d: round(counts[d] / n, 6) for d in range(1, 10)}


# ---------------------------------------------------------------------------
# Statistical tests
# ---------------------------------------------------------------------------


def chi_squared_test(observed: Dict[int, float], n: int) -> Dict[str, float]:
    """Chi-squared goodness-of-fit test against Benford's Law.

    Parameters
    ----------
    observed : dict
        Mapping digit (1-9) -> observed proportion.
    n : int
        Total number of observations.

    Returns
    -------
    dict with 'statistic' and 'p_value'.
    """
    if n <= 0:
        return {"statistic": 0.0, "p_value": 1.0}

    obs_counts = np.array([observed.get(d, 0.0) * n for d in range(1, 10)])
    exp_counts = np.array([BENFORD_EXPECTED[d] * n for d in range(1, 10)])
    # Normalize expected to match observed sum (scipy requires exact agreement)
    obs_sum = obs_counts.sum()
    exp_sum = exp_counts.sum()
    if exp_sum > 0:
        exp_counts = exp_counts * (obs_sum / exp_sum)

    # Guard against zero expected (cannot happen for Benford, but be safe)
    if np.any(exp_counts == 0):
        return {"statistic": 0.0, "p_value": 1.0}

    stat, p_val = stats.chisquare(obs_counts, f_exp=exp_counts)
    return {
        "statistic": round(float(stat), 6),
        "p_value": round(float(p_val), 6),
    }


def mean_absolute_deviation(observed: Dict[int, float]) -> Dict[str, Union[float, str]]:
    """Mean Absolute Deviation (MAD) conformity test.

    MAD = (1/9) * sum(|P_obs(d) - P_benford(d)|) for d = 1..9

    Returns
    -------
    dict with 'mad' (float) and 'classification' (str).
    """
    mad = sum(abs(observed.get(d, 0.0) - BENFORD_EXPECTED[d]) for d in range(1, 10)) / 9
    mad = round(mad, 6)

    if mad < 0.006:
        classification = "Close conformity"
    elif mad < 0.012:
        classification = "Acceptable conformity"
    elif mad < 0.015:
        classification = "Marginally acceptable"
    else:
        classification = "Nonconformity"

    return {"mad": mad, "classification": classification}


def ks_test(digits: List[int]) -> Dict[str, float]:
    """Kolmogorov-Smirnov test for first-digit data vs Benford CDF.

    Builds discrete CDFs manually and computes the KS statistic as the
    maximum absolute difference between the empirical and theoretical CDFs
    evaluated at every digit boundary.

    Returns
    -------
    dict with 'statistic' and 'p_value'.
    """
    if not digits:
        return {"statistic": 0.0, "p_value": 1.0}

    n = len(digits)

    # Count occurrences per digit
    counts = {d: 0 for d in range(1, 10)}
    for d in digits:
        if d in counts:
            counts[d] += 1

    # Build empirical CDF and Benford CDF at each digit point 1..9
    empirical_cdf = 0.0
    benford_cdf = 0.0
    max_diff = 0.0

    for d in range(1, 10):
        empirical_cdf += counts[d] / n
        benford_cdf += BENFORD_EXPECTED[d]
        diff = abs(empirical_cdf - benford_cdf)
        if diff > max_diff:
            max_diff = diff

    # Approximate p-value using the Kolmogorov distribution.
    # The KS distribution parameter is sqrt(n) * D.
    sqrt_n_d = math.sqrt(n) * max_diff
    # Use scipy's kstwo distribution (two-sided KS) for the p-value
    p_val = float(stats.kstwo.sf(max_diff, n))

    return {
        "statistic": round(max_diff, 6),
        "p_value": round(min(max(p_val, 0.0), 1.0), 6),
    }


# ---------------------------------------------------------------------------
# Deviation measures
# ---------------------------------------------------------------------------


def euclidean_deviation(observed: Dict[int, float]) -> float:
    """Euclidean distance (delta_B) between observed and Benford distributions.

    delta_B = sqrt(sum((P_obs(d) - P_benford(d))^2))
    """
    delta_b = math.sqrt(
        sum((observed.get(d, 0.0) - BENFORD_EXPECTED[d]) ** 2 for d in range(1, 10))
    )
    return round(delta_b, 6)


def per_digit_deviation(observed: Dict[int, float]) -> Dict[int, float]:
    """Per-digit deviation: epsilon(d) = P_obs(d) - P_benford(d) for d=1..9."""
    return {
        d: round(observed.get(d, 0.0) - BENFORD_EXPECTED[d], 6) for d in range(1, 10)
    }


# ---------------------------------------------------------------------------
# First-two-digit analysis
# ---------------------------------------------------------------------------


def first_two_digit_analysis(digits: List[int]) -> Dict:
    """Analyse first-two-digit distribution (values 10-99).

    Parameters
    ----------
    digits : list[int]
        Each element should be in 10..99.

    Returns
    -------
    dict with 'chi_squared' (statistic, p_value), 'top_deviations' (list),
    and 'n'.
    """
    if not digits:
        return {
            "chi_squared": {"statistic": 0.0, "p_value": 1.0},
            "top_deviations": [],
            "n": 0,
        }

    n = len(digits)

    # Count occurrences
    counts = {dd: 0 for dd in range(10, 100)}
    for dd in digits:
        if dd in counts:
            counts[dd] += 1

    obs_counts = np.array([counts[dd] for dd in range(10, 100)])
    exp_counts = np.array([BENFORD_FIRST_TWO[dd] * n for dd in range(10, 100)])

    # Chi-squared test (df = 90 - 1 = 89)
    stat, p_val = stats.chisquare(obs_counts, f_exp=exp_counts)

    # Compute deviations for each two-digit combo
    deviations = []
    for dd in range(10, 100):
        obs_prop = counts[dd] / n if n > 0 else 0.0
        exp_prop = BENFORD_FIRST_TWO[dd]
        dev = obs_prop - exp_prop
        deviations.append(
            {
                "digits": dd,
                "observed": round(obs_prop, 6),
                "expected": round(exp_prop, 6),
                "deviation": round(dev, 6),
            }
        )

    # Sort by absolute deviation descending, take top 10
    deviations.sort(key=lambda x: abs(x["deviation"]), reverse=True)
    top_deviations = deviations[:10]

    return {
        "chi_squared": {
            "statistic": round(float(stat), 6),
            "p_value": round(float(p_val), 6),
        },
        "top_deviations": top_deviations,
        "n": n,
    }


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


def classify_result(
    chi_sq_result: Dict[str, float],
    mad_result: Dict[str, Union[float, str]],
    delta_b: float,
) -> str:
    """Classify overall conformity to Benford's Law.

    Uses delta_B (Euclidean deviation) as the primary metric, per the
    Riner (2026) framework.  Thresholds:
        delta_B < 0.025  ->  CONFORMS   (near-exact Benford conformance)
        delta_B < 0.10   ->  MARGINAL   (moderate deviation)
        delta_B >= 0.10  ->  DEVIATES   (significant deviation)
    """
    if delta_b < 0.025:
        return "CONFORMS"
    elif delta_b < 0.10:
        return "MARGINAL"
    else:
        return "DEVIATES"


# ---------------------------------------------------------------------------
# Master analysis function
# ---------------------------------------------------------------------------


def run_full_analysis(values: List[Union[int, float]]) -> Dict:
    """Run the complete Benford's Law analysis on a list of raw numbers.

    Parameters
    ----------
    values : list
        Raw positive numbers. Non-positive, NaN, and Inf values are filtered.

    Returns
    -------
    dict with all analysis results.
    """
    # Filter to positive, finite values
    positive_values = []
    for v in values:
        try:
            fv = float(v)
            if fv > 0 and math.isfinite(fv):
                positive_values.append(fv)
        except (TypeError, ValueError):
            continue

    n = len(positive_values)

    # Edge case: insufficient data
    if n == 0:
        return {
            "n": 0,
            "n_original": len(values),
            "n_filtered": 0,
            "verdict": "MARGINAL",
            "message": "No positive values to analyse.",
            "first_digit": {},
            "first_two_digit": None,
        }

    # Extract first digits
    fd_list = [first_digit(v) for v in positive_values]
    fd_list = [d for d in fd_list if d is not None]

    if not fd_list:
        return {
            "n": n,
            "n_original": len(values),
            "n_filtered": n,
            "verdict": "MARGINAL",
            "message": "Could not extract any first digits.",
            "first_digit": {},
            "first_two_digit": None,
        }

    n_digits = len(fd_list)

    # Observed distribution
    obs = observed_distribution(fd_list)

    # Statistical tests
    chi_sq = chi_squared_test(obs, n_digits)
    mad = mean_absolute_deviation(obs)
    ks = ks_test(fd_list)
    delta_b = euclidean_deviation(obs)
    per_digit = per_digit_deviation(obs)

    # Classification
    verdict = classify_result(chi_sq, mad, delta_b)

    # First-two-digit analysis (only for larger samples)
    ftd_result = None
    if n_digits > 500:
        ftd_list = [first_two_digits(v) for v in positive_values]
        ftd_list = [d for d in ftd_list if d is not None]
        if ftd_list:
            ftd_result = first_two_digit_analysis(ftd_list)

    return {
        "n": n_digits,
        "n_original": len(values),
        "n_filtered": n,
        "verdict": verdict,
        "first_digit": {
            "observed_distribution": obs,
            "expected_distribution": {d: round(BENFORD_EXPECTED[d], 6) for d in range(1, 10)},
            "chi_squared": chi_sq,
            "mad": mad,
            "ks_test": ks,
            "euclidean_deviation": delta_b,
            "per_digit_deviation": per_digit,
        },
        "first_two_digit": ftd_result,
    }


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import random

    # ------------------------------------------------------------------
    # Test 1: Fibonacci numbers (classic Benford conformer)
    # ------------------------------------------------------------------
    print("=" * 60)
    print("Test 1: First 1000 Fibonacci numbers")
    print("=" * 60)

    fibs = []
    a, b = 1, 1
    for _ in range(1000):
        fibs.append(a)
        a, b = b, a + b

    result = run_full_analysis(fibs)
    print(f"  n            = {result['n']}")
    print(f"  verdict      = {result['verdict']}")
    print(f"  chi-sq stat  = {result['first_digit']['chi_squared']['statistic']}")
    print(f"  chi-sq p     = {result['first_digit']['chi_squared']['p_value']}")
    print(f"  MAD          = {result['first_digit']['mad']['mad']}")
    print(f"  MAD class    = {result['first_digit']['mad']['classification']}")
    print(f"  delta_b      = {result['first_digit']['euclidean_deviation']}")
    print(f"  KS stat      = {result['first_digit']['ks_test']['statistic']}")
    print(f"  KS p         = {result['first_digit']['ks_test']['p_value']}")
    print()

    assert result["verdict"] == "CONFORMS", (
        f"Fibonacci should CONFORM, got {result['verdict']}"
    )
    print("  PASSED: Fibonacci numbers conform to Benford's Law.\n")

    # ------------------------------------------------------------------
    # Test 2: Uniform random digits (should deviate)
    # ------------------------------------------------------------------
    print("=" * 60)
    print("Test 2: Uniform random digits 1-9 (1000 samples)")
    print("=" * 60)

    random.seed(42)
    # Generate numbers whose first digits are uniformly distributed 1-9.
    # Easiest approach: pick digits uniformly at random, then create numbers
    # that start with those digits.
    uniform_values = [random.randint(1, 9) * 10 ** random.randint(0, 5) for _ in range(1000)]

    result2 = run_full_analysis(uniform_values)
    print(f"  n            = {result2['n']}")
    print(f"  verdict      = {result2['verdict']}")
    print(f"  chi-sq stat  = {result2['first_digit']['chi_squared']['statistic']}")
    print(f"  chi-sq p     = {result2['first_digit']['chi_squared']['p_value']}")
    print(f"  MAD          = {result2['first_digit']['mad']['mad']}")
    print(f"  MAD class    = {result2['first_digit']['mad']['classification']}")
    print(f"  delta_b      = {result2['first_digit']['euclidean_deviation']}")
    print()

    assert result2["verdict"] == "DEVIATES", (
        f"Uniform digits should DEVIATE, got {result2['verdict']}"
    )
    print("  PASSED: Uniform random digits deviate from Benford's Law.\n")

    # ------------------------------------------------------------------
    # Edge-case tests
    # ------------------------------------------------------------------
    print("=" * 60)
    print("Test 3: Edge cases")
    print("=" * 60)

    # Empty list
    r = run_full_analysis([])
    assert r["n"] == 0
    print("  PASSED: empty list handled.")

    # All zeros
    r = run_full_analysis([0, 0, 0])
    assert r["n"] == 0
    print("  PASSED: all-zeros handled.")

    # Single value
    r = run_full_analysis([42])
    assert r["n"] == 1
    print("  PASSED: single value handled.")

    # Negative values filtered
    r = run_full_analysis([-5, -10, 100, 200])
    assert r["n"] == 2
    print("  PASSED: negative values filtered.")

    print("\nAll tests passed.")
