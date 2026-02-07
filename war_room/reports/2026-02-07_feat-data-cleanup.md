# Worker Report: feat/data-cleanup

**Date:** 2026-02-07
**Branch:** feat/data-cleanup

## What Changed and Why

The summary.json only tracked 4 tests (the quantum statistics validation tests) while 69 completed test results existed in results/individual/. The test_queue.json was similarly out of sync with only 4 entries.

### Changes:
1. **Rebuilt `results/summary.json`** — now includes all 69 completed datasets with their verdicts, delta_B, MAD, chi-squared p-values, KS p-values, and category information.
2. **Rebuilt `data/test_queue.json`** — now tracks all 85 tests (69 complete + 16 insufficient_data).
3. **Attempted to run 16 missing tests** from the CLAUDE.md master list — all 16 failed due to insufficient data points in their fetchers (<50 valid values).

## What Succeeded

- Full rebuild of summary.json with all 69 datasets
- Proper categorization and verdict tracking
- Queue now accurately reflects all test statuses

## What Failed

All 16 missing tests failed due to insufficient data in their fetchers:

| Test ID | Data Points | Reason |
|---------|------------|--------|
| galaxy_redshifts | 42 | Near threshold |
| asteroid_diameters | 40 | Near threshold |
| volcano_eruption_vei | 0 | No valid data (discrete 0-8 scale, zeros filtered) |
| river_lengths | 0 | Empty fetcher |
| lake_areas | 0 | Empty fetcher |
| island_areas | 4 | Barely any data |
| waterfall_heights | 41 | Near threshold |
| crater_diameters | 0 | Empty fetcher |
| genome_sizes | 0 | Empty fetcher |
| animal_lifespans | 0 | Empty fetcher |
| cell_counts_human | 40 | Near threshold |
| particle_masses | 38 | Near threshold |
| physical_constants_si | 40 | Near threshold |
| building_heights | 0 | Empty fetcher |
| airport_passengers | 0 | Empty fetcher |
| mersenne_primes | 14 | Only 14 known Mersenne primes exist |

To fix these, the data_fetchers.py would need expanded datasets. Several are near the 50-point threshold and could be fixed by adding ~10 more values each.

## Interesting Findings

### Verdict Distribution (69 tests):
- **CONFORMS (4):** Fibonacci numbers, Powers of 2, Bose-Einstein, ZIP code populations
- **MARGINAL (25):** Most real-world datasets fall here
- **DEVIATES (15):** Including S&P 500 prices (delta_B=0.796!), mountain heights (0.919), city populations (0.238)
- **EXPLAINED (23):** Deviations with scientific explanations (earthquake magnitudes, primes, etc.)
- **HAS MASS (2):** Fermi-Dirac and Maxwell-Boltzmann (paper's predicted deviation pattern)

### Notable Results:
- **Mountain heights** have the strongest deviation (delta_B=0.919) — extremely narrow range
- **Fibonacci & Powers of 2** are nearly perfect (delta_B ~0.003) — mathematically proven
- **S&P 500 prices** deviate strongly (delta_B=0.796) — narrow trading range
- **Bose-Einstein** conforms beautifully (delta_B=0.006) — validates the paper's central claim

## Questions for the War Room

1. Should we lower the minimum data threshold from 50 to 30 for some tests? Several fetchers are close (38-42 points).
2. The data fetchers for 8 tests return 0 values — these need real data to be added to scripts/data_fetchers.py.
3. Mersenne primes (n=14) will never reach 50 — should we flag this as a permanent skip?
