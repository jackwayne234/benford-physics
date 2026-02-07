# Benford's Law Testing Lab — CLAUDE.md

## Project Overview

This project systematically tests **Benford's Law** against real-world datasets across astrophysics, geology, biology, economics, mathematics, and more. It uses multiple Claude Code sub-agents working in parallel, writing results to a shared JSON store and a live dashboard website.

The analytical framework is based on the paper:
> **"Complete Monotonicity and Benford's Law: Deriving Quantum Statistics from the Significant Digit Distribution"** — C. Riner (2026)
>
> Key insight: a distribution satisfies Benford's law exactly iff it is completely monotonic (all coefficients non-negative in its exponential series). Deviations are governed by Dirichlet series factors computable via Fourier decomposition.

---

## Architecture

```
benfords-lab/
├── CLAUDE.md                  ← You are here (agent instructions)
├── data/
│   ├── test_queue.json        ← Master list of all tests
│   └── sources/               ← Cached raw data per test
├── results/
│   ├── summary.json           ← Aggregated results (website reads this)
│   └── individual/            ← Per-test detailed JSON results
├── scripts/
│   ├── benford_core.py        ← Core Benford analysis engine
│   ├── data_fetchers.py       ← Dataset acquisition functions
│   └── run_test.py            ← Single-test runner (agents call this)
├── website/
│   └── index.html             ← Dashboard (serves from results/)
└── Benford_Quantum_Statistics_Paper.pdf  ← Reference paper
```

---

## Agent Orchestration

### How to Spawn Agents

Use the **Task tool** to spawn 3–5 parallel sub-agents. Each agent independently:

1. Reads `data/test_queue.json`
2. Claims the next test with `"status": "queued"` (atomically sets it to `"in_progress"`)
3. Fetches or generates the dataset
4. Runs the full Benford analysis pipeline
5. Writes results to `results/individual/<test_id>.json`
6. Updates `results/summary.json` (append or update entry)
7. Sets the test status to `"complete"` (or `"error"`)
8. Moves to the next queued test

### Spawn Pattern

```
Agent 1: "You are Benford Agent 1. Run tests from data/test_queue.json starting
from the TOP of the queue. Use scripts/run_test.py for each test. Follow all
analysis protocols in CLAUDE.md."

Agent 2: "You are Benford Agent 2. Run tests from data/test_queue.json starting
from the BOTTOM of the queue. Use scripts/run_test.py for each test. Follow all
analysis protocols in CLAUDE.md."

Agent 3: "You are Benford Agent 3. Run tests from data/test_queue.json starting
from the MIDDLE of the queue. Use scripts/run_test.py for each test. Follow all
analysis protocols in CLAUDE.md."
```

### Concurrency Rules

- **File locking**: Use `fcntl.flock()` when reading/writing shared JSON files (test_queue.json, summary.json)
- **Atomic claims**: Read → find first "queued" → set "in_progress" with agent ID → write back, all under lock
- **Graceful failures**: If data fetch fails, set status to `"error"` with a message, move on
- **Minimum data**: Require ≥50 data points for a valid test; flag <100 as `"low_confidence"`

---

## Analysis Protocol

### Step 1: Data Collection

For each test, the agent must:

- Fetch data from the specified source (public API, Wikipedia tables, mathematical generation, web scrape)
- Extract numerical values — these must be **positive real numbers** (Benford's law applies to magnitudes)
- Strip units, convert to float, discard zeros and negatives
- Cache raw data to `data/sources/<test_id>.json` with metadata:
  ```json
  {
    "test_id": "black_hole_masses",
    "source_url": "https://...",
    "fetch_timestamp": "2026-02-06T...",
    "raw_count": 312,
    "valid_count": 298,
    "values": [6.5e9, 4.1e6, ...],
    "notes": "Masses in solar masses from the McConnell & Ma catalog"
  }
  ```

### Step 2: First-Digit Extraction

Extract the **first significant digit** d₁ ∈ {1,...,9} from each value:
```python
import math
def first_digit(x):
    """Extract first significant digit of a positive number."""
    if x <= 0:
        return None
    return int(str(f"{x:.15e}")[0])
```

Also extract **first two digits** d₁d₂ ∈ {10,...,99} if n > 500 (for finer resolution analysis).

### Step 3: Statistical Tests

Run ALL of the following:

#### Chi-Squared Test
Compare observed digit counts vs. Benford expected counts.
```
Expected count for digit d: E(d) = n × log₁₀(1 + 1/d)
χ² = Σ (O(d) - E(d))² / E(d)  for d = 1..9
df = 8, report p-value
```

#### Mean Absolute Deviation (MAD)
```
MAD = (1/9) × Σ |P_obs(d) - P_benford(d)|  for d = 1..9
```
Classification thresholds (per Nigrini):
- MAD < 0.006 → **Close conformity**
- 0.006 ≤ MAD < 0.012 → **Acceptable conformity**
- 0.012 ≤ MAD < 0.015 → **Marginally acceptable**
- MAD ≥ 0.015 → **Nonconformity**

#### Kolmogorov-Smirnov Test
Compare empirical CDF of first digits against Benford CDF. Report D-statistic and p-value.

#### Euclidean (L²) Deviation — δ_B
This is the primary metric from the Riner (2026) framework:
```
δ_B = √( Σ [P_obs(d) - log₁₀(1 + 1/d)]² )  for d = 1..9
```
- δ_B = 0 → exact conformance
- δ_B < 0.02 → strong conformance
- δ_B 0.02–0.05 → moderate conformance
- δ_B > 0.05 → significant deviation

#### Per-Digit Deviation — ε(d)
```
ε(d) = P_obs(d) - log₁₀(1 + 1/d)
```
Signed measure for each digit. Report all 9 values. Look for:
- Systematic patterns (oscillatory → possible alternating-sign structure)
- Excess at d=1 (possible fabrication or rounding bias)
- Deficit at d=1 with excess at d=9 (possible truncation or selection effects)

#### First-Two-Digit Analysis (if n > 500)
Extend to d₁d₂ ∈ {10,...,99} using generalized Benford:
```
P(d₁d₂) = log₁₀(1 + 1/(d₁d₂))
```

### Step 4: Result Classification

Rate each test using this rubric:

| Verdict | Symbol | Criteria |
|---------|--------|----------|
| **CONFORMS** | ✅ | p_χ² > 0.05 AND MAD < 0.012 AND δ_B < 0.03 |
| **MARGINAL** | ⚠️ | p_χ² > 0.01 OR (MAD 0.012–0.015) OR (δ_B 0.03–0.06) |
| **DEVIATES** | ❌ | p_χ² < 0.01 AND MAD > 0.015 AND δ_B > 0.06 |
| **INTERESTING** | 🔬 | Deviates in a scientifically notable pattern — explain why |

For datasets that DEVIATE, the agent should attempt to explain **why**:
- Is the data range too narrow (< 2 orders of magnitude)?
- Is there a selection effect or reporting bias?
- Is the underlying distribution uniform, normal, or otherwise non-multiplicative?
- Does the deviation pattern suggest alternating-sign structure (Fermi-Dirac–like)?

### Step 5: Write Results

Each completed test produces `results/individual/<test_id>.json`:

```json
{
  "test_id": "black_hole_masses",
  "display_name": "Supermassive Black Hole Masses",
  "category": "Astrophysics",
  "description": "Observed masses of supermassive black holes in solar mass units",
  "status": "complete",
  "data_points": 298,
  "source": "McConnell & Ma (2013) catalog + updates",
  "source_url": "https://...",
  "digit_distribution": {
    "1": 0.312, "2": 0.168, "3": 0.119, "4": 0.098,
    "5": 0.081, "6": 0.064, "7": 0.057, "8": 0.053, "9": 0.048
  },
  "expected_distribution": {
    "1": 0.301, "2": 0.176, "3": 0.125, "4": 0.097,
    "5": 0.079, "6": 0.067, "7": 0.058, "8": 0.051, "9": 0.046
  },
  "per_digit_deviation": {
    "1": 0.011, "2": -0.008, "3": -0.006, "4": 0.001,
    "5": 0.002, "6": -0.003, "7": -0.001, "8": 0.002, "9": 0.002
  },
  "chi_squared": {"statistic": 3.41, "p_value": 0.906},
  "mad": 0.004,
  "mad_classification": "Close conformity",
  "ks_test": {"statistic": 0.038, "p_value": 0.812},
  "delta_b": 0.016,
  "verdict": "CONFORMS",
  "confidence": "high",
  "notes": "Strong Benford conformance. Black hole masses span ~6 orders of magnitude, a key predictor of conformance.",
  "interesting_findings": "",
  "agent_id": "agent_1",
  "timestamp": "2026-02-06T14:30:00Z"
}
```

Then update `results/summary.json` (under file lock).

---

## Master Test Queue

### Astrophysics & Cosmology
| test_id | Display Name | Why It's Interesting |
|---------|-------------|---------------------|
| `black_hole_masses` | Supermassive Black Hole Masses | Span ~6 orders of magnitude; multiplicative growth |
| `star_distances` | Distances to Nearest Stars | Geometric distribution in 3D space |
| `galaxy_redshifts` | Galaxy Redshift Values | Expansion-driven; spans large range |
| `exoplanet_orbital_periods` | Exoplanet Orbital Periods | Kepler's laws → power-law relationship |
| `exoplanet_radii` | Exoplanet Radii | Detection bias may distort digits |
| `quasar_luminosities` | Quasar Luminosities | Extreme range; AGN physics |
| `cmb_anisotropy` | CMB Temperature Anisotropies | Quantum fluctuations from early universe |
| `gravitational_wave_strain` | Gravitational Wave Amplitudes | LIGO/Virgo detections |
| `pulsar_periods` | Pulsar Rotation Periods | Bimodal: millisecond vs. regular |
| `asteroid_diameters` | Asteroid Diameters | Power-law size distribution |
| `comet_orbital_periods` | Comet Orbital Periods | Wide range |
| `supernova_luminosities` | Type Ia Supernova Luminosities | Standard candles — narrow range |
| `solar_flare_energies` | Solar Flare Energies | Power-law distributed events |
| `gamma_ray_burst_durations` | Gamma-Ray Burst Durations | Bimodal (short vs. long) |

### Earth Sciences & Geology
| test_id | Display Name | Why It's Interesting |
|---------|-------------|---------------------|
| `earthquake_magnitudes` | Earthquake Magnitudes | Already log scale |
| `earthquake_energies` | Earthquake Energies (Joules) | Huge range |
| `volcano_eruption_vei` | Volcanic Explosivity Index | Discrete 0–8 scale |
| `river_lengths` | World River Lengths | 2–3 orders of magnitude |
| `mountain_heights` | Mountain Peak Elevations | Narrow range |
| `lake_areas` | Lake Surface Areas | Many orders of magnitude |
| `island_areas` | Island Areas | Power-law expected |
| `waterfall_heights` | Waterfall Heights | Moderate range |
| `ocean_depths` | Ocean Depth Measurements | Narrow range |
| `tectonic_plate_velocities` | Tectonic Plate Velocities | Very narrow |
| `crater_diameters` | Impact Crater Diameters | Power-law |
| `mineral_densities` | Mineral Densities | Narrow range |

### Biology & Life Sciences
| test_id | Display Name | Why It's Interesting |
|---------|-------------|---------------------|
| `animal_populations` | Species Population Estimates | ~10 orders of magnitude |
| `genome_sizes` | Genome Sizes Across Species | C-value paradox; huge range |
| `protein_molecular_weights` | Human Proteome Molecular Weights | Well-characterized |
| `tree_heights` | Maximum Tree Heights by Species | Moderate range |
| `animal_lifespans` | Maximum Lifespans by Species | 1 day to 500+ years |
| `cell_counts_human` | Cell Counts in Human Body | Huge range |
| `viral_genome_lengths` | Virus Genome Lengths | ~2k–2M nucleotides |
| `bird_migration_distances` | Bird Migration Distances | Wide range |
| `insect_species_per_family` | Species per Insect Family | Highly skewed |
| `bacterial_doubling_times` | Bacterial Doubling Times | Growth rates |

### Physics & Chemistry
| test_id | Display Name | Why It's Interesting |
|---------|-------------|---------------------|
| `atomic_weights` | Atomic Weights | Only 118 points — low n |
| `radioactive_half_lives` | Radioactive Half-Lives | ~50 orders of magnitude |
| `speed_of_sound` | Speed of Sound in Materials | Moderate range |
| `element_melting_points` | Element Melting Points | Narrow range |
| `element_boiling_points` | Element Boiling Points | Similar |
| `particle_masses` | Known Particle Masses | eV to GeV |
| `physical_constants_si` | Fundamental Physical Constants | Classic test |
| `nuclear_half_lives` | Nuclear Half-Lives | Confirmed by Ni & Ren (2008) |
| `hadron_widths` | Hadron Full Widths | Shao & Ma (2009) |
| `atomic_spectra_wavelengths` | Atomic Spectral Lines | Ralchenko & Pain (2024) |

### Demographics & Geography
| test_id | Display Name | Why It's Interesting |
|---------|-------------|---------------------|
| `country_populations` | Country Populations | Classic Benford test |
| `city_populations` | World City Populations | Strong conformance expected |
| `country_gdp` | GDP by Country | Financial data |
| `country_areas` | Country Land Areas | Moderate range |
| `building_heights` | Tallest Buildings | Narrow range |
| `airport_passengers` | Airport Passenger Counts | Large range |
| `us_county_populations` | US County Populations | ~3100 data points |
| `zip_code_populations` | US ZIP Code Populations | ~42,000 data points |
| `election_vote_counts` | Election Vote Counts | Fraud detection context |

### Mathematics & Number Theory
| test_id | Display Name | Why It's Interesting |
|---------|-------------|---------------------|
| `fibonacci_numbers` | Fibonacci Numbers (first 1000) | Proven to satisfy Benford exactly |
| `powers_of_2` | Powers of 2 | Equidistributed mod 1 → exact Benford |
| `prime_numbers` | First 10,000 Primes | Do NOT satisfy Benford |
| `factorials` | Factorials (1! to 200!) | Stirling → conformance |
| `collatz_lengths` | Collatz Sequence Lengths | Chaotic system |
| `twin_prime_gaps` | Twin Prime Gap Sizes | Number theory |
| `riemann_zeta_zeros` | Imaginary Parts of Zeta Zeros | Deep connection to ζ(s) in framework |
| `catalan_numbers` | Catalan Numbers | Combinatorial growth |
| `mersenne_primes` | Known Mersenne Primes | Very few but huge range |

### Finance & Economics
| test_id | Display Name | Why It's Interesting |
|---------|-------------|---------------------|
| `sp500_prices` | S&P 500 Closing Prices | Financial time series |
| `crypto_market_caps` | Cryptocurrency Market Caps | Volatile; wide range |
| `fortune500_revenue` | Fortune 500 Revenues | Well-characterized |
| `bitcoin_transaction_values` | Bitcoin Transaction Values | Decentralized |

### Health & Environmental
| test_id | Display Name | Why It's Interesting |
|---------|-------------|---------------------|
| `covid_cases_by_country` | COVID-19 Case Counts | Reporting anomalies? |
| `disease_incidence_rates` | Disease Incidence Rates | Epidemiological data |
| `carbon_emissions` | CO₂ Emissions by Country | Environmental policy |
| `species_counts_by_genus` | Species per Genus | Taxonomic distribution |

### Quantum Statistics — Validation Tests (Run First!)
| test_id | Display Name | Why It's Interesting |
|---------|-------------|---------------------|
| `bose_einstein_numerical` | Bose-Einstein Distribution | **Must conform exactly** — validates framework |
| `fermi_dirac_numerical` | Fermi-Dirac Distribution | **Must show periodic deviation** — validates Dirichlet factor |
| `maxwell_boltzmann_numerical` | Maxwell-Boltzmann Distribution | **Approximate conformance** |
| `planck_radiation_spectrum` | Planck Black-Body Spectrum | Bosonic (photon) distribution |

---

## Adding New Tests

Anyone can add tests to `data/test_queue.json`:
```json
{
  "test_id": "unique_snake_case_id",
  "display_name": "Human Readable Name",
  "category": "Category Name",
  "description": "What and why",
  "data_source": "url | 'generated:method_name' | 'api:endpoint'",
  "status": "queued",
  "priority": 1,
  "added_by": "username_or_anonymous",
  "added_date": "2026-02-06"
}
```

---

## Dependencies

```bash
pip install numpy scipy pandas matplotlib requests beautifulsoup4
```

---

## Important Notes

- The website reads `results/summary.json` — keep it updated after every completed test
- Run quantum statistics validation tests FIRST to confirm the engine works against known results
- For generated mathematical sequences, use Python's native arbitrary-precision `int` to avoid floating-point artifacts
- Minimum 50 data points required; flag <100 as low confidence
- If a dataset has fewer than 50 valid points, mark `"insufficient_data"` and skip
