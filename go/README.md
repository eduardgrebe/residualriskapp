# Risk Days Bootstrap Calculator - Go Implementation

High-performance Go implementation of the risk days bootstrap calculation for HIV transfusion transmission risk estimation.

## Overview

This Go implementation provides a significant performance improvement over the Python version through:
- Native compiled code execution
- Efficient goroutine-based parallelism
- Optimized numerical integration
- Lower memory overhead

Expected speedup: **10-50x faster** than the Python implementation.

## Building

### Prerequisites

- Go 1.26 or later
- Make (optional, for using Makefile commands)

### Build Commands

```bash
# Download dependencies
make deps

# Build the binary
make build

# Install to ~/.local/bin (no sudo)
make install

# Build for all platforms
make build-all
```

Manual build:
```bash
cd go
go mod download
go build -o bin/riskdays_go main.go
```

## Usage

### Command Line

The binary accepts JSON input via stdin or a file argument:

```bash
# From stdin
echo '{"k": 0.5, "doubling_time": 0.85, ...}' | ./bin/riskdays_go

# From file
./bin/riskdays_go input.json
```

### Input JSON Format

Baseline (window-period) model, with *k* drawn from an Inverse Gamma — a complete, runnable
payload (every field not shown takes its default):

```json
{
  "k": 0.000673,
  "doubling_time": 0.8542,
  "doubling_time_norm_sd": 0.0554,
  "lod50": 2.73,
  "lod50_sd": 0.193,
  "lod95_lod50_ratio": 4.51,
  "volume_transfused": 20,
  "volume_transfused_min": 15,
  "volume_transfused_max": 30,
  "pool_size": 16,
  "retests": 1,
  "k_invgamma_alpha": 2.0,
  "k_invgamma_beta": 0.002019,
  "n_bs": 10000,
  "seed": 126887,
  "threads": 7,
  "point_estimate": "median"
}
```

PrEP-breakthrough model — set `prep_mode` and supply the PrEP scalars (see
[PrEP-Breakthrough Parameters](#prep-breakthrough-parameters)):

```json
{
  "k": 0.000673,
  "doubling_time": 0.8542,
  "doubling_time_norm_sd": 0.0554,
  "lod50": 2.73,
  "lod50_sd": 0.193,
  "lod95_lod50_ratio": 4.51,
  "volume_transfused": 200,
  "volume_transfused_min": 100,
  "volume_transfused_max": 340,
  "pool_size": 16,
  "retests": 1,
  "k_invgamma_alpha": 2.0,
  "k_invgamma_beta": 0.002019,
  "prep_mode": true,
  "set_point": 336,
  "eclipse": 7.0,
  "a": 0.7,
  "b": 0.6,
  "offset": 1.0,
  "drug_effect": 1.0,
  "ser_min": 28.7,
  "ser_max": 250,
  "ser_alpha": 50.49434,
  "ser_beta": 1.15062,
  "n_bs": 10000,
  "seed": 126887,
  "threads": 7,
  "point_estimate": "median"
}
```

### Output JSON Format

Standard output:
```json
{
  "point_estimate": 1.75,
  "credible_interval": [0.85, 3.12],
  "range": [0.45, 5.67],
  "simulations": [1.2, 1.5, 1.8, ...]
}
```

Progress updates (stderr):
```json
{"type": "progress", "completed": 1000, "total": 10000, "percent": 0.1}
{"type": "progress", "completed": 2000, "total": 10000, "percent": 0.2}
```

Errors (stderr):
```json
{"type": "error", "message": "error description"}
```

## Python Integration

The webapp calls the Go binary via `residualrisk/_go.py` (part of the `residualrisk` package). The `use_go=True` flag is set by default in the Streamlit interface. The Python implementation remains available as a fallback.

```python
from residualrisk import risk_days_bs

# Default for the webapp: uses Go binary
result = risk_days_bs(k, doubling_time, ..., use_go=True)

# Fallback: pure Python (slow — avoid for n_bs > 1000)
result = risk_days_bs(k, doubling_time, ..., use_go=False)
```

The binary is located by `residualrisk.find_go_binary()`, which searches (in order): the `RESIDUALRISK_GO_BINARY` env var, `<repo>/go/bin/riskdays_go`, the platform binary bundled in the wheel (`residualrisk/_bin/riskdays_go-<goos>-<goarch>`, version-checked), `~/.local/bin/riskdays_go`, then `$PATH`.

## Parameters

### Required Parameters

- `k`: Infectivity parameter (float)
- `doubling_time`: Viral doubling time in days (float)
- `doubling_time_norm_sd`: Standard deviation for doubling time (float)
- `lod50`: Limit of detection 50% (copies/mL) (float)
- `lod50_sd`: Standard deviation for LOD50 (float)
- `lod95_lod50_ratio`: Ratio between LOD95 and LOD50 (float)
- `volume_transfused`: Average volume transfused in mL (float)
- `volume_transfused_min`: Minimum volume (float)
- `volume_transfused_max`: Maximum volume (float)
- `pool_size`: Minipool size for testing (int)
- `retests`: Number of retests (int)

### K Distribution (exactly one required)

The infectivity parameter *k* is drawn from **one** of four distributions. `Validate()`
**rejects a payload that specifies more than one**, and rejects a partially-specified one — a
distribution counts as specified if *any* of its fields is present, so `k_gamma_shape` without
`k_gamma_scale` is an error rather than silently falling through to another distribution.

- `k_posterior_sample`: Array of posterior draws for *k*, resampled with replacement (float[])
- `k_gamma_shape` **and** `k_gamma_scale`: Gamma distribution (float) — *legacy, deprecated*
- `k_invgamma_alpha` **and** `k_invgamma_beta`: Inverse Gamma, shape α and scale β (float).
  This is what the app uses by default. (The Python API's `k_invgamma_mode` alternative is
  converted to β before serialising, so the binary only ever sees α+β.)
- `k_lnmix_w`, `k_lnmix_mu1`, `k_lnmix_sigma1`, `k_lnmix_mu2`, `k_lnmix_sigma2`: two-component
  lognormal mixture — weight *w* on component 1, then each component's log-scale μ and σ (float).
  All five are required together.

### Optional Parameters (with defaults)

- `c0`: Initial viral concentration (default: 0.00025)
- `copies_per_virion`: RNA copies per virion, χ (default: 2)
- `alpha`: Significance level (default: 0.05)
- `z`: Z-score for test sensitivity (default: 1.6449)
- `limits`: Integration domain of the RDE integral, in days — `[lo, hi]` (default: `[-100, 500]`).
  Must be finite with `lo < hi`. The integrand is ≈0 well inside the default, which is a
  deliberately generous, numerically safe choice; widen it only for exotic parameterisations
  (e.g. a PrEP `ser_max` beyond 500 days, which the domain would otherwise silently truncate).
- `n_bs`: Number of simulations (default: 10000)
- `seed`: Random seed (default: 126887)
- `threads`: Number of parallel workers (default: CPU cores - 1)
- `point_estimate`: Method for point estimate - "primary parameters", "median", "mean", or "mode" (default: "primary parameters")
- `mode_precision`: Decimal precision for mode calculation (default: 2)
- `return_params`: Include the per-iteration sampled parameters in `simulations` (default: false)

### PrEP-Breakthrough Parameters

Set `prep_mode: true` to run the PrEP-breakthrough model (eclipse → exponential growth →
oscillating plateau, with NAT *and* serology detection) instead of the baseline window-period
model. The scalars below are then **required**: unlike the parameters above, they are
deliberately *not* defaulted from a `0` sentinel — a zero cannot be told apart from an omitted
field, and defaulting silently corrupted both a legitimate `0` (e.g. `a = 0`, a flat plateau)
and an invalid one (e.g. `drug_effect = 0`). `Validate()` rejects them if missing or degenerate.
The sampling *ranges* do have defaults, because their "unset" sentinel is the zero array, which
is distinct from a scalar 0.

- `prep_mode`: Enable the PrEP-breakthrough model (bool, default: false)
- `set_point`: Breakthrough plateau viral load, **RNA copies/mL** (float). Note the units: the
  model's concentration *C* is in virions/mL, so the engine divides by `copies_per_virion`.
- `eclipse`: Eclipse duration — delay before viral outgrowth begins, days (float)
- `a`: Plateau oscillation amplitude, as a fraction of the set-point (float). Must be `<= offset`,
  or the plateau viral load would go negative.
- `b`: Plateau oscillation angular frequency, rad/day (float)
- `offset`: Plateau oscillation centre, as a multiple of the set-point (float; **must be > 0**).
  `tcrit` targets this central level — `eclipse + doubling_time * log2((offset*set_point/χ)/C0)` —
  so the trajectory is continuous at the growth→plateau crossover. It is *exactly* a set-point
  multiplier: `(set_point, o, a)` is the same model as `(set_point*o, 1, a/o)`. Use 1 and vary
  `set_point` instead; the app does not expose it.
- `drug_effect`: Antiretroviral transmissibility-reduction factor, in `(0, 1]` (float; 1.0 = no
  reduction). A linear multiplier on the per-time infection probability.
- `ser_min`: Seroconversion onset — earliest serological detectability, days (float)
- `ser_max`: Seroconversion cutoff — serology certain beyond this, days (float; must exceed `ser_min`)
- `ser_alpha`: Weibull **scale** of the seroconversion-delay distribution, days (float)
- `ser_beta`: Weibull **shape** of the seroconversion-delay distribution (float)

Per-iteration sampling ranges (`[lo, hi]`, sampled `Uniform(lo, hi)`; omit to hold the scalar fixed):

- `set_point_dist_uniform` (default: `[19.1, 2265]`)
- `eclipse_dist_uniform` (default: `[4.0, 10.0]`)
- `a_dist_uniform`, `b_dist_uniform` (no default — `a` and `b` are held fixed unless given; the
  upper bound of `a_dist_uniform` must not exceed `offset`)
- `drug_effect_dist_uniform` (no default — held fixed unless given; must lie within `(0, 1]`)

## Performance

Typical performance on Apple M1 (8 cores):

- 10,000 simulations: ~5-10 seconds
- 25,000 simulations: ~12-25 seconds
- 100,000 simulations: ~50-100 seconds

Compare to Python (single core): 10,000 simulations can take 5-15 minutes.

## Testing

```bash
# Run all tests
make test

# Run with coverage
make coverage
```

## Architecture

```
main.go              - CLI interface, JSON I/O
riskdays/
  models.go          - Input/output data structures
  riskdays.go        - Main bootstrap orchestration
  integration.go     - Numerical integration (quad)
  probability.go     - Probability calculations
  helpers.go         - Utility functions
  random.go          - Random sampling
```

## Dependencies

- `gonum.org/v1/gonum` - Scientific computing (stats, integration, distributions)

## License

See parent project license.
