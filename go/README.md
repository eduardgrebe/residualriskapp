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

- Go 1.21 or later
- Make (optional, for using Makefile commands)

### Build Commands

```bash
# Download dependencies
make deps

# Build the binary
make build

# Install to /usr/local/bin
sudo make install

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

```json
{
  "k": 0.5,
  "doubling_time": 0.8542,
  "doubling_time_norm_sd": 0.0663,
  "lod50": 2.73,
  "lod50_sd": 0.193,
  "lod95_lod50_ratio": 4.51,
  "volume_transfused": 20,
  "volume_transfused_min": 15,
  "volume_transfused_max": 30,
  "pool_size": 16,
  "retests": 1,
  "c0": 0.00025,
  "copies_per_virion": 2,
  "alpha": 0.05,
  "z": 1.6449,
  "k_posterior_sample": [0.45, 0.52, 0.48, ...],
  "n_bs": 10000,
  "seed": 126887,
  "threads": 7,
  "point_estimate": "primary parameters",
  "mode_precision": 2
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

The webapp calls the Go binary via `app/residualrisk_go.py`. The `use_go=True` flag is set by default in the Streamlit interface. The Python implementation remains available as a fallback.

```python
from residualrisk import risk_days_bs

# Default for the webapp: uses Go binary
result = risk_days_bs(k, doubling_time, ..., use_go=True)

# Fallback: pure Python (slow — avoid for n_bs > 1000)
result = risk_days_bs(k, doubling_time, ..., use_go=False)
```

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

### K Distribution (one required)

- `k_posterior_sample`: Array of posterior samples for k (float[])
- OR `k_gamma_shape` AND `k_gamma_scale`: Gamma distribution parameters (float)

### Optional Parameters (with defaults)

- `c0`: Initial viral concentration (default: 0.00025)
- `copies_per_virion`: RNA copies per virion (default: 2)
- `alpha`: Significance level (default: 0.05)
- `z`: Z-score for test sensitivity (default: 1.6449)
- `n_bs`: Number of simulations (default: 10000)
- `seed`: Random seed (default: 126887)
- `threads`: Number of parallel workers (default: CPU cores - 1)
- `point_estimate`: Method for point estimate - "primary parameters", "median", "mean", or "mode" (default: "primary parameters")
- `mode_precision`: Decimal precision for mode calculation (default: 2)

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
