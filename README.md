# Residual Risk Estimation Tool

A simulation tool for estimating the residual risk of HIV transmission through blood transfusion during early infection when screening tests may not detect infection, either
* during the pre-NAT window period, or
* as a result of PrEP breakthrough infection

## Overview

This tool implements a Monte Carlo simulation model to estimate the infectious window period and residual transmission risk for HIV in blood donations. The primary window period model accounts for:

- Viral kinetics during early infection (exponential growth from initial concentration)
- Nucleic acid testing (NAT) sensitivity characteristics and limit of detection
- Minipool testing and retesting protocols
- Uncertainty in infectivity parameters using Bayesian posterior distributions from Belov et al. (2023)
- Variable transfusion volumes

The tool provides both a Python implementation and a high-performance Go implementation for large-scale parameter space exploration.

A Streamlit-based interactive web interface is provided that can make use of either the Python or Go implementations.

## Features

- **Interactive Web Interface**: Streamlit-based UI for parameter exploration and visualization
- **High-Performance Computation**: Go implementation provides 10-50x speedup over the Python fallback and is required for practical use
- **Flexible Parameterization**: Supports various NAT assays, pooling strategies, and viral kinetics models
- **Flexible k Input Distribution**: k can be sampled from posterior draws (human, animal, human-weighted) or a parametric Inverse Gamma distribution with user-specified α and β (or α and mode)
- **Credible Interval Estimation**: Bootstrap-based credible intervals for risk estimates

## Requirements

### For Using the Web Interface

- **Python**: 3.14.x
- **Package Manager**: [uv](https://github.com/astral-sh/uv) (recommended) or pip
- **Go**: 1.26+ (required to build the Go binary; see below)
- **Make**: (optional, for using build shortcuts)

### For Development

- **Python**: 3.14.x
- **uv**: For dependency management
- **Go**: 1.26+
- **Make**: (optional, for using build shortcuts)

## Installation

### Quick Start (Web Interface)

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd residualriskapp
   ```

2. **Install dependencies using uv** (recommended):
   ```bash
   # Install uv if you don't have it
   curl -LsSf https://astral.sh/uv/install.sh | sh

   # Create virtual environment and install dependencies
   uv sync

   # Activate the environment
   source .venv/bin/activate
   ```

   **Or using pip**:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   pip install -e .
   ```

3. **Run the web interface**:
   ```bash
   streamlit run app.py
   ```

4. **Open your browser** to http://localhost:8501

### Building the Go Implementation

The webapp defaults to the Go implementation. Without the binary, it falls back to the Python implementation, which is 10-50x slower and impractical for the simulation counts used in normal operation.

1. **Install Go** 1.26.0 or later from https://go.dev/dl/

2. **Build the Go binary** (from the repo root):
   ```bash
   bash scripts/build_go.sh
   ```
   Or directly:
   ```bash
   cd go
   make deps    # Download Go dependencies
   make build   # Build binary to go/bin/riskdays_go
   ```

3. **Verify the build**:
   ```bash
   ls -lh go/bin/riskdays_go
   ```

The Python code automatically detects the binary at `<repo>/go/bin/riskdays_go`. To point at a binary installed elsewhere, set `RESIDUALRISK_GO_BINARY=/absolute/path/to/riskdays_go`.

### Docker Deployment (Recommended for Production)

For production deployments or if you want to avoid installing dependencies locally, use Docker:

**Quick Start:**
```bash
# Using Docker Compose (easiest - builds and runs)
docker-compose up -d

# Or build and run separately
./docker/build.sh --load      # Auto-detects your architecture
./docker/run.sh
```

The container binds only to `127.0.0.1:8501` and does not handle SSL. For production, place a reverse proxy (nginx, Caddy, Traefik, etc.) in front of the container to handle SSL termination. A reference nginx configuration is provided at `docker/nginx/conf.d/app.conf`.

**Multi-Architecture Support:**

```bash
# Local build (auto-detects your architecture)
./docker/build.sh --load

# Build and push multi-arch to registry
./docker/build.sh --registry docker.io/username --push
```

The Docker image includes:
- Python 3.14 with all dependencies (managed by uv)
- Pre-compiled Go binary for high-performance computation
- Streamlit configured for production use
- Multi-architecture support (AMD64/ARM64)
- Non-root user for security
- Health checks and proper signal handling

**See [docker/README.md](docker/README.md) for complete Docker documentation including:**
- Multi-architecture builds
- Production deployment with reverse proxy
- Resource management and security
- Troubleshooting

Also available: [Quick Reference Guide](docker/QUICKREF.md)

## Usage

### Web Interface

Launch the Streamlit app:

```bash
streamlit run app.py
```

The interface allows you to:
- Adjust all model parameters interactively
- Select infectivity parameter (k) distribution: posterior samples (human, animal, human-weighted exponential-decay) or Inverse Gamma with user-specified α and β
- Choose the k point estimate summary (mode, median, or mean) used for the IWP point estimate
- Choose between Python and Go computation engines
- Visualize results with credible intervals
- Export results and simulation data

### Command Line (Go)

The Go implementation accepts JSON input:

```bash
# From file
./go/bin/riskdays_go input.json

# From stdin
echo '{
  "k": 0.5,
  "doubling_time": 0.85,
  "lod50": 2.73,
  "pool_size": 16,
  "n_bs": 10000
}' | ./go/bin/riskdays_go
```

See `go/README.md` for detailed documentation of the JSON schema and parameters.

### Python API

`residualrisk` is a proper installable package. Install it into the environment of any downstream analysis with `uv pip install -e /path/to/residualriskapp` (editable) or pin to a git tag/SHA for reproducibility.

```python
import residualrisk as rr

# Bootstrap risk-day equivalents (IWP) — using a posterior sample for k
rd_pe, rd_cri, rd_range, rdests = rr.risk_days_bs(
    k=0.013,
    doubling_time=20.5 / 24,
    doubling_time_norm_sd=1.33 / 24,
    lod50=2.73,
    lod50_sd=0.193,
    lod95_lod50_ratio=12.33 / 2.73,
    volume_transfused=20,
    volume_transfused_range=(15, 30),
    pool_size=16,
    retests=1,
    k_posterior_sample=k_samples,    # numpy array of posterior draws
    n_bs=10000,
    use_go=True,                      # use Go acceleration (10-50x faster)
)

# Alternative: sample k from an Inverse Gamma distribution (α=2, β=0.002019)
# k_pe can be the mode (β/(α+1)), median, or mean (β/(α-1)) of the distribution
rd_pe, rd_cri, rd_range, rdests = rr.risk_days_bs(
    k=0.002019 / 3,                  # mode = β/(α+1) = 0.002019/3
    doubling_time=20.5 / 24,
    doubling_time_norm_sd=1.33 / 24,
    lod50=2.73,
    lod50_sd=0.193,
    lod95_lod50_ratio=12.33 / 2.73,
    volume_transfused=20,
    volume_transfused_range=(15, 30),
    pool_size=16,
    retests=1,
    k_invgamma_alpha=2.0,
    k_invgamma_beta=0.002019,
    n_bs=10000,
    use_go=True,
)
print(f"RDEs point estimate: {rd_pe:.2f} days")
print(f"95% CrI: [{rd_cri[0]:.2f}, {rd_cri[1]:.2f}]")

# Combine RDEs with incidence to get residual risk
rr_pe, rr_cri, rr_sd = rr.residual_risk_rd(
    iwp_pe=rd_pe,
    iwp_bs=rdests,
    incidence=2.5 / 1e5,              # per person-year
    incidence_norm_sd=0.5 / 1e5,
    per=1e6,                          # report per 1 million transfusions
)
print(f"Residual risk: {rr_pe:.3f} per million (95% CrI {rr_cri[0]:.3f}–{rr_cri[1]:.3f})")

# Record provenance alongside outputs
print(f"residualrisk version: {rr.__version__}")
```

**Public API surface** (see `residualrisk/__init__.py`): `risk_days_bs`, `iwp_from_lookback_data`, `residual_risk_rd`, `get_cpu_core_count`, `mode_rounded`, `mode_kde`, `sample_invgamma`, `find_go_binary`, `__version__`.

### R integration (reticulate)

```r
library(reticulate)
use_virtualenv("/path/to/residualriskapp/.venv", required = TRUE)
rr <- import("residualrisk")
bs <- rr$risk_days_bs(...)           # returns a Python tuple; index with [[1]], [[2]], ...
```

## Dependencies

### Python (Core)

- `streamlit` - Web application framework
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `scipy` - Statistical functions and numerical integration
- `matplotlib` - Static plotting
- `seaborn` - Statistical visualization
- `plotly` - Interactive plotting
- `watchdog` - File monitoring for Streamlit

### Go

- `gonum.org/v1/gonum` - Scientific computing (statistics, integration, distributions)

## Model Description

The model estimates the **infectious window period (IWP)**: the time interval during which a donation contains infectious virus but falls below the NAT detection threshold. The IWP can be multiplied with HIV incidence to obtain the residual risk of HIV transfusion transmission.

### Key Steps

1. **Viral Growth**: Concentration increases exponentially: C(t) = C₀ × 2^(t/doubling_time)
2. **Detection Probability**: Based on LOD characteristics and pooling/retesting protocol
3. **Infectivity Probability**: P(infection) = 1 - exp(-k × n_copies), where n_copies depends on viral load and transfusion volume. k is sampled each bootstrap iteration from the chosen input distribution: a posterior sample array or a parametric Inverse Gamma(α, β).

The choice of input distribution for *k* is discussed in detail in the companion
analysis repository. See [`residualrisk_analysis/exploration/K_PARAM_INPUTDIST.md`](../residualrisk_analysis/exploration/K_PARAM_INPUTDIST.md)
for a systematic comparison of candidate distributions (inverse gamma, lognormal,
log-logistic, lognormal mixture) and recommendations with full justification.
4. **Bootstrap Simulation**: Samples parameter uncertainty to generate IWP distribution
5. **Risk Estimation**: Point estimate and credible intervals from bootstrap distribution

### Point Estimate Options

- **Primary parameters**: Uses mode/mean parameter values (default)
- **Median**: Median of bootstrap distribution
- **Mean**: Mean of bootstrap distribution
- **Mode**: Mode of bootstrap distribution

## Performance

### Go Implementation (default)
- Multi-threaded (uses N-1 CPU cores by default)
- 10,000 simulations: ~5-10 seconds
- 100,000 simulations: ~50-100 seconds
- **10-50x faster** than the Python fallback

### Python Implementation (fallback only)
- Single-threaded
- 10,000 simulations: ~5-15 minutes
- Available as a fallback if the Go binary is not built; not suitable for normal use

## Project Structure

```
residualriskapp/
├── app.py                   # Streamlit web application (imports the residualrisk package)
├── residualrisk/            # Installable Python package (core calculation engine)
│   ├── __init__.py          # Public API surface
│   ├── core.py              # Core calculation engine (bootstrap, integration, IWP)
│   └── _go.py               # Wrapper around the Go binary
├── static/                  # Pre-computed posterior distributions (Parquet)
│   ├── k_param_human.parquet
│   ├── k_param_animal.parquet
│   └── ...
├── go/                      # High-performance Go implementation
│   ├── main.go
│   ├── riskdays/            # Core Go package
│   └── README.md            # Go-specific documentation
├── scripts/
│   └── build_go.sh          # One-command Go binary build
├── tests/                   # Python test suite (targets `residualrisk.core`)
├── docker/                  # Docker build and deployment scripts
├── pyproject.toml           # Python project configuration (hatchling build backend)
└── .venv/                   # Virtual environment (created by uv sync)
```

## License

Copyright (C) 2025 Vitalant and Eduard Grebe Consulting
Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

## Author

**Eduard Grebe**
- Email: egrebe@vitalant.org
- Email: eduard@grebe.consulting
- Institution: Vitalant Research Institute

## Citation

If you use this tool in your research, please cite:

[Citation information to be added]

## References

The model makes use of established methodology and novel approaches for HIV transfusion transmission risk estimation. Key literature informing the model includes:

- Fiebig, E.W., et al. (2003). Dynamics of HIV viremia and antibody seroconversion in plasma donors: implications for diagnosis and staging of primary HIV infection. *AIDS*, 17(13):1871-1879. doi:[10.1097/00002030-200309050-00005](https://doi.org/10.1097/00002030-200309050-00005).

- Weusten J., et al. (2011) Refinement of a viral transmission risk model for blood donations in seroconversion window phase screened by nucleic acid testing in different pool sizes and repeat test algorithms. *Transfusion*, 51(1):203-15. doi:[10.1111/j.1537-2995.2010.02804.x](https://doi.org/10.1111/j.1537-2995.2010.02804.x).

- Grebe E., et al. (2020) HIV incidence in US first-time blood donors and transfusion risk with a 12-month deferral for men who have sex with men. *Blood*, 136(11):1359-1367. doi:[10.1182/blood.2020007003](https://doi.org/10.1182/blood.2020007003).

- Belov A., et al. (2023) Modeling the Risk of HIV Transfusion Transmission. *J Acquir Immune Defic Syndr*, 92(2):173-179. doi:[10.1097/QAI.0000000000003115](https://doi.org/10.1097/QAI.0000000000003115).

## Support and Contributions

For questions, bug reports, or feature requests:
- Email: egrebe@vitalant.org
- GitHub Issues: [Add once repo becomes public]

## Development

### Running Tests

```bash
# Python
pytest tests/

# Go
cd go
make test
```

The Python test suite covers the core calculation functions in `residualrisk/core.py`. It imports via `from residualrisk import core as rr` so the package must be installed in the environment first — `uv sync` (or `uv pip install -e .`) handles this. Tests that exercise the bootstrap simulation (both Python and Go implementations) require the Go binary to be built first — see "Building the Go Implementation" above.

### Adding Dependencies

```bash
# Python
uv add package-name

# Go
cd go
go get package-name
go mod tidy
```

## Acknowledgments

Developed at Vitalant Research Institute for blood safety research.

Authors:
- Eduard Grebe

Contributors:
- Brian Custer
- Michael P. Busch

