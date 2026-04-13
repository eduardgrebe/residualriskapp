# Claude Code Instructions for residualriskapp

## Project Overview

**Residual HIV Transfusion Transmission Risk Estimation Tool**

A production-ready Streamlit web application and simulation framework for estimating the residual risk of HIV transmission through blood transfusion during the pre-NAT window period (or PrEP breakthrough infection).

- **License**: GNU Affero General Public License v3.0
- **Author**: Eduard Grebe <egrebe@vitalant.org>
- **Institution**: Vitalant Research Institute
- **Python Version**: >=3.14, <3.15
- **Package Manager**: uv

## Repository Structure

```
residualriskapp/
├── README.md              # User-facing documentation
├── AGENTS.md              # This file — canonical project instructions for agents
├── LICENSE                # GNU AGPL v3.0
├── app.py                 # Streamlit web application (entry point)
├── residualrisk.py        # Core calculation engine
├── residualrisk_go.py     # Python wrapper for Go binary
├── residualrisk_prep.py   # PrEP model (incomplete, not yet wired into app)
├── main.py                # CLI entry point (placeholder)
├── go/                    # High-performance Go reimplementation (10-50x faster)
│   ├── main.go            # CLI interface with JSON I/O
│   └── riskdays/          # Core Go package
│       ├── models.go      # Data structures
│       ├── riskdays.go    # Bootstrap orchestration
│       ├── integration.go # Numerical integration
│       ├── probability.go # Probability calculations
│       ├── helpers.go     # Utility functions
│       ├── random.go      # Random sampling
│       └── riskdays_test.go
├── static/                # Pre-computed Bayesian posterior distributions (Parquet)
├── tests/                 # Python test suite
├── docker/                # Docker build and deployment scripts
├── .streamlit/            # Streamlit configuration
├── pyproject.toml         # Python project config (managed by uv)
└── uv.lock                # Locked dependency versions (do not edit manually)
```

## Core Application Files

- **`app.py`** — Streamlit web UI
  - Parameter input interface with expandable sections
  - Supports three RDE estimation methods: Mechanistic model, Lookback data, Mechanistic model with PrEP
  - Real-time calculation and result visualization
  - Entry point: `streamlit run app.py`

- **`residualrisk.py`** — Main calculation engine
  - Viral concentration dynamics
  - Infectivity probability calculations
  - Bootstrap simulation methods
  - Integration with Go acceleration via `residualrisk_go.py`

- **`residualrisk_go.py`** — Go binary wrapper
  - JSON-based communication with Go binary
  - Automatic fallback to Python if Go binary unavailable
  - Progress monitoring

- **`residualrisk_prep.py`** — PrEP model (incomplete)
  - PrEP viral dynamics and Monte Carlo simulation
  - Not yet invoked by the webapp (UI stubs exist in `app.py`)

## Technical Architecture

### Risk Estimation Model

1. **Viral Dynamics**: Exponential growth from initial concentration (C0) with doubling time
2. **Test Sensitivity**: Incorporates LOD (limit of detection) with uncertainty
3. **Infectivity**: Probabilistic model using posterior-sampled k parameter
4. **Bootstrap Simulation**: Monte Carlo sampling of parameter uncertainty
5. **Window Period Calculation**: Numerical integration to find infectious window period

### Key Parameters

**Viral Growth**:
- `C0` — Initial viral concentration (default: 0.00025 copies/mL)
- `doubling_time` — Viral doubling time in days
- `doubling_time_norm_sd` — Uncertainty in doubling time

**Test Characteristics**:
- `lod50` — 50% limit of detection (copies/mL)
- `lod50_sd` — Standard deviation of LOD50
- `lod95_lod50_ratio` — Ratio between 95% and 50% LOD
- `z` — Z-score for test sensitivity (default: 1.6449)
- `pool_size` — Minipool size for NAT testing
- `retests` — Number of retests performed

**Transmission**:
- `k` — Infectivity parameter (sampled from posterior)
- `volume_transfused` — Volume of blood transfused (mL)
- `volume_transfused_min/max` — Uncertainty range
- `copies_per_virion` — RNA copies per virion (default: 2)

**Simulation**:
- `n_bs` — Number of bootstrap simulations (default: 10,000)
- `seed` — Random seed for reproducibility
- `threads` — Parallel workers (default: CPU cores - 1)

### Static Data (`static/`)

Pre-computed posterior parameter distributions in Parquet format:
- `k_param_human.parquet` — Human infectivity parameter posterior
- `k_param_animal.parquet` — Animal model infectivity posterior
- `k_param_expdecay.parquet` — Exponential decay posterior
- `newnewdist_k_param.parquet` — Latest distribution samples
- `weibull_min_k_param.parquet` — Weibull minimum posterior
- `iwp_estimates_expdecay.parquet` — Infectious window period estimates

**Do not modify these files** — they are large pre-computed Bayesian posteriors (up to 8MB). Regeneration requires rerunning the upstream Bayesian analyses.

## PrEP Model Status

- `residualrisk_prep.py` contains a Python PrEP breakthrough risk model — **not yet wired into the webapp**
- `app.py` has parameter UI stubs for PrEP inputs under the "Mechanistic model with PrEP" dropdown
- Go implementation of PrEP model **does not exist yet**
- Plan: release initial version without PrEP, then add it alongside a manuscript

## Development Workflow

### Environment Setup

```bash
# Install uv if needed: https://github.com/astral-sh/uv
uv sync                 # Install/update dependencies
source .venv/bin/activate

uv add package-name     # Add a new dependency
```

### Running the Application

```bash
streamlit run app.py    # → http://localhost:8501
```

### Building the Go Implementation

The webapp defaults to the Go binary. Without it, it falls back to Python (10-50x slower, impractical for normal use).

```bash
cd go
make deps      # Download dependencies
make build     # Build binary to go/bin/riskdays_go
make test      # Run Go tests
```

The Python code automatically detects the binary at `go/bin/riskdays_go`.

### Testing

```bash
# Python
pytest tests/

# Go
cd go && make test
```

### Code Style

- All new Python files must include the AGPL v3.0 license header (copy from an existing file)
- Use type hints where practical
- Document complex calculations with references to literature/methodology
- Follow existing naming conventions

### Dependencies

**Core Python Stack** (see `pyproject.toml` for pinned versions):
- `streamlit` — Web application framework
- `pandas` — Data manipulation
- `numpy` — Numerical computing
- `scipy` — Scientific computing (stats, integration)
- `pyarrow` — Parquet file I/O
- `matplotlib`, `seaborn`, `plotly` — Visualization
- `watchdog` — File watching (Streamlit hot reload)

**Dev Dependencies**:
- `ruff` — Linter/formatter
- `pytest` — Test runner

**Go Dependencies**:
- `gonum.org/v1/gonum` — Scientific computing library

## Git Workflow

### Remote Access Limitation

**Agents do not have access to the SSH key** required to interact with the remote. Agents **cannot** execute git commands that communicate with the remote, including `git fetch`, `git pull`, `git push`. For any task involving remote interaction, provide the user with the exact commands to run. For example:

```bash
git fetch origin main
git rebase origin/main
git push --force-with-lease origin <branch-name>
```

Agents **can** safely run local-only git commands (`git status`, `git log`, `git diff`, `git add`, `git commit`, local `git rebase`).

## Common Tasks

### Adding a New Parameter

1. Update calculation functions in `residualrisk.py`
2. Add UI controls in `app.py`
3. Update Go implementation in `go/riskdays/` if needed
4. Update JSON schema in `go/riskdays/models.go`
5. Document in this file and in `README.md`

### Modifying the Model

1. Update mathematical functions in `residualrisk.py`
2. Verify numerical stability with test cases
3. Update Go implementation for consistency
4. Document methodology changes
5. Consider impact on existing posterior distributions in `static/`

### Updating Dependencies

```bash
uv add package-name@version          # Python
cd go && go get package@version && go mod tidy   # Go
```

## Important Notes

### License Compliance

This is AGPL v3.0 licensed software. Any modifications must:
- Include the license header in new files
- Maintain license notice when distributing
- Provide source code access for network users (AGPL network copyleft)

### Numerical Precision

- Calculations involve very small probabilities and concentrations
- Use appropriate numerical methods (log space where needed)
- Integration tolerances are calibrated for epidemiological accuracy

### Performance Considerations

- Python implementation: Single-core, suitable for interactive use only
- Go implementation: Multi-core, 10-50x faster, required for practical use
- For `n_bs > 25,000`, strongly prefer Go implementation

## References

- Fiebig, E.W., et al. (2003). Dynamics of HIV viremia and antibody seroconversion in plasma donors. *AIDS*, 17(13):1871-1879. doi:[10.1097/00002030-200309050-00005](https://doi.org/10.1097/00002030-200309050-00005).
- Weusten J., et al. (2011) Refinement of a viral transmission risk model for blood donations in seroconversion window phase screened by nucleic acid testing. *Transfusion*, 51(1):203-15. doi:[10.1111/j.1537-2995.2010.02804.x](https://doi.org/10.1111/j.1537-2995.2010.02804.x).
- Grebe E., et al. (2020) HIV incidence in US first-time blood donors and transfusion risk with a 12-month deferral for men who have sex with men. *Blood*, 136(11):1359-1367. doi:[10.1182/blood.2020007003](https://doi.org/10.1182/blood.2020007003).
- Belov A., et al. (2023) Modeling the Risk of HIV Transfusion Transmission. *J Acquir Immune Defic Syndr*, 92(2):173-179. doi:[10.1097/QAI.0000000000003115](https://doi.org/10.1097/QAI.0000000000003115).

