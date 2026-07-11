# Residual Risk Estimation Tool

A simulation tool for estimating the residual risk of HIV transmission through blood transfusion during early infection when screening tests may not detect infection, either
* during the pre-NAT window period, or
* as a result of PrEP breakthrough infection

## Overview

A simulation tool for estimating the residual risk of HIV transmission through blood transfusion during the pre-NAT window period (or due to PrEP breakthrough infection). Accounts for viral kinetics during early infection, NAT sensitivity, pooling strategy, and infectivity using the Belov et al. (2023) model. Provides a Streamlit-based web interface, a Python API and a high-performance Go implementation of the core simulation.

The tool provides both a Python and a high-performance Go implementation, wrapped in a Python package for ease of use.

A Streamlit-based interactive web interface is provided that can make use of either the Python or Go implementations.

## Features

- **Interactive Web Interface**: Streamlit-based UI for parameter exploration and visualization
- **High-Performance Computation**: Go implementation provides 10-50x speedup over the Python fallback and is required for practical use
- **PrEP-Breakthrough Modeling**: Optional oral (oPrEP) and injectable (iPrEP) breakthrough-infection risk components, layered additively on the baseline window-period risk
- **Canned NAT-Assay Presets**: Built-in published 50%/95% limits of detection (HIV-1 Group M) for seven blood-screening NAT assays, selectable by name (or enter LoDs manually)
- **Flexible Parameterization**: Supports various NAT assays, pooling strategies, and viral kinetics models
- **Flexible k Input Distribution**: k can be sampled from posterior draws (human, animal, human-weighted) or a parametric distribution: Inverse Gamma(α, β) or a two-component lognormal mixture (90% human + 10% animal by default)
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

### Installing the Python library (with bundled Go acceleration)

The `residualrisk` library is published to this project's **Codeberg PyPI registry**, and its
wheels **bundle the pre-compiled Go accelerator for every platform** (Linux / macOS / Windows,
x86-64 + arm64). So `pip install` gives you the fast Go engine with **no Go toolchain and no
separate build** — and if a platform's binary is ever unusable, the library transparently falls
back to the pure-Python engine.

> **Currently pre-release:** only pre-release builds are on the registry so far, so a bare
> `residualrisk` won't resolve until a stable `1.1.0` is published. Allow pre-releases — pip: add
> `--pre`; uv: add `--prerelease=allow` (works with `uv pip install` and `uv sync`) — or pin an
> exact version (e.g. `residualrisk==1.1.0b3`).

**pip** — the registry hosts only `residualrisk`, so add it *alongside* PyPI (dependencies such as
numpy/scipy still come from PyPI). Use `--extra-index-url`, **not** `--index-url` (the latter
replaces PyPI entirely and breaks dependency resolution):

```bash
pip install --extra-index-url https://codeberg.org/api/packages/eduardgrebe/pypi/simple/ residualrisk
```

**uv** (quick):

```bash
uv pip install --extra-index-url https://codeberg.org/api/packages/eduardgrebe/pypi/simple/ residualrisk
```

**uv, in a project's `pyproject.toml`** — declare the dependency, define the registry as a named
index, and pin `residualrisk` to it (everything else still resolves from PyPI):

```toml
[project]
dependencies = ["residualrisk"]

[[tool.uv.index]]
name = "codeberg"
url = "https://codeberg.org/api/packages/eduardgrebe/pypi/simple/"

[tool.uv.sources]
residualrisk = { index = "codeberg" }
```

Then `uv sync` (or `uv add residualrisk`) resolves it from Codeberg. If the registry is ever made
private, supply a Codeberg token via the index's auth (uv `UV_INDEX_CODEBERG_USERNAME` /
`UV_INDEX_CODEBERG_PASSWORD` env vars, or a `~/.netrc` entry for `codeberg.org`).

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
- Select infectivity parameter (k) distribution: posterior samples (human, animal, human-weighted exponential-decay), Inverse Gamma with user-specified α and β, or a two-component lognormal mixture with adjustable mixing weight
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
  "k": 0.000673,
  "doubling_time": 0.8542,
  "doubling_time_norm_sd": 0.2813,
  "lod50": 2.73, "lod50_sd": 0.53, "lod95_lod50_ratio": 3.5,
  "volume_transfused": 200, "volume_transfused_min": 100, "volume_transfused_max": 340,
  "pool_size": 16, "retests": 1,
  "k_invgamma_alpha": 2.0, "k_invgamma_beta": 0.002019,
  "n_bs": 10000
}' | ./go/bin/riskdays_go
```

See `go/README.md` for detailed documentation of the JSON schema and parameters.

> **Note:** the Go CLI takes explicit LoDs (`lod50` / `lod50_sd` / `lod95_lod50_ratio`).
> The canned-NAT-assay shortcut is Python-side only — use `risk_days_bs(assay="…")`
> (which resolves the preset to LoDs before calling Go); there is no `"assay"` JSON field.

### Python API

`residualrisk` is a proper installable package. Install it into the environment of any downstream analysis with `uv pip install -e /path/to/residualriskapp` (editable) or pin to a git tag/SHA for reproducibility.

```python
import residualrisk as rr

# Bootstrap risk-day equivalents (IWP) — using a posterior sample for k
rd_pe, rd_cri, rd_range, rdests, _ = rr.risk_days_bs(
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
    k_posterior_sample=k_samples,  # numpy array of posterior draws
    n_bs=10000,
    use_go=True,  # use Go acceleration (10-50x faster)
)

# Shortcut: use a canned NAT assay's published LoDs instead of the explicit
# lod50 / lod50_sd / lod95_lod50_ratio triplet (mutually exclusive with them).
rr.list_assays()  # {'ultrio': 'Procleix Ultrio (Tigris)', ...}
rd_pe, rd_cri, rd_range, rdests, _ = rr.risk_days_bs(
    k=0.013,
    doubling_time=20.5 / 24,
    doubling_time_norm_sd=1.33 / 24,
    assay="ultrio_elite",  # in place of lod50/lod50_sd/lod95_lod50_ratio
    volume_transfused=20,
    volume_transfused_range=(15, 30),
    pool_size=16,
    retests=1,
    k_posterior_sample=k_samples,
    n_bs=10000,
    use_go=True,
)
# Inspect a preset's numbers (+ provenance) without running a simulation:
lod = rr.lods_for_assay("ultrio_elite")
# lod.lod50, lod.lod50_sd, lod.lod95_lod50_ratio, lod.cp_per_iu, lod.iu_std

# Alternative: sample k from an Inverse Gamma distribution (α=2, β=0.002019)
# k_pe can be the mode (β/(α+1)), median, or mean (β/(α-1)) of the distribution
rd_pe, rd_cri, rd_range, rdests, _ = rr.risk_days_bs(
    k=0.002019 / 3,  # mode = β/(α+1) = 0.002019/3
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

# Alternative: sample k from a lognormal mixture (Recommendation B: 90% human + 10% animal)
rd_pe, rd_cri, rd_range, rdests, _ = rr.risk_days_bs(
    k=0.000649,  # approximate mixture mode
    doubling_time=20.5 / 24,
    doubling_time_norm_sd=1.33 / 24,
    lod50=2.73,
    lod50_sd=0.193,
    lod95_lod50_ratio=12.33 / 2.73,
    volume_transfused=20,
    volume_transfused_range=(15, 30),
    pool_size=16,
    retests=1,
    k_lnmix_w=0.90,
    k_lnmix_mu1=-7.2403,
    k_lnmix_sigma1=0.3241,
    k_lnmix_mu2=-3.7423,
    k_lnmix_sigma2=0.5258,
    n_bs=10000,
    use_go=True,
)

# Combine RDEs with incidence to get residual risk
rr_pe, rr_cri, rr_sd = rr.residual_risk_rd(
    iwp_pe=rd_pe,
    iwp_bs=rdests,
    incidence=2.5 / 1e5,  # per person-year
    incidence_norm_sd=0.5 / 1e5,
    per=1e6,  # report per 1 million transfusions
)
print(
    f"Residual risk: {rr_pe:.3f} per million (95% CrI {rr_cri[0]:.3f}–{rr_cri[1]:.3f})"
)

# Record provenance alongside outputs
print(f"residualrisk version: {rr.__version__}")
```

**Public API surface** (see `residualrisk/__init__.py`): `risk_days_bs`, `risk_days_prep_bs`, `iwp_from_lookback_data`, `residual_risk_rd`, `total_residual_risk_rd`, `get_cpu_core_count`, `mode_rounded`, `mode_kde`, `sample_invgamma`, `sample_lnmix`, `NAT_ASSAYS`, `lods_for_assay`, `list_assays`, `AssayLoD`, `find_go_binary`, `mode_hsm_go`, `mode_kde_go`, `risk_days_prep_bs_go`, `__version__`.

> **Canned assays — provisional Bio-Manguinhos SD.** `NAT_ASSAYS` carries published 50%/95% LoDs (HIV-1 Group M, copies/mL) for the supported assays. One caveat: the `biomanguinhos` (Brazilian NAT Platform, Bio-Manguinhos) entry has **no published LoD50 confidence interval** — Rocha et al. (2018) report point estimates only — so its `lod50_sd` is an *assumed* relative SD of **13%** (6.08 IU/mL ≈ 3.527 copies/mL). A value of **4.95 IU/mL (RSE 10.6%) was used in prior analyses.** This is provisional pending the per-dilution hit-rate table; see the note above `NAT_ASSAYS` in `residualrisk/assays.py`.

#### LoD50 relative standard error (RSE) by assay

RSE = `lod50_sd / lod50` (the coefficient of variation of the 50% LoD; invariant under the IU/mL→copies/mL conversion). For every assay except Bio-Manguinhos the SD derives from a 95% CI of the 50% LoD; Bio-Manguinhos uses an *assumed* RSE (see above). Underlying LoD/CI data are compiled in the companion analysis `residualrisk_analysis/assays/ASSAYS.qmd`.

| Assay | RSE on LoD50 | Source |
|---|---|---|
| `cobas_taqscreen_mpxv2` | 4.72% | Probit fit of Roche insert reactivity data (95% CI) |
| `cobas_mpx` | 6.04% | Roche cobas MPX CE/IVD insert (95% CI of 50% LoD) |
| `cobas_taqscreen_mpx` | 7.00% | Probit fit of Roche insert reactivity data (95% CI) |
| `ultrio_plus` | 7.07% | Grifols Procleix Ultrio Plus insert (95% CI of 50% LoD) |
| `ultrio` | 7.28% | Grifols Procleix Ultrio insert, dHIV-1 (95% CI of 50% LoD) |
| `ultrio_elite` | 7.55% | Grifols Procleix Ultrio Elite insert (95% CI of 50% LoD) |
| `biomanguinhos` | 13.00% (assumed, provisional) | Assumed RSE — Rocha et al. (2018) report point LoDs only (no CI) - ballparked from 24 replicates/dilution |

### R integration (reticulate)

Install the `residualrisk` Python package in a virtual environment if not installed:

```bash
uv venv
uv pip install --extra-index-url https://codeberg.org/api/packages/eduardgrebe/pypi/simple/ residualrisk
```

Launch R and install the `reticulate` package if not installed:

```r
install.packages("reticulate")
```

Import the Python package using `reticulate`:

```r
library(reticulate)
use_virtualenv("./.venv", required = TRUE) # Assumes venv is in the working directory
rr <- import("residualrisk")
bs <- rr$risk_days_bs(...)           # returns a Python tuple; index with [[1]], [[2]], ...
```

## Dependencies

### Python (Core)

- `streamlit` - Web application framework
- `polars` - Data manipulation
- `numpy` - Numerical computing
- `scipy` - Statistical functions and numerical integration
- `plotly` - Interactive plotting
- `watchdog` - File monitoring for Streamlit

### Go

- `gonum.org/v1/gonum` - Scientific computing (statistics, integration, distributions)

## Model Description

The model estimates the **infectious window period (IWP)**: the time interval during which a donation contains infectious virus but falls below the NAT detection threshold. The IWP can be multiplied with HIV incidence to obtain the residual risk of HIV transfusion transmission.

### Key Steps

1. **Viral Growth**: Concentration increases exponentially: C(t) = C₀ × 2^(t/doubling_time)
2. **Detection Probability**: Based on LOD characteristics and pooling/retesting protocol
3. **Infectivity Probability**: P(infection) = 1 - exp(-k × n_copies), where n_copies depends on viral load and transfusion volume. k is sampled each bootstrap iteration from the chosen input distribution: a posterior sample array, a parametric Inverse Gamma(α, β), or a two-component lognormal mixture.

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

### Scope and limitations

The tool estimates **risk-day-equivalents (RDEs) / the infectious window period**
and applies a **pre-computed incidence** to obtain residual risk. It deliberately
does **not** model population-level PrEP-use prevalence, donor self-deferral or
discard, or stratified (sex-, route-, or donor-type-specific) incidence — those
depend on a specific blood operator and time period.

**One product per run.** RDEs are estimated for a single transfused product at a
time, in both the baseline and PrEP models. Model each product *separately* by
entering that product's transfused plasma volume and range and running the tool
once per product — e.g. red-cell units (~20 mL residual plasma), fresh frozen
plasma (~200 mL), or platelets (using your own plasma-volume estimate). The tool
does **not** automate simultaneous multi-product estimation; this is an explicit
design choice, and combining products into a population figure is left to the user.

For a sophisticated analysis that needs those layers, build them in code on top of
the [Python API](#python-api) (which wraps the Go engine) — i.e. call the RDE
functions and assemble the population model around their outputs. For the **web
interface**, fold the disaggregated inputs (PrEP-use proportions, stratified
incidences, self-deferral/discard, etc.) into a single **effective incidence** for
each scenario (e.g. oral vs injectable PrEP) and supply that as the incidence
input. (See the project `TODO.md` → "Scope & validation" for the rationale and the
full list of population-level components that live outside the tool.)

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
├── app.py                   # Multipage entry point / st.navigation router
├── estimator.py             # Main estimator page (imports the residualrisk package)
├── pages/                   # Secondary Streamlit pages (Documentation, Credits)
├── docs/                    # In-app docs (theory.md, theory_prep.md, assays.md) + figures
├── residualrisk/            # Installable Python package (core calculation engine)
│   ├── __init__.py          # Public API surface
│   ├── core.py              # Baseline calculation engine (bootstrap, integration, IWP)
│   ├── prep.py              # PrEP-breakthrough model
│   ├── assays.py            # Canned NAT-assay LoD presets
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

## License and usage terms

Source code and text copyright © 2025–2026 Vitalant, with components © Eduard Grebe Consulting. All code is released under the GNU Affero General Public License v3.0 (AGPL).
Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>

You are free to use, copy and host instances of the app, for both noncommercial and commercial applications, as long as the creators are credited and the terms of the AGPL are complied with.

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

## Author

**Eduard Grebe**
- Email: egrebe@vitalant.org
- Email: eduard@grebe.consulting
- Institution: Vitalant Research Institute

Contributors:
- **Brian Custer** (Vitalant Research Institute) — conceptualization, supervision, oversight, guidance, and financial support
- **Vivian I. Avelino-Silva** — conceptualization
- **Marjorie D. Bravo** — collaboration and data curation
- **Michael P. Busch** (Vitalant Research Institute) — conceptualization and guidance
- **Artur Belov** (U.S. Food and Drug Administration) — infectivity model development and *k*-parameter posterior distributions (human and animal data)

## Citation

If you use this tool in your research, please cite it as:

> Grebe, E. (2026). *Residual HIV Transfusion Transmission Risk Estimation Tool* (Version 1.1.0) [Software]. Vitalant Research Institute. https://codeberg.org/eduardgrebe/residualriskapp

BibTeX:

```bibtex
@software{grebe2026rr,
  author       = {Grebe, Eduard},
  title        = {{Residual HIV Transfusion Transmission Risk Estimation Tool}},
  year         = {2026},
  version      = {1.1.0},
  url          = {https://codeberg.org/eduardgrebe/residualriskapp},
  organization = {Vitalant Research Institute},
  license      = {AGPL-3.0-or-later}
}
```

Machine-readable citation metadata is provided in [`CITATION.cff`](CITATION.cff). A DOI will be assigned for a future public release; please update your citation accordingly once one is available.

## References

The model makes use of established methodology and novel approaches for HIV transfusion transmission risk estimation. Key literature informing the model includes:

- Fiebig, E.W., et al. (2003). Dynamics of HIV viremia and antibody seroconversion in plasma donors: implications for diagnosis and staging of primary HIV infection. *AIDS*, 17(13):1871-1879. doi:[10.1097/00002030-200309050-00005](https://doi.org/10.1097/00002030-200309050-00005).

- Weusten, J.J.A.M., van Drimmelen, H.A.J., Lelie, N.P. (2002) Mathematic modeling of the risk of HBV, HCV, and HIV transmission by window-phase donations not detected by NAT. *Transfusion*, 42(5):537-548. doi:[10.1046/j.1537-2995.2002.00099.x](https://doi.org/10.1046/j.1537-2995.2002.00099.x).

- Weusten J., et al. (2011) Refinement of a viral transmission risk model for blood donations in seroconversion window phase screened by nucleic acid testing in different pool sizes and repeat test algorithms. *Transfusion*, 51(1):203-215. doi:[10.1111/j.1537-2995.2010.02804.x](https://doi.org/10.1111/j.1537-2995.2010.02804.x).

- Grebe E., et al. (2020) HIV incidence in US first-time blood donors and transfusion risk with a 12-month deferral for men who have sex with men. *Blood*, 136(11):1359-1367. doi:[10.1182/blood.2020007003](https://doi.org/10.1182/blood.2020007003).

- Grebe E, Avelino-Silva VI, Bravo MD, Busch MP, Custer B. (2025) Development of a risk assessment model of HIV transfusion transmission associated with undisclosed use of pre-exposure prophylaxis (PrEP) by blood donors. [ISBT Abstract PA28-L04; 35th Regional Congress of the ISBT, Milan, Italy.] *Vox Sanguinis*, 120(Suppl. 1):110.

- Belov A., et al. (2023) Modeling the Risk of HIV Transfusion Transmission. *J Acquir Immune Defic Syndr*, 92(2):173-179. doi:[10.1097/QAI.0000000000003115](https://doi.org/10.1097/QAI.0000000000003115).

## Support and Contributions

For questions, bug reports, or feature requests:
- Email: egrebe@vitalant.org or eduard@grebe.consulting
- Issue tracker: https://codeberg.org/eduardgrebe/residualriskapp/issues

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

## Funding and acknowledgments

Developed at [Vitalant Research Institute](https://research.vitalant.org) for blood safety research. Development of this tool was primarily sponsored by Vitalant Research Institute, with additional support from [Eduard Grebe Consulting](https://grebe.consulting).
