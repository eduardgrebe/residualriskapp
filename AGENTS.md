# Agent Instructions for residualriskapp

## Task Tracking

**`TODO.md` in the repository root is the canonical task tracker for this project.**

Before starting any work, read `TODO.md` to understand open tasks, in-progress work, and recently completed items. After completing a task, move it from the **Open** section to the **Completed** section (or update its checkbox). Add new tasks to the **Open** section as they are identified.

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
├── app.py                 # Multipage entry point / st.navigation router
├── estimator.py           # Main estimator page (Streamlit web UI)
├── pages/                 # Secondary Streamlit pages
│   ├── 1_Documentation.py
│   └── 2_Credits.py
├── residualrisk/          # Installable Python package (core calculation engine)
│   ├── __init__.py        # Public API surface (re-exports from core and _go)
│   ├── core.py            # Calculation engine (formerly residualrisk.py)
│   └── _go.py             # Go binary wrapper (formerly residualrisk_go.py)
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
├── scripts/
│   └── build_go.sh        # One-command wrapper for `make -C go build`
├── static/                # Pre-computed Bayesian posterior distributions (Parquet)
├── tests/                 # Python test suite
├── docker/                # Docker build and deployment scripts
├── .streamlit/            # Streamlit configuration
├── pyproject.toml         # Python project config (managed by uv, hatchling build backend)
└── uv.lock                # Locked dependency versions (do not edit manually)
```

### Public Python API

`residualrisk` is a proper installable package. `uv sync` (or `uv pip install -e .`) installs it into `.venv` so `import residualrisk` resolves from anywhere. The public surface — defined in `residualrisk/__init__.py` — is:

- `risk_days_bs`, `risk_days_prep_bs`, `iwp_from_lookback_data`, `residual_risk_rd` — top-level estimation functions
  - `risk_days_bs` accepts `integration_method` (`"gauss-legendre"` default, or `"quad"`). The default is a fixed 1000-point Gauss-Legendre rule matching the Go backend (robust; immune to the adaptive-quad "missed peak" failure on compact-support integrands). `"quad"` selects scipy adaptive Gauss-Kronrod and is **Python-path only** (`use_go=False`) — provided for reproducing prior analyses computed with quad. `use_go=True` with `integration_method="quad"` raises `ValueError`.
  - `risk_days_prep_bs` accepts the same `integration_method` kwarg with identical semantics. The default `"gauss-legendre"` is **more than future-proofing here**: the PrEP integrand has *compact support* (exactly zero before the eclipse phase and after the serology cutoff), so adaptive `"quad"` can silently miss a narrow active window and return ~0 — Gauss-Legendre fixes this. Production serology defaults give a wide window where the two methods agree (~6 sig figs); `"quad"` is retained for reproducing prior PrEP analyses (Python-path only; `use_go=True` + `"quad"` raises).
  - `risk_days_prep_bs` also takes `a_dist_uniform` / `b_dist_uniform` (default `None`). The sinusoidal set-point oscillation params `a` (amplitude) and `b` (frequency) are **held fixed** at their scalar values unless a `(lo, hi)` tuple is given, in which case they are sampled `Uniform(lo, hi)` per bootstrap iteration. `offset` is never varied. `a` — and the upper bound of `a_dist_uniform` — must be `<= offset` (enforced in Python and Go `Validate()`, and capped in the UI), since `a > offset` would drive the plateau viral load negative. Prior published analyses kept `a`/`b` fixed (a=0.7, b=0.6); the UI exposes this via a "Vary sinusoidal oscillation parameters (a, b)" checkbox (off by default; default ranges a∈(0.5,0.9), b∈(0.4,0.8)).
  - `risk_days_prep_bs` also takes `drug_effect` (default `1.0`) and `drug_effect_dist_uniform` (default `None`). `drug_effect` is an antiretroviral transmissibility-reduction factor in `(0, 1]` (1.0 = no reduction) for PrEP-breakthrough infectivity, applied as a **linear multiplier on the per-time infection probability inside the integrand** (`_drug_effect` in `prep.py` / `DrugEffectFactor` in `prep.go`) — held fixed at the scalar unless `drug_effect_dist_uniform=(lo, hi)` (with `0 < lo ≤ hi ≤ 1`) is given, then sampled `Uniform(lo, hi)` per bootstrap iteration. Default 1.0 leaves results bit-for-bit unchanged. Because the factor is constant in `t` it factors out of the RDE integral (numerically identical to scaling the RDE), but `_drug_effect`/`DrugEffectFactor` deliberately take `t` as a **placeholder for a future time-varying drug effect** (e.g. long-acting-injectable wash-out — see the deferred PK/PD modelling task in `TODO.md`); that placement is the only correct one once it depends on `t`. Set **independently per scenario** in the UI (oPrEP / iPrEP): a point input (default 1.0) plus a bootstrap range slider (default `(1.0, 1.0)` → fixed; e.g. `(0.5, 1.0)` reproduces the prior analysis). Enforced (`0 < drug_effect ≤ 1`, range `⊂ (0,1]`) in Python and Go `Validate()`. Matches prior `rr_prep_v3.py` (`drug_effect ~ Uniform(0.5, 1.0)`, median ≈ 0.75 = "25% reduction").
- `get_cpu_core_count`, `mode_rounded` — utility helpers used by the UI
- `mode_kde` — estimate the mode of a positive posterior via KDE on the log scale (pure-Python, slow on large posteriors; kept as fallback)
- `mode_kde_go` — fast Go-backed KDE mode estimation via the `--kde-mode` subcommand; `cap=40_000, n_grid=5_000` by default (< 0.1% error, ~0.9s for all three posteriors); used by `estimator.py` at load time with Python fallback to hardcoded values
- `sample_invgamma` — sample from an Inverse Gamma distribution; supports `alpha`+`beta` or `alpha`+`mode` parameterisations
- `sample_lnmix` — sample from a two-component lognormal mixture; parameters: `n, w, mu1, sigma1, mu2, sigma2, seed=None`
- `find_go_binary` — locator for the Go binary (honors `$RESIDUALRISK_GO_BINARY` env var)
- `__version__` — package version

Downstream analyses (e.g. R scripts via `reticulate`) should call these rather than reaching into `residualrisk.core` or `residualrisk._go`. Test code may import `residualrisk.core` directly to exercise private `_`-prefixed functions.

## Core Application Files

- **`app.py`** — Multipage entry point / router
  - Thin `st.navigation` router: defines the page list and the shared page config (title, favicon)
  - Holds `APP_VERSION`; renders the shared sidebar footer (VRI logo + centred app/library version caption) on every page
  - Sets explicit nav labels: **Estimator** (default page), **Documentation**, **Credits**
  - Entry point: `streamlit run app.py`

- **`estimator.py`** — Streamlit web UI (the Estimator page)
  - Parameter input interface with expandable sections
  - Supports three RDE estimation methods: Mechanistic model, Lookback data, Mechanistic model with PrEP
  - Real-time calculation and result visualization
  - Imports via the public API: `import residualrisk as rr`

- **`residualrisk/core.py`** — Main calculation engine
  - Viral concentration dynamics
  - Infectivity probability calculations
  - Bootstrap simulation methods
  - Integration with Go acceleration via `residualrisk/_go.py`

- **`residualrisk/_go.py`** — Go binary wrapper
  - JSON-based communication with Go binary
  - Automatic fallback to Python if Go binary unavailable
  - Progress monitoring
  - `find_go_binary()` search order: `$RESIDUALRISK_GO_BINARY` env var → `<repo>/go/bin/riskdays_go` → `/usr/local/bin/riskdays_go` → `$PATH`
  - `mode_kde_go()` — KDE mode via `riskdays_go --kde-mode`; pre-caps data in Python to minimise JSON payload

## Technical Architecture

### Risk Estimation Model

1. **Viral Dynamics**: Exponential growth from initial concentration (C0) with doubling time
2. **Test Sensitivity**: Incorporates LOD (limit of detection) with uncertainty
3. **Infectivity**: Probabilistic model using k sampled each bootstrap iteration from the chosen input distribution — either a posterior sample array (human, animal, or human-weighted exponential-decay) or a parametric Inverse Gamma distribution (α, β specified by the user).
4. **Bootstrap Simulation**: Monte Carlo sampling of parameter uncertainty
5. **Window Period Calculation**: Numerical integration to find infectious window period

**Design scope (single product per run).** The tool estimates RDEs / the
infectious window period for **one transfused product at a time**, in both the
baseline and PrEP models, then applies a **pre-computed incidence** to obtain
residual risk. Model each product separately by setting its transfused plasma
volume + range (`volume_transfused` / `volume_transfused_range`) and running once
per product (e.g. red cells ~20 mL residual plasma, FFP ~200 mL, platelets with
their own plasma-volume estimate). The tool deliberately does **not** automate
multi-product estimation, nor the population-level "Layer 2" aggregation
(PrEP-use prevalence, self-deferral / discard, stratified incidence) — that is the
user's, built on the Python API or folded into an "effective incidence" per
scenario. See `TODO.md` → "Scope & validation" for the full rationale.

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
- `k` — Infectivity parameter point estimate (used for IWP point estimate only; bootstrap samples from the chosen distribution)
- `k_posterior_sample` — Array of posterior draws for k (used when sampling from a posterior)
- `k_invgamma_alpha` — Shape parameter α for Inverse Gamma k distribution (omit or `None` for posterior-sample paths)
- `k_invgamma_beta` — Scale parameter β for Inverse Gamma k distribution (omit or `None` for posterior-sample paths)
- `volume_transfused` — Per-product transfused **plasma** volume (mL); single-product-per-run, so set this per product and run once per product (see Design scope above)
- `volume_transfused_min/max` — Uncertainty range for the plasma volume
- `copies_per_virion` — RNA copies per virion (default: 2)

**Simulation**:
- `n_bs` — Number of bootstrap simulations (default: 10,000)
- `seed` — Random seed for reproducibility
- `threads` — Parallel workers (default: CPU cores - 1)

### Static Data (`static/`)

Pre-computed posterior parameter distributions in Parquet format:
- `k_param_human.parquet` — Human infectivity parameter posterior
- `k_param_animal.parquet` — Animal model infectivity posterior
- `k_param_expdecay.parquet` — Human-weighted exponential-decay posterior

**Do not modify these files** — they are pre-computed Bayesian posteriors. Regeneration requires rerunning the upstream Bayesian analyses.

### Input Distribution for *k* — See Companion Analysis Repo

The choice of input parameter distribution for *k* (as opposed to using a raw
posterior sample directly) is documented in the companion analysis repository
`residualrisk_analysis`, which lives alongside this repo on the same host at
`../residualrisk_analysis/`. The relevant document is:

**`residualrisk_analysis/exploration/K_PARAM_INPUTDIST.md`**

It covers:
- Detailed characterisation of the human and animal posteriors
- Systematic evaluation of candidate parametric distributions (lognormal,
  inverse gamma, log-logistic, Burr XII, lognormal mixture)
- Quantile tables, survival function comparisons, and six diagnostic figures
- Two formal recommendations with scipy parameterisations and scientific
  justification:
  - **Recommendation A**: Inverse Gamma(α=2, β=0.002019) — smooth unimodal,
    power-law tail, mode at human posterior mode (0.000673). Note: α=2 is a
    deliberate conservative choice — it is far heavier-tailed than a best-fit
    InvGamma to the human posterior (MLE ≈ α=9.5) and encodes substantial
    additional uncertainty beyond what the Belov data alone support.
  - **Recommendation B**: 90% LN(human) + 10% LN(animal) mixture — best
    preserves human posterior bulk while giving explicit 10% weight to the
    animal-derived transmissibility range
- Guidance on sensitivity analysis

**Implementation status:**
- **Inverse Gamma**: fully implemented in both Python (`residualrisk/core.py`,
  `sample_invgamma()`) and Go (`go/riskdays/random.go`, `GenerateInvGamma()`),
  with UI wiring in `estimator.py`. Supports α+β or α+mode parameterisations.
  KDE modes of the three posteriors are pre-computed at load time via
  `mode_kde_go()` (Go KDE subprocess, ~0.9s total, < 0.1% error) cached by
  `@st.cache_data`, with hardcoded fallback if Go binary is unavailable.
- **Lognormal mixture**: fully implemented in Python (`residualrisk/core.py`,
  `sample_lnmix()`), Go (`go/riskdays/random.go`, `GenerateLogNormalMixture()`),
  bridge (`residualrisk/_go.py`), and UI (`estimator.py`). Parameters: `k_lnmix_w`,
  `k_lnmix_mu1`, `k_lnmix_sigma1`, `k_lnmix_mu2`, `k_lnmix_sigma2`. Default
  values (w=0.90, μ₁=−7.2403, σ₁=0.3241, μ₂=−3.7423, σ₂=0.5258) implement
  Recommendation B. UI provides a mixing-weight slider with optional advanced
  component-parameter editing; PE options are mode/median (numerical) and mean
  (analytic).

Agents modifying the *k* parameter handling, adding new posterior files to
`static/`, or implementing a custom input distribution for *k* should consult
this document first.

## PrEP Model Status

- `estimator.py` has UI stubs for the "Mechanistic model with PrEP" dropdown (currently shows a "not yet available" message and halts)
- No PrEP model source is present in the repo at this time; the prior `residualrisk_prep.py` was removed during the package restructure
- Go implementation of PrEP model **does not exist yet**
- Plan: release initial version without PrEP, then add a PrEP module (likely at `residualrisk/prep.py`) alongside a manuscript

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
bash scripts/build_go.sh          # one-command wrapper (idempotent)

# Or directly:
cd go
make deps      # Download dependencies
make build     # Build binary to go/bin/riskdays_go
make test      # Run Go tests
```

The Python code auto-detects the binary at `<repo>/go/bin/riskdays_go`. Override with `RESIDUALRISK_GO_BINARY=/absolute/path/to/riskdays_go` when running from a different install layout.

### Testing

```bash
# Python — tests import `from residualrisk import core as rr`
# and require the package to be installed (uv sync does this).
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

### Signing and Remote Access Limitation

**Agents must never run `git commit`, `git fetch`, `git pull`, `git push`, or any other git operation that requires SSH keys for commit signing or remote authentication.** The user has commit signing configured (GPG/SSH) and agents do not have access to those keys.

For any git operation that creates commits or touches the remote, provide the user with the exact commands to run. For example:

```bash
git add -A
git commit -m "your message"
git push origin <branch-name>
git tag -a vX.Y.Z -m "message" && git push origin vX.Y.Z
```

Agents **can** safely run read-only local git commands: `git status`, `git log`, `git diff`, `git branch`, `git show`.

## Versioning

This project uses **three independent version numbers** that can evolve separately:

| Version | Location | Tracks | Bump when |
|---|---|---|---|
| **App version** | `app.py` → `APP_VERSION` | Streamlit web application | UI, UX, or app-level changes |
| **Library version** | `residualrisk/__init__.py` → `__version__` | `residualrisk` Python package API | Calculation logic or public API changes |
| **Go version** | `go/riskdays/version.go` → `Version` | Go binary (`riskdays_go`) | Go implementation changes |

`pyproject.toml` reads its version dynamically from `residualrisk/__init__.py` via hatchling, so the installable package version always matches the library version. **Do not add a hardcoded `version =` field to `[project]` in `pyproject.toml`.**

The Go version is the single source of truth for the binary: it is embedded in the JSON output (`"version"` field) and printed by `riskdays_go --version`.

Both Python versions are displayed together in the app sidebar (`App vX.Y.Z · Library vX.Y.Z`).

**Git tags** track the **app version** — that's what users interact with.

### Releasing a new version

1. Decide which version(s) to bump (app, library, Go, or any combination).
2. Edit the relevant file(s):
   - App: `app.py` → `APP_VERSION`
   - Library: `residualrisk/__init__.py` → `__version__`
   - Go: `go/riskdays/version.go` → `Version`
3. Rebuild the Go binary if the Go version changed: `bash scripts/build_go.sh`
4. Commit the change.
5. Tag the commit with the new app version and push:
   ```bash
   git tag -a vX.Y.Z -m "Release vX.Y.Z"
   git push origin vX.Y.Z
   ```

## Common Tasks

### Adding a New Parameter

1. Update calculation functions in `residualrisk/core.py`
2. If it belongs on the public API, re-export it from `residualrisk/__init__.py` and add to `__all__`
3. Add UI controls in `estimator.py`
4. Update Go implementation in `go/riskdays/` if needed
5. Update JSON schema in `go/riskdays/models.go` and wire it through `residualrisk/_go.py`
6. Document in this file and in `README.md`

### Modifying the Model

1. Update mathematical functions in `residualrisk/core.py`
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

