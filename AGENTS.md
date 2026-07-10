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
│   ├── 1_Documentation.py # Tabbed docs (Baseline model & methods / NAT assay parameters / PrEP model)
│   └── 2_Credits.py
├── docs/                  # Markdown documentation rendered in-app, plus figures
│   ├── theory.md          # Baseline model technical documentation
│   ├── theory_prep.md     # PrEP-breakthrough model technical documentation
│   ├── assays.md          # NAT assay parameters (LoD sources, WHO IS, IU→copies)
│   ├── credits.md
│   └── figures/           # Generated figures (+ make_*_figures.py generators)
├── residualrisk/          # Installable Python package (core calculation engine)
│   ├── __init__.py        # Public API surface (re-exports from core, prep, assays, _go)
│   ├── core.py            # Baseline calculation engine (formerly residualrisk.py)
│   ├── prep.py            # PrEP-breakthrough model (viral dynamics, serology, bootstrap)
│   ├── assays.py          # Canned NAT-assay LoD presets (single source of truth)
│   └── _go.py             # Go binary wrapper (formerly residualrisk_go.py)
├── go/                    # High-performance Go reimplementation (10-50x faster)
│   ├── main.go            # CLI interface with JSON I/O
│   └── riskdays/          # Core Go package
│       ├── models.go      # Data structures
│       ├── riskdays.go    # Bootstrap orchestration
│       ├── integration.go # Numerical integration
│       ├── probability.go # Probability calculations
│       ├── prep.go        # PrEP viral dynamics, detection & infectivity
│       ├── prep_models.go # PrEP data structures
│       ├── prep_integration.go # PrEP integrator (compact-support Gauss-Legendre)
│       ├── kde.go         # KDE-log mode estimation (--kde-mode)
│       ├── hsm.go         # Half-sample-mode estimator
│       ├── helpers.go     # Utility functions
│       ├── random.go      # Random sampling
│       ├── version.go     # Go binary version (single source of truth)
│       └── *_test.go      # Go tests (riskdays, prep, kde)
├── scripts/
│   ├── build_go.sh        # One-command wrapper for `make -C go build`
│   └── run_tests.sh       # Go + Python test runner (sandbox-safe `fast` mode)
├── static/                # Pre-computed Bayesian posterior distributions (Parquet)
├── tests/                 # Python test suite
├── docker/                # Docker build and deployment scripts
├── .streamlit/            # Streamlit configuration
├── TODO.md                # Canonical task tracker (see "Task Tracking")
├── pyproject.toml         # Python project config (managed by uv, hatchling build backend)
└── uv.lock                # Locked dependency versions (do not edit manually)
```

### Public Python API

`residualrisk` is a proper installable package. `uv sync` (or `uv pip install -e .`) installs it into `.venv` so `import residualrisk` resolves from anywhere. The public surface — defined in `residualrisk/__init__.py` — is:

- `risk_days_bs`, `risk_days_prep_bs`, `iwp_from_lookback_data`, `residual_risk_rd`, `total_residual_risk_rd` — top-level estimation functions
  - `risk_days_bs` accepts `integration_method` (`"gauss-legendre"` default, or `"quad"`). The default is a fixed 1000-point Gauss-Legendre rule matching the Go backend (robust; immune to the adaptive-quad "missed peak" failure on compact-support integrands). `"quad"` selects scipy adaptive Gauss-Kronrod and is **Python-path only** (`use_go=False`) — provided for reproducing prior analyses computed with quad. `use_go=True` with `integration_method="quad"` raises `ValueError`.
  - `risk_days_bs` also accepts `assay` — a canned-NAT-assay slug (see `NAT_ASSAYS` / `lods_for_assay` below) whose published 50%/95% LoDs are used in place of the explicit `lod50`/`lod50_sd`/`lod95_lod50_ratio`. It is **mutually exclusive** with that triplet: pass exactly one of the two (passing both, or neither, raises `ValueError`), mirroring the mutually-exclusive *k* input modes.
  - `risk_days_prep_bs` accepts the same `integration_method` kwarg with identical semantics. The default `"gauss-legendre"` is **more than future-proofing here**: the PrEP integrand has *compact support* (exactly zero before the eclipse phase and after the serology cutoff), so adaptive `"quad"` can silently miss a narrow active window and return ~0 — Gauss-Legendre fixes this. Production serology defaults give a wide window where the two methods agree (~6 sig figs); `"quad"` is retained for reproducing prior PrEP analyses (Python-path only; `use_go=True` + `"quad"` raises).
  - `risk_days_prep_bs` also takes `a_dist_uniform` / `b_dist_uniform` (default `None`). The sinusoidal set-point oscillation params `a` (amplitude) and `b` (frequency) are **held fixed** at their scalar values unless a `(lo, hi)` tuple is given, in which case they are sampled `Uniform(lo, hi)` per bootstrap iteration. `offset` is never varied. `a` — and the upper bound of `a_dist_uniform` — must be `<= offset` (enforced in Python and Go `Validate()`; in the UI the sampled *a*-range slider is capped at offset, and a scalar `a > offset` is caught on Run and shown as an error rather than crashing), since `a > offset` would drive the plateau viral load negative. Prior published analyses kept `a`/`b` fixed (a=0.7, b=0.6); the UI exposes this via a "Vary sinusoidal oscillation parameters (a, b)" checkbox (off by default; default ranges a∈(0.5,0.9), b∈(0.4,0.8)).
  - `risk_days_prep_bs` also takes `drug_effect` (default `1.0`) and `drug_effect_dist_uniform` (default `None`). `drug_effect` is an antiretroviral transmissibility-reduction factor in `(0, 1]` (1.0 = no reduction) for PrEP-breakthrough infectivity, applied as a **linear multiplier on the per-time infection probability inside the integrand** (`_drug_effect` in `prep.py` / `DrugEffectFactor` in `prep.go`) — held fixed at the scalar unless `drug_effect_dist_uniform=(lo, hi)` (with `0 < lo ≤ hi ≤ 1`) is given, then sampled `Uniform(lo, hi)` per bootstrap iteration. Default 1.0 leaves results bit-for-bit unchanged. Because the factor is constant in `t` it factors out of the RDE integral (numerically identical to scaling the RDE), but `_drug_effect`/`DrugEffectFactor` deliberately take `t` as a **placeholder for a future time-varying drug effect** (e.g. long-acting-injectable wash-out — see the deferred PK/PD modelling task in `TODO.md`); that placement is the only correct one once it depends on `t`. Set **independently per scenario** in the UI (oPrEP / iPrEP): a point input (default 1.0) plus a bootstrap range slider (default `(1.0, 1.0)` → fixed; e.g. `(0.5, 1.0)` reproduces the prior analysis). Enforced (`0 < drug_effect ≤ 1`, range `⊂ (0,1]`) in Python and Go `Validate()`. Matches prior `rr_prep_v3.py` (`drug_effect ~ Uniform(0.5, 1.0)`, median ≈ 0.75 = "25% reduction").
  - Both `risk_days_bs` and `risk_days_prep_bs` return `(pe, cri, range, rdests, sim_df)` (a **5-tuple**); `sim_df` is `None` unless `return_sim_df=True`, in which case the per-iteration frame carries a `backend` column (`'go'`/`'python'`) tagged at dispatch (`core._append_backend`) — recording the engine that *actually* ran, even after a silent Go→Python fallback (now `logging.warning`-ed, not swallowed).
- `get_cpu_core_count`, `mode_rounded` — utility helpers used by the UI
- `mode_kde` — estimate the mode of a positive posterior via KDE on the log scale (pure-Python, slow on large posteriors; kept as fallback)
- `mode_kde_go` — fast Go-backed KDE mode estimation via the `--kde-mode` subcommand; defaults `cap=None, n_grid=1_000_000` (the values `estimator.py` passes at load time, cached via `@st.cache_data`; < 0.1% error vs the pure-Python mode), with a hardcoded fallback if the Go binary is unavailable
- `sample_invgamma` — sample from an Inverse Gamma distribution; supports `alpha`+`beta` or `alpha`+`mode` parameterisations
- `sample_lnmix` — sample from a two-component lognormal mixture; parameters: `n, w, mu1, sigma1, mu2, sigma2, seed=None`
- `NAT_ASSAYS`, `lods_for_assay`, `list_assays`, `AssayLoD` — canned NAT-assay limit-of-detection presets (HIV-1 Group M, copies/mL), defined in `residualrisk/assays.py` as the **single source of truth** (consumed by both the API and `estimator.py`). `NAT_ASSAYS` is a dict keyed by **slug** (`ultrio`, `ultrio_plus`, `ultrio_elite`, `cobas_taqscreen_mpx`, `cobas_taqscreen_mpxv2`, `cobas_mpx`, `biomanguinhos`); each entry carries `display_name`, `lod50`, `lod50_sd`, `lod95`, and the informational `cp_per_iu` / `iu_std` provenance fields (the upstream IU/mL→copies/mL factor and its WHO IS — **not** constant across assays). `lods_for_assay(slug)` returns an `AssayLoD` namedtuple (adds the derived `lod95_lod50_ratio`) and raises `ValueError` on an unknown slug; `list_assays()` returns `{slug: display_name}` for menus. Prefer `risk_days_bs(assay=…)` (or these helpers) over transcribing LoD numbers. **Provisional SD:** `biomanguinhos` (Brazilian NAT Platform, Bio-Manguinhos) has no published LoD50 CI — Rocha et al. (2018) report point 50%/95% LoDs only — so its `lod50_sd` is an *assumed* relative SD of **13%** (6.08 IU/mL ≈ 3.527 cp/mL); **4.95 IU/mL (RSE 10.6%) was the value used in prior analyses.** Revisit if the per-dilution hit-rate table becomes available; see the note above `NAT_ASSAYS` in `residualrisk/assays.py`.
- `find_go_binary` — locator for the Go binary (honors `$RESIDUALRISK_GO_BINARY` env var)
- `mode_hsm_go` — half-sample-mode estimate via the Go binary (`--hsm-mode`)
- `risk_days_prep_bs_go` — direct Go-backed PrEP bootstrap; the `use_go=True` path of `risk_days_prep_bs` dispatches here
- `__version__` — package version

Downstream analyses (e.g. R scripts via `reticulate`) should call these rather than reaching into `residualrisk.core` or `residualrisk._go`. Test code may import `residualrisk.core` directly to exercise private `_`-prefixed functions.

#### LoD50 relative standard error (RSE) by assay

RSE = `lod50_sd / lod50` (the coefficient of variation of the 50% LoD; invariant under the IU/mL→copies/mL conversion). For every assay except Bio-Manguinhos the SD derives from a 95% CI of the 50% LoD; Bio-Manguinhos uses an *assumed* RSE (see the provisional-SD note above and in `residualrisk/assays.py`). Underlying LoD/CI data are compiled in the companion analysis `residualrisk_analysis/assays/ASSAYS.qmd`.

| Assay | RSE on LoD50 | Source |
|---|---|---|
| `cobas_taqscreen_mpxv2` | 4.72% | Probit fit of Roche insert reactivity data (95% CI) |
| `cobas_mpx` | 6.04% | Roche cobas MPX CE/IVD insert (95% CI of 50% LoD) |
| `cobas_taqscreen_mpx` | 7.00% | Probit fit of Roche insert reactivity data (95% CI) |
| `ultrio_plus` | 7.07% | Grifols Procleix Ultrio Plus insert (95% CI of 50% LoD) |
| `ultrio` | 7.28% | Grifols Procleix Ultrio insert, dHIV-1 (95% CI of 50% LoD) |
| `ultrio_elite` | 7.55% | Grifols Procleix Ultrio Elite insert (95% CI of 50% LoD) |
| `biomanguinhos` | 13.00% (assumed, provisional) | Assumed RSE — Rocha et al. (2018) report point LoDs only (no CI) |

## Core Application Files

- **`app.py`** — Multipage entry point / router
  - Thin `st.navigation` router: defines the page list and the shared page config (title, favicon)
  - Holds `APP_VERSION`; renders the shared sidebar footer (VRI logo + centred app/library version caption) on every page
  - Sets explicit nav labels: **Estimator** (default page), **Documentation**, **Credits**
  - Entry point: `streamlit run app.py`

- **`estimator.py`** — Streamlit web UI (the Estimator page)
  - Parameter input interface with expandable sections
  - Two RDE estimation methods (selectbox): **Lookback data** and **Mechanistic model**
  - Optional **oral- and/or injectable-PrEP** breakthrough-infection RDE components (oPrEP/iPrEP checkboxes) layered on top, with per-component and additive total residual risk
  - NAT-assay dropdown (canned LoD presets from `residualrisk/assays.py`) or manual LoD entry
  - Real-time calculation and result visualization
  - Imports via the public API: `import residualrisk as rr`

- **`residualrisk/core.py`** — Baseline calculation engine
  - Viral concentration dynamics
  - Infectivity probability calculations
  - Bootstrap simulation methods (`risk_days_bs`, `residual_risk_rd`, `total_residual_risk_rd`)
  - Integration with Go acceleration via `residualrisk/_go.py`

- **`residualrisk/prep.py`** — PrEP-breakthrough model
  - Breakthrough viral dynamics (eclipse → growth → oscillating plateau; analytic `tcrit`)
  - **Units:** `set_point` / `set_point_dist_uniform` are a clinical breakthrough viral load in **RNA copies/mL**. The model's concentration `C` is in **virions/mL** (`k` is calibrated per RNA copy, Belov 2023), so `_find_tcrit` / `_vl_postbt` (and Go `FindTcrit` / `VLPostBT`) divide the set-point by `copies_per_virion` (χ=2). Fixed on branch `fix-prep-setpoint-units` (2026-07-08): a copies/mL set-point was previously used directly as `C` (virions/mL), running the plateau 2× high and under-stating the RDE by ≈20–45%.
  - NAT + serology (Weibull) detection; optional drug-effect transmissibility reduction
  - `risk_days_prep_bs` bootstrap (compact-support Gauss-Legendre integrator), Go-accelerated

- **`residualrisk/assays.py`** — Canned NAT-assay LoD presets
  - `NAT_ASSAYS` table (single source of truth) + `lods_for_assay` / `list_assays` helpers
  - Consumed by both the public API (`risk_days_bs(assay=…)`) and `estimator.py`; see the Public Python API section above

- **`residualrisk/_go.py`** — Go binary wrapper
  - JSON-based communication with Go binary
  - Automatic fallback to Python if Go binary unavailable
  - Progress monitoring
  - `find_go_binary()` search order: `$RESIDUALRISK_GO_BINARY` env var → `<repo>/go/bin/riskdays_go` → the wheel-bundled `residualrisk/_bin/riskdays_go-<goos>-<goarch>` (validated against `__version__`; built by `hatch_build.py`) → `~/.local/bin/riskdays_go` (sudo-free) → `$PATH`
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
  `mode_kde_go()` (Go KDE subprocess, < 0.1% error vs the pure-Python mode) cached
  by `@st.cache_data`, with hardcoded fallback if Go binary is unavailable.
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

The PrEP-breakthrough model is **fully implemented** across the stack:

- **Python** — `residualrisk/prep.py` (`risk_days_prep_bs`), exported from the public API. Supports the same *k* input distributions as the baseline, optional sinusoidal `a`/`b` uncertainty, and the optional `drug_effect` transmissibility-reduction factor (see the Public Python API section for the full kwarg set).
- **Go** — `go/riskdays/prep*.go`, bridged via `risk_days_prep_bs_go` in `_go.py`, with Python↔Go parity tests.
- **UI** — `estimator.py` exposes oral-PrEP (oPrEP) and injectable-PrEP (iPrEP) breakthrough risk via **checkboxes** (not a separate RDE-method dropdown option), each with its own parameter expander, plus per-component and additive total residual risk (`total_residual_risk_rd`).
- **Documentation** — `docs/theory_prep.md`, rendered in the Documentation page's **"PrEP model"** tab; **awaiting EG review** (see `TODO.md`).

Open follow-ups live in `TODO.md` — notably independent oral/injectable `drug_effect` draws for the total-risk credible interval, and the deferred PK/PD drug-concentration extension.

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

1. Decide which version(s) to bump (app, library, Go, or any combination). If a
   version currently carries an in-development `.dev0` suffix (see "In-development
   versions" below), drop it to the clean release version now — the `verify-version`
   gate rejects a tag whose `APP_VERSION` still has `.dev0`.
2. Edit the relevant file(s):
   - App: `app.py` → `APP_VERSION`
   - Library: `residualrisk/__init__.py` → `__version__`
   - Go: `go/riskdays/version.go` → `Version`
3. Rebuild the Go binary if the Go version changed: `bash scripts/build_go.sh`
4. Commit the change.
5. Tag the commit with the new app version and push (signed tag — this project's maintainer signs tags):
   ```bash
   git tag -s vX.Y.Z -m "Release vX.Y.Z"
   git push origin vX.Y.Z
   ```
   On a `v*` tag push, the `docker-publish.yml` `verify-version` job asserts that
   `APP_VERSION` in `app.py` equals the tag (minus the `v`) before any image is
   published — so the tag and `APP_VERSION` must match (e.g. tag `v1.1.0a7` ⇔
   `APP_VERSION = "1.1.0a7"`). Pre-release tags (PEP 440, e.g. `v1.1.0a7`) publish
   only their full-version image and never move the `latest` / `X.Y` tags.

### In-development versions (between releases)

Between tagged releases, set the version(s) under active change to the **next intended
release plus a PEP 440 dev segment** — e.g. `1.1.0b3.dev0` (use the canonical dotted
form, *not* `1.1.0b3dev0`). This marks `main` as unreleased work heading toward that
version:

- **Ordering is correct:** `1.1.0b2 < 1.1.0b3.dev0 < 1.1.0b3`. The `.dev0` target is
  only an intention — if the next release turns out different (e.g. the stable `1.1.0`),
  bump straight to it; ordering still holds.
- **Never tagged:** CI/Docker fire only on `v*` tags, so a `.dev0` version publishes
  nothing. The app sidebar shows e.g. `App v1.1.0b3.dev0`, flagging a running instance
  (including the live deployment) as an unreleased build rather than a tagged release.
- **Drop `.dev0` before tagging** (Releasing step 1): the `verify-version` gate rejects
  a tag whose `APP_VERSION` still carries `.dev0` — the safety net working as intended.
- Apply to whichever of the three versions have unreleased changes; the convention here
  is to **re-sync all three** to the same value when a change spans them (as with the
  `.dev0` bump itself).

## Distribution — the bundled-binary wheel

The `residualrisk` wheel **bundles the pre-compiled Go accelerator for every platform**, so
`pip install residualrisk` runs Go-accelerated with no Go toolchain and no separate build.

- **Build hook** (`hatch_build.py`, registered in `pyproject.toml`): at wheel-build time it
  cross-compiles `riskdays_go` for all six targets — `linux`/`darwin`/`windows` × `amd64`/`arm64` —
  static (`CGO_ENABLED=0`) and stripped, into `residualrisk/_bin/riskdays_go-<goos>-<goarch>[.exe]`,
  and `force_include`s them (the dir is git-ignored). Go cross-compiles all six from one machine (no
  QEMU/matrix). If `go` is absent the hook **skips** bundling → a pure-Python wheel. The result is a
  single **`py3-none-any` "fat" wheel** (~6 MB); `find_go_binary()` selects and version-checks the
  right binary at runtime (see its search order above). Per-platform wheels are possible (build one
  target + set `build_data["tag"]`), but the fat wheel is what we ship.
- **Publishing** (`.forgejo/workflows/release.yml`, on Codeberg's **hosted** runner — no self-hosted
  instance needed, unlike the Docker image): on a `v*` tag it verifies `tag == __version__`, runs
  `uv build`, publishes to the **Codeberg PyPI registry**
  (`https://codeberg.org/api/packages/eduardgrebe/pypi`), and creates a Codeberg Release with the
  artifacts. Requires Forgejo Actions enabled + a repo secret **`PACKAGE_TOKEN`** (Forgejo token,
  `write:package` scope); confirm the `runs-on` label matches the runner (`docker` on Codeberg).
- **Installing** (see `README.md` for pip / uv / `pyproject.toml`): use `--extra-index-url`, **not**
  `--index-url`, so dependencies still resolve from PyPI — e.g. `pip install --extra-index-url
  https://codeberg.org/api/packages/eduardgrebe/pypi/simple/ residualrisk`.
- **Cross-OS smoke test** (`.github/workflows/smoke-test.yml`, kept on GitHub even as the mirror,
  because Codeberg's hosted runners are Linux-only): builds the wheel on Linux, then *runs* the
  cross-compiled macOS/Windows binaries on real macOS/Windows runners via
  `scripts/smoke_bundled_binary.py` (`--version` + a small computation). **Gap:** `windows-arm64`
  has no hosted runner, so it ships untested (Windows/ARM runs the amd64 build via emulation; the
  pure-Python engine is the fallback).

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

