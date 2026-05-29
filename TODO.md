# TODO

## Open

### SESSION STATE (2026-05-29)

On `feature_prep_model`; `main` is fully merged in (merge commit `31118a6`),
working tree clean, all three versions at `1.1.0.dev0`. Everything below the
line — the PrEP model build-out, all code-review findings (H1, H2, M1, M2, M3,
L1, L2, L4, L5), the integration-robustness / truncnorm-positivity /
analytic-`tcrit` work, and the marker-based test runner — is **committed**.
`bash scripts/run_tests.sh fast` is robust (marker-based) on both branches.

**UNCOMMITTED (this session):** optional uniform uncertainty for the sinusoidal
`a`/`b` — resolves the `a`/`b`/`offset` design decision (see Completed) — plus
its tests and docs.

One item remains before the PrEP release, a **scientific task yours to drive**.
(Agents can't commit/push — provide exact git commands.)

### Remaining before the PrEP release

- [ ] **End-to-end validation of PrEP + baseline results against prior analyses.**
      Run the tool with input parameters from previous published baseline analyses
      (Grebe et al. 2020, ISBT 2024) and PrEP analyses (ISBT 2025) and confirm that
      RDE estimates are comparable to previously reported values (or shifted
      proportional to known changes). If results differ, understand and document
      why. Known shifts to expect: the truncnorm-positivity fix lowers RDE ≈8%
      (baseline) / ≈25% (PrEP); plus the k-distribution choice, GL vs quad, and
      analytic vs grid `tcrit`. Applies to both Python and Go. Unit/integration
      tests pass, but scientific validity requires reproducing known results.
      The full `bash scripts/run_tests.sh` suite (incl. the `ProcessPoolExecutor`
      tests — M2 and the varied-`a`/`b` `TestPrepPythonGoAgreementVariedAB`)
      passes **outside the sandbox** as of 2026-05-29; what remains here is the
      *scientific* reproduction of prior published results, not test-green.

---

## Completed

### PrEP sinusoidal `a`/`b` uncertainty (2026-05-29) — UNCOMMITTED

- [x] **Optional uniform uncertainty for the sinusoidal `a` (amplitude) and `b`
      (frequency).** Resolves the `a`/`b`/`offset` bootstrapping design decision:
      `a` and `b` are sampled `Uniform(lo, hi)` per bootstrap iteration when
      `a_dist_uniform` / `b_dist_uniform` are given; default `None` = fixed at the
      scalar (back-compatible, draws no RNG → reproducible). `offset` stays fixed.
      `a` and the upper bound of `a_dist_uniform` must be `<= offset` (enforced in
      Python and Go `Validate()`, capped in the UI) since `a > offset` drives the
      plateau viral load negative. Confirmed against the analysis repo
      (`residualrisk_analysis` `rr_prep*.py`) that prior published analyses always
      held `a`/`b`/`offset` fixed and varied `set_point`/`eclipse` uniform — no
      prior template for varying `a`/`b`. Six layers wired: `prep.py`,
      `go/riskdays/models.go`, `go/riskdays/riskdays.go`, `go/main.go`
      (per-iteration `a`/`b` emitted as binary output columns → correct Go-path
      `sim_df`), `residualrisk/_go.py`, `app.py` ("Vary sinusoidal oscillation
      parameters (a, b)" checkbox, off by default, ranges a∈(0.5,0.9) / b∈(0.4,0.8)).
      Tests: Go-path fixed-vs-varied `sim_df` (`test_go_prep_fixed_ab_default`,
      `test_go_prep_varied_ab`), `a > offset` validation (`test_scalar_a_exceeds_offset_raises`,
      `test_a_dist_upper_exceeds_offset_raises`), and `TestPrepPythonGoAgreementVariedAB`
      (Python↔Go equivalency with `a`/`b` varied — median ~1.4%, CrI bounds ≲12% at
      n_bs=2000; primary-parameters PE still agrees to ~1e-9 since the PE uses the
      scalar `a`/`b`). Docs updated in `AGENTS.md`.

### Code-review findings (2026-05-26 review) — all resolved

- [x] **H1** — Lookback branch resets `prep_oral_run`/`prep_inj_run` so stale PrEP
      results don't leak into the lookback display / additive total. *(4e82576)*
- [x] **H2** — clamp viral load `max(0, Cv)` in Python `_vl_postbt` + Go `VLPostBT`
      (`a > offset` no longer crashes Python / yields negative Go probability). *(4e82576)*
- [x] **M1** — Go `SetDefaults` PrEP serology + uniform-range defaults aligned to
      Python production values. *(4e82576)*
- [x] **M2** — real Python↔Go PrEP cross-validation (`TestPrepPythonGoAgreement`):
      primary-params PE within 1e-9 + median/CrI agreement at n_bs=2000 (bumped
      from 500 so the InvGamma(α=2) heavy-tail upper CrI is stable). *(1f1e6b3, calibrated in 31118a6)*
- [x] **M3** — production-serology Go golden `TestRiskDaysPrep_GoldenValue_Production`
      (≈3.0919 vs Simpson truth); narrow-window 1.0086 golden comment corrected. *(1f1e6b3)*
- [x] **L1** — `point_estimate="mode"` is KDE-log on both sides: PrEP Python now uses
      `_kde_mode_log` (was `mode_rounded`), matching Go `KDEModeLog` and baseline Python.
- [x] **L2 + L4** — analytic `_find_tcrit`; deleted grid `_vl_postbt_vec`. Fixes the
      empty-argmin crash, ~4× faster Python PrEP, machine-precision Python≡Go. *(93f4e34)*
- [x] **L5** — Go `version.go` bumped to 1.1.0.dev0; `integration.go` "adaptive quad"
      comment corrected to fixed Gauss-Legendre; `app.py:1062` PrEP "(Python-only)"
      comment fixed; `test_prep_go_parity.py` given the full AGPL header.

### PrEP model build-out (2026-05-20/21)

- [x] Design: oPrEP/iPrEP as independent additive scenarios; PrEP breakthrough risk
      layered on top of baseline window-period risk.
- [x] `app.py` wiring — oPrEP/iPrEP checkboxes, parameter expanders, side-by-side
      results, per-component + additive residual risk, combined download.
- [x] `risk_days_prep_bs()` with `return_sim_df` (PrEP-specific columns); prep.py cleanup.
- [x] Integration tests — `tests/test_prep_bootstrap.py`.
- [x] Full Go port — `prep*.go`, `risk_days_prep_bs_go()` bridge, `use_go` dispatcher,
      cross-validation tests.

### truncnorm positivity fix (2026-05-29, `fix_truncnorm_positivity` → `main`)

- [x] **`truncnorm.rvs(0, np.inf, mean, sd)` truncated at the mean, not 0.**
      scipy's `a`/`b` bounds are in *standard deviations from* `loc`, so `a=0`
      maps to the lower bound `loc` — discarding the lower half of the
      distribution and inflating sampled means by ≈`0.8*sd`. Affected
      `_risk_days_bs_python` (`doubling_time`, `lod50`) and `residual_risk_rd`
      (`incidence`). Fixed via a shared `_sample_positive_normal(mean, sd, n)`
      helper (`a = -mean/sd`, truncating at 0, matching Go's
      `GenerateTruncatedNormal`). Baseline RDE median **−8.3%** (1.355 → 1.243),
      now matching Go to **0.6%**; prior results were biased conservative (high).
      Main-lineage library version 0.9.5 → 0.9.6 (subsumed by `1.1.0.dev0` on
      this branch). Regression guard: `TestSamplePositiveNormal`. *Surfaced by*
      the PrEP M2 Python↔Go agreement test; pre-existing since the original code.
      **PrEP (`prep.py`) had the same two sites — fixed in the merge commit
      (`31118a6`) using the merged `_sample_positive_normal` helper (+ `_vl_postbt`
      restructured so the growth exponential is only evaluated for `t ≤ tcrit`,
      avoiding a discarded overflow now that small `doubling_time` is sampleable).**

### Test runner: marker-based sandbox filter + dedup + thread tuning (2026-05-29, `feature_prep_model`, commit 01226b8)

- [x] **`scripts/run_tests.sh fast` now reliable.** The old `fast` mode used a
      fragile `-k` name-substring filter (`not (Python or Agreement or Bootstrap
      or agree_with_python)`) that was wrong both ways: it **dropped** the new
      sandbox-safe `TestPrepIntegrationMethod` (the whole `test_prep_bootstrap.py`
      file matched "Bootstrap") and **leaked** 6 `ProcessPoolExecutor` tests
      (`test_prep_k_distributions.py::TestPrepBsKDistributions`) that then failed
      in the sandbox. Replaced with an explicit `@pytest.mark.multiprocessing`
      marker on the pool tests, registered in `pyproject.toml`; `fast` now uses
      `-m "not multiprocessing"`. (Backported to `main` in commits 64a0910/6ba5668
      + dedup/threads in 8483702, so the merge was conflict-free on the test files.)
- [x] **Deduplicated within-file test names.** Renamed the ~20 colliding
      method occurrences across parallel backend/dist classes to be globally
      unique (e.g. `test_sanity_checks` → `…_python` / `…_go`; `test_lnmix` →
      `…_sampler` / `…_bootstrap`; `test_different_seeds_differ` → `…_theory` /
      `…_go`). No within-file duplicate test names remain.
- [x] **Thread tuning.** Added guarded `TEST_THREADS = max(1,
      get_cpu_core_count() - 1)` to the large-`n_bs` Python-pool suites
      (`test_residualrisk.py`, `test_invgamma_parity.py`, `test_lnmix_parity.py`)
      so they finish quickly outside the sandbox; tiny prep suites stay at
      `threads=1` (macOS uses `spawn`, so many workers for trivial `n_bs` would
      add re-import overhead). Worker count never affects results (sampling is in
      the main process; reduction is order-independent). Also removed two
      pre-existing unused imports (`scipy.stats`, `math`) flagged by ruff.

### PrEP integrator fix (2026-05-28, `feature_prep_model`, commit 82bab58)

- [x] **Gauss-Legendre default for the PrEP integrand.** `_risk_days_prep` and
      `risk_days_prep_bs` now take `integration_method` (`"gauss-legendre"`
      default | `"quad"`), mirroring the baseline. GL via the shared
      `_integrate_gauss_legendre` from `core.py`; `"quad"` preserves the
      historical `limit=500`; invalid value raises `ValueError`. Threaded through
      the bootstrap (`partial(_risk_days_prep, integration_method=…)`) and the
      primary-parameters PE call; `use_go=True` + `"quad"` raises (Go is always
      GL — not passed to the bridge). This is a genuine **correctness** fix: the
      PrEP integrand has compact support, so at a narrow serology window quad
      returns ~0 while GL recovers the true value:
      - narrow (α=9.1, β=5.2): GL **1.00865** (truth 1.00864) vs quad **5e-18**.
      - production (α=50.49434, β=1.15062): GL 3.09185 vs quad 3.09187 → **no
        production shift** (rel 6e-6).
      5 new sandbox-safe tests in `TestPrepIntegrationMethod`
      (`tests/test_prep_bootstrap.py`). AGENTS.md public-API note added (and
      `risk_days_prep_bs` listed).

### Baseline integration robustness (`fix_integration_quadrature` → `main`, merged 2026-05-28)

- [x] **Integration robustness.** Baseline `_risk_days` now defaults to a fixed
      1000-point Gauss-Legendre rule (`integration_method="gauss-legendre"`),
      matching the Go backend to machine precision; `integration_method="quad"`
      (scipy adaptive Gauss-Kronrod) remains selectable on the Python path for
      reproducing prior analyses (`use_go=True` + `quad` raises `ValueError`).
      Threaded through `_risk_days`, `_risk_days_bs_python`, `risk_days_bs`. Added
      an overflow guard in `_concentration` (caps `t/doubling_time` at 700).
      Corrected the misleading "adaptive quad" comment in
      `go/riskdays/integration.go`. 8 new sandbox-safe tests in
      `TestIntegrationMethod`. *Context:* scipy adaptive `quad` silently returns
      ~0 on **compact-support** integrands when its initial Gauss-Kronrod nodes
      miss the active window; the baseline integrand has noncompact
      (exponential-tail) support so quad was robust there — GL was adopted for
      Python↔Go parity and future-proofing. (The PrEP integrand does have compact
      support → handled by the PrEP integrator fix above.)

### PrEP model

- [x] k-distribution parity — extract shared `_sample_k()` helper in `core.py`; add InvGamma + LN-mixture kwargs to `risk_days_prep_bs()`; 9 unit tests in `tests/test_prep_k_distributions.py` (2026-05-21)
- [x] PrEP UI widgets — eclipse, oPrEP/iPrEP set points + ranges, oPrEP/iPrEP seroconversion Weibull params (2026-05-20)
- [x] Initial `residualrisk/prep.py` module — viral dynamics, infectivity, non-detection, bootstrap (2026-05-20)

### Standard mechanistic model

- [x] Plot histogram does not render after a successful **Mechanistic model** run — confirmed working (2026-05-20)
- [x] Port KDE mode estimation (`_kde_mode_log`) to Go (`go/riskdays/kde.go`)
- [x] Implement Inverse Gamma sampling in Go (`go/riskdays/random.go`)
- [x] Native InvGamma sampling in Python backend — `k_invgamma_alpha` / `k_invgamma_beta` / `k_invgamma_mode` kwargs in `core.py` and `_go.py`
- [x] Rename `sample_invgamma` params from `a`/`scale` to `alpha`/`beta` throughout
- [x] Parity test suite — `tests/test_invgamma_parity.py` (32 tests: theoretical moments, KS, Go/Python agreement)
- [x] InvGamma UI wiring in `app.py` — parameter inputs, PE selectbox (mode/median/mean), mode-source radio, cached KDE modes
- [x] Dynamic alpha constraints — min lowered to 0.01; mean PE disabled when α ≤ 1; help text explains tail weight
- [x] Documentation updates — `AGENTS.md` and `README.md` updated for InvGamma, public API surface, k distribution status
- [x] Lognormal mixture — Go: `GenerateLogNormalMixture` in `random.go`; struct fields + `Validate()` in `models.go`; dispatch in `riskdays.go`; 9 Go tests
- [x] Lognormal mixture — Python: `sample_lnmix()` in `core.py`, `_risk_days_bs_python()` dispatch, `risk_days_bs()` public kwargs; exported from `__init__.py`
- [x] Lognormal mixture — Bridge: lnmix kwargs + JSON serialisation + sim_df regeneration in `_go.py`
- [x] Lognormal mixture — UI: mixing-weight slider, advanced component-param editing, PE selectbox (mode/median/mean), derived statistics caption; all scoping variables set in non-lnmix paths
- [x] Lognormal mixture — Parity tests: `tests/test_lnmix_parity.py` (32 tests: theoretical stats, component isolation, KS, Go-only sanity)
- [x] Version bumps: library 0.9.3, app 0.9.3, Go binary 0.9.3
- [x] Go KDE mode via `--kde-mode` subcommand in `go/main.go`; `mode_kde_go()` in `residualrisk/_go.py` exported from `__init__.py`; `load_data()` in `app.py` uses Go KDE (~0.9s for all 3 posteriors, 30× faster than Python KDE) with hardcoded fallback if Go binary absent

---

## Deferred

### Design considerations for future distributions

If additional distributions are needed beyond InvGamma and LN-mixture, consider
refactoring from individual kwargs to a general dispatch API:
`k_distribution` (string) + `k_dist_params` (dict). Not needed for two parametric
distributions but would be cleaner for N > 3.

---

## Notes

- `k_gamma_shape` / `k_gamma_scale` in Python and Go are legacy — kept for backward compat,
  deprecated in comments. They sample from plain Gamma, not Inverse Gamma. Do not extend.
- Python and Go use independent RNGs (NumPy PCG-64 vs Gonum xoshiro256++): same seed gives
  reproducible results *within* each implementation but not *across* them. This is expected
  and documented in `_go.py`.
- The companion analysis repo (`../residualrisk_analysis/exploration/K_PARAM_INPUTDIST.md`)
  documents the rationale for choosing InvGamma(α=2, β=0.002019) and the 90/10 mixture as
  the two primary recommendations.
- Lognormal mixture default parameters:
  - Component 1 (human): w=0.90, μ₁=−7.2403, σ₁=0.3241
  - Component 2 (animal): w=0.10, μ₂=−3.7423, σ₂=0.5258
  - Mixture mode ≈ 0.000649, median ≈ 0.000750, mean ≈ 0.003389
