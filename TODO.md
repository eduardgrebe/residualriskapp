# TODO

## Open

### SESSION STATE (2026-05-28)

`main`'s baseline integration fix is merged into `feature_prep_model`; all three
versions at `1.1.0.dev0`. H1/H2/M1 committed. **The PrEP integrator fix is now
DONE but UNCOMMITTED** (see Completed → "PrEP integrator fix"); 3 files changed:
`residualrisk/prep.py`, `tests/test_prep_bootstrap.py`, `AGENTS.md`. Remaining
before the PrEP release: M2/M3 (next item), the `a`/`b`/`offset` design decision,
and end-to-end validation. (Agents can't commit/push — give the user the exact
git commands.)

### ⚠️ M2/M3 — PrEP Go golden + real Python↔Go parity test

**Status:** OPEN. Now unblocked by the PrEP integrator fix (GL on both sides).

- **M3 (Go test, `go/riskdays/prep_test.go`):** keep the 1.0086 golden but
  correct its comment — it was validated vs Python *Simpson* / Go, NOT vs Python
  `quad` (which returns 5e-18 at those narrow-window params). Add a second golden
  at production serology (≈3.0919, where all methods agree).
- **M2 (Python↔Go parity, `tests/test_prep_go_parity.py`):** add a real
  `risk_days_prep_bs(use_go=False)` vs `use_go=True` comparison (median / CrI
  within tolerance). NOTE: Python and Go agree only to **~1e-5** for PrEP, not
  machine precision, because Python still computes `tcrit` on a 0.1-day grid
  (`_vl_postbt_vec`) while Go uses the analytic `tcrit` (the open L2/L4 item).
  Use a tolerance like rel=1e-3, or do L2/L4 first for machine-precision parity.
  Requires ProcessPoolExecutor → run outside the sandbox.

### PrEP model — `feature_prep_model` branch

- [x] **Design decision: oPrEP vs iPrEP scenario handling.** *(Resolved 2026-05-21.)*
      `risk_days_prep_bs()` remains a single-scenario function.
      Populations don't overlap → risk is additive.

- [x] **Design decision: PrEP as additive layer on baseline risk.** *(Resolved 2026-05-21.)*
      PrEP breakthrough risk layers *on top of* baseline window-period risk.
      UI design:
      - Remove "Mechanistic model with PrEP" from the RDE method dropdown.
      - Add two independent checkboxes below the dropdown:
        ☐ Include oral PrEP breakthrough risk
        ☐ Include injectable PrEP breakthrough risk
      - When checked, corresponding PrEP parameter expander appears.
      - Run button always runs baseline first, then PrEP bootstraps if checked.
      - Results: baseline RDE (always), oPrEP RDE (if checked), iPrEP RDE
        (if checked), displayed side-by-side.
      - Residual risk: baseline incidence input (exists), plus separate PrEP
        breakthrough incidence inputs (prepopulated with defaults). Total
        residual risk = sum of all components.

- [x] **Wire PrEP into app.py.** *(Done 2026-05-21.)*
      Removed "Mechanistic model with PrEP" from dropdown. Added two independent
      checkboxes (Include oral PrEP / Include injectable PrEP) below the dropdown,
      only shown when mechanistic model is selected. Shared PrEP params expander
      (eclipse, a/b/offset) + per-modality expanders (set point, seroconversion
      Weibull). Button dispatch runs oPrEP and/or iPrEP bootstraps after baseline.
      Incidence section shows PrEP breakthrough incidence inputs when checkboxes ticked.
      Results show baseline + oPrEP + iPrEP RDEs; residual risk section shows
      per-component breakdown + additive total. Download button combines all scenario
      DataFrames with a `scenario` column.

- [x] **Add `return_sim_df` support to `risk_days_prep_bs()`.** *(Done 2026-05-21.)*
      Returns a per-iteration Polars DataFrame including PrEP-specific columns
      (set_point, eclipse, ser_*, a, b, offset) alongside standard columns.

- [x] **Code cleanup in `prep.py`:** *(Done 2026-05-21.)*
      Removed dead `_vl_noarv()`. Moved all inline imports to top-level.
      Integration limits `(-100, 500)` kept — sufficient for 265-day dynamics.

- [ ] **Design decision: bootstrapping `a`, `b`, `offset`.** The sinusoidal
      set-point oscillation parameters are currently fixed across all bootstrap
      iterations. Decide whether to add uncertainty distributions (e.g. uniform
      ranges) or keep them fixed. Scientific decision.

- [x] **Integration tests for `risk_days_prep_bs()`.** *(Done 2026-05-21.)*
      `tests/test_prep_bootstrap.py` — 39 tests across 7 classes: result
      structure (8), sim_df validation (8), reproducibility (2), point-estimate
      methods (3), k-distribution paths (6), PrEP parameter effects (4),
      edge cases (3), Go dispatch (5). Python-path tests require
      ProcessPoolExecutor (fail in macOS sandbox, pass outside).

- [x] **Go implementation of PrEP model.** *(Done 2026-05-21.)* Ported PrEP
      breakthrough infection model from Python to Go. All five phases complete:

  - [x] **Phase 1: Go PrEP functions** — `prep.go` (SinVaried, FindTcrit
        (analytic), VLPostBT, ProbInfectiousPrep, ProbNondetectionSerology,
        ProbNondetectionPrep, ProbInfectiousNondetectionPrep),
        `prep_models.go` (PrepInnerParams), `prep_integration.go`
        (RiskDaysPrep), `models.go` (PrEP fields, SetDefaults, Validate),
        `riskdays.go` (PrEP bootstrap path), `main.go` (PrEP binary output).
  - [x] **Phase 2: Go tests** — 30 tests in `prep_test.go` (unit +
        integration, golden value cross-validated against Python).
  - [x] **Phase 3: Python wrapper** — `risk_days_prep_bs_go()` in `_go.py`.
  - [x] **Phase 4: Dispatcher + app.py wiring** — `use_go` param in
        `risk_days_prep_bs()`; `app.py` passes `use_go=use_go_acceleration`.
  - [x] **Phase 5: Cross-validation** — `tests/test_prep_go_parity.py`
        (5 tests, all pass).

- [ ] **End-to-end validation of PrEP + baseline results against prior analyses.**
      Run the tool with input parameters from previous published baseline analyses
      (Grebe et al. 2020, ISBT 2024) and PrEP analyses (ISBT 2025) and confirm that 
      RDE estimates are comparable to previously reported values (or shifted 
      proportional to the major shift in k distribution used). If results
      differ, understand and document why (e.g. different k distribution,
      quadrature method, analytic vs grid tcrit). This applies to both the
      Python and Go implementations. The code passes unit/integration tests, but
      scientific validity requires reproducing known results with known inputs.

### Code review findings (2026-05-26)

From a review of the 14 PrEP commits in range `6a7f233..722ab36` (~3280 LOC).
Priority order, highest first.

#### High severity

- [x] **H1 — Lookback path doesn't reset `prep_oral_run` / `prep_inj_run`.**
      *(Done 2026-05-28 — Lookback branch in `app.py` now sets both
      `prep_oral_run` and `prep_inj_run` to `False` on a successful run.
      UNCOMMITTED.)*
      `app.py:1184–1219` (the `elif rde_method == "Lookback data":` branch in
      the button handler) sets `sims_run`, `rde_method_run`, `bs`, `samp`,
      `sim_df` but leaves the PrEP session keys untouched. After a
      Mechanistic+PrEP run followed by a Lookback run, the results section
      (`app.py:1333, 1351`) still renders the stale oPrEP / iPrEP RDE
      histograms, and the residual-risk section (`app.py:1408, 1440`) adds
      stale PrEP components into `_rr_total_pe` — so "Total residual risk
      (additive)" reports the wrong number.
      Fix: in the Lookback branch, set
      `st.session_state["prep_oral_run"] = False` and
      `st.session_state["prep_inj_run"] = False` (and ideally clear
      `iwp_pe_prep_*` / `samp_prep_*` / `sim_df_prep_*`).

- [x] **H2 — Python PrEP integrand crashes with `TypeError` when `a > offset`.**
      *(Done 2026-05-28 — chose fix option (b): clamp `max(0.0, Cv)` in Python
      `_vl_postbt` and Go `VLPostBT`. Physically correct (VL can't be negative),
      keeps full UI range, and fixes the Go negative-probability bug too.
      UNCOMMITTED.)*
      `_prob_nondetection_prep` in `residualrisk/prep.py:147–158` falls through
      both `if Cc == 0.0` and `elif Cc > 0.0` and returns `None` when `Cc < 0`.
      `Cc < 0` whenever `offset + a·sin(b·(t−tcrit)) < 0`, i.e. `a > offset`.
      The UI permits `prep_a ∈ [0, 2]` and `prep_offset ∈ [0, 2]`
      (`app.py:757–783`), so this is reachable. Verified by direct call:
      at `a=2.0, offset=1.0, t=32.24`, `_vl_postbt → −336`, the integrand
      raises `TypeError: unsupported operand type(s) for *: 'float' and 'NoneType'`.
      Go does not crash but `ProbInfectiousCopies(negative_n_copies, k) =
      1 − exp(positive) < 0` produces a negative probability.
      Fix options (cleanest first): (a) constrain UI so `a ≤ offset`;
      (b) clamp `Cv = max(Cv, 0)` in `_vl_postbt` / `VLPostBT`;
      (c) make `_prob_nondetection_prep` explicit on the negative branch.

#### Medium severity

- [x] **M1 — Go `SetDefaults` PrEP serology defaults don't match Python.**
      *(Done 2026-05-28 — serology defaults aligned to production
      `ser_min=28.7, ser_max=250, ser_alpha=50.49434, ser_beta=1.15062`; added
      `SetPointDistUniform=(19.1,2265)` and `EclipseDistUniform=(4.0,10.0)`
      defaults so a direct CLI caller matches Python. Go tests still pass.
      UNCOMMITTED.)*
      `go/riskdays/models.go:153–186` uses `SerMin=10, SerMax=500,
      SerAlpha=9.1, SerBeta=5.2`. Python `risk_days_prep_bs` (and the
      `risk_days_prep_bs_go` bridge, and the app UI) uses
      `ser_min=28.7, ser_max=250, ser_alpha=50.49434, ser_beta=1.15062`.
      Latent in production (bridge always passes explicit values from its own
      Python-matching defaults) but a direct CLI caller that omits these
      gets silently wrong defaults. Also: `SetPointDistUniform` and
      `EclipseDistUniform` have no defaults at all → direct caller omitting
      them gets degenerate (zero-width) sampling.
      Fix: align `SetDefaults` PrEP block with Python production defaults;
      add defaults (or `Validate()` errors) for the uniform-range fields.

- [ ] **M2 — `test_prep_go_parity.py` doesn't actually cross-validate Python vs Go.**
      *(NOT STARTED — 2026-05-28. Plan: add `TestPrepPythonGoAgreement` comparing
      `risk_days_prep_bs(use_go=False)` vs `use_go=True` on median/CrI within
      tolerance at PRODUCTION serology params, modeled on baseline
      `TestPythonGoAgreement`. At production params Python≡Go already holds
      (~3.09), so this test will PASS. A narrow-window parity case would expose
      the integration bug below and will only pass AFTER the quad→Gauss-Legendre
      fix — add it as part of that fix, not before. NOTE: Python PrEP path needs
      ProcessPoolExecutor → can't run inside the macOS sandbox; verify outside.)*
      Despite the filename and docstring ("Cross-validate Python and Go PrEP
      bootstrap results"), none of the five tests compares a Python-computed
      number to a Go-computed number — they are all Go-only sanity /
      reproducibility / dispatch checks. The baseline suite has a real
      `TestPythonGoAgreement` in `tests/test_residualrisk.py` (medians / PE
      within tolerance); the PrEP suite has no equivalent.
      Combined with M3 below, the Go PrEP implementation has never been
      numerically validated against the Python reference at production params.
      Fix: add `TestPrepPythonGoAgreement` modeled on the baseline (PE /
      median / quantiles within tolerance over a few hundred bootstraps at
      production serology defaults).

- [ ] **M3 — Go PrEP golden values use non-production serology params.**
      *(IN PROGRESS — 2026-05-28. Reference values computed (below); Go test not
      yet edited. `RiskDaysPrep` single-call at the standard non-serology params
      with serology varied:*
      *- OLD/test serology (α=9.1, β=5.2): Simpson truth = 1.008635, Go = 1.0086,
        Python `quad` = **5.0e-18** (quad fails — see integration bug above).*
      *- PRODUCTION serology (α=50.49434, β=1.15062): Simpson truth = 3.091868,
        Go = 3.09188, Python `quad` = 3.091867 (all agree).*
      *Plan: keep the 1.0086 golden with a CORRECTED comment ("validated vs
      Simpson/Go truth, NOT vs Python `quad`, which mis-integrates this narrow
      window") and add a second golden at production params (≈3.0919). Do this
      together with / after the integration fix.)*
      `go/riskdays/prep_test.go:25–48` `defaultPrepParams` uses
      `SerMin=10, SerMax=500, SerAlpha=9.1, SerBeta=5.2` (same numbers as M1).
      The "cross-validated against Python" golden values
      (`TestRiskDaysPrep_GoldenValue → 1.0086`,
      `TestProbInfectiousNondetectionPrep_CrossValidate → 4.153e-2`,
      `TestProbNondetectionSerology_InWindow → 0.9565`) are therefore
      validated at serology params the production code never uses. The
      production params (`α=50.49, β=1.15`) give a much wider, slower-decaying
      active integration window than the test params (`β=5.2` cuts sharply
      around t≈15–25 days). `quad.Fixed(integrand, -100, 500, 1000, nil, 0)`
      is only demonstrated to agree with Python Simpson within 1% in the
      narrow-window case.
      Fix: keep the existing golden value with an updated comment, and add a
      second golden value computed at production serology params.

#### Low severity

- [ ] **L1 — `point_estimate="mode"` differs between Python and Go.**
      `riskDaysBSPrep case "mode"` uses `KDEModeLog(rdests, 1_000_000, 0,
      input.Threads)` (`go/riskdays/riskdays.go:349–350`); Python
      `risk_days_prep_bs` uses `mode_rounded(rdests,
      precision=mode_precision)` (`residualrisk/prep.py:493`).
      `mode_precision` is sent to Go but ignored. PrEP faithfully mirrors
      baseline Go (`riskdays.go:187–188`) — so it is consistent with the
      pre-existing baseline pattern, but the cross-impl PE discrepancy is real
      for any user who picks "mode" and toggles Go/Python.
      Decide on one mode algorithm for both sides, or document that "mode"
      PE is implementation-dependent.

- [ ] **L2 — Python `_vl_postbt_vec` can raise on empty `argmin`.**
      `residualrisk/prep.py:43–54`:
      `idx = np.array([np.where(concentration > set_point)]).min()` raises
      `ValueError: zero-size array to reduction operation minimum` if
      exponential growth never exceeds `set_point` within
      `np.arange(0, 265, 0.1)`. In practice the sampled
      `(set_point, doubling_time, eclipse)` ranges keep `tcrit < 265`, so this
      is unreachable through the UI today, but the failure mode is silent.
      Go's analytic `FindTcrit` handles arbitrarily large `tcrit` gracefully.
      Either switch Python to the analytic formula (which is also L4 — a
      large perf win) or guard with `if not idx_arr.size: raise ValueError(...)`.

- [ ] **L4 — Python PrEP path needlessly slow (`_vl_postbt_vec` per-eval).**
      `_prob_infectious_prep` and `_prob_nondetection_prep` each rebuild a
      2650-element grid + sine array (`_vl_postbt_vec` with
      `np.arange(0, 265, 0.1)`) on every integrand evaluation just to extract
      `tcrit`. `_prob_infectious_nondetection_prep` calls both → tcrit is
      recomputed twice per evaluation, then discarded; the `conc_attenuated`
      (`tmp`) return value is never read anywhere. Replacing with the
      analytic `tcrit = eclipse + dt*log2(set_point/C0)` (as Go does) would
      cut Python PrEP integration time by a large factor and let
      `_vl_postbt_vec` be deleted. Also fixes L2.

- [ ] **L5 — Misc nits.**
      - `go/riskdays/version.go` not bumped despite adding the entire PrEP
        path to the Go binary (per `CLAUDE.md` versioning rules).
      - `app.py:1062` comment says "PrEP bootstrap runs (Python-only)" but
        the calls pass `use_go=use_go_acceleration`. Stale.
      - `tests/test_prep_go_parity.py` carries the abbreviated AGPL header
        instead of the full header used elsewhere (per `CLAUDE.md`
        "All new Python files must include the AGPL v3.0 license header").
      - `go/riskdays/integration.go:39` comment still says
        "adaptive quadrature / Equivalent to scipy.integrate.quad" but
        `quad.Fixed` is fixed-order Gauss-Legendre. Pre-existing, unrelated
        to PrEP, but worth fixing while in the area.

### Pre-existing on `main` (file against `main`, not this branch)



---

## Completed

### Test runner: marker-based sandbox filter + dedup + thread tuning (2026-05-29, `feature_prep_model`) — UNCOMMITTED

- [x] **`scripts/run_tests.sh fast` now reliable.** The old `fast` mode used a
      fragile `-k` name-substring filter (`not (Python or Agreement or Bootstrap
      or agree_with_python)`) that was wrong both ways: it **dropped** the new
      sandbox-safe `TestPrepIntegrationMethod` (the whole `test_prep_bootstrap.py`
      file matched "Bootstrap") and **leaked** 6 `ProcessPoolExecutor` tests
      (`test_prep_k_distributions.py::TestPrepBsKDistributions`) that then failed
      in the sandbox. Replaced with an explicit `@pytest.mark.multiprocessing`
      marker on the 69 pool tests (10 class-level + method-level on mixed
      classes), registered in `pyproject.toml`; `fast` now uses
      `-m "not multiprocessing"`. Deterministic partition: **159 safe / 69 pool /
      228 total**; the safe set now includes `TestPrepIntegrationMethod` (5) and
      `TestIntegrationMethod` (8). `bash scripts/run_tests.sh fast` is green in
      the sandbox (Go PASSED, Python 159 passed / 69 deselected, 0 PermissionError).
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

### PrEP integrator fix (2026-05-28, `feature_prep_model`) — UNCOMMITTED

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
      (`tests/test_prep_bootstrap.py`); also removed a pre-existing unused
      `import math` there. AGENTS.md public-API note added (and
      `risk_days_prep_bs` listed). No version bump (already `1.1.0.dev0`; no Go
      change — Go already uses GL). Python↔Go agree to ~1e-5 for PrEP (Python
      still uses grid `tcrit`; see L2/L4). *Files:* `residualrisk/prep.py`,
      `tests/test_prep_bootstrap.py`, `AGENTS.md`. *Verified in-sandbox:* 5 new +
      8 baseline integration tests pass, ruff clean, Go tests pass; the 69
      full-suite failures are all the known `ProcessPoolExecutor` SemLock
      `PermissionError` (run `pytest` outside the sandbox to confirm those).

### Baseline integration robustness (merged from `main` 2026-05-28)

- [x] **Integration robustness (`fix_integration_quadrature` → `main` → merged into `feature_prep_model`).**
      Baseline `_risk_days` now defaults to a fixed 1000-point Gauss-Legendre
      rule (`integration_method="gauss-legendre"`), matching the Go backend to
      machine precision; `integration_method="quad"` (scipy adaptive
      Gauss-Kronrod) remains selectable on the Python path for reproducing
      prior analyses (`use_go=True` + `quad` raises `ValueError`). Threaded
      through `_risk_days`, `_risk_days_bs_python`, `risk_days_bs`. Added an
      overflow guard in `_concentration` (caps `t/doubling_time` at 700) so the
      GL rule, which always samples near the upper limit, can't raise
      `OverflowError` at small `doubling_time`. Corrected the misleading
      "adaptive quad" comment in `go/riskdays/integration.go`. 8 new
      sandbox-safe tests in `TestIntegrationMethod`. *Context:* scipy adaptive
      `quad` silently returns ~0 on **compact-support** integrands when its
      initial Gauss-Kronrod nodes miss the active window; the baseline integrand
      has noncompact (exponential-tail) support so quad was actually robust
      there — GL was adopted for Python↔Go parity and future-proofing. **The
      PrEP integrand HAS compact support (eclipse + serology cutoffs) and is
      where quad genuinely fails → applying the same GL default to
      `_risk_days_prep` is the remaining open item on this branch (see Open).**
      *NOTE:* the `ProcessPoolExecutor` bootstrap tests can't run in the sandbox
      (SemLock `PermissionError`); run `pytest` outside the sandbox to confirm.

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
