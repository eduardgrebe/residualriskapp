# TODO

## Open

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

- [ ] **Integration tests for `risk_days_prep_bs()`.** Full bootstrap tests
      exercising each k-distribution path. Unit tests for `_sample_k()` already
      exist in `tests/test_prep_k_distributions.py`.

### Deferred (post-initial-release)

- [ ] **Go acceleration for PrEP model.** No PrEP code exists in `go/`. Without
      it the PrEP bootstrap is ~10-50× slower than the standard model. Porting
      requires: viral dynamics, serology non-detection, integration → Go.

### Pre-existing on `main` (file against `main`, not this branch)



---

## Completed

### PrEP model

- [x] k-distribution parity — extract shared `_sample_k()` helper in `core.py`; add InvGamma + LN-mixture kwargs to `risk_days_prep_bs()`; 9 unit tests in `tests/test_prep_k_distributions.py` (2026-05-21)
- [x] PrEP UI widgets — eclipse, oPrEP/iPrEP set points + ranges, oPrEP/iPrEP seroconversion Weibull params (2026-05-20)
- [x] Initial `residualrisk/prep.py` module — viral dynamics, infectivity, non-detection, bootstrap (2026-05-20)

### Standard mechanistic model

- [x] Plot histogram does not render after a successful **Mechanistic model** run — confirmed working (2026-05-20)
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
