# TODO

## Completed

- [x] Port KDE mode estimation (`_kde_mode_log`) to Go (`go/riskdays/kde.go`)
- [x] Implement Inverse Gamma sampling in Go (`go/riskdays/random.go`)
- [x] Native InvGamma sampling in Python backend — `k_invgamma_alpha` / `k_invgamma_beta` / `k_invgamma_mode` kwargs in `core.py` and `_go.py`
- [x] Rename `sample_invgamma` params from `a`/`scale` to `alpha`/`beta` throughout
- [x] Parity test suite — `tests/test_invgamma_parity.py` (32 tests: theoretical moments, KS, Go/Python agreement)

---

## Immediate priority: InvGamma UI wiring (`app.py`)

The InvGamma backend is complete; the UI currently sets `k_param = None` and tries
to call `rr.mode_kde(k_param)`, which fails for the InvGamma code path.

### What to do

- [ ] Extend the "Transmissibility model" selectbox to include **"Inverse Gamma distribution"**
  - Keep existing posterior-sample options unchanged
- [ ] When "Inverse Gamma" is selected, show parameter inputs (in the existing advanced options expander or a new one):
  - `alpha` (shape) — default 2.0
  - `beta` (scale) — default 0.002019
  - Display read-only computed mode = `beta / (alpha + 1)` for user reference
- [ ] Fix `k_param_pe` logic for InvGamma: compute mode analytically (`beta / (alpha + 1)`) rather than from a KDE on `k_param`
- [ ] In the simulation call, pass `k_invgamma_alpha` and `k_invgamma_beta` to `rr.risk_days_bs()` (and do NOT pass `k_posterior_sample`)
- [ ] Display results the same way as the posterior-sample paths (the output of `risk_days_bs` is identical in structure)
- [ ] Smoke-test: `streamlit run app.py` → select InvGamma → run → verify results display correctly → download sims CSV and check that k column contains positive values

---

## Deferred: Additional k input distributions

Hold until InvGamma UI is complete and confirmed working.

### Design decision required first

The InvGamma backend uses individual kwargs (`k_invgamma_alpha`, `k_invgamma_beta`).
Before adding Lognormal and Mixture, decide whether to:

- **Option A (individual kwargs):** add `k_lognormal_mu`, `k_lognormal_sigma`,
  `k_mixture_w`, `k_mixture_mu1`, `k_mixture_sigma1`, `k_mixture_mu2`, `k_mixture_sigma2` —
  consistent with the current InvGamma approach, no breaking change.
- **Option B (general dispatch):** refactor to `k_distribution` (string) +
  `k_dist_params` (dict) — cleaner for N distributions, requires updating all call sites.

### Target distributions

| # | Distribution | scipy params | Mode | Median |
|---|-------------|-------------|------|--------|
| B | Lognormal Mixture (90/10) | 90% `lognorm(s=0.3241, scale=exp(-7.2403))` + 10% `lognorm(s=0.5258, scale=exp(-3.7423))` | ≈0.000649 | ≈0.000750 |
| C | Lognormal | `lognorm(s=1.0, scale=exp(-6.3038))` | 0.000673 | 0.001829 |

### Tasks (pending design decision)

**Python `core.py`:**
- [ ] Add Lognormal sampling branch to `_risk_days_bs_python()` and `risk_days_bs()`
- [ ] Add LognormalMixture sampling branch (Bernoulli component selection)
- [ ] Add `_k_distribution_mode()` helper for PE calculation (invgamma, lognormal, mixture)

**Go `go/riskdays/`:**
- [ ] `random.go` — `GenerateLogNormal(mu, sigma float64, n int) []float64`
- [ ] `random.go` — `GenerateLogNormalMixture(w, mu1, sigma1, mu2, sigma2 float64, n int) []float64`
- [ ] `models.go` — add struct fields for lognormal and mixture params; update `Validate()`
- [ ] `riskdays.go` — add dispatch branches for lognormal and mixture
- [ ] `riskdays_test.go` — tests for `GenerateLogNormal`, `GenerateLogNormalMixture`, and integration tests
- [ ] Rebuild binary: `bash scripts/build_go.sh`

**`_go.py` bridge:**
- [ ] Add JSON-building branches for lognormal and mixture
- [ ] Update `sim_df` regeneration for lognormal and mixture

**UI `app.py`:**
- [ ] Add "Lognormal mixture" and "Lognormal" options to k input selectbox
- [ ] Show relevant parameter inputs per selection
- [ ] Wire PE calculation for each new distribution

**Tests:**
- [ ] `test_lognormal_parity.py` — same structure as `test_invgamma_parity.py`:
  theoretical stats, bootstrap k samples, Go/Python quantile agreement
- [ ] Edge cases: `n_bs=1`, mixture `w=0` / `w=1` (degenerate), very small scale params

---

## Documentation (do after InvGamma UI complete)

- [ ] `AGENTS.md` — document InvGamma params in Technical Architecture section
- [ ] `README.md` — add InvGamma to key parameters list
- [ ] Add to the above when additional distributions are implemented

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
