# TODO

## In progress

### Technical documentation — `docs/theory.md` (base mechanistic + lookback models)

Scientific-style technical documentation for the **base** model (on `main`; PrEP to
follow on its branch). Sources: the 2020 notebooks (`residualrisk_analysis/notebooks/
Residual Risk.ipynb`, `IWP from Operational Data.ipynb`) and the published methodology
(`literature/Grebe et al. 2020 Blood RR Appendix.pdf`) — reuse that text where possible,
updated for the improvements now in the package/Go (GL integration default,
truncnorm-positivity fix, k input-distribution options, KDE-log mode, lookback
Jeffreys-posterior samples, Go acceleration).

**Initial draft complete (2026-05-29)** — `docs/theory.md` (§§1–9) plus three
generated figures (`docs/figures/`, with `make_theory_figures.py` generator).
Decisions used: GitHub-flavoured `$…$`/`$$…$$` KaTeX (renders in Streamlit too);
comprehensive methods chapter (§4); k input distributions documented in full (§5.2);
figures generated. **Awaiting EG review.**

- [x] Confirm format/scope decisions (rendering, detail depth, k-dist in full, figures)
- [x] §1 Introduction & scope · §2 Notation table
- [x] §3 Mechanistic IWP model (viral dynamics; dose-response; NAT non-detection; IWP
      integral; worst-case note)
- [x] §4 Numerical methods (1000-pt Gauss-Legendre vs quad; limits/overflow; Go + RNG
      independence)
- [x] §5 Uncertainty analysis (truncnorm-positivity fix; k input distributions in full;
      point-estimate options incl. KDE-log mode; credible intervals)
- [x] §6 IWP → residual risk (incidence × IWP / 365.25; 1-in-x; incidence uncertainty)
- [x] §7 Lookback / operational-data model (Poisson/χ² estimator + exact CI;
      Jeffreys-Gamma posterior; CI vs CrI distinction)
- [x] §8 Default parameters & worked example · §9 References
- [x] Three explanatory figures (dose-response; IWP construction; k input distributions)
- [x] Wired into the app: `pages/1_Documentation.py` renders `theory.md` as a
      collapsible accordion — each top-level `## ` section in an `st.expander`
      (titles form the TOC, §1 open by default), with figure references split out to
      `st.image` (st.markdown does not load relative local images) and text/LaTeX
      rendered via `st.markdown`
- [x] Visual check that LaTeX renders in `st.markdown` — **confirmed (2026-05-29)**
      via rendered screenshots: display/inline KaTeX and `$…$` inside the
      notation/parameter/k-distribution **tables** all typeset correctly in
      Streamlit 1.58 (the one KaTeX-in-Streamlit nuance I could not verify headless).
      Also fixed a cosmetic figure defect (an overlapping rotated label in fig3) and
      regenerated the figure.
- [ ] **EG review pass** (verify equations, prose, citations; confirm worked-example
      numbers)
- [ ] _(later)_ PrEP model documentation — on the PrEP branch

## Completed

- [x] **truncnorm positivity fix (2026-05-29, `fix_truncnorm_positivity`).**
      `stats.truncnorm.rvs(0, np.inf, mean, sd)` truncates at the **mean**, not
      at 0 — scipy's `a`/`b` bounds are in *standard deviations from* `loc`, so
      `a=0` maps to the lower bound `loc`. This discarded the lower half of the
      intended distribution and inflated sampled means by ≈`0.8*sd`. Affected
      `_risk_days_bs_python` (`doubling_time`, `lod50`) and `residual_risk_rd`
      (`incidence`). Fixed via a shared `_sample_positive_normal(mean, sd, n)`
      helper (`a = -mean/sd`, truncating at 0, matching Go's
      `GenerateTruncatedNormal`). Effect at standard baseline params: RDE median
      **−8.3%** (1.355 → 1.243), now matching Go to **0.6%** (was an ~8% gap);
      sampled doubling_time 0.898 → 0.854 (nominal), lod50 2.88 → 2.74 (nominal).
      Current (buggy) results were biased conservative (high). Library version
      0.9.5 → 0.9.6 (calculation change; no Go change — Go was already correct).
      Regression guard: `TestSamplePositiveNormal`. *Surfaced by* the new PrEP M2
      Python↔Go
      agreement test on `feature_prep_model` (the first numerical Python↔Go PrEP
      *distribution* comparison); the bug is pre-existing since the original code.
      **PrEP (`prep.py`) has the same two sites — fixed separately on
      `feature_prep_model` after this merges to `main`.**

- [x] **Integration robustness (2026-05-28, `fix_integration_quadrature`).**
      Baseline `_risk_days` now defaults to a fixed 1000-point Gauss-Legendre
      rule (`integration_method="gauss-legendre"`), matching the Go backend to
      machine precision; `integration_method="quad"` (scipy adaptive
      Gauss-Kronrod) remains selectable on the Python path for reproducing
      prior analyses (`use_go=True` + `quad` raises `ValueError`). Threaded
      through `_risk_days`, `_risk_days_bs_python`, `risk_days_bs`. Added an
      overflow guard in `_concentration` (caps `t/doubling_time` at 700) so the
      GL rule, which always samples near the upper limit, can't raise
      `OverflowError` at small `doubling_time`. Corrected the misleading
      "adaptive quad" comment in `go/riskdays/integration.go`. Versions bumped
      to 0.9.5 (library + Go). 8 new sandbox-safe tests in
      `TestIntegrationMethod`. *Context:* scipy adaptive `quad` silently returns
      ~0 on compact-support integrands when its initial Gauss-Kronrod nodes miss
      the active window; the baseline integrand has noncompact (exponential-tail)
      support so quad was actually robust there, but GL is adopted for
      Python↔Go parity and future-proofing. The PrEP integrand DOES have compact
      support (eclipse + serology cutoffs) and is where quad genuinely fails —
      apply the same GL fix on `feature_prep_model`.
      *NOTE:* the 15 `ProcessPoolExecutor` bootstrap tests can't run in the
      sandbox (SemLock `PermissionError`); run `pytest` outside the sandbox to
      confirm. No hardcoded absolute-value goldens depend on the method, and the
      quad→GL shift is ~5e-8 (within all tolerances), so they're expected green.
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
