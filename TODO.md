# TODO

## Open

### SESSION STATE (2026-05-29)

On `feature_prep_model`; `main` is fully merged in (merge commit `31118a6`),
working tree clean, all three versions at `1.1.0.dev0`. Everything below the
line — the PrEP model build-out, all code-review findings (H1, H2, M1, M2, M3,
L1, L2, L4, L5), the integration-robustness / truncnorm-positivity /
analytic-`tcrit` work, the marker-based test runner, and the sinusoidal
`a`/`b` uniform-uncertainty feature (commit `40a6cb2`), and the PrEP drug-effect
(transmissibility-reduction) parameter (commit `05b9e76`) — is **committed**.
`bash scripts/run_tests.sh fast` is robust (marker-based) on both branches.
(Note: the drug-effect `@multiprocessing` tests still need an outside-sandbox run;
in-sandbox we verified ruff, Go tests, 33 sandbox-safe prep tests, and the Go
bridge.)

This session ("prep-model-cleanup") then crystallised the **scope, validation, and
documentation** of the tool: the Layer 1 (tool) vs Layer 2 (Python API / effective
incidence) boundary, the one-product-per-run design, the assay-defaults task, and
the PK/PD deferred objective — see "## Scope & validation", "Assay defaults &
calibration", and Deferred. (`README.md` + `AGENTS.md` updated to match.)

The end-to-end scientific validation — the last PrEP-release gate — is now
**discharged** (2026-05-29); see Completed → "End-to-end PrEP validation".
(Agents can't commit/push — provide exact git commands.)

The base-model **technical documentation** (`docs/theory.md` + figures, rendered as
a collapsible accordion on the app's Documentation page) was written on `main` and
is merged in here; it still awaits EG review. **PrEP-model documentation** is now the
active documentation task on this branch — see In progress → "PrEP model
documentation".

### Remaining before the PrEP release

- _Nothing release-blocking outstanding._ The end-to-end PrEP validation gate is
  discharged (see Completed). The still-open items below (assay defaults, in-app
  guidance) and the Deferred items are quality/scope work, not release gates.
- [ ] _(Optional, non-blocking)_ **Baseline-only reproduction (Grebe et al. 2020
      / ISBT 2024).** The validation reproduced the **PrEP** pipeline against ISBT
      2025; the baseline window-period numbers were not separately reproduced
      (same engine, covered by unit tests + the truncnorm-fix tests). A short
      baseline replication in `../residualriskapp_validation` would close this too.

### Assay defaults & calibration

- [ ] **Re-check and populate default NAT assay parameters (LoD50 / LoD95).**
      Re-derive the tool's default assay parameters — LoD50, LoD95 (and the
      implied `lod95_lod50_ratio`) — from manufacturer **package inserts**.
      Assays to cover (at least):
    - Roche cobas TaqScreen **MPX**
    - Roche cobas TaqScreen **MPX v2.0**
    - Grifols Procleix **Ultrio Plus**
    - Grifols Procleix **Ultrio Elite**

      Specifics to get right:
    - **Revisit the IU/mL → copies/mL conversion** (the canonical home for the
      IU↔copies/mL LOD-conversion question raised in the `rr_prep_v3` review).
      Refer to the existing research already performed and applied in
      `~/dev/residual_risk/roche_rr_latam/` — that repo also holds the Roche
      package insert and the correct LoD numbers.
    - **HIV type / group specificity matters:** LoD differs by HIV-1 vs HIV-2 and
      by HIV-1 group (M vs O). **Provide defaults *only* for the most common HIV-1
      group (Group M).** Do **not** ship defaults for HIV-2 or the less common
      HIV-1 groups.
    - **WHO International Standard version:** analytical-sensitivity LoD values are
      calibrated against a specific WHO standard; the IU values (and hence the
      IU→copies conversion) depend on which standard/version was used. Record, per
      assay/insert, which WHO standard version the quoted LoD is calibrated to, and
      keep that consistent across assays and the conversion factor.

---

## Scope & validation

### Scope: what the tool models (Layer 1) vs what it does not (Layer 2)

**Decision (2026-05-29).** The tool — the `residualrisk` Python package, its Go
engine, and the Streamlit webapp — covers **Layer 1**: the mechanistic
window-period model that estimates **risk-day-equivalents (RDEs)**, plus a final
step that applies a **pre-computed incidence** to an RDE distribution to obtain
residual risk (`residual_risk_rd` = incidence × RDE / 365.25 × per).

It deliberately does **not** implement **Layer 2** — the population residual-risk
aggregation (donor-stratum PrEP-use prevalence, self-deferral / discard, sex- and
route-stratified incidence, product mix; full component list in the reference
below). Layer 2 is too complex and too operator-/time-specific to reflect in the
webapp for now; some aspects may be added in future (deferred — see Deferred →
"Population PrEP-use & stratified-incidence modelling in the webapp").

**Recommended approach for a sophisticated modelling exercise.** Build Layer 2 in
code using the **Python API** (which wraps the fast Go engine): call
`risk_days_prep_bs` for the RDE distributions and assemble the population layers
around them. For the **webapp**, fold the disaggregated inputs — PrEP-use
population proportions, stratified incidences, self-deferral / discard, etc. —
into a single **"effective incidence"** for each of oPrEP and iPrEP and supply
those as the fully-baked incidence inputs. The webapp's job is RDEs + applying
that pre-computed incidence; the disaggregation happens upstream.

### One product per run — by design

The tool estimates RDEs for **one transfused product at a time**, in both the
baseline and PrEP models. The user models each product separately by entering that
product's **transfused plasma volume and volume range** (`volume_transfused` /
`volume_transfused_range`) and running the tool once per product — e.g. red-cell
units (~20 mL residual plasma), fresh frozen plasma (~200 mL), or platelets (using
their own plasma-volume estimate). This is an **explicit design choice**: the tool
will **not** automate simultaneous multi-product RDE estimation — the user drives
that, and any cross-product mix-weighting is part of Layer 2 (above).

- [ ] **Documentation must make the one-product-per-run design explicit.** Done in
      `README.md` ✓ and `AGENTS.md` ✓; **remaining: in-app guidance / help text**
      in `app.py` (make clear that `volume_transfused` is the per-product
      transfused plasma volume, that the user runs the tool once per product —
      RBC / FFP / platelets / … — and that multi-product estimation is
      intentionally not automated).

### Validation

- [ ] **Replicate the `rr_prep_v3` analysis (with refinements) as an ad-hoc
      script driving the Python package.** Build Layer 2 on top of the package's
      RDE outputs (Python API → Go engine) and reproduce the prior published
      residual-risk numbers. This both (a) validates the mechanistic engine
      end-to-end against a known result and (b) is the concrete vehicle for the
      end-to-end validation release gate (Open → "Remaining before the PrEP
      release"). Expect known, explainable shifts vs the notebook: the
      truncnorm-positivity fix, Gauss-Legendre vs quad, analytic vs grid `tcrit`,
      and the k-distribution choice. Belongs in the analysis repo
      (`residualrisk_analysis`), not the tool.

### Layer 2 reference — population residual-risk aggregation

The components that translate per-window-period RDEs into population residual
risk, as implemented ad hoc in `compute_risks` (`rr_prep_v3.py`) and **not** in
the tool. Recorded here so we stay on top of it as the tool develops. (The
mechanistic RDE engine itself — `sin_varied`, `vl_postbt`, `prob_*`, `risk_days`
— is fully ported; the tool even improves on the notebook by fixing its
`truncnorm(0, inf, mean, sd)` positivity bug, defaulting to Gauss-Legendre, and
adding optional `a`/`b` + `drug_effect` variation.)

- **Donor strata & counts** — first-time vs repeat donors × male/female, with
  operator-specific annual donation counts.
- **PrEP-use prevalence by stratum** — oral/injectable × sex × FT/RD, with a
  repeat-donor:first-time ratio (≈1/13) and sex ratios (oral F:M ≈1/12,
  injectable ≈1/6).
- **Self-deferral rate** (`sd_rate`) — fraction of PrEP-using donors who
  self-defer and do not donate (prior: `Uniform(0.1, 0.7)`).
- **Disclosure / discard rate** (`disc_rate`) — fraction of non-self-deferred
  PrEP-user donations disclosed → discarded (prior: `Uniform(0.5, 0.75)`).
  Combined: at-risk donations = usage × (1 − `sd_rate`) × (1 − `disc_rate`).
- **Sex- and route-specific incidence** — separate male/female oral & injectable
  incidence, with an injectable:oral ratio (≈0.33); the tool/webapp take a single
  oPrEP and a single iPrEP incidence (hence the "effective incidence" guidance
  above).
- **Cross-product mix-weighting** — combining per-product residual risks (RBC,
  FFP, platelets, …) into a single population figure weighted by each product's
  share of transfusions. Note: the *per-product* RDE is by design a separate tool
  run (different plasma volume — see "One product per run" above); only the
  population mix-weighting is Layer 2.

The prior published RBC/FFP residual-risk numbers were produced by `compute_risks`
fed from the mechanistic RDE, so the validation task above reproduces them by
re-running that aggregation on the *tool's* RDE outputs — which also confirms the
tool's RDE is the correct input to the published pipeline.

---

## In progress

### Technical documentation — base model (`docs/theory.md`)

Scientific-style technical documentation for the **base** model (written on `main`,
merged here). Sources: the 2020 notebooks (`residualrisk_analysis/notebooks/
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

### PrEP model documentation (this branch)

Companion to the base-model documentation — the PrEP breakthrough model
(`residualrisk/prep.py` / `go/riskdays/prep*.go`), for the same scientific audience.

**Initial draft complete (2026-05-31)** — `docs/theory_prep.md` (§§1–12, self-contained,
cross-references `theory.md`) plus `docs/figures/make_prep_figures.py` and four generated
figures (pipeline, viral dynamics, NAT/serology detection windows, RDE construction).
AI-written from the code (`prep.py` / `prep*.go`), the analysis notebook (`rr_prep_v3.py`),
the serology-Weibull derivation (`weibull_serology_nondetection.R`, Seed 2021 data), and
the ISBT 2025 deck (`residualriskapp_validation`). **Structure decision:** written as a
*separate* file — the base `theory.md` states PrEP is "documented separately", and this
keeps the under-review base doc clean; it can be merged into `theory.md` later if one
document is preferred. **Awaiting EG review.**

- [x] Decide structure (separate `docs/theory_prep.md`; keep-vs-merge deferred to EG)
- [x] §1 Overview & scope (undisclosed PrEP; two-layer structure; residual-risk relationship)
- [x] §2 Notation (PrEP-specific symbols) · §3 Breakthrough viral dynamics (eclipse → growth →
      oscillating plateau; analytic `tcrit`; oscillation rationale; `a ≤ offset`)
- [x] §4 Detection on PrEP (NAT as base §3.3; serology Weibull + its Seed-2021 derivation;
      combined screen) · §5 Infectivity & the drug-effect factor
- [x] §6 PrEP RDE integral (compact support; GL necessity) · §7 oPrEP vs iPrEP (differences
      table; additive combination; structural caveat re injectable wash-out)
- [x] §8 Uncertainty/bootstrap · §9 Layer-2 population aggregation (documented; outside the tool)
- [x] §10 Default parameters & worked results (RDE + residual-risk tables; sensitivity) ·
      §11 Assumptions & limitations · §12 References
- [x] Figures — `docs/figures/make_prep_figures.py` (faithful to `prep.py`): pipeline,
      viral dynamics, NAT/serology detection, RDE construction
- [ ] **Wire into the app** — render `theory_prep.md` (a second Documentation page, or append
      it to the existing accordion page); decide alongside the keep-separate-vs-merge question
- [ ] **Resolve drafting flags** (marked "Reviewer note" / "to confirm" in the draft): the
      copies-vs-virions convention for `set_point`; the serology-derivation eclipse (6 d) vs
      model eclipse (7 d); injectable set-point provenance; incidence-input sources; and the
      Custer 2020 / Eshleman 2023 citations (DOIs unverified)
- [ ] **EG review pass (PrEP)** — verify equations, prose, citations, and worked numbers

---

## Completed

### End-to-end PrEP validation — release gate discharged (2026-05-29)

- [x] **End-to-end validation of the PrEP pipeline against the prior ad-hoc
      analysis and the published results.** Realised as a standalone repo,
      `../residualriskapp_validation` (paired Jupytext notebook
      `notebooks/validation.py` + executed `.ipynb`), which editable-installs this
      package and drives the **Go engine**. Findings:
    - **Reproduces `rr_prep_v3` end-to-end.** The package matches the frozen
      mechanistic RDE outputs (`rd_*_bs.npy`) for all four products
      ({RBC, FFP} × {oral, injectable}) up to a single, uniform **≈ −6% downward
      shift**, *fully attributed* to the truncated-normal positivity fix (buggy-
      truncnorm reconstruction matches the `.npy` to <0.3%; GL-vs-quad and
      analytic-vs-grid `tcrit` negligible; k handling identical via gamma fit).
    - **Matches the published ISBT 2025 base case.** Layer-2 aggregation
      (`compute_risks`, ported into the validation repo) on the original RDEs
      reproduces the presented figures (`presentations/`): RBC **1-in-110.8M** vs
      published **110** (95% CrI 24–1,437 vs 24–1,405); Plasma **1-in-74.6M** vs
      **75** (17–968 vs 17–942); increase-over-baseline +7.7% / +7.0% vs +7.7% /
      +6.9%. The current package gives ~6% *safer* figures (RBC 1-in-119M, Plasma
      1-in-79M) — the truncnorm-fix shift, expressed in the published metric.
    - **Dominant sensitivity = the k-distribution choice** (animal-derived vs
      human-anchored default InvGamma), which moves the estimate far more than the
      truncnorm fix — corroborated by the published tornado slides.
    - Note: the **baseline-only** comparison (Grebe 2020 / ISBT 2024) was not
      separately reproduced (optional follow-up — see Open). The full
      `bash scripts/run_tests.sh` suite (incl. the `ProcessPoolExecutor` tests)
      passes **outside the sandbox**.

### PrEP drug-effect (transmissibility reduction) parameter (2026-05-29, commit `05b9e76`)

- [x] **Optional antiretroviral drug-effect (transmissibility-reduction) factor
      for PrEP.** `drug_effect ∈ (0, 1]` (1.0 = no reduction) — a linear
      multiplier on the per-time infection probability, applied *inside the
      integrand* (`_drug_effect` in `prep.py` / `DrugEffectFactor` in `prep.go`),
      held fixed at the scalar unless `drug_effect_dist_uniform=(lo, hi)` is given
      (then sampled `Uniform(lo, hi)` per bootstrap iteration). Matches the prior
      analysis (`rr_prep_v3.py`: `Uniform(0.5, 1.0)`, median 0.75) — but because
      the factor is constant in `t` it factors out of the RDE integral, so the
      in-integrand placement is numerically identical to scaling the RDE while
      being the only correct placement should it ever become time-varying.
      `_drug_effect` / `DrugEffectFactor` take `t` as a deliberate placeholder for
      that future time-varying (long-acting-injectable wash-out) extension — see
      the deferred PK/PD task. Full-stack like the `a`/`b` feature: `prep.py`, Go
      (`prep.go` / `prep_models.go` / `models.go` / `riskdays.go` / `main.go` +
      `_go.py` bridge, with a per-iteration `drug_effect` binary column → correct
      Go-path `sim_df`), and `app.py` (per-scenario point input default 1.0 +
      bootstrap range slider default (1.0, 1.0); on by default, value 1.0 → no
      shift). `0 < drug_effect ≤ 1` and the range `⊂ (0, 1]` are enforced in
      Python and Go `Validate()` and capped in the UI. Default 1.0 leaves results
      bit-for-bit unchanged (verified: production PrEP RDE 3.09188 = the existing
      golden value; Python≡Go PE). Tests: direct-integration linearity +
      validation (`TestPrepIntegrationMethod`), full-bootstrap element-wise
      scaling + backward-compat + varied range (`TestPrepBsDrugEffect`), Go-path
      fixed/varied/linear (`TestPrepGoParity`), and
      `TestPrepPythonGoAgreementDrugEffect` (Python↔Go parity with `drug_effect`
      varied). Docs in `AGENTS.md`.

### PrEP sinusoidal `a`/`b` uncertainty (2026-05-29, commit `40a6cb2`)

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

### PK/PD-driven drug-concentration modelling for PrEP (longer-term)

A substantial modelling extension that needs careful design — a longer-term
objective, not an immediate task. Today the drug effect enters only as a
constant transmissibility-reduction scalar (`_drug_effect` in `prep.py` /
`DrugEffectFactor` in `prep.go` — a deliberate placeholder that takes `t` but
currently ignores it). A more faithful model would represent the drug
concentration itself over time and let the downstream quantities follow from it:

- **Explicit drug-concentration dynamics** for oral and injectable PrEP —
  fluctuating around a steady state under (sub-optimal) adherence for oral PrEP,
  and decaying after the last dose for long-acting injectables.
- **Couple viral load, detectability and transmissibility to concentration.**
  Theoretically these move *inversely* to drug concentration: as concentration
  falls the protective effect weakens, so breakthrough viral load rises — which
  in turn raises both NAT/serology detectability and transmission probability.
  They should therefore be derived together as functions of concentration rather
  than set independently. Promoting `_drug_effect` to a genuine function of `t`
  (e.g. an exponential wash-out from the last-injection time) is the narrow first
  step; the fuller version drives viral load from concentration as well.
- **Injectable PrEP is the key motivation.** Breakthrough infections on
  long-acting injectables typically occur as the drug washes out — precisely the
  regime a constant factor handles worst. This also surfaces a structural point:
  the current sigmoidal / oscillating-plateau viral-load model is most
  appropriate for **sub-optimal-adherence oral PrEP** and less so for
  **drug-washout injectable PrEP**, so the VL model *structure* itself may need
  rethinking for the injectable case.
- **Drug-specific potency.** Agents differ materially (e.g. the more effective
  lenacapavir vs. cabotegravir), so parameters — and possibly the functional
  form — should be drug-specific.
- **Prerequisite:** review the PK/PD literature for the relevant PrEP drugs to
  ground both the concentration dynamics and the concentration→effect
  relationships.

### Population PrEP-use & stratified-incidence modelling in the webapp

Deferred possible future work: surface some of Layer 2 — population PrEP-use
prevalence, self-deferral / discard, sex- and route-stratified incidence —
directly in the webapp, instead of requiring users to pre-compute an "effective
incidence." See **## Scope & validation** for the decision (Layer 2 stays out of
the webapp for now), the recommended Python-API approach, and the full Layer 2
component reference. (The IU/mL ↔ copies/mL LOD conversion surfaced in the same
review is tracked under Open → "Assay defaults & calibration".)

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
- *Idea (not a hard task):* consider **Mermaid** diagrams to make the algorithm/architecture
  documentation clearer — e.g. a flowchart of the IWP integrand (viral dynamics →
  dose-response × non-detection → integral → bootstrap → residual risk) and a component
  diagram (Streamlit app ↔ `residualrisk` package ↔ Go engine). Caveat: Mermaid renders
  natively in fenced `mermaid` code blocks on GitHub, but `st.markdown` does **not** render
  it, so the accordion Documentation page would need a helper (e.g. `streamlit-mermaid`) or
  the diagrams pre-rendered to SVG/PNG like the existing figures. Weigh that before adopting.
