# TODO: Add Additional k Input Parameter Distributions

**Created:** 2026-04-27
**Status:** Not started

Based on the analysis in [`../residualrisk_analysis/exploration/K_PARAM_INPUTDIST.md`](../residualrisk_analysis/exploration/K_PARAM_INPUTDIST.md).

---

## Target Distributions

| # | Distribution | scipy params | Mode | Median |
|---|-------------|-------------|------|--------|
| A | Inverse Gamma | `invgamma(a=2.0, scale=0.002019)` | 0.000673 | 0.001203 |
| B | Lognormal Mixture (90/10) | 90% `lognorm(s=0.3241, scale=exp(-7.2403))` + 10% `lognorm(s=0.5258, scale=exp(-3.7423))` | ≈0.000649 | ≈0.000750 |
| C | Lognormal | `lognorm(s=1.0, scale=exp(-6.3038))` | 0.000673 | 0.001829 |

---

## 1. Python Core — `residualrisk/core.py`

### 1.1 Replace `k_gamma_shape`/`k_gamma_scale` with general distribution dispatch

Currently the only parametric option is `np.random.gamma` (plain Gamma, not Inverse Gamma as intended). Replace with a `k_distribution` string + `k_dist_params` dict.

- [ ] Add `_sample_k(k_distribution, k_dist_params, n_bs, *, k_posterior_sample=None, seed=None) -> np.ndarray` — internal helper that dispatches to the correct sampler:
  - `"posterior"` — `np.random.choice(k_posterior_sample, n_bs)` (backward compat)
  - `"invgamma"` — `scipy.stats.invgamma(a=params["a"], scale=params["scale"]).rvs(n_bs)`
  - `"lognormal"` — `scipy.stats.lognorm(s=params["s"], scale=params["scale"]).rvs(n_bs)`
  - `"lognormal_mixture"` — Bernoulli(1−w) to pick component, then draw from the appropriate lognormal
- [ ] Update `_risk_days_bs_python()` signature: replace `k_gamma_shape`/`k_gamma_scale` with `k_distribution=None` and `k_dist_params=None`
- [ ] Wire `_sample_k()` call in place of the current `if/elif/else` k-sampling block
- [ ] **Backward compat:** if `k_posterior_sample` is passed without `k_distribution`, default to `"posterior"` mode. If old `k_gamma_shape`/`k_gamma_scale` are passed, auto-convert to `k_distribution="gamma"` with appropriate params and emit a deprecation warning.
- [ ] **Primary parameters PE mode:** when `point_estimate="primary parameters"`, the PE integration needs the distribution's *mode* as the single `k` value (not a random draw). Add `_k_distribution_mode(k_distribution, k_dist_params, k_posterior_sample=None) -> float` helper.

### 1.2 Update `risk_days_bs()` (public wrapper)

- [ ] Mirror the same signature changes (replace gamma params with k_distribution/k_dist_params)
- [ ] Pass through to both Python and Go paths

---

## 2. Public API — `residualrisk/__init__.py`

- [ ] Consider exporting distribution parameter constants as a convenience:

```python
K_DIST_INVGAMMA = {"a": 2.0, "scale": 0.002019}
K_DIST_LOGNORMAL = {"s": 1.0, "scale": 0.001829}  # scale = exp(μ) = exp(-6.3038)
K_DIST_LOGNORMAL_MIXTURE = {
    "w": 0.90,
    "mu1": -7.2403, "sigma1": 0.3241,
    "mu2": -3.7423, "sigma2": 0.5258,
}
```

- [ ] Add to `__all__` if exported

---

## 3. Go Implementation — `go/riskdays/`

### 3.1 `models.go`

- [ ] Add new fields to `RiskDaysInput`:
  ```go
  KDistribution string             `json:"k_distribution,omitempty"`
  KDistParams   map[string]float64 `json:"k_dist_params,omitempty"`
  ```
- [ ] Update `SetDefaults()` — default `KDistribution` to `"posterior"` if empty and `KPosteriorSample` is set
- [ ] Update `Validate()`:
  - Require valid `KDistribution` + `KDistParams` combination
  - Valid values: `"posterior"`, `"invgamma"`, `"lognormal"`, `"lognormal_mixture"`
  - Keep backward compat: empty `KDistribution` + populated `KPosteriorSample` → `"posterior"`
- [ ] Deprecate `KGammaShape`/`KGammaScale` (keep in struct for now, add deprecation comment)

### 3.2 `random.go`

- [ ] Add `GenerateInvGamma(alpha, beta float64, n int) []float64`
  - Gonum has `distuv.InverseGamma` — verify availability
  - Fallback: sample from `Gamma(alpha, 1/beta)` and return `1/x`
- [ ] Add `GenerateLogNormal(mu, sigma float64, n int) []float64`
  - Use `distuv.LogNormal{Mu: mu, Sigma: sigma, Src: rg.rng}`
- [ ] Add `GenerateLogNormalMixture(w, mu1, sigma1, mu2, sigma2 float64, n int) []float64`
  - Use `rand.Float64() < w` to pick component, then draw

### 3.3 `riskdays.go`

- [ ] **Dispatch k sampling in `RiskDaysBS()`** — replace the current `if/else` on `KPosteriorSample` vs `KGammaShape`/`KGammaScale` with a `switch input.KDistribution` block:
  ```go
  switch input.KDistribution {
  case "posterior":
      ks = rng.BootstrapChoice(input.KPosteriorSample, input.NBS)
  case "invgamma":
      ks = rng.GenerateInvGamma(input.KDistParams["a"], input.KDistParams["scale"], input.NBS)
  case "lognormal":
      ks = rng.GenerateLogNormal(input.KDistParams["mu"], input.KDistParams["sigma"], input.NBS)
  case "lognormal_mixture":
      ks = rng.GenerateLogNormalMixture(
          input.KDistParams["w"],
          input.KDistParams["mu1"], input.KDistParams["sigma1"],
          input.KDistParams["mu2"], input.KDistParams["sigma2"],
          input.NBS,
      )
  }
  ```
- [ ] **Primary parameters PE mode:** compute the distribution *mode* and use it as the single `K` for the integration:
  - invgamma: `beta / (alpha + 1)` → `input.KDistParams["scale"] / (input.KDistParams["a"] + 1)`
  - lognormal: `exp(mu - sigma²)` → `math.Exp(input.KDistParams["mu"] - input.KDistParams["sigma"]*input.KDistParams["sigma"])`
  - lognormal_mixture: use the human component's mode for simplicity (≈0.000649), or compute KDE mode numerically
  - posterior: same as before (use `input.K` which is the point estimate passed in)

### 3.4 `riskdays_test.go`

- [ ] Add `TestGenerateInvGamma` — verify shape: mean of 10k samples ≈ β/(α−1) when α > 1
- [ ] Add `TestGenerateLogNormal` — verify log(samples) are normal
- [ ] Add `TestGenerateLogNormalMixture` — verify mixing ratio approximately correct (~90% from component 1)
- [ ] Add integration test: `TestRiskDaysBS_InvGamma` — run full BS with invgamma, check outputs
- [ ] Add integration test: `TestRiskDaysBS_LogNormal` — run full BS with lognormal
- [ ] Add integration test: `TestRiskDaysBS_LogNormalMixture` — run full BS with mixture

---

## 4. Go/Python Bridge — `residualrisk/_go.py`

- [ ] Build JSON input with `k_distribution` and `k_dist_params` instead of `k_gamma_shape`/`k_gamma_scale`:
  ```python
  "k_distribution": "invgamma",
  "k_dist_params": {"a": 2.0, "scale": 0.002019},
  ```
- [ ] **Update `sim_df` regeneration** (the Python-side re-sampling for the downloadable DataFrame):
  - Call `_sample_k(k_distribution, k_dist_params, n_bs, k_posterior_sample=k_posterior_sample, seed=seed)` to regenerate `ks`
  - This ensures the sim_df's `k` column matches what Go used (same distribution + seed)
  - Note: Python and Go RNGs differ, so values won't match exactly. Document this.
- [ ] Add backward compat: if `k_gamma_shape`/`k_gamma_scale` are passed without `k_distribution`, auto-convert to `k_distribution="gamma"`

---

## 5. UI — `app.py`

### 5.1 Replace/extend transmissibility model dropdown

- [ ] Replace the three-option "Transmissibility model" selectbox with a more flexible UI:
  - **New selectbox: "k input method"** with options:
    1. "Posterior sample (human)" — existing
    2. "Posterior sample (animal)" — existing
    3. "Posterior sample (human-weighted / expdecay)" — existing
    4. "Inverse Gamma distribution" — new
    5. "Lognormal mixture (90% human / 10% animal)" — new
    6. "Lognormal distribution" — new
  - Keep the existing "Transmissibility parameter: posterior..." selectbox (mean/median/mode) but only show it for posterior options
  - For parametric options, compute the PE from the distribution mode directly

### 5.2 Add distribution parameter inputs

- [ ] When "Inverse Gamma" is selected, show:
  - `a` (shape) — default 2.0
  - `scale` (β) — default 0.002019
- [ ] When "Lognormal" is selected, show:
  - `s` (σ) — default 1.0
  - `mu` (log-scale mean) — default −6.3038
- [ ] When "Lognormal mixture" is selected, show:
  - `w` (weight of human component) — default 0.90
  - Component 1 (human): `mu1` (−7.2403), `sigma1` (0.3241)
  - Component 2 (animal): `mu2` (−3.7423), `sigma2` (0.5258)

### 5.3 Wire into simulation call

- [ ] In the "Run simulations" button handler, build `k_distribution` and `k_dist_params` based on UI state
- [ ] Pass them to `rr.risk_days_bs()` instead of (or in addition to) `k_posterior_sample`
- [ ] Update `k_param_pe` logic — for parametric distributions, compute mode analytically

---

## 6. Python Tests — `tests/test_residualrisk.py`

- [ ] **`test_sample_k_posterior`** — verify backward compat: passing `k_posterior_sample` without `k_distribution` works
- [ ] **`test_sample_k_invgamma`** — sample 10k values, check median ≈ 0.0012, mode ≈ 0.00067
- [ ] **`test_sample_k_lognormal`** — sample 10k values, check log(samples) normal-ish
- [ ] **`test_sample_k_lognormal_mixture`** — sample 10k values, verify ~90% from human component
- [ ] **`test_risk_days_bs_invgamma_python`** — end-to-end with invgamma, n_bs=500, verify results are reasonable
- [ ] **`test_risk_days_bs_lognormal_python`** — end-to-end with lognormal
- [ ] **`test_risk_days_bs_mixture_python`** — end-to-end with mixture
- [ ] **`test_risk_days_bs_invgamma_go`** — Go path (skip if binary not found)
- [ ] **`test_risk_days_bs_lognormal_go`** — Go path
- [ ] **`test_risk_days_bs_mixture_go`** — Go path
- [ ] **`test_k_distribution_cross_impl`** — Python vs Go consistency (within 15% tolerance, matching existing tests)
- [ ] **Edge cases:**
  - `n_bs = 1` with each distribution
  - Mixture `w = 0` and `w = 1` (degenerate, should match single lognormal)
  - Very small `scale` parameters

---

## 7. Documentation

- [ ] **Update `AGENTS.md`** — document the new `k_distribution` / `k_dist_params` API in the Technical Architecture section
- [ ] **Update `README.md`** — add brief description of the new k input options under key parameters
- [ ] **Add docstrings** to:
  - `_sample_k()` and `_k_distribution_mode()` in `core.py`
  - `GenerateInvGamma`, `GenerateLogNormal`, `GenerateLogNormalMixture` in `random.go`
  - New fields in `models.go`

---

## 8. Cross-Cutting Concerns

- [ ] **Seed reproducibility:**
  - Verify Python and Go produce identical `ks` for the same distribution + seed
  - Current RNGs: Python uses NumPy PCG-64; Go uses `golang.org/x/exp/rand` (xoshiro256++)
  - These will NOT produce identical samples. Document this clearly in `_go.py` and test tolerances.
  - The important thing is that running twice with the same seed on the same implementation gives identical results.

- [ ] **Build Go binary:** after all Go changes, run `bash scripts/build_go.sh` and verify:
  ```bash
  echo '{"k_distribution": "invgamma", "k_dist_params": {"a": 2.0, "scale": 0.002019}, "n_bs": 100, "threads": 1}' | go/bin/riskdays_go
  ```

- [ ] **Full app smoke test:** `streamlit run app.py` → test each k input method → verify results display correctly → download sims CSV and check k column.

---

## Suggested Order of Work

| Step | Component | Rationale |
|------|-----------|-----------|
| 1 | Python core (`core.py`) | Foundation — everything else depends on the sampling logic |
| 2 | Python tests | Validate sampling correctness before porting to Go |
| 3 | Go `random.go` + `models.go` | Add new sampling functions + data structures |
| 4 | Go `riskdays.go` | Wiring |
| 5 | Go tests | Validate Go side |
| 6 | `_go.py` bridge | JSON schema changes |
| 7 | `app.py` UI | New dropdown + params |
| 8 | End-to-end smoke test | Full integration check |
| 9 | Documentation | `AGENTS.md`, `README.md`, docstrings |

---

## Notes

- The existing `k_gamma_shape`/`k_gamma_scale` path in both Python and Go appears to be a leftover — it samples from a plain Gamma distribution, not Inverse Gamma. The analysis document does not recommend a Gamma distribution. Deprecate gracefully rather than removing (in case external callers exist).
- The three existing posterior sample options (animal, human, human-weighted) should remain available alongside the new distributions.
- The analysis document recommends InvGamma(α=2) and the 90/10 mixture as the two primary recommendations. The lognormal(σ=1.0) is a simpler alternative. UI defaults should reflect this prioritization.
