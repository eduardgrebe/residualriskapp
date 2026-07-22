# Changelog

Versions follow [PEP 440](https://peps.python.org/pep-0440/). Git tags track the **app**
version (`app.py` → `APP_VERSION`); the `residualrisk` Python library and the Go binary carry
their own independent version numbers (`residualrisk/__init__.py` → `__version__` and
`go/riskdays/version.go` → `Version`).

## Versioning note — why the first public line is 1.1.0

The **1.1.0** line — which adds the PrEP-breakthrough-infection model — is the **first public
release line**. Earlier `0.1.x` / `0.9.x` tags were internal/preview builds of the baseline
(pre-PrEP) app. A `1.0.0` release was anticipated for that baseline but was **never cut**: when
the PrEP feature matured, its branch (`feature_prep_model`) was mainlined onto `main` and the
project moved straight to `1.1.0`. So `1.1.0` does **not** imply a public `1.0.0` predecessor.

## [Unreleased] — 1.1.0 (pre-release)

The PrEP-breakthrough model, mainlined onto `main` (2026-07). Working toward the first stable
`1.1.0`; currently at beta pre-release `1.1.0b4`. Highlights:

- Oral- and injectable-PrEP breakthrough-infection RDE components, with per-component and
  additive total residual risk.
- Canned NAT-assay limit-of-detection presets (single source of truth) alongside manual entry.
- Inverse-Gamma and lognormal-mixture input distributions for the infectivity parameter *k*.
- Go-accelerated PrEP bootstrap, with Python↔Go parity tests.
- Theme-aware sidebar branding, and assorted UI, documentation, and numerical-robustness fixes.
- **Beta 1 (`1.1.0b1`)** — pre-release-review hardening: guarded PrEP oscillation inputs; engine
  input validation across Python and Go (LoDs, doubling time, PrEP scalars); and accurate backend
  reporting (fallback logged, total-risk CrI label gated, `backend` column on simulation outputs).
- **Beta 2 (`1.1.0b2`)** — the reported RDE point estimate now defaults to the bootstrap **mode**
  (always within the credible interval) instead of the plug-in "primary parameters" value, which for
  right-skewed distributions could fall in the far right tail — even above the upper CrI; choosing
  "primary parameters" now shows a tail caveat. UI-only.
- **Beta 3 (`1.1.0b3`)** — **distribution & packaging.** The `residualrisk` wheel now bundles the
  pre-compiled Go accelerator for every platform (Linux/macOS/Windows × amd64/arm64), so
  `pip install residualrisk` runs Go-accelerated with **no Go toolchain** and no separate build, and
  the library is published to the project's **Codeberg PyPI registry** on each release. The app now
  shows a **pre-release warning banner** under the title on non-stable builds. Under the hood:
  tag-triggered release CI (Codeberg package + GitHub Docker image, with pre-release tags flagged as
  such on the Releases page) and a cross-OS bundled-binary smoke test.

- **Beta 4 (`1.1.0b4`)** — **model and API correctness.** Two changes alter
  behaviour for callers who reach beyond the app's defaults; neither moves a previously published
  estimate.
  - **PrEP: the viral-load trajectory is now continuous at the growth→plateau crossover.** `tcrit`
    solved for the time growth reaches the bare set-point, while the plateau *begins* at
    `offset × set_point` — so for any `offset ≠ 1` the modelled viral load jumped instantaneously
    by exactly a factor of `offset`. `tcrit` now targets the plateau's central level. **Bit-for-bit
    unchanged at the default `offset = 1`**, which every shipped and published result used. The
    **Offset control has been removed from the UI**: after the fix it is exactly a set-point
    multiplier — `(set_point, o, a)` ≡ `(set_point·o, 1, a/o)` — so it added no expressive power
    while being a first-order lever on the answer (~5× across its old range, since the breakthrough
    plateau sits near the pooled-NAT detection threshold and most of the PrEP risk accrues there).
    Vary the **set point** instead: clinical units, its own bootstrap range, and it drives `tcrit`
    correctly. The `offset` argument remains on the Python API.
  - **The four *k* input distributions are now mutually exclusive.** They were a silent priority
    cascade (posterior > gamma > invgamma > lnmix), so specifying two ran the higher-priority one
    with no warning — a stale `k_posterior_sample` in a reused parameter dict quietly turned an
    "InvGamma sensitivity analysis" back into the posterior. Specifying more than one, or a partial
    one, now raises. The app was never affected.
  - **Go acceleration is no longer reported when the binary cannot actually run.** Every candidate
    binary is smoke-tested (`--version`); a present-but-unrunnable one (wrong architecture, corrupt,
    missing a shared library) is skipped instead of being announced as available and then silently
    degrading to the 10–50× slower Python engine at the first calculation.
  - The RDE **integration domain** (`limits`, default `[-100, 500]` days) is exposed on both
    backends — previously hardcoded in Go and undrivable from Python, so a PrEP `ser_max` beyond
    500 days was silently truncated with no way to widen it.
  - Fixes: the baseline point estimate silently ignored a custom `z`; the Go simulation output did
    not record `limits`; a nondeterministic `ModeRounded` tie-break; the Go KDE `cap` subsampled
    *with* replacement where Python does so *without*.
  - Testing and CI: a bootstrap golden-regression test (per backend, plus a cross-backend agreement
    check), a PrEP per-iteration recompute invariant, and a `docker build` check on pull requests and
    `main` — the Docker image was previously built only on release tags, so a broken Dockerfile
    surfaced at release time.

- **Beta 5 (`1.1.0b5`)** — **Dependency bump release.** Chore release with updated dependencies to ensure final release will use current versions of all packages.

_The detailed record lives in the git tag history (`v1.1.0a1` …) and `TODO.md`._
