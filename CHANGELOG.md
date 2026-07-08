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
`1.1.0`; currently at beta pre-release `1.1.0b2`. Highlights:

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

_The detailed record lives in the git tag history (`v1.1.0a1` …) and `TODO.md`._
