# Structural Review: `estimator.py`

**File:** `estimator.py` (1885 lines)
**Date:** 2026-07-22
**Scope:** Overall structure, design quality, and footguns. No changes made.

---

## Overall Assessment

The file is **Streamlit-idiomatic** — one long procedural script that reruns top-to-bottom on every interaction. That's the expected pattern for Streamlit, not a flaw. The code is **well-commented** in the tricky spots (Go/Python backend alignment, seed callback timing, offset exposure rationale) and consistently uses the public API (`import residualrisk as rr`) rather than reaching into internals. Separation between UI and computation engine is clean.

At 1885 lines this file is doing a lot, and the procedural style creates specific fragility patterns.

---

## What's Good

- **Single source of truth** for assays (`rr.NAT_ASSAYS`), Go binary discovery, and KDE modes — no duplicated data.
- **`_prep_k_kwargs` dict** avoids repeating 10 k-distribution kwargs across the oral/injectable PrEP calls.
- **`_rr_row` helper** for the results table — clean, readable.
- **Session state keys for PrEP** are pre-initialized in a loop (lines 84–97) — good pattern.
- **`used_go` tracking** with the comment explaining why you gate on binary availability, not just the user's selection.
- **Stale-state cleanup** when switching to lookback (clearing PrEP flags).
- **`_prerelease_notice`** is a clean, well-documented utility.

---

## Footguns (ordered by severity)

### 1. Conditional variable definitions — NameError landmines

Many variables are defined only inside `if is_mechanistic_ui:` or `if rde_method == "Lookback data":` blocks but used in the run-button handler. This is *currently safe* because the handler checks the same conditions, but it's fragile:

| Variable | Defined in | Used in | Risk |
|---|---|---|---|
| `n_threads` | `if is_mechanistic_ui:` block | run handler | Safe today, breaks if a new code path references it without the guard |
| `use_go_acceleration` | `if is_mechanistic_ui:` block | run handler (lnmixture k_pe) | Same |
| `lod50`, `lod50_sd`, `lod95_lod50_ratio` | `if is_mechanistic_ui:` → `with model_param_container:` | run handler | Same |
| `pool_size`, `retests` | Same | Same | Same |
| `k_param`, `k_param_dist`, `k_param_pe` | Same | Same | Same |
| `inc_perpy`, `inc_perpy_sd` | `if calculate_rr:` inside incidence container | Display section under `if calculate_rr:` | Safe today, but NOT pre-initialized to `None` like the PrEP incidence vars are |

The pattern works because Streamlit reruns everything and the guards match. But one misplaced reference — say, someone adds a residual-risk display that references `inc_perpy` without checking `calculate_rr` — and you get a `NameError` at runtime.

**Mitigation:** Pre-initialize everything at the top (like you already do for `inc_prep_oral_perpy = None` etc.), or factor the variable-producing blocks into functions with explicit returns.

### 2. CrI/range recomputed from quantiles, ignoring library-returned values

The run handler unpacks all 5 returns from `risk_days_bs` and stores them:

```python
( st.session_state["iwp_pe"],
  st.session_state["iwp_cri"],   # ← set but never read
  st.session_state["iwp_range"], # ← set but never read
  st.session_state["bs"],
  st.session_state["sim_df"],
) = rr.risk_days_bs(...)
```

But the display code ignores `iwp_cri` and `iwp_range`, recomputing from the raw sample:

```python
iwp_cri = (
    st.session_state["samp"]["iwp"].quantile(alpha / 2),
    st.session_state["samp"]["iwp"].quantile(1 - alpha / 2),
)
iwp_range = (
    st.session_state["samp"]["iwp"].min(),
    st.session_state["samp"]["iwp"].max(),
)
```

Today these agree (both are quantiles of the same bootstrap array). But if the library ever changes its CrI method (e.g., HPD interval, BCa), the display silently diverges from what the library computed. The stored `iwp_cri` and `iwp_range` are set but never consumed.

The same recomputation happens for oral and injectable PrEP results too — six redundant quantile calls per display.

### 3. `k_pe` fallback chain silently defaults instead of raising

In the InvGamma branch:

```python
else:
    k_pe = k_invgamma_beta / (k_invgamma_alpha + 1)  # fallback to mode
```

If `k_invgamma_pe_choice` is somehow `None` or an unexpected string, the code silently falls through to the mode. Same pattern for lnmixture. A `raise ValueError(f"Unexpected PE choice: ...")` would surface the bug immediately instead of producing a plausible-but-wrong number.

### 4. `convert_for_download` returns `None` for unknown format

```python
def convert_for_download(df, file_format="csv"):
    if file_format == "csv":
        ...
    elif file_format == "parquet":
        ...
    else:
        return None
```

This `None` flows into `st.download_button(data=None)`, giving the user an empty/corrupt download with no error. Should `raise ValueError`.

### 5. Session state initialization is inconsistent

Some keys are pre-initialized at the top (`iwp_pe`, `sims_run`, all the PrEP keys), but others are not:

| Key | Pre-initialized? | Set where? |
|---|---|---|
| `iwp_pe` | ✅ Yes | Top of file |
| `iwp_cri` | ❌ No | Run handler only |
| `iwp_range` | ❌ No | Run handler only |
| `bs` | ❌ No | Run handler only |
| `sim_df` | ❌ No | Run handler only |
| `sim_df_prep_oral` | ✅ Yes | Top of file (in loop) |
| `used_go` | ❌ No | Run handler only; read with `.get("used_go", True)` |

All uninitialized keys are guarded by `if st.session_state["sims_run"]:`, so they don't crash today. But the inconsistency makes it easy to add a new read path that forgets the guard. The `.get("used_go", True)` default is a silent assumption that could mislabel a caption after a lookback run + residual risk calculation (though it's guarded by `len(total_components) > 1` in practice).

### 6. Lookback IDI parsing: silent precedence, no validation

```python
if uploaded_idi is not None:
    idis = ...  # CSV wins silently
elif idi_text.strip():
    idis = ...  # text ignored if CSV also present
```

If the user uploads a CSV AND types in the text area, the CSV wins with no indication. Also, no validation that parsed IDI values are positive or plausible — a negative IDI or a value of 50000 would flow straight into `iwp_from_lookback_data`.

### 7. `samp` vs `sim_df` fallback creates type ambiguity

```python
st.session_state["samp"] = pl.DataFrame({"iwp": st.session_state["bs"]})
if st.session_state["sim_df"] is None:
    st.session_state["sim_df"] = st.session_state["samp"]
```

After this, `sim_df` can be either a rich per-iteration frame (with `backend` column, all sampled parameters) or a one-column frame. Downstream code that expects the rich schema would silently get the simple one.

### 8. Stale comment at the top

```python
# Expects streamlit to be run from the root of the repository
# streamlit run app/app.py
```

The actual entry point is `streamlit run app.py` (confirmed by `app.py` itself and AGENTS.md). The `app/app.py` path was from a previous directory layout.

---

## Minor Nits

- **`import re`** appears at the top of the file AND inside the lookback UI block AND inside the lookback button handler. Redundant.
- **`import math as _math`** imported twice (lnmixture UI block + run handler). Harmless but noisy.
- **`time.sleep(0.3)`** ×3 (baseline + oral + inj progress bars) adds ~1s to every run for cosmetic purposes. Worth knowing it compounds.
- **`plot_histogram`** parameterises `x` but hardcodes `labels={"iwp": ...}` — the label mapping silently misses if called with a different column.
- **`_prerelease_notice` regexes** (`r"a\d"`, `r"b\d"`) are loose — they'd match any version containing those patterns. Fine for PEP 440 but slightly fragile.
- **No type hints** on most functions, despite `_prerelease_notice` having them and AGENTS.md recommending them. Inconsistent.
- **Dead code**: commented-out parameter defaults at the top (lines 42–57) and the logo footer at the bottom. The footer has an explanatory comment; the top block doesn't.
- **`volume_range_default`** is computed from `volume_pe` but only used as the slider's initial `value=`. Streamlit ignores `value=` on reruns, so changing `volume_pe` after first render doesn't update the range slider. Correct behavior, but a comment would prevent confusion.
- **`alpha=0` edge case** not explicitly handled. `sig_level` becomes 100, quantiles become `min()` and `max()`, reflecting the full range. The library would receive `alpha=0.00` — if its CrI computation has boundary-sensitive logic, that could be an issue.

---

## Summary

The design is **sound for a Streamlit app of this complexity**. The main structural risk is the conditional-variable-definition pattern — it works because the guards align today, but it's a maintenance hazard as the file grows. The most actionable items are:

1. Stop recomputing CrI/range from quantiles when the library already returned them.
2. Make the `k_pe` fallback chains raise on unexpected paths.
3. Make session state initialization consistent.
4. Fix the stale `app/app.py` comment.
5. Add validation and feedback to the lookback IDI parsing.
