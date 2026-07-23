# Residual HIV Transfusion Transmission Risk Estimator
# Copyright (C) 2025-2026 Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Expects streamlit to be run from the root of the repository
# streamlit run app.py with venv activated, or
# uv run streamlit run app.py
import io
import random
import re
import statistics
import time
from pathlib import Path

import numpy as np
import polars as pl
import plotly.express as px
import scipy.stats as stats
import streamlit as st

import residualrisk as rr
import residualrisk.prep as rrprep

# Static assets (pre-computed posteriors + Vitalant Research Institute branding)
_STATIC_DIR = Path(__file__).parent / "static"

# Set default values
# this keeps resetting to this value, so I am going to get rid of it
# seed = 126887

n_cpu = rr.get_cpu_core_count()

# alpha = 0.05

C0 = 0.00025
# doubling_time = 0.8542 # 20.5/24 Fiebig et al. (2003)
# doubling_time_sd = (23.4/24 - 18.2/24) / (2 * 1.96)
# IUs_per_copy = 1.72 # The WHO uses this conversion factor
# lod50 = 4.7 / IUs_per_copy
# lod50_sd = (5.3 - 4.0) / (2 * 1.96) / IUs_per_copy
# lod95_lod50_ratio = 21.2 / 4.7
# pool_size = 16
# retests = 1

# Canned NAT-assay LoD presets now live in the installable package
# (residualrisk/assays.py) as the single source of truth; the API and this UI
# consume the same table. NAT_ASSAYS is keyed by slug, each entry carrying a
# display_name plus lod50 / lod50_sd / lod95 / cp_per_iu / iu_std. See the
# residualrisk.assays module docstring for the derivation and conversion notes.
NAT_ASSAYS = rr.NAT_ASSAYS
MANUAL_LOD_OPTION = "Enter limits of detection"

if "seed" not in st.session_state:
    st.session_state["seed"] = random.randint(1, 999999)

if "sims_run" not in st.session_state:
    st.session_state["sims_run"] = False

if "rde_method_run" not in st.session_state:
    st.session_state["rde_method_run"] = None

if "iwp_pe_lookback" not in st.session_state:
    st.session_state["iwp_pe_lookback"] = None

if "iwp_ci_lookback" not in st.session_state:
    st.session_state["iwp_ci_lookback"] = None

if "iwp_pe" not in st.session_state:
    st.session_state["iwp_pe"] = None

if "iwp_pe_last" not in st.session_state:
    st.session_state["iwp_pe_last"] = None

# PrEP session state
for _k in (
    "iwp_pe_prep_oral",
    "samp_prep_oral",
    "sim_df_prep_oral",
    "iwp_pe_prep_inj",
    "samp_prep_inj",
    "sim_df_prep_inj",
):
    if _k not in st.session_state:
        st.session_state[_k] = None

for _k in ("prep_oral_run", "prep_inj_run"):
    if _k not in st.session_state:
        st.session_state[_k] = False


@st.cache_data
def load_data():
    # Use Path to ensure files are loaded relative to this script, not cwd
    static_dir = _STATIC_DIR

    k_animal = (
        pl
        .read_parquet(static_dir / "k_param_animal.parquet", columns=["k"])
        .get_column("k")
        .to_numpy()
    )
    k_human = (
        pl
        .read_parquet(static_dir / "k_param_human.parquet", columns=["k"])
        .get_column("k")
        .to_numpy()
    )
    k_expdecay = (
        pl
        .read_parquet(static_dir / "k_param_expdecay.parquet", columns=["k"])
        .get_column("k")
        .to_numpy()
    )

    # KDE modes via Go binary (~1.5s total, 30× faster than Python KDE).
    # Falls back to hardcoded values if Go binary is unavailable.
    _go_bin = rr.find_go_binary()
    if _go_bin is not None:
        k_human_mode = rr.mode_kde_go(k_human, cap=None, n_grid=100_000)
        k_animal_mode = rr.mode_kde_go(k_animal, cap=None, n_grid=100_000)
        k_expdecay_mode = rr.mode_kde_go(k_expdecay, cap=None, n_grid=100_000)
    else:
        # Hardcoded fallback (computed with Python KDE on full posteriors).
        # TODO: remove once Go binary is always available in deployment.
        # k_human_mode = rr.mode_kde(k_human)   # ~14s — too slow for startup
        # k_animal_mode = rr.mode_kde(k_animal)  # ~13s
        # k_expdecay_mode = rr.mode_kde(k_expdecay)  # ~18s
        k_human_mode = 0.0006716945195203189
        k_animal_mode = 0.020860882277014912
        k_expdecay_mode = 0.0005854393416702409

    return k_animal, k_human, k_expdecay, k_human_mode, k_animal_mode, k_expdecay_mode


@st.cache_data
def convert_for_download(df, file_format="csv"):
    if file_format == "csv":
        return df.write_csv().encode("utf-8")
    elif file_format == "parquet":
        buf = io.BytesIO()
        df.write_parquet(buf)
        return buf.getvalue()
    else:
        return None


def plot_histogram(
    data, x="iwp", colour=None, histnorm=None, title="Distribution of simulated RDEs"
):
    fig = px.histogram(
        data,
        x=x,
        color=colour,
        histnorm=histnorm,
        barmode="overlay",
        labels={"iwp": "risk day equivalents (RDEs)"},
        title=title,
    )
    return fig


def _prerelease_notice(version: str) -> str | None:
    """Pre-release warning for the app banner, keyed off the PEP 440 version
    string; ``None`` for a stable (``X.Y.Z``) release. ``.dev`` is checked first
    so a dev build of a beta (e.g. ``1.1.0b3.dev1``) reads "unstable test build",
    not "beta"."""
    if ".dev" in version:
        kind = "an unstable test build"
    elif re.search(r"a\d", version):
        kind = "an alpha release"
    elif re.search(r"b\d", version):
        kind = "a beta release"
    elif re.search(r"rc\d", version):
        kind = "a release candidate"
    else:
        return None
    return f"This is {kind} ({version}). Use with caution."


header_container = st.container()

header_container.write("""
# Residual HIV Transfusion Transmission Risk Estimator
Tool for estimating the residual risk of HIV transfusion transmission with NAT screening.
""")

# Pre-release banner (unstable / alpha / beta / rc) — shown for non-stable builds only.
_prerelease_msg = _prerelease_notice(rr.__version__)
if _prerelease_msg:
    header_container.warning(_prerelease_msg)

if "k_human" not in st.session_state:
    (
        st.session_state["k_animal"],
        st.session_state["k_human"],
        st.session_state["k_expdecay"],
        st.session_state["k_human_mode"],
        st.session_state["k_animal_mode"],
        st.session_state["k_expdecay_mode"],
    ) = load_data()

rde_method = st.selectbox(
    "RDE estimation method",
    options=[
        "Lookback data",
        "Mechanistic model",
    ],
    index=1,
    help="Risk day equivalents (RDEs) are equivalent to the infectious window "
    "period (IWP). Lookback data: estimates the IWP directly from "
    "lookback investigation data. Mechanistic model: simulates the "
    "IWP from viral dynamics and assay parameters. ",
)

is_mechanistic_ui = rde_method == "Mechanistic model"

include_prep_oral = False
include_prep_inj = False
if is_mechanistic_ui:
    include_prep_oral = st.checkbox(
        "Include oral PrEP breakthrough risk",
        value=False,
        help="Add a separate oral PrEP (oPrEP) breakthrough infection RDE "
        "estimate on top of the baseline window-period risk.",
    )
    include_prep_inj = st.checkbox(
        "Include injectable PrEP breakthrough risk",
        value=False,
        help="Add a separate injectable PrEP (iPrEP) breakthrough infection "
        "RDE estimate on top of the baseline window-period risk.",
    )

st.sidebar.write("Number of CPU cores: ", n_cpu)


def _generate_random_seed() -> None:
    # Set the new seed in the button's on_click callback — callbacks run BEFORE the
    # widgets are re-rendered, so the number_input below shows the new seed on the
    # same click. Assigning it in an `if button:` block *after* the input renders
    # only took effect on the next rerun (the old "needs two clicks" behaviour).
    st.session_state["seed"] = random.randint(1, 999999)


# The number_input is bound to session_state["seed"] via key= (seeded above), so
# the callback and manual entry share one source of truth.
st.sidebar.number_input(
    "Specify a seed value:",
    min_value=1,
    max_value=999999,
    step=1,
    help="Random seed for the Monte Carlo draws. Fix it to reproduce a run exactly, or click Generate random seed.",
    key="seed",
)
st.sidebar.button("Generate random seed", on_click=_generate_random_seed)

if is_mechanistic_ui:
    implementation = st.sidebar.selectbox(
        "Simulation implementation",
        options=["Go", "Python"],
        index=0,  # Go is default
        help="Compute backend. Go is the fast multi-core default; Python is a slower single-core reference/fallback.",
    )
    use_go_acceleration = implementation == "Go"
    if use_go_acceleration:
        if rr.find_go_binary() is None:
            st.sidebar.warning(
                "Go binary not found. Simulations fall back to the Python "
                "implementation, which is significantly slower — and mode point "
                "estimates use a coarser KDE grid (5 000 vs the Go path's 100 000), "
                "which differs slightly (~0.1%)."
            )

sim_param_container = st.expander(
    "Simulation settings", expanded=True, icon=":material/menu_open:"
)
if is_mechanistic_ui:
    trans_param_container = st.expander(
        "Transmissibility parameters", expanded=True, icon=":material/menu_open:"
    )
    model_param_container = st.expander(
        "RDE model parameters", expanded=True, icon=":material/menu_open:"
    )
if rde_method == "Lookback data":
    lookback_param_container = st.expander(
        "Lookback data parameters", expanded=True, icon=":material/menu_open:"
    )
if include_prep_oral or include_prep_inj:
    prep_shared_container = st.expander(
        "PrEP shared parameters", expanded=True, icon=":material/menu_open:"
    )
if include_prep_oral:
    prep_oral_container = st.expander(
        "Oral PrEP (oPrEP) parameters", expanded=True, icon=":material/menu_open:"
    )
if include_prep_inj:
    prep_inj_container = st.expander(
        "Injectable PrEP (iPrEP) parameters", expanded=True, icon=":material/menu_open:"
    )
incidence_param_container = st.expander(
    "Incidence parameters", expanded=True, icon=":material/menu_open:"
)

output_container = st.container()

with sim_param_container:
    col1, col2 = st.columns(2)

    if is_mechanistic_ui:
        point_estimate = col1.selectbox(
            "Select method for point estimate of RDEs",
            options=["primary parameters", "median", "mode", "mean"],
            index=2,
            help="How the reported RDE point estimate is computed: the mode (default), "
            "median or mean of the bootstrap distribution, or 'primary parameters' — a "
            "single evaluation at the input point values. The bootstrap summaries always "
            "lie within the credible interval; 'primary parameters' need not.",
        )
        if point_estimate == "primary parameters":
            col1.warning(
                "'Primary parameters' evaluates the RDE once at the input point values. "
                "Because the RDE distribution is right-skewed, this estimate can fall in "
                "the far right tail — occasionally above the upper credible limit. The "
                "mode or median summarise the distribution more representatively.",
                icon="⚠️",
            )
    else:
        point_estimate = None

    plot_norm = col1.selectbox(
        "Normalise histogram",
        options=[None, "probability density"],
        index=1,
        help="Y-axis of the RDE histogram: raw counts, or a probability density that integrates to 1.",
    )

    alpha = col1.number_input(
        "Significance level (𝛼)",
        min_value=0.00,
        max_value=0.20,
        value=0.05,
        step=0.01,
        help="Significance level for the reported credible interval; it covers the central (1 - α) of the bootstrap RDE distribution (0.05 → 95%).",
    )

    sig_level = round((1 - alpha) * 100)

    n_sims = col2.select_slider(
        "Select number of simulations",
        options=[
            1_000,
            5_000,
            10_000,
            25_000,
            50_000,
            100_000,
            250_000,
            500_000,
            1_000_000,
        ],
        value=100_000,
        help="Number of simulations to obtain uncertainty distribution of RDEs",
    )

    if is_mechanistic_ui:
        n_threads = col2.slider(
            "Select number of CPU cores to use",
            min_value=1,
            max_value=n_cpu,
            value=max(1, n_cpu),
            step=1,
            help="Number of CPU cores to use for parallel simulation runs",
        )

if is_mechanistic_ui:
    with trans_param_container:
        col1, col2 = st.columns(2)

        k_param_distribution_choice = col1.selectbox(
            "Select transmissibility parameter distribution to sample from",
            options=[
                "Belov human posterior",
                "Belov animal posterior",
                "Human-weighted exponential decay distribution",
                "Inverse Gamma distribution",
                "Lognormal mixture distribution",
            ],
            index=3,
            help="Distribution for the infectivity parameter k, sampled each bootstrap iteration. See the Documentation page for the options.",
        )
        match k_param_distribution_choice:
            case "Belov human posterior":
                k_param_dist = "human"
            case "Belov animal posterior":
                k_param_dist = "animal"
            case "Human-weighted exponential decay distribution":
                k_param_dist = "human_weighted"
            case "Inverse Gamma distribution":
                k_param_dist = "invgamma"
            case "Lognormal mixture distribution":
                k_param_dist = "lnmixture"
            case _:
                k_param_dist = None  # This shouldn't happen

        # PE selectbox for non-InvGamma, non-lnmixture paths.
        # For InvGamma and lnmixture the PE selectbox is deferred until after
        # parameters are defined so its options can depend on them.
        if k_param_dist not in ("invgamma", "lnmixture"):
            k_invgamma_pe_choice = None
            k_param_pe = col2.selectbox(
                "Transmissibility parameter point estimate: posterior...",
                options=["mode", "median", "mean"],
                index=0,
                help=(
                    "Which summary statistic of the posterior distribution to use "
                    "as the k point estimate when computing the IWP point estimate. "
                    "Does not affect bootstrap sampling."
                ),
            )
        else:
            k_param_pe = None

        if k_param_dist == "human":
            k_param = st.session_state["k_human"]
        elif k_param_dist == "animal":
            k_param = st.session_state["k_animal"]
        elif k_param_dist == "human_weighted":
            k_param = st.session_state["k_expdecay"]
        else:
            k_param = None

        # InvGamma parameter inputs — shown only when InvGamma is selected
        if k_param_dist == "invgamma":
            st.divider()
            ig_col1, ig_col2 = st.columns([1, 2])

            k_invgamma_alpha = ig_col1.number_input(
                "α (shape)",
                min_value=0.01,
                max_value=20.0,
                value=2.0,
                step=0.05,
                format="%.2f",
                help=(
                    "Shape parameter of the Inverse Gamma distribution. "
                    "Decrease α for a heavier right tail (more weight on large k values); "
                    "increase α to concentrate the distribution more tightly around the mode. "
                    "Recommended value: 2 (power-law tail, infinite variance by design). "
                    "α > 1 is required for a finite mean; α > 2 for finite variance. "
                    "The 'mean' point estimate option is disabled when α ≤ 1."
                ),
            )

            # PE selectbox placed in col2 of the top row; rendered there even though
            # defined here — Streamlit column objects accept widgets at any point.
            # "mean" is excluded when α ≤ 1 because the mean is infinite.
            _ig_pe_options = (
                ["mode", "median", "mean"]
                if k_invgamma_alpha > 1.0
                else ["mode", "median"]
            )
            k_invgamma_pe_choice = col2.selectbox(
                "Transmissibility parameter point estimate: distribution...",
                options=_ig_pe_options,
                index=0,
                help=(
                    "Which summary statistic of the Inverse Gamma distribution to use "
                    "as the k point estimate when computing the IWP point estimate. "
                    "Does not affect bootstrap sampling. "
                    "'Mean' is only available when α > 1."
                ),
            )

            ig_param_by = ig_col2.radio(
                "Parameterise by",
                options=["Mode (recommended)", "β (scale)"],
                index=0,
                horizontal=True,
            )

            if ig_param_by == "Mode (recommended)":
                k_human_mode_val = st.session_state["k_human_mode"]
                mode_col, custom_col = st.columns([2, 1])
                ig_mode_source = mode_col.radio(
                    "Mode value",
                    options=[
                        f"Human posterior ({k_human_mode_val:.6f})",
                        "Custom",
                    ],
                    index=0,
                )
                if "Human" in ig_mode_source:
                    k_invgamma_mode = k_human_mode_val
                else:
                    k_invgamma_mode = custom_col.number_input(
                        "Custom mode",
                        min_value=1e-7,
                        max_value=1.0,
                        value=float(k_human_mode_val),
                        format="%.6f",
                        step=0.000001,
                        help="Mode of the Inverse Gamma distribution.",
                    )
                k_invgamma_beta = k_invgamma_mode * (k_invgamma_alpha + 1)
                st.caption(
                    f"β = mode × (α + 1) = {k_invgamma_mode:.6f} × "
                    f"{k_invgamma_alpha + 1:.2f} = {k_invgamma_beta:.6f}"
                )

            else:  # "β (scale)"
                ig_beta_col, _ = st.columns(2)
                k_invgamma_beta = ig_beta_col.number_input(
                    "β (scale)",
                    min_value=1e-8,
                    max_value=1.0,
                    value=0.002019,
                    format="%.6f",
                    step=0.000001,
                    help="Scale parameter of the Inverse Gamma distribution.",
                )
                k_invgamma_mode = k_invgamma_beta / (k_invgamma_alpha + 1)
                st.caption(
                    f"mode = β / (α + 1) = {k_invgamma_beta:.6f} / "
                    f"{k_invgamma_alpha + 1:.2f} = {k_invgamma_mode:.6f}"
                )

        else:
            k_invgamma_alpha = None
            k_invgamma_beta = None

        # Lognormal mixture parameter inputs — shown only when lnmixture is selected
        if k_param_dist == "lnmixture":
            st.divider()

            # Default mixture parameters (Recommendation B from K_PARAM_INPUTDIST.md)
            _LN_W_DEF = 0.90
            _LN_MU1_DEF = -7.2403
            _LN_S1_DEF = 0.3241
            _LN_MU2_DEF = -3.7423
            _LN_S2_DEF = 0.5258

            lnmix_col1, lnmix_col2 = st.columns([1, 2])

            k_lnmix_w = lnmix_col1.slider(
                "Mixing weight (human component)",
                min_value=0.0,
                max_value=1.0,
                value=_LN_W_DEF,
                step=0.01,
                format="%.2f",
                help=(
                    "Weight placed on the human posterior component (component 1). "
                    "Remainder (1 − w) goes to the animal posterior component (component 2). "
                    "Recommended default: 0.90 (90% human, 10% animal)."
                ),
            )

            # PE selectbox placed in col2 of the top row
            k_lnmix_pe_choice = col2.selectbox(
                "Transmissibility parameter point estimate: distribution...",
                options=["mode", "median", "mean"],
                index=0,
                help=(
                    "Which summary statistic of the lognormal mixture to use as the k "
                    "point estimate when computing the IWP point estimate. "
                    "Does not affect bootstrap sampling. "
                    "'Mean' is analytic; 'mode' and 'median' are computed numerically."
                ),
            )

            # Advanced: edit component parameters
            lnmix_advanced = lnmix_col2.checkbox(
                "Advanced: edit component parameters",
                value=False,
                help=(
                    "Edit the log-scale mean (μ) and log-scale standard deviation (σ) "
                    "of each mixture component. Defaults are the MLE fits to the human "
                    "and animal k posteriors from the companion analysis."
                ),
            )

            if lnmix_advanced:
                adv_col1, adv_col2 = st.columns(2)
                k_lnmix_mu1 = adv_col1.number_input(
                    "μ₁ (human, log-scale mean)",
                    value=_LN_MU1_DEF,
                    format="%.4f",
                    step=0.01,
                    help="Log-scale mean for component 1 (human). Default: −7.2403.",
                )
                k_lnmix_sigma1 = adv_col1.number_input(
                    "σ₁ (human, log-scale SD)",
                    min_value=0.001,
                    value=_LN_S1_DEF,
                    format="%.4f",
                    step=0.01,
                    help="Log-scale SD for component 1 (human). Default: 0.3241.",
                )
                k_lnmix_mu2 = adv_col2.number_input(
                    "μ₂ (animal, log-scale mean)",
                    value=_LN_MU2_DEF,
                    format="%.4f",
                    step=0.01,
                    help="Log-scale mean for component 2 (animal). Default: −3.7423.",
                )
                k_lnmix_sigma2 = adv_col2.number_input(
                    "σ₂ (animal, log-scale SD)",
                    min_value=0.001,
                    value=_LN_S2_DEF,
                    format="%.4f",
                    step=0.01,
                    help="Log-scale SD for component 2 (animal). Default: 0.5258.",
                )
            else:
                k_lnmix_mu1 = _LN_MU1_DEF
                k_lnmix_sigma1 = _LN_S1_DEF
                k_lnmix_mu2 = _LN_MU2_DEF
                k_lnmix_sigma2 = _LN_S2_DEF

            # Derived statistics display
            import math as _math

            _lnmix_comp1_median = _math.exp(k_lnmix_mu1)
            _lnmix_comp2_median = _math.exp(k_lnmix_mu2)
            _lnmix_mean = k_lnmix_w * _math.exp(k_lnmix_mu1 + k_lnmix_sigma1**2 / 2) + (
                1 - k_lnmix_w
            ) * _math.exp(k_lnmix_mu2 + k_lnmix_sigma2**2 / 2)
            st.caption(
                f"Component 1 median: {_lnmix_comp1_median:.6f} &nbsp;|&nbsp; "
                f"Component 2 median: {_lnmix_comp2_median:.6f} &nbsp;|&nbsp; "
                f"Mixture mean: {_lnmix_mean:.6f}"
            )

        else:
            k_lnmix_w = None
            k_lnmix_mu1 = None
            k_lnmix_sigma1 = None
            k_lnmix_mu2 = None
            k_lnmix_sigma2 = None
            k_lnmix_pe_choice = None

    with model_param_container:
        col1, col2 = st.columns(2)

        volume_pe = col1.number_input(
            "Average volume transfused (mL)",
            min_value=1,
            max_value=500,
            value=20,
            step=1,
            help="Point-estimate residual plasma volume of the transfused product (~20 mL red cells, ~200 mL plasma). Model one product per run.",
        )
        volume_range_default = (round(0.75 * volume_pe), round(1.5 * volume_pe))
        volume_range = col1.slider(
            "Range of volumes transfused (mL)",
            min_value=1,
            max_value=500,
            value=volume_range_default,
            help="Plausible range for the transfused plasma volume; sampled uniformly each bootstrap iteration.",
        )

        doubling_time_hours = col1.number_input(
            "Viral outgrowth doubling time (hours)",
            min_value=12.0,
            max_value=48.0,
            value=20.5,
            step=0.25,
            help="Doubling time of viral concentration during early ramp-up (default 20.5 h; Fiebig et al. 2003).",
        )
        doubling_time = doubling_time_hours / 24
        doubling_time_hours_sd = col1.number_input(
            "Viral outgrowth doubling time SD",
            min_value=0.0,
            max_value=10.0,
            value=1.33,
            step=0.01,
            help="Standard deviation (hours) of the doubling time; drawn from a truncated normal each bootstrap iteration.",
        )
        doubling_time_sd = doubling_time_hours_sd / 24

        id_nat = col1.checkbox(
            "Individual donation NAT screening",
            value=False,
            help="Screen each donation individually (pool size 1) instead of in a minipool. Disables the minipool inputs below.",
        )

        if not id_nat:
            pool_size = col2.number_input(
                "Minipool size",
                min_value=2,
                max_value=96,
                value=16,
                step=1,
                help="Donations combined per minipool NAT test; larger pools dilute a positive donation, raising the effective limit of detection.",
            )
            retests = col2.number_input(
                "Number of retests (pool resolution)",
                min_value=0,
                max_value=5,
                value=1,
                step=1,
                help="Number of individual-donation retests performed when resolving a reactive minipool.",
            )
        else:
            pool_size = 1
            retests = 0

        nat_assay_options = list(NAT_ASSAYS) + [MANUAL_LOD_OPTION]
        nat_assay = col2.selectbox(
            "Select NAT assay",
            options=nat_assay_options,
            index=nat_assay_options.index("ultrio_elite"),
            format_func=lambda key: (
                NAT_ASSAYS[key]["display_name"] if key in NAT_ASSAYS else key
            ),
            help="Select a NAT assay to use its published 50%/95% limits of "
            "detection (copies/mL, HIV-1 Group M), or choose "
            f"'{MANUAL_LOD_OPTION}' to enter values manually.",
        )

        if nat_assay == MANUAL_LOD_OPTION:
            # Manual entry, pre-populated with the cobas MPX defaults.
            # lod50 must be strictly positive: lod95_lod50_ratio divides by it.
            lod50 = col2.number_input(
                "NAT assay 50% LoD (copies/mL)",
                min_value=0.01,
                max_value=500.0,
                value=NAT_ASSAYS["cobas_mpx"]["lod50"],
                step=0.01,
                help="Viral concentration detected 50% of the time (copies/mL, HIV-1 Group M).",
            )
            lod50_sd = col2.number_input(
                "NAT assay 50% LoD SD (copies/mL)",
                min_value=0.0,
                max_value=500.0,
                value=NAT_ASSAYS["cobas_mpx"]["lod50_sd"],
                step=0.001,
                format="%.4f",
                help="Standard deviation of the 50% LoD (copies/mL); drawn from a truncated normal each bootstrap iteration.",
            )
            lod95 = col2.number_input(
                "NAT assay 95% LoD (copies/mL)",
                min_value=0.0,
                max_value=500.0,
                value=NAT_ASSAYS["cobas_mpx"]["lod95"],
                step=0.01,
                help="Viral concentration detected 95% of the time (copies/mL); with the 50% LoD it sets the slope of the detection curve.",
            )
            conversion_note = (
                "Limits of detection entered directly in copies/mL "
                "(no IU/mL → copies/mL conversion applied)."
            )
        else:
            # Canned assay: populate from the lookup table and show (read-only).
            lod50 = NAT_ASSAYS[nat_assay]["lod50"]
            lod50_sd = NAT_ASSAYS[nat_assay]["lod50_sd"]
            lod95 = NAT_ASSAYS[nat_assay]["lod95"]
            cp_per_iu = NAT_ASSAYS[nat_assay]["cp_per_iu"]
            iu_std = NAT_ASSAYS[nat_assay]["iu_std"]
            col2.number_input(
                "NAT assay 50% LoD (copies/mL)",
                value=lod50,
                step=0.01,
                disabled=True,
            )
            col2.number_input(
                "NAT assay 50% LoD SD (copies/mL)",
                value=lod50_sd,
                step=0.001,
                format="%.4f",
                disabled=True,
            )
            col2.number_input(
                "NAT assay 95% LoD (copies/mL)",
                value=lod95,
                step=0.01,
                disabled=True,
            )
            conversion_note = (
                f"IU/mL → copies/mL conversion factor in use: {cp_per_iu} cp/IU "
                f"(HIV-1 Group M, {iu_std})."
            )

        st.caption(conversion_note)
        st.caption("95% LoD : 50% LoD ratio will be fixed for simulations.")
        # fix_lod95_lod50_ratio = col2.checkbox(
        #     "Fix 95% LoD:50% LoD ratio",
        #     value = True
        # )
        lod95_lod50_ratio = lod95 / lod50

if rde_method == "Lookback data":
    with lookback_param_container:
        import re

        col1, col2 = st.columns(2)

        n_transmissions_lb = col1.number_input(
            "Number of confirmed transfusion transmissions",
            min_value=0,
            max_value=10000,
            value=0,
            step=1,
            help="Confirmed HIV transmissions from prior donations identified "
            "through lookback investigation.",
        )
        neg_diag_delay = col1.number_input(
            "Negative test diagnostic delay (days)",
            min_value=0.0,
            max_value=60.0,
            value=5.0,
            step=0.5,
            help="Diagnostic delay of the most sensitive test applied at the "
            "prior (negative) donation.",
        )
        pos_diag_delay = col1.number_input(
            "Positive test diagnostic delay (days)",
            min_value=0.0,
            max_value=60.0,
            value=10.0,
            step=0.5,
            help="Diagnostic delay of the least sensitive positive test at the "
            "seroconversion donation.",
        )

        col2.write("**Inter-donation intervals (IDIs)**")
        uploaded_idi = col2.file_uploader(
            "Upload CSV (single column, no header)",
            type=["csv"],
            help="One IDI value per row in days.",
        )
        idi_text = col2.text_area(
            "Or enter IDI values (days), one per line or comma-separated",
            value="",
            height=150,
            placeholder="105\n98\n120\n...",
        )

if include_prep_oral or include_prep_inj:
    with prep_shared_container:
        col1, col2 = st.columns(2)

        eclipse = col1.number_input(
            "Eclipse period (days)",
            min_value=1,
            max_value=30,
            value=7,
            step=1,
            help="Duration of the eclipse phase before viral RNA becomes detectable.",
        )
        eclipse_range = col1.slider(
            "Eclipse period range (days)",
            min_value=1,
            max_value=20,
            value=(4, 10),
            step=1,
            help="Uncertainty range for the eclipse period (sampled uniformly).",
        )

        col2.write("**Sinusoidal set-point oscillation parameters**")
        # The plateau oscillates as set_point * (offset + a*sin(b*t)), with the offset
        # fixed at the library default of 1.0 — so the plateau centres on the set
        # point and the amplitude a is a plain fraction of it. The offset is
        # deliberately NOT exposed: it is exactly redundant with the set point
        # (offset=o with amplitude a reproduces set_point*o with amplitude a/o) while
        # being a first-order lever on the answer, so it was a trap — the set point is
        # the parameter to vary, in interpretable clinical units and with its own
        # bootstrap range. Library callers can still pass offset= to risk_days_prep_bs.
        prep_a = col2.number_input(
            "Amplitude (a)",
            min_value=0.0,
            max_value=1.0,
            value=0.7,
            step=0.05,
            format="%.2f",
            help="Amplitude of the sinusoidal oscillation around the viral set point, "
            "as a fraction of it: the plateau swings between (1 - a) and (1 + a) times "
            "the set point. Capped at 1 — a larger amplitude would drive the plateau "
            "viral load negative.",
        )
        prep_b = col2.number_input(
            "Frequency (b)",
            min_value=0.0,
            max_value=5.0,
            value=0.6,
            step=0.05,
            format="%.2f",
            help="Frequency of sinusoidal oscillation.",
        )

        prep_vary_sin = col2.checkbox(
            "Vary sinusoidal oscillation parameters (a, b)",
            value=False,
            help="When enabled, the amplitude (a) and frequency (b) are sampled "
            "uniformly across bootstrap iterations instead of held fixed at the "
            "values above.",
        )
        if prep_vary_sin:
            prep_a_dist_uniform = col2.slider(
                "Amplitude (a) range",
                min_value=0.0,
                max_value=1.0,
                value=(0.5, 0.9),
                step=0.05,
                format="%.2f",
                help="Uniform sampling range for the amplitude a. Capped at 1 — a "
                "larger amplitude would drive the plateau viral load negative.",
            )
            prep_b_dist_uniform = col2.slider(
                "Frequency (b) range",
                min_value=0.0,
                max_value=5.0,
                value=(0.4, 0.8),
                step=0.05,
                format="%.2f",
                help="Uniform sampling range for the frequency b.",
            )
        else:
            prep_a_dist_uniform = None
            prep_b_dist_uniform = None

if include_prep_oral:
    with prep_oral_container:
        col1, col2 = st.columns(2)

        vl_setpoint_oral = col1.number_input(
            "Viral load set point (c/mL)",
            min_value=1,
            max_value=5000,
            value=340,
            step=10,
            help="Median oPrEP breakthrough viral load set point.",
        )
        vl_setpoint_range_oral = col1.slider(
            "Viral load set point range (c/mL)",
            min_value=1,
            max_value=5000,
            value=(10, 2270),
            step=10,
            help="Uncertainty range for the oPrEP viral load set point (sampled uniformly).",
        )
        drug_effect_oral = col1.number_input(
            "Drug effect (transmissibility factor)",
            min_value=0.05,
            max_value=1.0,
            value=1.0,
            step=0.05,
            format="%.2f",
            help=(
                "Multiplicative transmissibility-reduction factor from antiretroviral "
                "drug presence in oPrEP breakthrough infection (1.0 = no reduction). "
                "Used for the primary-parameters point estimate."
            ),
        )
        drug_effect_range_oral = col1.slider(
            "Drug effect range (bootstrap)",
            min_value=0.05,
            max_value=1.0,
            value=(0.75, 1.0),
            step=0.05,
            help=(
                "Uncertainty range for the oPrEP drug effect, sampled uniformly per "
                "bootstrap iteration. Set to (1.0, 1.0) for no reduction; e.g. "
                "(0.5, 1.0) reproduces the prior analysis."
            ),
        )
        # (1.0, 1.0) → no range; hold fixed at the point value (pass None).
        prep_drug_effect_dist_oral = (
            None if drug_effect_range_oral[0] >= 1.0 else drug_effect_range_oral
        )

        seroconversion_min_oral = col2.number_input(
            "Time to seroconversion min (days)",
            min_value=0,
            max_value=500,
            value=29,
            step=1,
            help="Minimum time from infection to seroconversion in oPrEP users.",
        )
        seroconversion_max_oral = col2.number_input(
            "Time to seroconversion max (days)",
            min_value=0,
            max_value=500,
            value=250,
            step=1,
            help="Maximum time from infection to seroconversion in oPrEP users.",
        )
        seroconversion_weibull_alpha_oral = col2.number_input(
            "Seroconversion Weibull scale (α)",
            min_value=0.0,
            max_value=500.0,
            value=50.49434,
            step=0.001,
            format="%.5f",
            help="Scale parameter of the Weibull-family seroconversion curve for oPrEP.",
        )
        seroconversion_weibull_beta_oral = col2.number_input(
            "Seroconversion Weibull shape (β)",
            min_value=0.0,
            max_value=500.0,
            value=1.15062,
            step=0.001,
            format="%.5f",
            help="Shape parameter of the Weibull-family seroconversion curve for oPrEP.",
        )

if include_prep_inj:
    with prep_inj_container:
        col1, col2 = st.columns(2)

        vl_setpoint_inj = col1.number_input(
            "Viral load set point (c/mL)",
            min_value=1,
            max_value=5000,
            value=30,
            step=10,
            help="Median iPrEP breakthrough viral load set point.",
        )
        vl_setpoint_range_inj = col1.slider(
            "Viral load set point range (c/mL)",
            min_value=1,
            max_value=5000,
            value=(10, 2500),
            step=10,
            help="Uncertainty range for the iPrEP viral load set point (sampled uniformly).",
        )
        drug_effect_inj = col1.number_input(
            "Drug effect (transmissibility factor)",
            min_value=0.05,
            max_value=1.0,
            value=1.0,
            step=0.05,
            format="%.2f",
            help=(
                "Multiplicative transmissibility-reduction factor from antiretroviral "
                "drug presence in iPrEP breakthrough infection (1.0 = no reduction). "
                "Used for the primary-parameters point estimate. NB: for long-acting "
                "injectables the effect is expected to wane as the drug washes out."
            ),
        )
        drug_effect_range_inj = col1.slider(
            "Drug effect range (bootstrap)",
            min_value=0.05,
            max_value=1.0,
            value=(0.75, 1.0),
            step=0.05,
            help=(
                "Uncertainty range for the iPrEP drug effect, sampled uniformly per "
                "bootstrap iteration. Set to (1.0, 1.0) for no reduction; e.g. "
                "(0.5, 1.0) reproduces the prior analysis."
            ),
        )
        # (1.0, 1.0) → no range; hold fixed at the point value (pass None).
        prep_drug_effect_dist_inj = (
            None if drug_effect_range_inj[0] >= 1.0 else drug_effect_range_inj
        )

        seroconversion_min_inj = col2.number_input(
            "Time to seroconversion min (days)",
            min_value=0,
            max_value=500,
            value=42,
            step=1,
            help="Minimum time from infection to seroconversion in iPrEP users.",
        )
        seroconversion_max_inj = col2.number_input(
            "Time to seroconversion max (days)",
            min_value=0,
            max_value=500,
            value=250,
            step=1,
            help="Maximum time from infection to seroconversion in iPrEP users.",
        )
        seroconversion_weibull_alpha_inj = col2.number_input(
            "Seroconversion Weibull scale (α)",
            min_value=0.0,
            max_value=500.0,
            value=90.88988,
            step=0.001,
            format="%.5f",
            help="Scale parameter of the Weibull-family seroconversion curve for iPrEP.",
        )
        seroconversion_weibull_beta_inj = col2.number_input(
            "Seroconversion Weibull shape (β)",
            min_value=0.0,
            max_value=500.0,
            value=3.048339,
            step=0.001,
            format="%.5f",
            help="Shape parameter of the Weibull-family seroconversion curve for iPrEP.",
        )

inc_prep_oral_perpy = None
inc_prep_oral_perpy_sd = None
inc_prep_inj_perpy = None
inc_prep_inj_perpy_sd = None

with incidence_param_container:
    calculate_rr = st.checkbox(
        "Calculate residual risk (incidence x RDEs)",
        value=False,
        help="Combine the RDEs with a donor-population HIV incidence to report residual transfusion-transmission risk (incidence × RDE).",
    )
    if calculate_rr:
        st.write("**Baseline (non-PrEP) incidence**")
        inc_per100k = st.number_input(
            "Incidence (/100,000 PY)",
            min_value=0.001,
            max_value=10000.000,
            value=2.500,
            step=0.5,
            help="Baseline (non-PrEP) HIV incidence in the donor population, per 100,000 person-years.",
        )
        inc_perpy = inc_per100k / 100000
        inc_perpd = inc_per100k / 100000 / 365.25
        inc_per100k_sd = st.number_input(
            "Incidence (/100,000 PY) SD",
            min_value=0.001,
            max_value=10000.000,
            value=inc_per100k * 0.2,
            step=0.01,
            help="Standard deviation of the incidence estimate (per 100,000 PY); drawn each bootstrap iteration.",
        )
        inc_perpy_sd = inc_per100k_sd / 100000
        inc_perpd_sd = inc_per100k_sd / 100000 / 365.25
        st.text(
            f"Relative standard error on incidence: {inc_per100k_sd / inc_per100k * 100:.1f}%"
        )

        if include_prep_oral:
            st.write("**Oral PrEP breakthrough infection incidence**")
            inc_prep_oral_per100k = st.number_input(
                "oPrEP breakthrough incidence (/100,000 PY)",
                min_value=0.001,
                max_value=10000.000,
                value=0.200,
                step=0.1,
                help="Incidence of HIV breakthrough infection among oral PrEP users.",
                key="inc_prep_oral",
            )
            inc_prep_oral_per100k_sd = st.number_input(
                "oPrEP breakthrough incidence SD (/100,000 PY)",
                min_value=0.001,
                max_value=10000.000,
                value=inc_prep_oral_per100k * 0.3,
                step=0.01,
                help="Standard deviation of oPrEP breakthrough incidence.",
                key="inc_prep_oral_sd",
            )
            inc_prep_oral_perpy = inc_prep_oral_per100k / 100000
            inc_prep_oral_perpy_sd = inc_prep_oral_per100k_sd / 100000
            st.text(
                f"Relative standard error: {inc_prep_oral_per100k_sd / inc_prep_oral_per100k * 100:.1f}%"
            )

        if include_prep_inj:
            st.write("**Injectable PrEP breakthrough infection incidence**")
            inc_prep_inj_per100k = st.number_input(
                "iPrEP breakthrough incidence (/100,000 PY)",
                min_value=0.001,
                max_value=10000.000,
                value=0.050,
                step=0.01,
                help="Incidence of HIV breakthrough infection among injectable PrEP users.",
                key="inc_prep_inj",
            )
            inc_prep_inj_per100k_sd = st.number_input(
                "iPrEP breakthrough incidence SD (/100,000 PY)",
                min_value=0.001,
                max_value=10000.000,
                value=inc_prep_inj_per100k * 0.3,
                step=0.001,
                help="Standard deviation of iPrEP breakthrough incidence.",
                key="inc_prep_inj_sd",
            )
            inc_prep_inj_perpy = inc_prep_inj_per100k / 100000
            inc_prep_inj_perpy_sd = inc_prep_inj_per100k_sd / 100000
            st.text(
                f"Relative standard error: {inc_prep_inj_per100k_sd / inc_prep_inj_per100k * 100:.1f}%"
            )

button_label = "Run simulations" if is_mechanistic_ui else "Calculate RDEs"
if st.sidebar.button(button_label):
    if rde_method == "Mechanistic model":
        progressbar = st.sidebar.progress(0, text="Running simulations...")
        try:
            if k_param_dist == "invgamma":
                if k_invgamma_pe_choice == "mode":
                    k_pe = k_invgamma_beta / (k_invgamma_alpha + 1)
                elif k_invgamma_pe_choice == "median":
                    k_pe = stats.invgamma.ppf(
                        0.5, a=k_invgamma_alpha, scale=k_invgamma_beta
                    )
                elif k_invgamma_pe_choice == "mean":
                    k_pe = k_invgamma_beta / (k_invgamma_alpha - 1)
                else:
                    k_pe = k_invgamma_beta / (k_invgamma_alpha + 1)  # fallback to mode
            elif k_param_dist == "lnmixture":
                import math as _math

                if k_lnmix_pe_choice == "mean":
                    k_pe = k_lnmix_w * _math.exp(
                        k_lnmix_mu1 + k_lnmix_sigma1**2 / 2
                    ) + (1 - k_lnmix_w) * _math.exp(k_lnmix_mu2 + k_lnmix_sigma2**2 / 2)
                else:
                    # Numerical mode or median from a large sample
                    # Use cached default-param values if parameters are at defaults to avoid delay
                    _lnmix_defaults = (0.90, -7.2403, 0.3241, -3.7423, 0.5258)
                    _lnmix_current = (
                        k_lnmix_w,
                        k_lnmix_mu1,
                        k_lnmix_sigma1,
                        k_lnmix_mu2,
                        k_lnmix_sigma2,
                    )
                    if _lnmix_current == _lnmix_defaults:
                        _lnmix_sample = st.session_state.get("k_lnmix_default_sample")
                        if _lnmix_sample is None:
                            _lnmix_sample = rr.sample_lnmix(
                                100_000, *_lnmix_defaults, seed=42
                            )
                            st.session_state["k_lnmix_default_sample"] = _lnmix_sample
                    else:
                        _lnmix_sample = rr.sample_lnmix(
                            100_000,
                            k_lnmix_w,
                            k_lnmix_mu1,
                            k_lnmix_sigma1,
                            k_lnmix_mu2,
                            k_lnmix_sigma2,
                            seed=42,
                        )
                    if k_lnmix_pe_choice == "median":
                        k_pe = float(np.median(_lnmix_sample))
                    else:  # mode
                        # Go path uses the Go binary (FFT, 100_000); the Python
                        # backend uses a pure-Python coarse grid (5_000). We do NOT
                        # route to Go when the user explicitly chose Python.
                        if use_go_acceleration:
                            k_pe = rr.mode_kde_go(
                                _lnmix_sample, cap=None, n_grid=100_000
                            )
                            if k_pe is None:  # Go unavailable → pure-Python fallback
                                k_pe = rr.mode_kde(
                                    _lnmix_sample, cap=None, n_grid=5_000
                                )
                        else:
                            k_pe = rr.mode_kde(_lnmix_sample, cap=None, n_grid=5_000)
            elif k_param_pe == "mode":
                _mode_key = {
                    "human": "k_human_mode",
                    "animal": "k_animal_mode",
                    "human_weighted": "k_expdecay_mode",
                }.get(k_param_dist)
                k_pe = st.session_state[_mode_key] if _mode_key else None
            elif k_param_pe == "mean":
                k_pe = statistics.mean(k_param)
            elif k_param_pe == "median":
                k_pe = statistics.median(k_param)
            else:
                k_pe = None  # should not happen
            (
                st.session_state["iwp_pe"],
                st.session_state["iwp_cri"],
                st.session_state["iwp_range"],
                st.session_state["bs"],
                st.session_state["sim_df"],
            ) = rr.risk_days_bs(
                k_pe,
                doubling_time,
                doubling_time_sd,
                lod50,
                lod50_sd,
                lod95_lod50_ratio,
                volume_pe,
                volume_range,
                pool_size,
                retests,
                k_posterior_sample=k_param,
                k_invgamma_alpha=k_invgamma_alpha,
                k_invgamma_beta=k_invgamma_beta,
                k_gamma_scale=None,
                k_gamma_shape=None,
                k_lnmix_w=k_lnmix_w,
                k_lnmix_mu1=k_lnmix_mu1,
                k_lnmix_sigma1=k_lnmix_sigma1,
                k_lnmix_mu2=k_lnmix_mu2,
                k_lnmix_sigma2=k_lnmix_sigma2,
                alpha=alpha,
                n_bs=n_sims,
                point_estimate=point_estimate,
                seed=st.session_state["seed"],
                threads=n_threads,
                progress=progressbar,
                return_sim_df=True,
                use_go=use_go_acceleration,
            )
            st.session_state["sims_run"] = True
            st.session_state["rde_method_run"] = "Mechanistic model"
            # Record the backend *actually* used, not merely requested: a requested
            # Go run silently falls back to Python when the binary is unavailable, and
            # the total-risk CrI is a valid joint interval only when components share
            # their per-iteration k / viral-dynamics / LOD / volume draws — which only
            # the Go backend does (Python draws in a different order, not aligned).
            # Gate on binary availability so the CrI caption isn't mislabelled "exact"
            # after a fallback. (Each sim_df also carries a ground-truth `backend`
            # column from the library; a rare present-but-crashes fallback is logged.)
            st.session_state["used_go"] = (
                use_go_acceleration and rr.find_go_binary() is not None
            )
            st.session_state["samp"] = pl.DataFrame({"iwp": st.session_state["bs"]})
            # Fallback: if sim_df is None (e.g., from Go implementation), use samp
            if st.session_state["sim_df"] is None:
                st.session_state["sim_df"] = st.session_state["samp"]
            progressbar.progress(1.0, text="Simulations complete!")
            # Brief pause to show completion, then clear progress bar
            time.sleep(0.3)
            progressbar.empty()

            # --- PrEP bootstrap runs (Go or Python, per the selected implementation) ---
            _prep_k_kwargs = dict(
                k_posterior_sample=k_param,
                k_invgamma_alpha=k_invgamma_alpha,
                k_invgamma_beta=k_invgamma_beta,
                k_gamma_scale=None,
                k_gamma_shape=None,
                k_lnmix_w=k_lnmix_w,
                k_lnmix_mu1=k_lnmix_mu1,
                k_lnmix_sigma1=k_lnmix_sigma1,
                k_lnmix_mu2=k_lnmix_mu2,
                k_lnmix_sigma2=k_lnmix_sigma2,
            )

            if include_prep_oral:
                prep_oral_bar = st.sidebar.progress(
                    0,
                    text=f"Running oral PrEP simulations ({'Go' if use_go_acceleration else 'Python'})...",
                )
                (
                    st.session_state["iwp_pe_prep_oral"],
                    st.session_state["iwp_cri_prep_oral"],
                    st.session_state["iwp_range_prep_oral"],
                    st.session_state["bs_prep_oral"],
                    st.session_state["sim_df_prep_oral"],
                ) = rrprep.risk_days_prep_bs(
                    k_pe,
                    doubling_time,
                    doubling_time_sd,
                    lod50,
                    lod50_sd,
                    lod95_lod50_ratio,
                    volume_pe,
                    volume_range,
                    pool_size,
                    retests,
                    set_point=vl_setpoint_oral,
                    set_point_dist_uniform=vl_setpoint_range_oral,
                    eclipse=eclipse,
                    eclipse_dist_uniform=eclipse_range,
                    a=prep_a,
                    b=prep_b,
                    a_dist_uniform=prep_a_dist_uniform,
                    b_dist_uniform=prep_b_dist_uniform,
                    drug_effect=drug_effect_oral,
                    drug_effect_dist_uniform=prep_drug_effect_dist_oral,
                    ser_min=seroconversion_min_oral,
                    ser_max=seroconversion_max_oral,
                    ser_alpha=seroconversion_weibull_alpha_oral,
                    ser_beta=seroconversion_weibull_beta_oral,
                    **_prep_k_kwargs,
                    alpha=alpha,
                    n_bs=n_sims,
                    point_estimate=point_estimate,
                    seed=st.session_state["seed"],
                    threads=n_threads,
                    progress=prep_oral_bar,
                    return_sim_df=True,
                    use_go=use_go_acceleration,
                )
                st.session_state["prep_oral_run"] = True
                st.session_state["samp_prep_oral"] = pl.DataFrame({
                    "iwp": st.session_state["bs_prep_oral"]
                })
                if st.session_state["sim_df_prep_oral"] is None:
                    st.session_state["sim_df_prep_oral"] = st.session_state[
                        "samp_prep_oral"
                    ]
                prep_oral_bar.progress(1.0, text="Oral PrEP simulations complete!")
                time.sleep(0.3)
                prep_oral_bar.empty()
            else:
                st.session_state["prep_oral_run"] = False

            if include_prep_inj:
                prep_inj_bar = st.sidebar.progress(
                    0,
                    text=f"Running injectable PrEP simulations ({'Go' if use_go_acceleration else 'Python'})...",
                )
                (
                    st.session_state["iwp_pe_prep_inj"],
                    st.session_state["iwp_cri_prep_inj"],
                    st.session_state["iwp_range_prep_inj"],
                    st.session_state["bs_prep_inj"],
                    st.session_state["sim_df_prep_inj"],
                ) = rrprep.risk_days_prep_bs(
                    k_pe,
                    doubling_time,
                    doubling_time_sd,
                    lod50,
                    lod50_sd,
                    lod95_lod50_ratio,
                    volume_pe,
                    volume_range,
                    pool_size,
                    retests,
                    set_point=vl_setpoint_inj,
                    set_point_dist_uniform=vl_setpoint_range_inj,
                    eclipse=eclipse,
                    eclipse_dist_uniform=eclipse_range,
                    a=prep_a,
                    b=prep_b,
                    a_dist_uniform=prep_a_dist_uniform,
                    b_dist_uniform=prep_b_dist_uniform,
                    drug_effect=drug_effect_inj,
                    drug_effect_dist_uniform=prep_drug_effect_dist_inj,
                    ser_min=seroconversion_min_inj,
                    ser_max=seroconversion_max_inj,
                    ser_alpha=seroconversion_weibull_alpha_inj,
                    ser_beta=seroconversion_weibull_beta_inj,
                    **_prep_k_kwargs,
                    alpha=alpha,
                    n_bs=n_sims,
                    point_estimate=point_estimate,
                    seed=st.session_state["seed"],
                    threads=n_threads,
                    progress=prep_inj_bar,
                    return_sim_df=True,
                    use_go=use_go_acceleration,
                )
                st.session_state["prep_inj_run"] = True
                st.session_state["samp_prep_inj"] = pl.DataFrame({
                    "iwp": st.session_state["bs_prep_inj"]
                })
                if st.session_state["sim_df_prep_inj"] is None:
                    st.session_state["sim_df_prep_inj"] = st.session_state[
                        "samp_prep_inj"
                    ]
                prep_inj_bar.progress(1.0, text="Injectable PrEP simulations complete!")
                time.sleep(0.3)
                prep_inj_bar.empty()
            else:
                st.session_state["prep_inj_run"] = False

        except ValueError as e:
            progressbar.empty()
            st.sidebar.error(f"Invalid parameters: {e}")
    elif rde_method == "Lookback data":
        import re

        idis = None
        try:
            if uploaded_idi is not None:
                idi_df = pl.read_csv(uploaded_idi, has_header=False)
                idis = idi_df.to_series(0).to_list()
            elif idi_text.strip():
                parts = re.split(r"[,\n\r\s]+", idi_text.strip())
                idis = [float(v) for v in parts if v.strip()]
            else:
                st.sidebar.error("Please enter IDI values or upload a CSV file.")
        except ValueError as e:
            st.sidebar.error(f"Could not parse IDI values: {e}")

        if idis is not None:
            try:
                iwp_pe_lb, iwp_ci_lb, iwp_samples_lb = rr.iwp_from_lookback_data(
                    n_transmissions=n_transmissions_lb,
                    intervals=idis,
                    negative_diagnostic_delay=neg_diag_delay,
                    positive_diagnostic_delay=pos_diag_delay,
                    alpha=alpha,
                    n_bs=n_sims,
                    seed=st.session_state["seed"],
                )
                st.session_state["iwp_pe_lookback"] = iwp_pe_lb
                st.session_state["iwp_ci_lookback"] = iwp_ci_lb
                st.session_state["bs"] = list(iwp_samples_lb)
                st.session_state["samp"] = pl.DataFrame({"iwp": iwp_samples_lb})
                st.session_state["sim_df"] = st.session_state["samp"]
                st.session_state["sims_run"] = True
                st.session_state["rde_method_run"] = "Lookback data"
                # PrEP breakthrough risk is not applicable to the lookback
                # method; clear any stale PrEP run state so its results and
                # residual-risk components don't leak into the lookback display.
                st.session_state["prep_oral_run"] = False
                st.session_state["prep_inj_run"] = False
            except ValueError as e:
                st.sidebar.error(f"Error: {e}")


# Show plot in app
output_container.write("""
### Outputs
""")

# Debug only
# print(rr.mode_rounded(st.session_state["k_human"], precision = 5))
# print(rr.mode_rounded(st.session_state["k_animal"], precision = 5))


if st.session_state["sims_run"]:
    download_format = st.sidebar.selectbox(
        "Simulation download format",
        options=[
            "csv",
            "parquet",
        ],
        index=1,
        help="File format for the downloadable simulation output: CSV (universal) or Parquet (compact, typed).",
    )
    # Build combined download: baseline + any PrEP scenarios
    _dl_frames = []
    if "sim_df" in st.session_state and st.session_state["sim_df"] is not None:
        _bl_df = st.session_state["sim_df"].with_columns(
            pl.lit("baseline").alias("scenario")
        )
        _dl_frames.append(_bl_df)
    else:
        _bl_df = st.session_state["samp"].with_columns(
            pl.lit("baseline").alias("scenario")
        )
        _dl_frames.append(_bl_df)

    if (
        st.session_state.get("prep_oral_run")
        and st.session_state["sim_df_prep_oral"] is not None
    ):
        _dl_frames.append(
            st.session_state["sim_df_prep_oral"].with_columns(
                pl.lit("oral_prep").alias("scenario")
            )
        )
    if (
        st.session_state.get("prep_inj_run")
        and st.session_state["sim_df_prep_inj"] is not None
    ):
        _dl_frames.append(
            st.session_state["sim_df_prep_inj"].with_columns(
                pl.lit("injectable_prep").alias("scenario")
            )
        )

    _combined_df = pl.concat(_dl_frames, how="diagonal_relaxed")
    res_dl = convert_for_download(_combined_df, file_format=download_format)
    st.sidebar.download_button(
        label="Download simulations",
        data=res_dl,
        file_name=f"iwp_simulations.{download_format}",
        mime="text/csv" if download_format == "csv" else "application/octet-stream",
        icon=":material/download:",
    )

if not st.session_state["sims_run"]:
    st.sidebar.write("Downloads will be available once an estimation has been run.")
else:
    st.sidebar.write("Outputs are from most recent estimation run.")

if not st.session_state["sims_run"]:
    output_container.info("Run an estimation to see results.")
else:
    # Warn if the displayed results are from a different estimation method
    if st.session_state["rde_method_run"] != rde_method:
        output_container.warning(
            "Displayed results are from the "
            f"**{st.session_state['rde_method_run']}** method. "
            "Run using the selected method to update."
        )

    # Use the point estimate returned directly by risk_days_bs
    method_match = st.session_state["rde_method_run"] == rde_method
    if method_match:
        if rde_method == "Mechanistic model":
            iwp_pe = st.session_state["iwp_pe"]
        elif rde_method == "Lookback data":
            iwp_pe = st.session_state["iwp_pe_lookback"]
        else:
            iwp_pe = None
        # Save PE so it can be shown when method is switched
        st.session_state["iwp_pe_last"] = iwp_pe
    else:
        iwp_pe = st.session_state["iwp_pe_last"]

    # Interval label: Bayesian CrI for mechanistic model, frequentist CI for lookback
    interval_label = (
        "CrI" if st.session_state["rde_method_run"] == "Mechanistic model" else "CI"
    )

    iwp_cri = (
        st.session_state["samp"]["iwp"].quantile(alpha / 2),
        st.session_state["samp"]["iwp"].quantile(1 - alpha / 2),
    )
    iwp_range = (
        st.session_state["samp"]["iwp"].min(),
        st.session_state["samp"]["iwp"].max(),
    )

    if iwp_pe is not None:
        output_container.write(
            f"**Baseline** RDEs PE: {iwp_pe:.2f} ({sig_level:.0f}% {interval_label}: "
            f"{iwp_cri[0]:.2f} to {iwp_cri[1]:.2f}; "
            f"Range: {iwp_range[0]:.2f} to {iwp_range[1]:.2f})"
        )

    fig = plot_histogram(
        st.session_state["samp"],
        histnorm=plot_norm,
        title="Distribution of simulated RDEs (baseline)",
    )
    output_container.plotly_chart(fig, width="stretch")

    # --- PrEP RDE results ---
    if st.session_state.get("prep_oral_run"):
        iwp_pe_oral = st.session_state.get("iwp_pe_prep_oral")
        samp_oral = st.session_state.get("samp_prep_oral")
        if iwp_pe_oral is not None and samp_oral is not None:
            oral_cri = (
                samp_oral["iwp"].quantile(alpha / 2),
                samp_oral["iwp"].quantile(1 - alpha / 2),
            )
            oral_range = (samp_oral["iwp"].min(), samp_oral["iwp"].max())
            output_container.write(
                f"**Oral PrEP** RDEs PE: {iwp_pe_oral:.2f} "
                f"({sig_level:.0f}% {interval_label}: "
                f"{oral_cri[0]:.2f} to {oral_cri[1]:.2f}; "
                f"Range: {oral_range[0]:.2f} to {oral_range[1]:.2f})"
            )
            fig_oral = plot_histogram(
                samp_oral,
                histnorm=plot_norm,
                title="Distribution of simulated RDEs (oral PrEP)",
            )
            output_container.plotly_chart(fig_oral, width="stretch")

    if st.session_state.get("prep_inj_run"):
        iwp_pe_inj = st.session_state.get("iwp_pe_prep_inj")
        samp_inj = st.session_state.get("samp_prep_inj")
        if iwp_pe_inj is not None and samp_inj is not None:
            inj_cri = (
                samp_inj["iwp"].quantile(alpha / 2),
                samp_inj["iwp"].quantile(1 - alpha / 2),
            )
            inj_range = (samp_inj["iwp"].min(), samp_inj["iwp"].max())
            output_container.write(
                f"**Injectable PrEP** RDEs PE: {iwp_pe_inj:.2f} "
                f"({sig_level:.0f}% {interval_label}: "
                f"{inj_cri[0]:.2f} to {inj_cri[1]:.2f}; "
                f"Range: {inj_range[0]:.2f} to {inj_range[1]:.2f})"
            )
            fig_inj = plot_histogram(
                samp_inj,
                histnorm=plot_norm,
                title="Distribution of simulated RDEs (injectable PrEP)",
            )
            output_container.plotly_chart(fig_inj, width="stretch")

    if calculate_rr:
        if iwp_pe is None or iwp_pe <= 0:
            output_container.warning(
                "Residual risk cannot be calculated: IWP point estimate is zero or "
                "undefined. Run the calculation first, or check your inputs."
            )
        else:
            rr_pe, rr_cri, rr_sd = rr.residual_risk_rd(
                iwp_pe,
                st.session_state["samp"]["iwp"],
                inc_perpy,
                inc_perpy_sd,
                per=1e6,
                seed=st.session_state["seed"],
                alpha=alpha,
                one_in_x=False,
            )
            rr_onein_pe, rr_onein_cri, rr_onein_sd = rr.residual_risk_rd(
                iwp_pe,
                st.session_state["samp"]["iwp"],
                inc_perpy,
                inc_perpy_sd,
                per=None,
                seed=st.session_state["seed"],
                alpha=alpha,
                one_in_x=True,
            )

            # Build one results table instead of a long list of lines. Each
            # scenario is a row; the point estimate is shown with its credible
            # interval in parentheses, per 10^6 transfusions and as "1 in N".
            def _rr_row(label, pe, cri, onein_pe, onein_cri):
                per_m = f"{pe:.5f} ({cri[0]:.5f} – {cri[1]:.5f})"
                freq = f"1 in {onein_pe:,.0f} (1 in {onein_cri[1]:,.0f} – 1 in {onein_cri[0]:,.0f})"
                return (label, per_m, freq)

            rr_rows = [_rr_row("Baseline", rr_pe, rr_cri, rr_onein_pe, rr_onein_cri)]

            # Per-population (iwp_pe, iwp_bs, incidence, incidence_sd) for the
            # joint total-risk CrI (see rr.total_residual_risk_rd).
            total_components = [
                (
                    iwp_pe,
                    st.session_state["samp"]["iwp"].to_numpy(),
                    inc_perpy,
                    inc_perpy_sd,
                )
            ]

            if (
                st.session_state.get("prep_oral_run")
                and st.session_state.get("samp_prep_oral") is not None
            ):
                iwp_pe_oral = st.session_state.get("iwp_pe_prep_oral")
                if (
                    iwp_pe_oral is not None
                    and iwp_pe_oral > 0
                    and inc_prep_oral_perpy is not None
                ):
                    rr_oral_pe, rr_oral_cri, _ = rr.residual_risk_rd(
                        iwp_pe_oral,
                        st.session_state["samp_prep_oral"]["iwp"],
                        inc_prep_oral_perpy,
                        inc_prep_oral_perpy_sd,
                        per=1e6,
                        seed=st.session_state["seed"],
                        alpha=alpha,
                        one_in_x=False,
                    )
                    rr_oral_onein_pe, rr_oral_onein_cri, _ = rr.residual_risk_rd(
                        iwp_pe_oral,
                        st.session_state["samp_prep_oral"]["iwp"],
                        inc_prep_oral_perpy,
                        inc_prep_oral_perpy_sd,
                        per=None,
                        seed=st.session_state["seed"],
                        alpha=alpha,
                        one_in_x=True,
                    )
                    rr_rows.append(
                        _rr_row(
                            "Oral PrEP",
                            rr_oral_pe,
                            rr_oral_cri,
                            rr_oral_onein_pe,
                            rr_oral_onein_cri,
                        )
                    )
                    total_components.append((
                        iwp_pe_oral,
                        st.session_state["samp_prep_oral"]["iwp"].to_numpy(),
                        inc_prep_oral_perpy,
                        inc_prep_oral_perpy_sd,
                    ))

            if (
                st.session_state.get("prep_inj_run")
                and st.session_state.get("samp_prep_inj") is not None
            ):
                iwp_pe_inj = st.session_state.get("iwp_pe_prep_inj")
                if (
                    iwp_pe_inj is not None
                    and iwp_pe_inj > 0
                    and inc_prep_inj_perpy is not None
                ):
                    rr_inj_pe, rr_inj_cri, _ = rr.residual_risk_rd(
                        iwp_pe_inj,
                        st.session_state["samp_prep_inj"]["iwp"],
                        inc_prep_inj_perpy,
                        inc_prep_inj_perpy_sd,
                        per=1e6,
                        seed=st.session_state["seed"],
                        alpha=alpha,
                        one_in_x=False,
                    )
                    rr_inj_onein_pe, rr_inj_onein_cri, _ = rr.residual_risk_rd(
                        iwp_pe_inj,
                        st.session_state["samp_prep_inj"]["iwp"],
                        inc_prep_inj_perpy,
                        inc_prep_inj_perpy_sd,
                        per=None,
                        seed=st.session_state["seed"],
                        alpha=alpha,
                        one_in_x=True,
                    )
                    rr_rows.append(
                        _rr_row(
                            "Injectable PrEP",
                            rr_inj_pe,
                            rr_inj_cri,
                            rr_inj_onein_pe,
                            rr_inj_onein_cri,
                        )
                    )
                    total_components.append((
                        iwp_pe_inj,
                        st.session_state["samp_prep_inj"]["iwp"].to_numpy(),
                        inc_prep_inj_perpy,
                        inc_prep_inj_perpy_sd,
                    ))

            # Total (additive) residual risk with a joint credible interval.
            # The component IWP bootstrap arrays share their per-iteration k /
            # viral-dynamics / LOD / volume draws (Go backend, common seed), so
            # summing per iteration and taking quantiles yields a valid joint
            # CrI; incidence is drawn independently per population.
            #
            # TODO: oral and injectable PrEP currently share their PrEP-specific
            # draws too (eclipse, a, b, AND drug_effect), because they run off the
            # same seed. drug_effect in particular should be specified AND drawn
            # independently per scenario — see TODO.md ("Independent PrEP
            # drug_effect for oral vs injectable"). That needs the inject-shared-
            # arrays refactor (shared k/dt/lod/volume, independent PrEP draws);
            # same-seed gives shared-everything, different-seed would break the
            # shared-parameter alignment that makes this total CrI valid.
            if len(total_components) > 1:
                t_pe, t_cri, t_onein_pe, t_onein_cri = rr.total_residual_risk_rd(
                    total_components,
                    per=1e6,
                    seed=st.session_state["seed"],
                    alpha=alpha,
                )
                rr_rows.append((
                    "**Total (additive)**",
                    f"**{t_pe:.5f} ({t_cri[0]:.5f} – {t_cri[1]:.5f})**",
                    f"**1 in {t_onein_pe:,.0f} (1 in {t_onein_cri[1]:,.0f} – 1 in {t_onein_cri[0]:,.0f})**",
                ))

            output_container.subheader("Residual risk of HIV transfusion transmission")
            _rr_table = (
                "| Scenario | Residual risk (per 10⁶ transfusions) | Residual risk (1 transmission in X transfusions) |\n"
                "|:--|--:|--:|\n"
            )
            for _label, _per_m, _freq in rr_rows:
                _rr_table += f"| {_label} | {_per_m} | {_freq} |\n"
            # Centre the results table as a whole within the content column. The
            # markdown table sizes to its content, so margin:auto centres the
            # block. Scoped to this page; the results table is its only table.
            output_container.markdown(
                "<style>[data-testid='stMarkdown'] table"
                " { margin-left: auto !important; margin-right: auto !important; }</style>",
                unsafe_allow_html=True,
            )
            output_container.markdown(_rr_table)

            _caption = f"Point estimates; ranges in parentheses are {sig_level:.0f}% {interval_label}s."
            if len(total_components) > 1:
                if st.session_state.get("used_go", True):
                    _caption += (
                        f" The additive total's {interval_label} is computed per iteration "
                        "(components share the sampled k, viral-dynamics, LOD and volume "
                        "draws; the populations' incidence rates are assumed independent)."
                    )
                else:
                    _caption += (
                        f" The additive total's {interval_label} is approximate: the Python "
                        "backend does not align the shared parameters across components, so "
                        "they are treated as independent (incidences assumed independent "
                        "regardless). Use the Go backend for the exact shared-parameter interval."
                    )
            output_container.caption(_caption)

# The shared sidebar footer (VRI logo + app/library version caption) is rendered
# by the app.py router so it appears on every page.

# # Vitalant Research Institute logo — centred footer at the bottom of the page.
# st.divider()
# _, _footer_logo, _ = st.columns([1, 2, 1])
# _footer_logo.image(str(_STATIC_DIR / "vri_logo.png"), width="stretch")
