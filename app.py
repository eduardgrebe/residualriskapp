# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025  Vitalant and Eduard Grebe Consulting
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
# streamlit run app/app.py
import io
import random
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

APP_VERSION = "1.1.0.dev0"

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


@st.cache_data
def load_data():
    # Use Path to ensure files are loaded relative to this script, not cwd
    static_dir = Path(__file__).parent / "static"

    k_animal = pl.read_parquet(static_dir / "k_param_animal.parquet", columns=["k"]).get_column("k").to_numpy()
    k_human = pl.read_parquet(static_dir / "k_param_human.parquet", columns=["k"]).get_column("k").to_numpy()
    k_expdecay = pl.read_parquet(static_dir / "k_param_expdecay.parquet", columns=["k"]).get_column("k").to_numpy()

    # KDE modes via Go binary (~1.5s total, 30× faster than Python KDE).
    # Falls back to hardcoded values if Go binary is unavailable.
    _go_bin = rr.find_go_binary()
    if _go_bin is not None:
        k_human_mode = rr.mode_kde_go(k_human, cap=None, n_grid=1_000_000)
        k_animal_mode = rr.mode_kde_go(k_animal, cap=None, n_grid=1_000_000)
        k_expdecay_mode = rr.mode_kde_go(k_expdecay, cap=None, n_grid=1_000_000)
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


def plot_histogram(data, x="iwp", colour=None, histnorm=None):
    fig = px.histogram(
        data,
        x=x,
        color=colour,
        histnorm=histnorm,
        barmode="overlay",
        labels={"iwp": "risk day equivalents (RDEs)"},
        title="Distribution of RDEs",
    )
    return fig


header_container = st.container()

header_container.write("""
# Residual HIV-TT Risk Estimator
Tool for estimating the residual risk of HIV transfusion transmission with NAT screening.
""")

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
        "Mechanistic model with PrEP",
    ],
    index=1,
    help="Risk day quivalents (RDEs) are equivalent to the infectious window "
    "period (IWP). Lookback data: estimates the IWP directly from "
    "lookback investigation data. Mechanistic model: simulates the "
    "IWP from viral dynamics and assay parameters. ",
)

is_mechanistic_ui = rde_method in ("Mechanistic model", "Mechanistic model with PrEP")

st.sidebar.write("Number of CPU cores: ", n_cpu)

st.session_state["seed"] = st.sidebar.number_input(
    "Specify a seed value:",
    min_value=1,
    max_value=999999,
    value=st.session_state["seed"],
    step=1,
    help="Placeholder help text",
)
if st.sidebar.button("Generate random seed"):
    st.session_state["seed"] = random.randint(1, 999999)

if is_mechanistic_ui:
    implementation = st.sidebar.selectbox(
        "Simulation implementation",
        options=["Go", "Python"],
        index=0,  # Go is default
        help="Placeholder help text",
    )
    use_go_acceleration = implementation == "Go"
    if use_go_acceleration:
        if rr.find_go_binary() is None:
            st.sidebar.warning(
                "Go binary not found. Simulations will fall back to the Python "
                "implementation, which is significantly slower."
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
if rde_method == "Mechanistic model with PrEP":
    prep_param_container = st.expander(
        "PrEP parameters", expanded=True, icon=":material/menu_open:"
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
            index=0,
            help="Placeholder help text",
        )
    else:
        point_estimate = None

    plot_norm = col1.selectbox(
        "Normalise histogram",
        options=[None, "probability density"],
        index=1,
        help="Placeholder help text",
    )

    alpha = col1.number_input(
        "Significance level (𝛼)",
        min_value=0.00,
        max_value=0.20,
        value=0.05,
        step=0.01,
        help="Placeholder help text",
    )

    sig_level = round((1 - alpha) * 100)

    n_sims = col2.select_slider(
        "Select number of simulations",
        options=[1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000],
        value=25000,
        help="Placeholder help text",
    )

    if is_mechanistic_ui:
        n_threads = col2.slider(
            "Select number of CPU cores to use",
            min_value=1,
            max_value=n_cpu,
            value=n_cpu,
            step=1,
            help="Placeholder help text",
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
            help="Placeholder help text",
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
                ["mode", "median", "mean"] if k_invgamma_alpha > 1.0 else ["mode", "median"]
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
            _LN_W_DEF   = 0.90
            _LN_MU1_DEF = -7.2403
            _LN_S1_DEF  = 0.3241
            _LN_MU2_DEF = -3.7423
            _LN_S2_DEF  = 0.5258

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
                k_lnmix_mu1    = _LN_MU1_DEF
                k_lnmix_sigma1 = _LN_S1_DEF
                k_lnmix_mu2    = _LN_MU2_DEF
                k_lnmix_sigma2 = _LN_S2_DEF

            # Derived statistics display
            import math as _math
            _lnmix_comp1_median = _math.exp(k_lnmix_mu1)
            _lnmix_comp2_median = _math.exp(k_lnmix_mu2)
            _lnmix_mean = (
                k_lnmix_w * _math.exp(k_lnmix_mu1 + k_lnmix_sigma1**2 / 2)
                + (1 - k_lnmix_w) * _math.exp(k_lnmix_mu2 + k_lnmix_sigma2**2 / 2)
            )
            st.caption(
                f"Component 1 median: {_lnmix_comp1_median:.6f} &nbsp;|&nbsp; "
                f"Component 2 median: {_lnmix_comp2_median:.6f} &nbsp;|&nbsp; "
                f"Mixture mean: {_lnmix_mean:.6f}"
            )

        else:
            k_lnmix_w      = None
            k_lnmix_mu1    = None
            k_lnmix_sigma1 = None
            k_lnmix_mu2    = None
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
            help="Placeholder help text",
        )
        volume_range_default = (round(0.75 * volume_pe), round(1.5 * volume_pe))
        volume_range = col1.slider(
            "Range of volumes transfused (mL)",
            min_value=1,
            max_value=500,
            value=volume_range_default,
            help="Placeholder help text",
        )

        doubling_time_hours = col1.number_input(
            "Viral outgrowth doubling time (hours)",
            min_value=12.0,
            max_value=48.0,
            value=20.5,
            step=0.25,
            help="Placeholder help text",
        )
        doubling_time = doubling_time_hours / 24
        doubling_time_hours_sd = col1.number_input(
            "Viral outgrowth doubling time SD",
            min_value=0.0,
            max_value=10.0,
            value=1.33,
            step=0.01,
            help="Placeholder help text",
        )
        doubling_time_sd = doubling_time_hours_sd / 24

        id_nat = col1.checkbox(
            "Individual donation NAT screening",
            value=False,
            help="Placeholder help text",
        )

        if not id_nat:
            pool_size = col2.number_input(
                "Minipool size",
                min_value=2,
                max_value=96,
                value=16,
                step=1,
                help="Placeholder help text",
            )
            retests = col2.number_input(
                "Number of retests (pool resolution)",
                min_value=0,
                max_value=5,
                value=1,
                step=1,
                help="Placeholder help text",
            )
        else:
            pool_size = 1
            retests = 0

        lod50 = col2.number_input(
            "NAT assay 50% LoD (copies/mL)",
            min_value=0.0,
            max_value=500.0,
            value=2.73,
            step=0.01,
            help="Placeholder help text",
        )
        lod50_sd = col2.number_input(
            "NAT assay 50% LoD SD (copies/mL)",
            min_value=0.0,
            max_value=500.0,
            value=0.193,
            step=0.001,
            help="Placeholder help text",
        )
        lod95 = col2.number_input(
            "NAT assay 95% LoD (copies/mL)",
            min_value=0.0,
            max_value=500.0,
            value=12.33,
            step=0.01,
            help="Placeholder help text",
        )
        st.text("95% LoD : 50% LoD ratio will be fixed for simulations.")
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

if rde_method == "Mechanistic model with PrEP":
    with prep_param_container:
        col1, col2 = st.columns(2)

        eclipse = col1.number_input(
            "Eclipse period (days)",
            min_value=1,
            max_value=30,
            value=7,
            step=1,
        )
        eclipse_range = col1.slider(
            "Eclipse period range (days)",
            min_value=1,
            max_value=20,
            value=(4, 10),
            step=1,
        )
        vl_setpoint_oral = col1.number_input(
            "oPrEP viral load setpoint (c/mL)",
            min_value=1,
            max_value=5000,
            value=340,
            step=10,
        )
        vl_setpoint_range_oral = col1.slider(
            "oPrEP viral load setpoint range (c/mL)",
            min_value=1,
            max_value=5000,
            value=(10, 2270),
            step=10,
        )
        vl_setpoint_inj = col1.number_input(
            "iPrEP viral load setpoint (c/mL)",
            min_value=1,
            max_value=5000,
            value=30,
            step=10,
        )
        vl_setpoint_range_inj = col1.slider(
            "iPrEP viral load setpoint range (c/mL)",
            min_value=1,
            max_value=5000,
            value=(10, 2500),
            step=10,
        )

        seroconversion_min_oral = col2.number_input(
            "oPrEP time to seroconversion min (days)",
            min_value=0,
            max_value=500,
            value=29,
            step=1,
        )
        seroconversion_max_oral = col2.number_input(
            "oPrEP time to seroconversion max (days)",
            min_value=0,
            max_value=500,
            value=250,
            step=1,
        )
        seroconversion_weibull_alpha_oral = col2.number_input(
            "oPrEP time to seroconversion Weibul shape (α)",
            min_value=0.0,
            max_value=500.0,
            value=50.49434,
            step=0.001,
        )
        seroconversion_weibull_beta_oral = col2.number_input(
            "oPrEP time to seroconversion Weibul scale (β)",
            min_value=0.0,
            max_value=500.0,
            value=1.15062,
            step=0.001,
        )

        seroconversion_min_inj = col2.number_input(
            "iPrEP time to seroconversion min (days)",
            min_value=0,
            max_value=500,
            value=42,
            step=1,
        )
        seroconversion_max_inj = col2.number_input(
            "iPrEP time to seroconversion max (days)",
            min_value=0,
            max_value=500,
            value=250,
            step=1,
        )
        seroconversion_weibull_alpha_inj = col2.number_input(
            "iPrEP time to seroconversion Weibul shape (α)",
            min_value=0.0,
            max_value=500.0,
            value=90.88988,
            step=0.001,
        )
        seroconversion_weibull_beta_inj = col2.number_input(
            "iPrEP time to seroconversion Weibul scale (β)",
            min_value=0.0,
            max_value=500.0,
            value=3.048339,
            step=0.001,
        )

with incidence_param_container:
    calculate_rr = st.checkbox(
        "Calculate residual risk (incidence x RDEs)",
        value=False,
        help="Placeholder help text",
    )
    if calculate_rr:
        inc_per100k = st.number_input(
            "Incidence (/100,000 PY)",
            min_value=0.001,
            max_value=10000.000,
            value=2.500,
            step=0.5,
            help="Placeholder help text",
        )
        inc_perpy = inc_per100k / 100000
        inc_perpd = inc_per100k / 100000 / 365.25
        inc_per100k_sd = st.number_input(
            "Incidence (/100,000 PY) SD",
            min_value=0.001,
            max_value=10000.000,
            value=inc_per100k * 0.2,
            step=0.01,
            help="Placeholder help text",
        )
        inc_perpy_sd = inc_per100k_sd / 100000
        inc_perpd_sd = inc_per100k_sd / 100000 / 365.25
        st.text(
            f"Relative standard error on incidence: {inc_per100k_sd / inc_per100k * 100:.1f}%"
        )

button_label = (
    "Run simulations" if is_mechanistic_ui else "Calculate RDEs"
)
if st.sidebar.button(button_label):
    if rde_method == "Mechanistic model":
        progressbar = st.sidebar.progress(0, text="Running simulations...")
        if k_param_dist == "invgamma":
            if k_invgamma_pe_choice == "mode":
                k_pe = k_invgamma_beta / (k_invgamma_alpha + 1)
            elif k_invgamma_pe_choice == "median":
                k_pe = stats.invgamma.ppf(0.5, a=k_invgamma_alpha, scale=k_invgamma_beta)
            elif k_invgamma_pe_choice == "mean":
                k_pe = k_invgamma_beta / (k_invgamma_alpha - 1)
            else:
                k_pe = k_invgamma_beta / (k_invgamma_alpha + 1)  # fallback to mode
        elif k_param_dist == "lnmixture":
            import math as _math
            if k_lnmix_pe_choice == "mean":
                k_pe = (
                    k_lnmix_w * _math.exp(k_lnmix_mu1 + k_lnmix_sigma1**2 / 2)
                    + (1 - k_lnmix_w) * _math.exp(k_lnmix_mu2 + k_lnmix_sigma2**2 / 2)
                )
            else:
                # Numerical mode or median from a large sample
                # Use cached default-param values if parameters are at defaults to avoid delay
                _lnmix_defaults = (0.90, -7.2403, 0.3241, -3.7423, 0.5258)
                _lnmix_current = (
                    k_lnmix_w, k_lnmix_mu1, k_lnmix_sigma1, k_lnmix_mu2, k_lnmix_sigma2
                )
                if _lnmix_current == _lnmix_defaults:
                    _lnmix_sample = st.session_state.get("k_lnmix_default_sample")
                    if _lnmix_sample is None:
                        _lnmix_sample = rr.sample_lnmix(100_000, *_lnmix_defaults, seed=42)
                        st.session_state["k_lnmix_default_sample"] = _lnmix_sample
                else:
                    _lnmix_sample = rr.sample_lnmix(
                        100_000, k_lnmix_w, k_lnmix_mu1, k_lnmix_sigma1,
                        k_lnmix_mu2, k_lnmix_sigma2, seed=42
                    )
                if k_lnmix_pe_choice == "median":
                    k_pe = float(np.median(_lnmix_sample))
                else:  # mode
                    k_pe = rr.mode_kde(_lnmix_sample, cap=1_000_000, n_grid=1_000_000)
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
        st.session_state["samp"] = pl.DataFrame({"iwp": st.session_state["bs"]})
        # Fallback: if sim_df is None (e.g., from Go implementation), use samp
        if st.session_state["sim_df"] is None:
            st.session_state["sim_df"] = st.session_state["samp"]
        progressbar.progress(1.0, text="Simulations complete!")
        # Brief pause to show completion, then clear progress bar
        time.sleep(0.3)
        progressbar.empty()

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
        help="Placeholder help text",
    )
    if "sim_df" in st.session_state and st.session_state["sim_df"] is not None:
        res_dl = convert_for_download(
            st.session_state["sim_df"], file_format=download_format
        )
    else:
        res_dl = convert_for_download(
            st.session_state["samp"], file_format=download_format
        )
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
            f"RDEs PE: {iwp_pe:.2f} ({sig_level:.0f}% {interval_label}: "
            f"{iwp_cri[0]:.2f} to {iwp_cri[1]:.2f}; "
            f"Range: {iwp_range[0]:.2f} to {iwp_range[1]:.2f})"
        )

    fig = plot_histogram(st.session_state["samp"], histnorm=plot_norm)
    output_container.plotly_chart(fig, width="stretch")

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
            output_container.write(
                f"RR PE: {rr_pe:.5f} /million transfusions ({sig_level:.0f}% {interval_label}: {rr_cri[0]:.5f} to {rr_cri[1]:.5f})"
            )
            output_container.write(
                f"RR PE: 1 transmission in {rr_onein_pe:,.0f} transfusions ({sig_level:.0f}% {interval_label}: {rr_onein_cri[1]:,.0f} to {rr_onein_cri[0]:,.0f})"
            )

st.sidebar.divider()
st.sidebar.caption(f"App v{APP_VERSION} · Library v{rr.__version__}")
