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
import random
import statistics
import time
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import scipy.stats as stats
import streamlit as st

import residualrisk as rr

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
iwp_pe_primpar_animal_default = 4.49
iwp_pe_primpar_human_default = 0.83
iwp_pe_primpar_expdecay_default = 1.81

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

if "iwp_pe_primpar" not in st.session_state:
    st.session_state["iwp_pe_primpar"] = None

if "iwp_pe_last" not in st.session_state:
    st.session_state["iwp_pe_last"] = None


@st.cache_data
def load_data():
    # Use Path to ensure files are loaded relative to this script, not cwd
    static_dir = Path(__file__).parent / "static"

    ests = pd.read_parquet(static_dir / "iwp_estimates_expdecay.parquet")
    # lowercase = lambda x: str(x).lower()
    # ests.rename(lowercase, axis='columns', inplace=True)
    # ests["id"] = ests.index
    # ests_long = ests.melt(id_vars=["id"], var_name="product", value_name="iwp")

    k_animal = np.array(pd.read_parquet(static_dir / "k_param_animal.parquet").k)
    k_human = np.array(pd.read_parquet(static_dir / "k_param_human.parquet").k)
    k_expdecay = np.array(pd.read_parquet(static_dir / "k_param_expdecay.parquet").k)

    return ests, k_animal, k_human, k_expdecay


@st.cache_data
def convert_for_download(df, file_format="csv"):
    if file_format == "csv":
        return df.to_csv().encode("utf-8")
    elif file_format == "parquet":
        return df.to_parquet()
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

if "samp" not in st.session_state:
    (
        st.session_state["samp"],
        st.session_state["k_animal"],
        st.session_state["k_human"],
        st.session_state["k_expdecay"],
    ) = load_data()
# defaults only, not after running simulations

rde_method = st.selectbox(
    "RDE estimation method",
    options=["Lookback data", "Mechanistic model", "Mechanistic model with PrEP (coming soon)"],
    index=1,
    help="Risk day quivalents (RDEs) are equivalent to the infectious window "
    "period (IWP). Lookback data: estimates the IWP directly from "
    "lookback investigation data. Mechanistic model: simulates the "
    "IWP from viral dynamics and assay parameters. ",
)

if rde_method == "Mechanistic model with PrEP (coming soon)":
    st.info("The PrEP model is not yet available. Please select another method.")
    st.stop()

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

if rde_method == "Mechanistic model":
    implementation = st.sidebar.selectbox(
        "Simulation implementation",
        options=["Go", "Python"],
        index=0,  # Go is default
        help="Placeholder help text",
    )
    use_go_acceleration = implementation == "Go"
    if use_go_acceleration:
        from residualrisk_go import find_go_binary
        if find_go_binary() is None:
            st.sidebar.warning(
                "Go binary not found. Simulations will fall back to the Python "
                "implementation, which is significantly slower."
            )

sim_param_container = st.expander(
    "Simulation settings", expanded=True, icon=":material/menu_open:"
)
if rde_method == "Mechanistic model":
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
incidence_param_container = st.expander(
    "Incidence parameters", expanded=True, icon=":material/menu_open:"
)

output_container = st.container()

with sim_param_container:
    col1, col2 = st.columns(2)

    if rde_method == "Mechanistic model":
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

    if rde_method == "Mechanistic model":
        n_threads = col2.slider(
            "Select number of CPU cores to use",
            min_value=1,
            max_value=n_cpu,
            value=n_cpu,
            step=1,
            help="Placeholder help text",
        )

if rde_method == "Mechanistic model":
    with trans_param_container:
        col1, col2 = st.columns(2)

        belov_model = col1.selectbox(
            "Select transmissibility model",
            options=["Belov animal model", "Belov human model", "Belov human-weighted"],
            index=2,
            help="Placeholder help text",
        )

        k_param_pe = col2.selectbox(
            "Transmissibility parameter: posterior...",
            options=["mean", "median", "mode"],
            index=1,
            help="Placeholder help text",
        )

        if belov_model == "Belov animal model":
            k_param = st.session_state["k_animal"]
        elif belov_model == "Belov human model":
            k_param = st.session_state["k_human"]
        elif belov_model == "Belov human-weighted":
            k_param = st.session_state["k_expdecay"]

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

button_label = "Run simulations" if rde_method == "Mechanistic model" else "Calculate RDEs"
if st.sidebar.button(button_label):
    if rde_method == "Mechanistic model":
        progressbar = st.sidebar.progress(0, text="Running simulations...")
        if k_param_pe == "mode":
            k_pe = rr.mode_rounded(
                k_param, precision=5
            )  # use mode rounded to 5 decimal places
        elif k_param_pe == "mean":
            k_pe = statistics.mean(k_param)
        elif k_param_pe == "median":
            k_pe = statistics.median(k_param)
        else:
            k_pe = None  # should error out but never happen
        (
            st.session_state["iwp_pe_primpar"],
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
            k_gamma_scale=None,
            k_gamma_shape=None,
            alpha=alpha,
            n_bs=n_sims,
            point_estimate="primary parameters",  # always run with primary parameters to get iwp_pe_primpar and store in session state -- calculate other methods in app
            seed=st.session_state["seed"],
            threads=n_threads,
            progress=progressbar,
            return_sim_df=True,
            use_go=use_go_acceleration,
        )
        st.session_state["sims_run"] = True
        st.session_state["rde_method_run"] = "Mechanistic model"
        st.session_state["samp"] = pd.DataFrame(st.session_state["bs"], columns=["iwp"])
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
                idi_df = pd.read_csv(uploaded_idi, header=None)
                idis = idi_df.iloc[:, 0].tolist()
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
                st.session_state["samp"] = pd.DataFrame(iwp_samples_lb, columns=["iwp"])
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
    st.sidebar.write(
        "Downloads will be available once an estimation has been run."
    )
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

    # Calculate point estimate based on whichever method was actually run
    method_match = st.session_state["rde_method_run"] == rde_method
    if method_match:
        if rde_method == "Mechanistic model":
            if point_estimate == "median":
                iwp_pe = st.session_state["samp"]["iwp"].median()
            elif point_estimate == "mean":
                iwp_pe = st.session_state["samp"]["iwp"].mean()
            elif point_estimate == "mode":
                iwp_pe = stats.mode(np.array(st.session_state["samp"]["iwp"]).round(2)).mode
            else:
                iwp_pe = st.session_state["iwp_pe_primpar"]
        elif rde_method == "Lookback data":
            iwp_pe = st.session_state["iwp_pe_lookback"]
        else:
            iwp_pe = None
        # Save PE so it can be shown when method is switched
        st.session_state["iwp_pe_last"] = iwp_pe
    else:
        iwp_pe = st.session_state["iwp_pe_last"]

    # Interval label: Bayesian CrI for mechanistic model, frequentist CI for lookback
    interval_label = "CrI" if st.session_state["rde_method_run"] == "Mechanistic model" else "CI"

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
