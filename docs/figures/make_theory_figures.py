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

"""Generate the figures embedded in ``docs/theory.md`` (base mechanistic model
and k input distributions).

Run from the repo root:

    uv run python docs/figures/make_theory_figures.py

Figures are written next to this script as ``fig1_*.png`` … ``fig3_*.png``.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
import scipy.stats as stats  # noqa: E402

from residualrisk import core as rr  # noqa: E402

OUT = Path(__file__).resolve().parent
STATIC = OUT.parent.parent / "static"

# Representative MP-NAT / pRBC scenario (Procleix Ultrio Plus, pool of 16).
SCEN = dict(
    C0=0.00025, doubling_time=0.8542, volume_transfused=20.0,
    copies_per_virion=2, pool_size=16, retests=1,
    lod50=2.73, lod95_lod50_ratio=12.33 / 2.73, z=1.6449,
)
K_HUMAN, K_ANIMAL = 0.000673, 0.020918  # human / animal posterior modes

plt.rcParams.update({
    "figure.dpi": 130, "font.size": 10,
    "axes.grid": True, "grid.alpha": 0.3, "legend.framealpha": 0.9,
})


def fig1_dose_response():
    n = np.logspace(0, 5, 600)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for k, lab, c in [(K_HUMAN, f"k = {K_HUMAN:.2g} (human mode)", "C0"),
                      (K_ANIMAL, f"k = {K_ANIMAL:.2g} (animal mode)", "C1")]:
        ax.plot(n, 1 - np.exp(-k * n), label=lab, color=c, lw=2)
    ax.set_xscale("log")
    ax.set_xlabel("RNA copies transfused, $n$")
    ax.set_ylabel(r"$P_\mathrm{infectious}(n) = 1 - e^{-kn}$")
    ax.set_ylim(0, 1.02)
    ax.set_title("HIV transfusion-transmission dose–response")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "fig1_dose_response.png")
    plt.close(fig)


def fig2_iwp_construction():
    t = np.linspace(-5, 35, 1200)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
    for ax, k, title in [(axes[0], K_HUMAN, "human-posterior mode  $k$"),
                         (axes[1], K_ANIMAL, "animal-posterior mode  $k$")]:
        p_inf = np.array([rr._prob_infectious(
            ti, SCEN["C0"], SCEN["doubling_time"], SCEN["volume_transfused"],
            k, SCEN["copies_per_virion"]) for ti in t])
        p_nd = np.array([rr._prob_nondetection(
            ti, SCEN["copies_per_virion"], SCEN["C0"], SCEN["doubling_time"],
            SCEN["pool_size"], SCEN["lod50"], SCEN["lod95_lod50_ratio"],
            SCEN["retests"], SCEN["z"]) for ti in t])
        joint = p_inf * p_nd
        iwp = rr._risk_days(
            SCEN["copies_per_virion"], SCEN["C0"], SCEN["doubling_time"],
            SCEN["volume_transfused"], k, SCEN["pool_size"], SCEN["lod50"],
            SCEN["lod95_lod50_ratio"], SCEN["retests"])
        ax.plot(t, p_inf, "--", color="C0", label="infectivity")
        ax.plot(t, p_nd, "--", color="C1", label="non-detection")
        ax.plot(t, joint, "-", color="C3", lw=2, label="joint (integrand)")
        ax.fill_between(t, 0, joint, color="C3", alpha=0.25)
        ax.set_title(f"{title}\nIWP $= {iwp:.2f}$ days")
        ax.set_xlabel("days since infection")
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=8, loc="upper left")
    axes[0].set_ylabel("probability")
    fig.suptitle("Infectious window period: overlap of infectivity and non-detection",
                 y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_iwp_construction.png", bbox_inches="tight")
    plt.close(fig)


def _kde_logx(data, grid):
    kde = stats.gaussian_kde(np.log(data), bw_method="silverman")
    return kde(np.log(grid)) / grid  # change of variables back to k-scale


def fig3_k_distributions():
    kh = pl.read_parquet(STATIC / "k_param_human.parquet")["k"].to_numpy()
    ka = pl.read_parquet(STATIC / "k_param_animal.parquet")["k"].to_numpy()
    ig = stats.invgamma(a=2.0, scale=0.002019)               # Recommendation A
    hc = stats.lognorm(s=0.3241, scale=np.exp(-7.2403))      # mixture human comp
    ac = stats.lognorm(s=0.5258, scale=np.exp(-3.7423))      # mixture animal comp

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))

    grid = np.logspace(-4, -0.5, 700)
    ax = axes[0]
    ax.plot(grid, _kde_logx(kh, grid), color="green", lw=2, label="human posterior")
    ax.plot(grid, _kde_logx(ka, grid), color="crimson", lw=2, label="animal posterior")
    ax.plot(grid, ig.pdf(grid), color="C0", ls="--", lw=2, label="InvGamma(α=2)  [Rec A]")
    ax.plot(grid, 0.9 * hc.pdf(grid) + 0.1 * ac.pdf(grid), color="C1", ls=":", lw=2,
            label="90/10 LN mixture  [Rec B]")
    ax.set_xscale("log")
    ax.set_xlabel("$k$")
    ax.set_ylabel("density")
    ax.set_title("Input distributions for $k$ (PDF, log-$x$)")
    ax.legend(fontsize=8)

    xs = np.logspace(-4, -0.3, 500)
    ax = axes[1]
    ax.plot(xs, ig.sf(xs), color="C0", ls="--", lw=2, label="InvGamma(α=2)")
    ax.plot(xs, 0.9 * hc.sf(xs) + 0.1 * ac.sf(xs), color="C1", ls=":", lw=2,
            label="90/10 mixture")
    for p in (np.percentile(ka, 5), np.median(ka), np.percentile(ka, 95)):
        ax.axvline(p, color="crimson", alpha=0.4, lw=1)  # animal P5/P50/P95 (see caption)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("$k$")
    ax.set_ylabel("$P(K > k)$")
    ax.set_title("Tail behaviour (survival function)")
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(OUT / "fig3_k_distributions.png")
    plt.close(fig)


if __name__ == "__main__":
    fig1_dose_response()
    fig2_iwp_construction()
    fig3_k_distributions()
    print(f"Figures written to {OUT}")
