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

"""Generate the figures embedded in ``docs/theory_prep.md`` (the PrEP-breakthrough
extension). Faithful to the implementation: every curve is computed by calling the
production functions in ``residualrisk.prep``.

Run from the repo root:

    uv run --with matplotlib python docs/figures/make_prep_figures.py

Figures are written next to this script as ``prep_fig0_*.png`` … ``prep_fig3_*.png``.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

from residualrisk import prep  # noqa: E402

OUT = Path(__file__).resolve().parent

# Production parameters (ISBT 2025 analysis); see docs/theory_prep.md §10.1.
SHARED = dict(
    C0=0.00025,
    doubling_time=0.8542,
    eclipse=7.0,
    a=0.7,
    b=0.6,
    offset=1.0,
    pool_size=16,
    retests=1,
    copies_per_virion=2.0,
    z=1.6449,
    lod50=4.7 / 1.72,
    lod95_lod50_ratio=21.2 / 4.7,
)
ORAL = dict(
    set_point=336.0, ser_min=28.7, ser_max=250.0, ser_alpha=50.49434, ser_beta=1.15062
)
INJ = dict(
    set_point=25.0, ser_min=42.0, ser_max=250.0, ser_alpha=90.88988, ser_beta=3.048339
)
K_ANIMAL = 0.024464  # animal-posterior median (the PrEP analysis k)
VOL_RBC = 20.0
VOL_FFP = 200.0

plt.rcParams.update(
    {
        "figure.dpi": 130,
        "font.size": 10,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "legend.framealpha": 0.9,
    }
)


def _vl_copies(t, scen):
    """Viral concentration in copies/mL (= chi * C) along the PrEP trajectory."""
    tcrit = prep._find_tcrit(
        SHARED["eclipse"], SHARED["C0"], SHARED["doubling_time"], scen["set_point"]
    )
    c = np.array(
        [
            prep._vl_postbt(
                ti,
                SHARED["eclipse"],
                SHARED["C0"],
                SHARED["doubling_time"],
                scen["set_point"],
                SHARED["a"],
                SHARED["b"],
                SHARED["offset"],
                tcrit,
            )
            for ti in t
        ]
    )
    return SHARED["copies_per_virion"] * c


def _p_nat(t, scen):
    return np.array(
        [
            prep._prob_nondetection_prep(
                ti,
                SHARED["copies_per_virion"],
                SHARED["C0"],
                SHARED["doubling_time"],
                SHARED["eclipse"],
                scen["set_point"],
                SHARED["a"],
                SHARED["b"],
                SHARED["offset"],
                SHARED["pool_size"],
                SHARED["lod50"],
                SHARED["lod95_lod50_ratio"],
                SHARED["retests"],
                SHARED["z"],
            )
            for ti in t
        ]
    )


def _p_ser(t, scen):
    return np.array(
        [
            prep._prob_nondetection_serology_prep(
                ti,
                scen["ser_min"],
                scen["ser_max"],
                scen["ser_alpha"],
                scen["ser_beta"],
            )
            for ti in t
        ]
    )


def _p_inf(t, scen, vol):
    return np.array(
        [
            prep._prob_infectious_prep(
                ti,
                SHARED["eclipse"],
                SHARED["C0"],
                SHARED["doubling_time"],
                scen["set_point"],
                SHARED["a"],
                SHARED["b"],
                SHARED["offset"],
                vol,
                K_ANIMAL,
                SHARED["copies_per_virion"],
                1.0,
            )
            for ti in t
        ]
    )


def fig0_pipeline():
    stages = [
        ("PrEP user\ndonates", "L2"),
        ("Non-\ndisclosure", "L2"),
        ("Breakthrough\ninfection", "L2"),
        ("Escapes NAT\n& serology", "L1"),
        ("Infectious\ncomponent", "L1"),
        ("Transfusion\ntransmission", "out"),
    ]
    colors = {"L2": "#cfe8f3", "L1": "#f6d9c9", "out": "#dcdcdc"}
    fig, ax = plt.subplots(figsize=(11, 2.9))
    ax.set_xlim(0, len(stages))
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.grid(False)
    w = 0.82
    for i, (label, layer) in enumerate(stages):
        x = i + 0.5
        ax.add_patch(
            FancyBboxPatch(
                (x - w / 2, 0.34),
                w,
                0.32,
                boxstyle="round,pad=0.02",
                fc=colors[layer],
                ec="black",
                lw=1.1,
            )
        )
        ax.text(x, 0.5, label, ha="center", va="center", fontsize=9)
        if i < len(stages) - 1:
            ax.annotate(
                "",
                xy=(i + 1 + 0.5 - w / 2, 0.5),
                xytext=(x + w / 2, 0.5),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color="0.3"),
            )
    # grouping labels
    ax.text(
        1.5,
        0.82,
        "Population layer (§9) — user-supplied",
        ha="center",
        fontsize=9,
        style="italic",
        color="#16607f",
    )
    ax.text(
        4.0,
        0.18,
        "Mechanistic RDE — the tool (§§2–8)",
        ha="center",
        fontsize=9,
        style="italic",
        color="#9c4a21",
    )
    fig.suptitle("Residual-risk pipeline for undisclosed PrEP use", y=1.0)
    fig.tight_layout()
    fig.savefig(OUT / "prep_fig0_pipeline.png", bbox_inches="tight")
    plt.close(fig)


def fig1_viral_dynamics():
    t = np.linspace(0, 60, 1400)
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    for scen, lab, c in [
        (ORAL, "oral PrEP (set-point 336)", "C0"),
        (INJ, "injectable PrEP (set-point 25)", "C1"),
    ]:
        copies = _vl_copies(t, scen)
        copies = np.where(
            t < SHARED["eclipse"], np.nan, copies
        )  # line starts at eclipse
        ax.plot(t, copies, color=c, lw=2, label=lab)
    thr = (
        SHARED["pool_size"] * SHARED["lod50"]
    )  # minipool 50%-detection threshold, copies/mL
    ax.axhline(
        thr,
        color="gray",
        ls=":",
        lw=1.5,
        label=f"minipool NAT $X_{{50}}$ ≈ {thr:.0f} copies/mL",
    )
    ax.axvline(SHARED["eclipse"], color="0.6", ls="--", lw=1)
    ax.text(
        SHARED["eclipse"], 4e3, " eclipse", fontsize=8, color="0.4", ha="left", va="top"
    )
    ax.set_yscale("log")
    ax.set_ylim(0.5, 6e3)
    ax.set_xlim(0, 60)
    ax.set_xlabel("days since infection, $t$")
    ax.set_ylabel("viral concentration (copies/mL)")
    ax.set_title("Post-breakthrough viral dynamics on PrEP")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "prep_fig1_viral_dynamics.png")
    plt.close(fig)


def fig2_detection():
    t = np.linspace(0, 260, 1600)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
    for scen, lab, c in [(ORAL, "oral", "C0"), (INJ, "injectable", "C1")]:
        axes[0].plot(t, _p_nat(t, scen), color=c, lw=2, label=lab)
        axes[1].plot(t, _p_ser(t, scen), color=c, lw=2, label=lab)
    axes[0].set_title("NAT non-detection  $P_{nd}^{NAT}(t)$")
    axes[1].set_title("Serological non-detection  $S(t)$")
    for ax in axes:
        ax.set_xlabel("days since infection")
        ax.set_ylim(-0.02, 1.04)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("probability")
    fig.suptitle(
        "Detection windows extend on PrEP: blunted viral load (NAT) and delayed seroconversion (serology)",
        y=1.02,
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(OUT / "prep_fig2_detection.png", bbox_inches="tight")
    plt.close(fig)


def fig3_iwp_construction():
    # Both panels are oral PrEP (set-point 336) at the nominal parameters, so the
    # RDEs (~5 / ~8 d) are representative of the reported medians; they differ only
    # in transfused volume (RBC 20 mL vs plasma 200 mL), illustrating the volume
    # effect. (The injectable nominal set-point of 25 sits at the extreme low end of
    # its sampled U(5, 2500) range and gives an unrepresentative ~55-day point
    # window, so it is not used here; see theory_prep.md §7, §10.)
    t = np.linspace(0, 80, 1400)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
    for ax, vol, title in [
        (axes[0], VOL_RBC, "oral PrEP, RBC (20 mL)"),
        (axes[1], VOL_FFP, "oral PrEP, plasma (200 mL)"),
    ]:
        p_inf, p_nat, p_ser = _p_inf(t, ORAL, vol), _p_nat(t, ORAL), _p_ser(t, ORAL)
        joint = p_inf * p_nat * p_ser
        rde = prep._risk_days_prep(
            SHARED["copies_per_virion"],
            SHARED["C0"],
            SHARED["doubling_time"],
            ORAL["set_point"],
            SHARED["eclipse"],
            SHARED["a"],
            SHARED["b"],
            SHARED["offset"],
            vol,
            K_ANIMAL,
            SHARED["pool_size"],
            SHARED["lod50"],
            SHARED["lod95_lod50_ratio"],
            SHARED["retests"],
            ORAL["ser_min"],
            ORAL["ser_max"],
            ORAL["ser_alpha"],
            ORAL["ser_beta"],
            SHARED["z"],
        )
        ax.plot(t, p_inf, "--", color="C0", lw=1.4, label="infectious")
        ax.plot(t, p_nat, "--", color="C1", lw=1.4, label="NAT non-detection")
        ax.plot(t, p_ser, "--", color="C2", lw=1.4, label="serology non-detection")
        ax.plot(t, joint, "-", color="C3", lw=2.2, label="joint (integrand)")
        ax.fill_between(t, 0, joint, color="C3", alpha=0.22)
        ax.set_title(f"{title}\nRDE $= {rde:.2f}$ days")
        ax.set_xlabel("days since infection")
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=8, loc="upper right")
    axes[0].set_ylabel("probability")
    fig.suptitle(
        "PrEP risk-day-equivalents: area under (infectious × NAT-undetected × serology-undetected)",
        y=1.02,
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(OUT / "prep_fig3_iwp_construction.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig0_pipeline()
    fig1_viral_dynamics()
    fig2_detection()
    fig3_iwp_construction()
    print(f"PrEP figures written to {OUT}")
