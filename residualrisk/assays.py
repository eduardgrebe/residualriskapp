# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025-2026  Vitalant and Eduard Grebe Consulting
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

"""Canned NAT-assay limit-of-detection presets (HIV-1 Group M).

Single source of truth for the published 50%/95% limits of detection of the
supported blood-screening NAT assays, in copies/mL. Both the Streamlit UI
(``estimator.py``) and the public API consume this table, so a caller can select
an assay by slug instead of transcribing the LoD numbers.

The 50% LoD SDs are derived from each manufacturer's 95% CI of the 50% LoD
(reported in IU/mL): the coefficient of variation
``CoV = (CI_hi - CI_lo) / 3.92 / PE`` is invariant under the (multiplicative)
IU<->copies conversion, so it is computed in IU/mL and applied to the copies/mL
point estimate (``lod50_sd = CoV * lod50``). The Procleix Ultrio value uses the
discriminatory dHIV-1 (Tigris) CI; cobas MPX uses the EDTA-plasma CI.

``cp_per_iu`` is the IU/mL -> copies/mL conversion factor applied (upstream) to
the manufacturer IU/mL LoDs to obtain the copies/mL values stored here. It is
NOT constant: HIV-1 Group M is calibrated against three WHO International
Standard generations, each with its own factor (0.6 for the 1st IS 97/656, 0.58
for the 2nd IS 97/650, 0.35 for the 3rd IS 10/152). ``iu_std`` records that
standard. ``cp_per_iu`` and ``iu_std`` are informational (surfaced in the UI);
they feed no calculation. See ``../residualrisk_analysis/assays/ASSAYS.qmd`` for
the full conversion rationale.
"""

from typing import NamedTuple

# slug -> entry. ``display_name`` drives the UI label; the remaining fields are
# documented in the module docstring above.
#
# PROVISIONAL SD -- "biomanguinhos" (Brazilian NAT Platform, Bio-Manguinhos): Rocha
# et al. (2018, Transfusion vol. 58; Table 2, p.865) publishes only point 50%/95%
# LoDs for HIV-1 -- 46.77 / 95.86 IU/mL vs WHO 2nd IS 97/650 -- from a PROBIT on
# just 24 replicates/dilution, with no CI, fiducial limits, or per-dilution
# hit-rate table. So unlike every other assay (whose SD comes from a manufacturer
# 95% CI), its LoD50 SD cannot be derived and must be assumed. Two values, as
# relative SD (RSE = SD / LoD50):
#   * 4.95 IU/mL = 2.88 copies/mL, RSE 10.6% -- originally ballparked; used in
#     prior analyses.
#   * 6.08 IU/mL = 3.527 copies/mL, RSE 13% -- USED HERE (lod50_sd below). The
#     24-replicate study warrants a larger relative SD than the well-powered
#     assays (RSE ~5-8%); anchoring on Procleix Ultrio (similar steep curve,
#     ~120 reps/dilution, RSE ~7.3%) and scaling RSE ~ sigma * sqrt(1/N) gives
#     ~13%.
# REVISIT if the authors supply the per-dilution hit-rate table: probit-refit it
# for a proper delta-method SE (as residualrisk_analysis/assays/ASSAYS.qmd does
# for the cobas TaqScreen assays, which likewise publish no 50% LoD).
NAT_ASSAYS = {
    "ultrio": {"display_name": "Procleix Ultrio (Tigris)", "lod50": 5.0, "lod50_sd": 0.364, "lod95": 12.2, "cp_per_iu": 0.6, "iu_std": "WHO 1st IS 97/656 (dHIV-1)"},
    "ultrio_plus": {"display_name": "Procleix Ultrio Plus (Tigris)", "lod50": 2.7, "lod50_sd": 0.191, "lod95": 12.3, "cp_per_iu": 0.58, "iu_std": "WHO 2nd IS 97/650"},
    "ultrio_elite": {"display_name": "Procleix Ultrio Elite (Panther)", "lod50": 3.1, "lod50_sd": 0.234, "lod95": 10.4, "cp_per_iu": 0.58, "iu_std": "WHO 2nd IS 97/650"},
    "cobas_taqscreen_mpx": {"display_name": "cobas TaqScreen MPX (s 201)", "lod50": 5.5, "lod50_sd": 0.385, "lod95": 29.4, "cp_per_iu": 0.6, "iu_std": "WHO 1st IS 97/656"},
    "cobas_taqscreen_mpxv2": {"display_name": "cobas TaqScreen MPX v2.0 (s 201)", "lod50": 5.3, "lod50_sd": 0.250, "lod95": 26.8, "cp_per_iu": 0.58, "iu_std": "WHO 2nd IS 97/650"},
    "cobas_mpx": {"display_name": "cobas MPX (5800/6800/8800)", "lod50": 1.3, "lod50_sd": 0.0785, "lod95": 9.0, "cp_per_iu": 0.35, "iu_std": "WHO 3rd IS 10/152"},
    "biomanguinhos": {"display_name": "Brazilian NAT Platform (Bio-Manguinhos)", "lod50": 27.13, "lod50_sd": 3.527, "lod95": 55.6, "cp_per_iu": 0.58, "iu_std": "WHO 2nd IS 97/650"},  # SD provisional (RSE 13%) -- see note above
}


class AssayLoD(NamedTuple):
    """Resolved limit-of-detection parameters for a single canned NAT assay.

    ``lod50``, ``lod50_sd`` and ``lod95_lod50_ratio`` are the three values
    ``risk_days_bs`` consumes; ``lod95``, ``cp_per_iu`` and ``iu_std`` are
    carried for display/provenance.
    """

    slug: str
    display_name: str
    lod50: float
    lod50_sd: float
    lod95: float
    lod95_lod50_ratio: float
    cp_per_iu: float
    iu_std: str


def lods_for_assay(assay):
    """Resolve a canned-assay slug to its limit-of-detection parameters.

    Parameters
    ----------
    assay : str
        One of the keys of ``NAT_ASSAYS`` (e.g. ``"ultrio_elite"``). See
        ``list_assays()`` for the slug -> display-name mapping.

    Returns
    -------
    AssayLoD
        Named tuple with ``lod50``, ``lod50_sd``, ``lod95``,
        ``lod95_lod50_ratio`` (= ``lod95 / lod50``), and the ``cp_per_iu`` /
        ``iu_std`` provenance fields.

    Raises
    ------
    ValueError
        If ``assay`` is not a known slug.
    """
    try:
        entry = NAT_ASSAYS[assay]
    except KeyError:
        raise ValueError(
            f"Unknown assay {assay!r}. Available: {', '.join(NAT_ASSAYS)}"
        ) from None
    return AssayLoD(
        slug=assay,
        display_name=entry["display_name"],
        lod50=entry["lod50"],
        lod50_sd=entry["lod50_sd"],
        lod95=entry["lod95"],
        lod95_lod50_ratio=entry["lod95"] / entry["lod50"],
        cp_per_iu=entry["cp_per_iu"],
        iu_std=entry["iu_std"],
    )


def list_assays():
    """Return ``{slug: display_name}`` for every canned NAT assay.

    Convenience for building menus (e.g. from R via reticulate).
    """
    return {slug: entry["display_name"] for slug, entry in NAT_ASSAYS.items()}
