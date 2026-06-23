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

"""
Provenance and integrity tests for the canned NAT assay limit-of-detection
table (``NAT_ASSAYS``) defined in ``estimator.py``.

The 50% LoD SD for each assay is *not* an arbitrary constant: it is derived
from the manufacturer's reported 95% confidence interval of the 50% LoD (in
IU/mL) via the coefficient of variation

    CoV = (CI_hi - CI_lo) / (2 * 1.96) / PE_iu

which is invariant under the (multiplicative) IU<->copies conversion, applied
to the copies/mL point estimate:

    lod50_sd = CoV * lod50  (copies/mL)

These tests re-derive each SD from the source CIs and assert it matches the
value baked into the table, pinning the arithmetic and documenting where the
numbers came from. A typo or a future "just tweak this number" edit will fail
here.

``estimator.py`` is a Streamlit script that executes ``st.*`` calls at import
time, so it cannot be imported in a plain pytest process. We therefore read the
two module-level constants out of its source with ``ast`` (no execution, no
Streamlit dependency).
"""

import ast
import math
from pathlib import Path

import pytest

ESTIMATOR_PATH = Path(__file__).resolve().parent.parent / "estimator.py"


def _load_constant(name):
    """Extract a module-level literal constant from estimator.py without
    importing (and thus executing) the Streamlit script."""
    source = ESTIMATOR_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(ESTIMATOR_PATH))
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == name:
                    return ast.literal_eval(node.value)
    raise AssertionError(f"{name} not found at module level in {ESTIMATOR_PATH}")


NAT_ASSAYS = _load_constant("NAT_ASSAYS")
MANUAL_LOD_OPTION = _load_constant("MANUAL_LOD_OPTION")

# Source data backing each SD: the manufacturer 50% LoD point estimate and 95%
# CI **in IU/mL** (HIV-1 Group M). For Procleix Ultrio the discriminatory
# dHIV-1 (Tigris) CI is used; for cobas MPX the EDTA-plasma CI is used. These
# are the figures from which the copies/mL SDs in NAT_ASSAYS were derived.
SOURCE_LOD50_CI_IU = {
    # assay name: (point_estimate, ci_low, ci_high)  -- IU/mL
    "Procleix Ultrio (Tigris)": (8.4, 7.2, 9.6),
    "Procleix Ultrio Plus (Tigris)": (4.7, 4.0, 5.3),
    "Procleix Ultrio Elite (Panther)": (5.4, 4.5, 6.1),
    "cobas TaqScreen MPX (s 201)": (9.1, 8.0, 10.5),
    "cobas TaqScreen MPX v2.0 (s 201)": (9.2, 8.4, 10.1),
    "cobas MPX (5800/6800/8800)": (3.8, 3.4, 4.3),
}

# Key used to pre-populate the manual-entry fields (see estimator.py).
MANUAL_DEFAULT_ASSAY = "cobas MPX (5800/6800/8800)"

# 95% normal-approximation divisor: a 95% CI spans +/- 1.96 SD.
_CI_DIVISOR = 2 * 1.96


@pytest.mark.parametrize("name", list(SOURCE_LOD50_CI_IU))
def test_lod50_sd_matches_cov_derivation(name):
    """Each table SD equals CoV(IU/mL CI) * lod50(copies/mL), to the precision
    at which it is stored (3-4 decimal places)."""
    pe_iu, lo_iu, hi_iu = SOURCE_LOD50_CI_IU[name]
    cov = (hi_iu - lo_iu) / _CI_DIVISOR / pe_iu

    lod50 = NAT_ASSAYS[name]["lod50"]
    derived_sd = cov * lod50
    stored_sd = NAT_ASSAYS[name]["lod50_sd"]

    # Stored values are rounded to 3 or 4 decimal places; accept either.
    candidates = [round(derived_sd, n) for n in (3, 4)]
    assert any(math.isclose(stored_sd, c, rel_tol=0, abs_tol=1e-9) for c in candidates), (
        f"{name}: stored lod50_sd={stored_sd} does not match CoV derivation "
        f"{derived_sd:.6f} (CoV={cov:.6f}, lod50={lod50}); "
        f"expected one of {candidates}"
    )


def test_source_and_table_cover_the_same_assays():
    """The provenance table and the canned table must stay in lockstep — no
    canned assay without a documented source CI, and vice versa."""
    assert set(SOURCE_LOD50_CI_IU) == set(NAT_ASSAYS)


@pytest.mark.parametrize("name", list(NAT_ASSAYS))
def test_lod_invariants(name):
    """Basic physical sanity for every canned assay."""
    entry = NAT_ASSAYS[name]
    lod50 = entry["lod50"]
    lod50_sd = entry["lod50_sd"]
    lod95 = entry["lod95"]

    assert lod50 > 0, f"{name}: lod50 must be positive"
    assert lod50_sd > 0, f"{name}: lod50_sd must be positive"
    assert lod95 > lod50, f"{name}: lod95 ({lod95}) must exceed lod50 ({lod50})"
    # 95:50 ratio should be a sane test-sensitivity slope (well above 1).
    assert 1.0 < lod95 / lod50 < 50.0, f"{name}: implausible lod95:lod50 ratio"


def test_manual_option_is_distinct_and_default_exists():
    """The manual-entry sentinel must not collide with an assay name, and the
    assay used to pre-populate the manual fields must exist."""
    assert isinstance(MANUAL_LOD_OPTION, str) and MANUAL_LOD_OPTION
    assert MANUAL_LOD_OPTION not in NAT_ASSAYS
    assert MANUAL_DEFAULT_ASSAY in NAT_ASSAYS
