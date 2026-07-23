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

"""Tests for the canned-assay public API: the ``lods_for_assay`` resolver and
the ``assay=`` selection path on ``risk_days_bs``.

Most tests here are pool-free: the resolver is a pure lookup, and the
``risk_days_bs`` validation (assay XOR explicit LoDs) raises *before* any
ProcessPoolExecutor is created. Only the end-to-end parity test runs a real
bootstrap, so it carries the ``multiprocessing`` marker like the rest of the
suite (excluded in sandboxed runs via ``-m 'not multiprocessing'``).
"""

import math

import numpy as np
import pytest

import residualrisk as rr

# Small synthetic k posterior so the bootstrap is self-contained and fast.
_RNG = np.random.default_rng(0)
K_POSTERIOR = _RNG.exponential(scale=0.01, size=500)

# Valid risk_days_bs inputs (values mirror tests/test_residualrisk.py DEFAULTS),
# minus the LoD specification — supplied per-test as assay= or an explicit triplet.
BASE_KWARGS = dict(
    k=0.013,
    doubling_time=20.5 / 24,
    doubling_time_norm_sd=1.33 / 24,
    volume_transfused=20,
    volume_transfused_range=(15, 30),
    pool_size=16,
    retests=1,
    k_posterior_sample=K_POSTERIOR,
    n_bs=200,
    seed=42,
    threads=1,
    use_go=False,
)


# ---------------------------------------------------------------------------
# lods_for_assay / list_assays (pure lookups)
# ---------------------------------------------------------------------------


def test_lods_for_assay_fields():
    a = rr.lods_for_assay("ultrio_elite")
    assert isinstance(a, rr.AssayLoD)
    assert a.slug == "ultrio_elite"
    assert a.display_name == "Procleix Ultrio Elite (Panther)"
    assert a.lod50 == 3.1
    assert a.lod50_sd == 0.234
    assert a.lod95 == 10.4
    assert math.isclose(a.lod95_lod50_ratio, 10.4 / 3.1)
    assert a.cp_per_iu == 0.58
    assert a.iu_std


def test_lods_for_assay_unknown_slug_raises():
    with pytest.raises(ValueError, match="Unknown assay"):
        rr.lods_for_assay("not_a_real_assay")


def test_list_assays_matches_registry():
    listing = rr.list_assays()
    assert set(listing) == set(rr.NAT_ASSAYS)
    assert listing["ultrio_elite"] == "Procleix Ultrio Elite (Panther)"


# ---------------------------------------------------------------------------
# risk_days_bs assay= validation (raises before any pool is created)
# ---------------------------------------------------------------------------


def test_assay_and_explicit_lods_together_raises():
    with pytest.raises(ValueError, match="not both"):
        rr.risk_days_bs(**BASE_KWARGS, assay="ultrio_elite", lod50=3.1)


def test_neither_assay_nor_lods_raises():
    with pytest.raises(ValueError, match="lod50"):
        rr.risk_days_bs(**BASE_KWARGS)  # no assay and no explicit LoD triplet


def test_risk_days_bs_unknown_assay_raises():
    with pytest.raises(ValueError, match="Unknown assay"):
        rr.risk_days_bs(**BASE_KWARGS, assay="not_a_real_assay")


# ---------------------------------------------------------------------------
# End-to-end parity: assay= is a pure lookup of the explicit LoD triplet
# ---------------------------------------------------------------------------


@pytest.mark.multiprocessing
def test_assay_path_matches_explicit_lods():
    """Same seed + the LoDs resolved from the assay slug => identical bootstrap.

    Compares the point estimate (order-independent) and the sorted bootstrap
    array (robust to ProcessPoolExecutor completion-order)."""
    a = rr.lods_for_assay("ultrio_elite")
    pe_assay, _, _, bs_assay, _ = rr.risk_days_bs(**BASE_KWARGS, assay="ultrio_elite")
    pe_expl, _, _, bs_expl, _ = rr.risk_days_bs(
        **BASE_KWARGS,
        lod50=a.lod50,
        lod50_sd=a.lod50_sd,
        lod95_lod50_ratio=a.lod95_lod50_ratio,
    )
    assert pe_assay == pe_expl
    np.testing.assert_array_equal(
        np.sort(np.asarray(bs_assay)), np.sort(np.asarray(bs_expl))
    )
