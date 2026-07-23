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

"""Tests for ``total_residual_risk_rd`` — the additive-total residual-risk
credible interval (baseline + oral-PrEP + injectable-PrEP).

Most tests use synthetic, pre-aligned IWP arrays and do NOT spawn the
ProcessPoolExecutor bootstrap, so they are sandbox-safe (no ``multiprocessing``
marker needed). One integration test exercises the Go backend to guard the
shared-parameter alignment that makes the summed-quantile interval a valid joint
CrI; it is skipped when the Go binary is unavailable.

Run from the repo root:
    pytest tests/test_total_residual_risk.py
"""

import numpy as np
import pytest

import residualrisk as rrpkg
from residualrisk import core as rr
from residualrisk import risk_days_bs, risk_days_prep_bs, total_residual_risk_rd

SEED = 4242


def _iwp(n, loc=8.0, scale=1.5, seed=1):
    """A positive synthetic IWP/RDE bootstrap sample of length ``n``."""
    rng = np.random.default_rng(seed)
    return np.abs(rng.normal(loc, scale, n)) + 0.05


# --------------------------------------------------------------------------- #
# Point estimate
# --------------------------------------------------------------------------- #


def test_pe_is_additive_sum_of_component_pes():
    n = 5000
    iwp = _iwp(n)
    comps = [
        (8.0, iwp, 0.0003, 0.00005),
        (12.0, iwp, 0.0010, 0.00020),
        (6.0, iwp, 0.0007, 0.00010),
    ]
    pe, cri, onein_pe, onein_cri = total_residual_risk_rd(comps, per=1e6, seed=SEED)
    expected_prob = sum(inc * iwp_pe / 365.25 for (iwp_pe, _, inc, _) in comps)
    assert pe == pytest.approx(expected_prob * 1e6, rel=1e-12)
    assert onein_pe == pytest.approx(1.0 / expected_prob, rel=1e-12)


# --------------------------------------------------------------------------- #
# Credible-interval shape and consistency
# --------------------------------------------------------------------------- #


def test_cri_ordered_positive_and_brackets_pe():
    n = 8000
    comps = [
        (8.0, _iwp(n, loc=8.0, seed=1), 0.0005, 0.00010),
        (10.0, _iwp(n, loc=10.0, seed=2), 0.0008, 0.00015),
    ]
    pe, cri, onein_pe, onein_cri = total_residual_risk_rd(comps, per=1e6, seed=SEED)
    lo, hi = cri
    assert 0 < lo < hi
    assert lo <= pe <= hi
    # 1-in-N CrI is the inverse: cri[0] is the smaller-N (higher-risk) end.
    assert 0 < onein_cri[0] < onein_cri[1]
    assert onein_cri[0] <= onein_pe <= onein_cri[1]
    # The two representations are consistent.
    assert onein_pe == pytest.approx(1e6 / pe, rel=1e-9)


def test_onein_cri_is_inverse_of_perm_cri():
    n = 6000
    comps = [(8.0, _iwp(n), 0.0006, 0.00012), (9.0, _iwp(n, seed=3), 0.0009, 0.00018)]
    pe, cri, onein_pe, onein_cri = total_residual_risk_rd(comps, per=1e6, seed=SEED)
    # per-million quantile (1 - a/2) maps to the smaller-N end (onein_cri[0]).
    # The two are inverses at the order-statistic level; they differ only by
    # quantile linear-interpolation done in 1/x- vs x-space, hence the loose tol.
    assert onein_cri[0] == pytest.approx(1e6 / cri[1], rel=1e-6)
    assert onein_cri[1] == pytest.approx(1e6 / cri[0], rel=1e-6)


# --------------------------------------------------------------------------- #
# Reproducibility
# --------------------------------------------------------------------------- #


def test_reproducible_same_seed():
    n = 3000
    comps = [(8.0, _iwp(n), 0.0005, 0.0001), (10.0, _iwp(n, seed=2), 0.0008, 0.00015)]
    a = total_residual_risk_rd(comps, seed=SEED)
    b = total_residual_risk_rd(comps, seed=SEED)
    assert a[0] == b[0]
    assert a[1] == b[1]
    assert a[2] == b[2]
    assert a[3] == b[3]


def test_different_seed_changes_cri_not_pe():
    n = 3000
    comps = [(8.0, _iwp(n), 0.0005, 0.0001)]
    a = total_residual_risk_rd(comps, seed=1)
    b = total_residual_risk_rd(comps, seed=2)
    assert a[0] == pytest.approx(b[0])  # PE uses point incidence -> seed-independent
    assert a[1] != b[1]  # CrI uses sampled incidence -> seed-dependent


# --------------------------------------------------------------------------- #
# Equivalence with the validated residual_risk_rd
# --------------------------------------------------------------------------- #


def test_single_component_matches_residual_risk_rd():
    """A one-component total must reproduce residual_risk_rd bit-for-bit:
    component 0 uses seed+0 == seed and the same _sample_positive_normal draws."""
    n = 6000
    iwp = _iwp(n)
    inc, sd, iwp_pe = 0.0007, 0.00012, 9.0
    pe, cri, onein_pe, onein_cri = total_residual_risk_rd(
        [(iwp_pe, iwp, inc, sd)], per=1e6, seed=SEED
    )
    r_pe, r_cri, _ = rr.residual_risk_rd(
        iwp_pe, iwp, inc, sd, per=1e6, seed=SEED, one_in_x=False
    )
    o_pe, o_cri, _ = rr.residual_risk_rd(
        iwp_pe, iwp, inc, sd, per=None, seed=SEED, one_in_x=True
    )
    assert pe == pytest.approx(r_pe, rel=1e-12)
    assert onein_pe == pytest.approx(o_pe, rel=1e-12)
    np.testing.assert_allclose(cri, r_cri, rtol=1e-9)
    np.testing.assert_allclose(onein_cri, o_cri, rtol=1e-9)


# --------------------------------------------------------------------------- #
# Dependence structure: independent incidence, shared (correlated) IWP
# --------------------------------------------------------------------------- #


def test_incidence_drawn_independently_across_components():
    """Constant IWP removes IWP variance, so the total variance is driven purely
    by incidence. Two identical components drawn from independent seeds give
    Var(T) ~= 2x a single component's (independent sum) — not ~4x, which a
    shared/comonotonic incidence would give."""
    n = 300_000
    c = 10.0
    iwp = np.full(n, c)
    inc, sd = 0.001, 0.0002
    *_, samp1 = total_residual_risk_rd(
        [(c, iwp, inc, sd)], per=1e6, seed=SEED, return_samps=True
    )
    *_, samp2 = total_residual_risk_rd(
        [(c, iwp, inc, sd), (c, iwp, inc, sd)], per=1e6, seed=SEED, return_samps=True
    )
    assert np.var(samp2) == pytest.approx(2 * np.var(samp1), rel=0.05)


def test_shared_iwp_induces_positive_correlation():
    """With VARYING shared IWP, two identical components are positively
    correlated through the shared IWP, so Var(T) must EXCEED the
    independent-sum value (2x a single component's variance). This is the
    property that makes the joint CrI wider than an independence approximation."""
    n = 300_000
    iwp = _iwp(n, loc=8.0, scale=3.0, seed=7)  # substantial IWP variance
    inc, sd = 0.001, 0.0002
    *_, samp1 = total_residual_risk_rd(
        [(8.0, iwp, inc, sd)], per=1e6, seed=SEED, return_samps=True
    )
    *_, samp2 = total_residual_risk_rd(
        [(8.0, iwp, inc, sd), (8.0, iwp, inc, sd)],
        per=1e6,
        seed=SEED,
        return_samps=True,
    )
    assert np.var(samp2) > 2.2 * np.var(samp1)


# --------------------------------------------------------------------------- #
# Input validation
# --------------------------------------------------------------------------- #


def test_empty_components_raises():
    with pytest.raises(ValueError):
        total_residual_risk_rd([])


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        total_residual_risk_rd([
            (8.0, np.ones(100), 0.001, 0.0001),
            (8.0, np.ones(99), 0.001, 0.0001),
        ])


def test_nonpositive_incidence_or_iwp_raises():
    with pytest.raises(ValueError):
        total_residual_risk_rd([(8.0, np.ones(100), -0.001, 0.0001)])
    with pytest.raises(ValueError):
        total_residual_risk_rd([(-1.0, np.ones(100), 0.001, 0.0001)])


def test_empty_iwp_bs_raises():
    """residual_risk_rd with an empty iwp_bs (or all-nonpositive products) has no
    samples to summarise — raise cleanly rather than surfacing an opaque IndexError
    from np.quantile([])."""
    with pytest.raises(ValueError):
        rr.residual_risk_rd(1.0, [], 0.001, 0.0001)


# --------------------------------------------------------------------------- #
# Integration: the shared-parameter alignment the CrI depends on
# --------------------------------------------------------------------------- #


@pytest.mark.skipif(rrpkg.find_go_binary() is None, reason="Go binary not available")
def test_go_aligns_shared_params_across_components():
    """The total CrI is valid only because baseline/oral/inj share their
    per-iteration k / doubling_time / lod50 / volume draws (Go backend, common
    seed). Guard that invariant — if the Go draw order ever changes so the shared
    params diverge, this test fails and the total CrI silently degrades to an
    independence approximation."""
    shared = dict(
        k=0.000673,
        doubling_time=20.5 / 24,
        doubling_time_norm_sd=1.33 / 24,
        lod50=2.73,
        lod50_sd=0.193,
        lod95_lod50_ratio=12.33 / 2.73,
        volume_transfused=20,
        volume_transfused_range=(15, 30),
        pool_size=16,
        retests=1,
        k_invgamma_alpha=2.0,
        k_invgamma_beta=0.002019,
        n_bs=1500,
        seed=987,
        threads=2,
        return_sim_df=True,
        use_go=True,
    )
    *_, base = risk_days_bs(**shared)
    *_, oral = risk_days_prep_bs(set_point=336, **shared)
    *_, inj = risk_days_prep_bs(set_point=600, **shared)
    for col in ["k", "doubling_time", "lod50", "volume_transfused"]:
        np.testing.assert_array_equal(base[col].to_numpy(), oral[col].to_numpy())
        np.testing.assert_array_equal(base[col].to_numpy(), inj[col].to_numpy())
