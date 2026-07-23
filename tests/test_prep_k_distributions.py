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

"""
Tests for k-distribution flexibility in the PrEP model
(residualrisk.prep.risk_days_prep_bs) and the shared _sample_k helper.

Verifies that all four k-distribution paths (posterior sample, legacy gamma,
inverse gamma, lognormal mixture) work correctly in both the shared helper
and the PrEP bootstrap function.
"""

import numpy as np
import pytest

from residualrisk.core import _sample_k
from residualrisk.prep import risk_days_prep_bs


# ---------------------------------------------------------------------------
# Shared test fixtures
# ---------------------------------------------------------------------------

COMMON_PREP_KWARGS = dict(
    k=0.001,
    doubling_time=1.1,
    doubling_time_norm_sd=0.14,
    lod50=36.1,
    lod50_sd=4.6,
    lod95_lod50_ratio=3.01,
    volume_transfused=300,
    volume_transfused_range=(250, 350),
    pool_size=16,
    retests=0,
    seed=42,
    n_bs=20,
    threads=1,
    point_estimate="median",
)


# ---------------------------------------------------------------------------
# _sample_k helper tests
# ---------------------------------------------------------------------------


class TestSampleK:
    """Tests for the shared _sample_k dispatch helper."""

    def test_posterior_sample_sampler(self):
        np.random.seed(0)
        posterior = np.array([0.001, 0.002, 0.003, 0.004, 0.005])
        ks = _sample_k(100, seed=0, k_posterior_sample=posterior)
        assert len(ks) == 100
        assert all(k in posterior for k in ks)

    def test_gamma_legacy(self):
        np.random.seed(0)
        ks = _sample_k(100, seed=0, k_gamma_shape=2.0, k_gamma_scale=0.001)
        assert len(ks) == 100
        assert all(k > 0 for k in ks)

    def test_invgamma_alpha_beta(self):
        np.random.seed(0)
        ks = _sample_k(1000, seed=0, k_invgamma_alpha=2.0, k_invgamma_beta=0.002019)
        assert len(ks) == 1000
        assert all(k > 0 for k in ks)
        # Theoretical mean of InvGamma(2, 0.002019) = beta / (alpha - 1)
        expected_mean = 0.002019 / (2.0 - 1)
        assert abs(np.mean(ks) - expected_mean) / expected_mean < 0.15

    def test_invgamma_alpha_mode(self):
        np.random.seed(0)
        mode = 0.000673
        alpha = 2.0
        ks = _sample_k(1000, seed=0, k_invgamma_alpha=alpha, k_invgamma_mode=mode)
        assert len(ks) == 1000
        # beta = mode * (alpha + 1) = 0.000673 * 3 = 0.002019
        expected_beta = mode * (alpha + 1)
        expected_mean = expected_beta / (alpha - 1)
        assert abs(np.mean(ks) - expected_mean) / expected_mean < 0.15

    def test_invgamma_missing_beta_and_mode(self):
        np.random.seed(0)
        with pytest.raises(ValueError, match="exactly one of k_invgamma_beta"):
            _sample_k(100, seed=0, k_invgamma_alpha=2.0)

    def test_invgamma_beta_and_mode_together_raise(self):
        """Both beta and mode is now an error, matching the public sample_invgamma()
        contract — _sample_k used to silently prefer beta."""
        np.random.seed(0)
        with pytest.raises(ValueError, match="exactly one of k_invgamma_beta"):
            _sample_k(
                100,
                seed=0,
                k_invgamma_alpha=2.0,
                k_invgamma_beta=0.002,
                k_invgamma_mode=0.000673,
            )

    def test_lnmix_sampler(self):
        np.random.seed(0)
        ks = _sample_k(
            1000,
            seed=0,
            k_lnmix_w=0.90,
            k_lnmix_mu1=-7.2403,
            k_lnmix_sigma1=0.3241,
            k_lnmix_mu2=-3.7423,
            k_lnmix_sigma2=0.5258,
        )
        assert len(ks) == 1000
        assert all(k > 0 for k in ks)

    def test_lnmix_missing_params(self):
        np.random.seed(0)
        with pytest.raises(ValueError, match="lnmix parameters"):
            _sample_k(100, seed=0, k_lnmix_w=0.90, k_lnmix_mu1=-7.0)

    def test_no_distribution_raises_sampler(self):
        np.random.seed(0)
        with pytest.raises(ValueError, match="A k-distribution must be specified"):
            _sample_k(100, seed=0)

    def test_multiple_k_modes_raise(self):
        """Two k-distributions is an error, not a silent priority cascade.

        This test replaces `test_posterior_takes_priority`, which asserted the old
        behaviour: posterior + invgamma silently ran the *posterior*. That is exactly
        the failure mode being fixed — a stale k_posterior_sample in a reused config
        dict turned an InvGamma sensitivity analysis back into the posterior with no
        warning. Specifying both must now raise.
        """
        np.random.seed(0)
        posterior = np.array([0.001, 0.001, 0.001])
        with pytest.raises(ValueError, match="Exactly one k-distribution"):
            _sample_k(
                50,
                seed=0,
                k_posterior_sample=posterior,
                k_invgamma_alpha=2.0,
                k_invgamma_beta=0.002,
            )

    def test_partial_gamma_spec_raises(self):
        """A half-specified mode used to fall through to whichever *other* mode was
        populated: k_gamma_shape alone + k_invgamma_* silently ran the InvGamma."""
        np.random.seed(0)
        with pytest.raises(ValueError, match="Exactly one k-distribution"):
            _sample_k(
                50,
                seed=0,
                k_gamma_shape=2.0,  # no k_gamma_scale
                k_invgamma_alpha=2.0,
                k_invgamma_beta=0.002,
            )


# ---------------------------------------------------------------------------
# PrEP bootstrap with each k-distribution
# ---------------------------------------------------------------------------


class TestPrepBsKDistributions:
    """Integration tests: risk_days_prep_bs with each k-distribution path."""

    @pytest.mark.multiprocessing
    def test_posterior_sample_bootstrap(self):
        posterior = np.random.default_rng(99).exponential(0.001, size=200)
        rd_pe, rd_cri, rd_range, rdests, sim_df = risk_days_prep_bs(
            **COMMON_PREP_KWARGS,
            k_posterior_sample=posterior,
        )
        assert rd_pe is not None
        assert len(rdests) == COMMON_PREP_KWARGS["n_bs"]
        assert all(r >= 0 for r in rdests)
        assert sim_df is None  # return_sim_df=False by default

    @pytest.mark.multiprocessing
    def test_invgamma(self):
        rd_pe, rd_cri, rd_range, rdests, sim_df = risk_days_prep_bs(
            **COMMON_PREP_KWARGS,
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
        )
        assert rd_pe is not None
        assert len(rdests) == COMMON_PREP_KWARGS["n_bs"]
        assert all(r >= 0 for r in rdests)
        assert sim_df is None

    @pytest.mark.multiprocessing
    def test_invgamma_mode(self):
        rd_pe, rd_cri, rd_range, rdests, sim_df = risk_days_prep_bs(
            **COMMON_PREP_KWARGS,
            k_invgamma_alpha=2.0,
            k_invgamma_mode=0.000673,
        )
        assert rd_pe is not None
        assert len(rdests) == COMMON_PREP_KWARGS["n_bs"]
        assert sim_df is None

    @pytest.mark.multiprocessing
    def test_lnmix_bootstrap(self):
        rd_pe, rd_cri, rd_range, rdests, sim_df = risk_days_prep_bs(
            **COMMON_PREP_KWARGS,
            k_lnmix_w=0.90,
            k_lnmix_mu1=-7.2403,
            k_lnmix_sigma1=0.3241,
            k_lnmix_mu2=-3.7423,
            k_lnmix_sigma2=0.5258,
        )
        assert rd_pe is not None
        assert len(rdests) == COMMON_PREP_KWARGS["n_bs"]
        assert all(r >= 0 for r in rdests)
        assert sim_df is None

    @pytest.mark.multiprocessing
    def test_legacy_gamma(self):
        rd_pe, rd_cri, rd_range, rdests, sim_df = risk_days_prep_bs(
            **COMMON_PREP_KWARGS,
            k_gamma_shape=2.0,
            k_gamma_scale=0.001,
        )
        assert rd_pe is not None
        assert len(rdests) == COMMON_PREP_KWARGS["n_bs"]
        assert sim_df is None

    def test_no_distribution_raises_bootstrap(self):
        with pytest.raises(ValueError, match="A k-distribution must be specified"):
            risk_days_prep_bs(**COMMON_PREP_KWARGS)

    def test_multiple_k_modes_raise_bootstrap(self):
        """Rejected pre-dispatch, so it raises in-process (not inside a worker) and
        raises identically whichever backend was requested."""
        with pytest.raises(ValueError, match="Exactly one k-distribution"):
            risk_days_prep_bs(
                **COMMON_PREP_KWARGS,
                k_posterior_sample=np.array([0.001, 0.002]),
                k_invgamma_alpha=2.0,
                k_invgamma_beta=0.002,
            )

    @pytest.mark.multiprocessing
    def test_return_sim_df(self):
        """When return_sim_df=True, a Polars DataFrame with expected columns is returned."""
        import polars as pl

        rd_pe, rd_cri, rd_range, rdests, sim_df = risk_days_prep_bs(
            **COMMON_PREP_KWARGS,
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
            return_sim_df=True,
        )
        assert isinstance(sim_df, pl.DataFrame)
        assert len(sim_df) == COMMON_PREP_KWARGS["n_bs"]
        assert "iwp" in sim_df.columns
        assert "set_point" in sim_df.columns
        assert "eclipse" in sim_df.columns
        assert "k" in sim_df.columns
