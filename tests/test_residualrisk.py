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
Unit and integration tests for residualrisk.py core calculation functions.

Bootstrap simulation tests (test_risk_days_bs_*) run both the Python and Go
implementations and check that results are self-consistent and agree with each
other within a reasonable tolerance.  Because Python and Go use independent
RNGs, exact numerical equality is not expected; the tolerance is set at 15%
on distributional summaries (median, CrI) for n_bs=500.

Run from the app/ directory:
    pytest tests/
"""

import math

import numpy as np
import pytest

from residualrisk import core as rr

# ---------------------------------------------------------------------------
# Shared parameter fixtures
# ---------------------------------------------------------------------------

# Default parameter values matching app.py UI defaults
DEFAULTS = dict(
    C0=0.00025,
    doubling_time=20.5 / 24,  # hours → days
    doubling_time_norm_sd=1.33 / 24,
    volume_transfused=20,  # mL
    volume_transfused_range=(15, 30),
    copies_per_virion=2,
    pool_size=16,
    retests=1,
    lod50=2.73,
    lod50_sd=0.193,
    lod95=12.33,
    z=1.6449,
    k=0.013,
)
DEFAULTS["lod95_lod50_ratio"] = DEFAULTS["lod95"] / DEFAULTS["lod50"]

# A small synthetic k posterior: 1000 draws centred around realistic values
_RNG = np.random.default_rng(seed=0)
K_POSTERIOR = np.concatenate(
    [
        _RNG.exponential(scale=0.003, size=333),
        _RNG.exponential(scale=0.008, size=334),
        _RNG.exponential(scale=0.020, size=333),
    ]
)

# Bootstrap settings used in all bootstrap tests
BS_KWARGS = dict(
    k=DEFAULTS["k"],
    doubling_time=DEFAULTS["doubling_time"],
    doubling_time_norm_sd=DEFAULTS["doubling_time_norm_sd"],
    lod50=DEFAULTS["lod50"],
    lod50_sd=DEFAULTS["lod50_sd"],
    lod95_lod50_ratio=DEFAULTS["lod95_lod50_ratio"],
    volume_transfused=DEFAULTS["volume_transfused"],
    volume_transfused_range=DEFAULTS["volume_transfused_range"],
    pool_size=DEFAULTS["pool_size"],
    retests=DEFAULTS["retests"],
    k_posterior_sample=K_POSTERIOR,
    n_bs=500,
    seed=42,
    threads=2,
)


# ---------------------------------------------------------------------------
# _concentration
# ---------------------------------------------------------------------------


class TestConcentration:
    def test_at_t0_equals_C0(self):
        assert rr._concentration(
            DEFAULTS["C0"], DEFAULTS["doubling_time"], 0
        ) == pytest.approx(DEFAULTS["C0"])

    def test_at_one_doubling_time_doubles(self):
        dt = DEFAULTS["doubling_time"]
        assert rr._concentration(DEFAULTS["C0"], dt, dt) == pytest.approx(
            2 * DEFAULTS["C0"]
        )

    def test_at_two_doubling_times_quadruples(self):
        dt = DEFAULTS["doubling_time"]
        assert rr._concentration(DEFAULTS["C0"], dt, 2 * dt) == pytest.approx(
            4 * DEFAULTS["C0"]
        )

    def test_increases_with_t(self):
        dt = DEFAULTS["doubling_time"]
        C0 = DEFAULTS["C0"]
        values = [rr._concentration(C0, dt, t) for t in [0, 1, 5, 10, 20]]
        assert all(a < b for a, b in zip(values, values[1:]))


# ---------------------------------------------------------------------------
# _prob_infectious_copies
# ---------------------------------------------------------------------------


class TestProbInfectiousCopies:
    def test_zero_copies_gives_zero(self):
        assert rr._prob_infectious_copies(0, 0.013) == pytest.approx(0.0)

    def test_standard_single_hit_formula(self):
        # 1 - exp(-k * n) for n=100, k=0.01 → 1 - exp(-1) ≈ 0.6321
        expected = 1.0 - math.exp(-0.01 * 100)
        assert rr._prob_infectious_copies(100, 0.01) == pytest.approx(
            expected, rel=1e-9
        )

    def test_large_copies_approaches_one(self):
        assert rr._prob_infectious_copies(1_000_000, 0.013) == pytest.approx(
            1.0, abs=1e-6
        )

    def test_result_bounded_between_zero_and_one(self):
        for n in [0, 1, 10, 100, 1000]:
            p = rr._prob_infectious_copies(n, 0.013)
            assert 0.0 <= p <= 1.0

    def test_increases_with_copies(self):
        ks = [rr._prob_infectious_copies(n, 0.013) for n in [0, 1, 10, 100, 1000]]
        assert all(a <= b for a, b in zip(ks, ks[1:]))

    def test_increases_with_k(self):
        ps = [rr._prob_infectious_copies(50, k) for k in [0.001, 0.01, 0.1]]
        assert all(a < b for a, b in zip(ps, ps[1:]))


# ---------------------------------------------------------------------------
# _prob_infectious_copies_wc (worst-case, threshold model)
# ---------------------------------------------------------------------------


class TestProbInfectiousCopiesWC:
    def test_below_threshold_is_zero(self):
        assert rr._prob_infectious_copies_wc(0) == 0.0
        assert rr._prob_infectious_copies_wc(1.9) == 0.0

    def test_at_threshold_is_one(self):
        assert rr._prob_infectious_copies_wc(2) == 1.0

    def test_above_threshold_is_one(self):
        assert rr._prob_infectious_copies_wc(1000) == 1.0


# ---------------------------------------------------------------------------
# _prob_pos_init
# ---------------------------------------------------------------------------


class TestProbPosInit:
    def test_at_pool_lod50_gives_half(self):
        # When C = pool_size * lod50, the pooled concentration equals lod50,
        # so the initial detection probability should be 0.5.
        C = DEFAULTS["pool_size"] * DEFAULTS["lod50"]
        p = rr._prob_pos_init(
            C,
            DEFAULTS["doubling_time"],
            DEFAULTS["pool_size"],
            DEFAULTS["lod50"],
            DEFAULTS["lod95_lod50_ratio"],
            DEFAULTS["z"],
        )
        assert p == pytest.approx(0.5, abs=1e-6)

    def test_increases_with_concentration(self):
        lod50 = DEFAULTS["lod50"]
        concentrations = [lod50, 10 * lod50, 100 * lod50, 1000 * lod50]
        probs = [
            rr._prob_pos_init(
                C,
                DEFAULTS["doubling_time"],
                DEFAULTS["pool_size"],
                lod50,
                DEFAULTS["lod95_lod50_ratio"],
                DEFAULTS["z"],
            )
            for C in concentrations
        ]
        assert all(a < b for a, b in zip(probs, probs[1:]))

    def test_invalid_pool_size_raises(self):
        with pytest.raises(Exception):
            rr._prob_pos_init(
                10.0,
                DEFAULTS["doubling_time"],
                0,
                DEFAULTS["lod50"],
                DEFAULTS["lod95_lod50_ratio"],
                DEFAULTS["z"],
            )

    def test_non_integer_pool_size_raises(self):
        with pytest.raises(Exception):
            rr._prob_pos_init(
                10.0,
                DEFAULTS["doubling_time"],
                1.5,
                DEFAULTS["lod50"],
                DEFAULTS["lod95_lod50_ratio"],
                DEFAULTS["z"],
            )


# ---------------------------------------------------------------------------
# _prob_neg_retest
# ---------------------------------------------------------------------------


class TestProbNegRetest:
    def test_zero_retests_returns_zero(self):
        result = rr._prob_neg_retest(
            100.0,
            DEFAULTS["doubling_time"],
            DEFAULTS["pool_size"],
            DEFAULTS["lod50"],
            DEFAULTS["lod95_lod50_ratio"],
            0,
            DEFAULTS["z"],
        )
        assert result == 0

    def test_high_concentration_gives_near_zero(self):
        # At very high viral load the retest should almost certainly detect,
        # so the probability of a negative retest should be near zero.
        result = rr._prob_neg_retest(
            1e8,
            DEFAULTS["doubling_time"],
            DEFAULTS["pool_size"],
            DEFAULTS["lod50"],
            DEFAULTS["lod95_lod50_ratio"],
            1,
            DEFAULTS["z"],
        )
        assert result == pytest.approx(0.0, abs=1e-6)

    def test_decreases_with_concentration(self):
        # Higher concentration → more likely to detect → lower prob of neg retest
        concentrations = [
            DEFAULTS["lod50"],
            10 * DEFAULTS["lod50"],
            1000 * DEFAULTS["lod50"],
        ]
        probs = [
            rr._prob_neg_retest(
                C,
                DEFAULTS["doubling_time"],
                DEFAULTS["pool_size"],
                DEFAULTS["lod50"],
                DEFAULTS["lod95_lod50_ratio"],
                1,
                DEFAULTS["z"],
            )
            for C in concentrations
        ]
        assert all(a > b for a, b in zip(probs, probs[1:]))

    def test_invalid_retests_raises(self):
        with pytest.raises(Exception):
            rr._prob_neg_retest(
                10.0,
                DEFAULTS["doubling_time"],
                DEFAULTS["pool_size"],
                DEFAULTS["lod50"],
                DEFAULTS["lod95_lod50_ratio"],
                -1,
                DEFAULTS["z"],
            )


# ---------------------------------------------------------------------------
# _prob_nondetection
# ---------------------------------------------------------------------------


class TestProbNondetection:
    def _call(self, t):
        return rr._prob_nondetection(
            t,
            DEFAULTS["copies_per_virion"],
            DEFAULTS["C0"],
            DEFAULTS["doubling_time"],
            DEFAULTS["pool_size"],
            DEFAULTS["lod50"],
            DEFAULTS["lod95_lod50_ratio"],
            DEFAULTS["retests"],
            DEFAULTS["z"],
        )

    def test_very_early_time_is_one(self):
        # Long before detectable viremia, probability of non-detection is ~1.
        assert self._call(-50) == pytest.approx(1.0, abs=1e-6)

    def test_very_late_time_is_zero(self):
        # Far into viremia, viral load is so high that non-detection is ~0.
        assert self._call(100) == pytest.approx(0.0, abs=1e-6)

    def test_decreases_over_time(self):
        times = [-20, -10, 0, 10, 20]
        probs = [self._call(t) for t in times]
        assert all(a >= b for a, b in zip(probs, probs[1:]))


# ---------------------------------------------------------------------------
# _risk_days  (deterministic point-estimate integral)
# ---------------------------------------------------------------------------


class TestRiskDays:
    def _call(self, **overrides):
        params = dict(
            copies_per_virion=DEFAULTS["copies_per_virion"],
            C0=DEFAULTS["C0"],
            doubling_time=DEFAULTS["doubling_time"],
            volume_transfused=DEFAULTS["volume_transfused"],
            k=DEFAULTS["k"],
            pool_size=DEFAULTS["pool_size"],
            lod50=DEFAULTS["lod50"],
            lod95_lod50_ratio=DEFAULTS["lod95_lod50_ratio"],
            retests=DEFAULTS["retests"],
            z=DEFAULTS["z"],
        )
        params.update(overrides)
        return rr._risk_days(**params)

    def test_golden_value_default_params(self):
        # Reference value computed with fixed default parameters.
        assert self._call() == pytest.approx(3.7207, rel=1e-3)

    def test_idnat_golden_value(self):
        # Individual-donation NAT (pool_size=1, retests=0) reduces risk days.
        assert self._call(pool_size=1, retests=0) == pytest.approx(0.9001, rel=1e-3)

    def test_higher_k_gives_more_risk_days(self):
        rd_low = self._call(k=0.005)
        rd_high = self._call(k=0.05)
        assert rd_high > rd_low

    def test_larger_volume_gives_more_risk_days(self):
        rd_small = self._call(volume_transfused=5)
        rd_large = self._call(volume_transfused=100)
        assert rd_large > rd_small

    def test_idnat_less_than_minipool(self):
        # Individual-donation NAT detects earlier than a 16-pool.
        rd_minipool = self._call(pool_size=16, retests=1)
        rd_idnat = self._call(pool_size=1, retests=0)
        assert rd_idnat < rd_minipool

    def test_result_is_positive(self):
        assert self._call() > 0


# ---------------------------------------------------------------------------
# risk_days_bs  (bootstrap simulations — Python and Go)
# ---------------------------------------------------------------------------


def _assert_bs_result_sane(result, n_bs):
    """Shared structural assertions for any bootstrap result."""
    rd_pe, rd_cri, rd_range, rdests, _ = result
    assert rd_pe > 0, "Point estimate must be positive"
    assert len(rdests) == n_bs, "Simulation count must equal n_bs"
    assert all(r > 0 for r in rdests), "All simulation values must be positive"
    assert rd_cri[0] < rd_cri[1], "CrI lower bound must be less than upper bound"
    assert rd_range[0] <= rd_range[1], "Range lower bound must not exceed upper bound"
    assert rd_range[0] <= rd_cri[0], "Range min must be ≤ CrI lower bound"
    assert rd_cri[1] <= rd_range[1], "CrI upper bound must be ≤ range max"


class TestRiskDaysBsPython:
    def test_returns_correct_structure(self):
        result = rr.risk_days_bs(**BS_KWARGS, use_go=False)
        assert len(result) == 5

    def test_sanity_checks(self):
        result = rr.risk_days_bs(**BS_KWARGS, use_go=False)
        _assert_bs_result_sane(result, BS_KWARGS["n_bs"])

    def test_point_estimate_matches_risk_days(self):
        # With point_estimate="primary parameters" the pe should equal _risk_days
        # evaluated at the primary (non-bootstrapped) parameter values.
        result = rr.risk_days_bs(
            **BS_KWARGS, use_go=False, point_estimate="primary parameters"
        )
        expected_pe = rr._risk_days(
            copies_per_virion=DEFAULTS["copies_per_virion"],
            C0=DEFAULTS["C0"],
            doubling_time=DEFAULTS["doubling_time"],
            volume_transfused=DEFAULTS["volume_transfused"],
            k=DEFAULTS["k"],
            pool_size=DEFAULTS["pool_size"],
            lod50=DEFAULTS["lod50"],
            lod95_lod50_ratio=DEFAULTS["lod95_lod50_ratio"],
            retests=DEFAULTS["retests"],
        )
        assert result[0] == pytest.approx(expected_pe, rel=1e-6)

    def test_reproducible_with_same_seed(self):
        r1 = rr.risk_days_bs(**BS_KWARGS, use_go=False)
        r2 = rr.risk_days_bs(**BS_KWARGS, use_go=False)
        assert r1[0] == r2[0]
        # as_completed returns futures in non-deterministic order, so compare
        # sorted values rather than the raw list.
        assert sorted(r1[3]) == sorted(r2[3])

    def test_different_seeds_give_different_results(self):
        r1 = rr.risk_days_bs(**{**BS_KWARGS, "seed": 42}, use_go=False)
        r2 = rr.risk_days_bs(**{**BS_KWARGS, "seed": 99999}, use_go=False)
        assert sorted(r1[3]) != sorted(r2[3])

    def test_returns_sim_df_when_requested(self):
        result = rr.risk_days_bs(**BS_KWARGS, use_go=False, return_sim_df=True)
        assert len(result) == 5
        sim_df = result[4]
        assert len(sim_df) == BS_KWARGS["n_bs"]
        assert "iwp" in sim_df.columns
        assert "k" in sim_df.columns

    def test_invalid_n_bs_raises(self):
        with pytest.raises(ValueError):
            rr.risk_days_bs(**{**BS_KWARGS, "n_bs": 0}, use_go=False)

    def test_missing_k_distribution_raises(self):
        kwargs = {k: v for k, v in BS_KWARGS.items() if k != "k_posterior_sample"}
        with pytest.raises(ValueError):
            rr.risk_days_bs(**kwargs, use_go=False)


class TestRiskDaysBsGo:
    def test_returns_correct_structure(self):
        result = rr.risk_days_bs(**BS_KWARGS, use_go=True)
        assert len(result) == 5

    def test_sanity_checks(self):
        result = rr.risk_days_bs(**BS_KWARGS, use_go=True)
        _assert_bs_result_sane(result, BS_KWARGS["n_bs"])

    def test_point_estimate_matches_risk_days(self):
        result = rr.risk_days_bs(
            **BS_KWARGS, use_go=True, point_estimate="primary parameters"
        )
        expected_pe = rr._risk_days(
            copies_per_virion=DEFAULTS["copies_per_virion"],
            C0=DEFAULTS["C0"],
            doubling_time=DEFAULTS["doubling_time"],
            volume_transfused=DEFAULTS["volume_transfused"],
            k=DEFAULTS["k"],
            pool_size=DEFAULTS["pool_size"],
            lod50=DEFAULTS["lod50"],
            lod95_lod50_ratio=DEFAULTS["lod95_lod50_ratio"],
            retests=DEFAULTS["retests"],
        )
        assert result[0] == pytest.approx(expected_pe, rel=1e-6)

    def test_returns_sim_df_when_requested(self):
        result = rr.risk_days_bs(**BS_KWARGS, use_go=True, return_sim_df=True)
        assert len(result) == 5
        sim_df = result[4]
        assert len(sim_df) == BS_KWARGS["n_bs"]
        assert "iwp" in sim_df.columns


class TestPythonGoAgreement:
    """
    Python and Go use independent RNGs so results will not be identical.
    We check that distributional summaries agree within 15%, which is a
    reasonable tolerance for n_bs=500.
    """

    def test_point_estimates_agree(self):
        # Both implementations compute the pe via the same deterministic call
        # to _risk_days, so the point estimates should be identical.
        py = rr.risk_days_bs(
            **BS_KWARGS, use_go=False, point_estimate="primary parameters"
        )
        go = rr.risk_days_bs(
            **BS_KWARGS, use_go=True, point_estimate="primary parameters"
        )
        assert py[0] == pytest.approx(go[0], rel=1e-6)

    def test_simulation_medians_agree_within_tolerance(self):
        py = rr.risk_days_bs(**BS_KWARGS, use_go=False)
        go = rr.risk_days_bs(**BS_KWARGS, use_go=True)
        py_median = np.median(py[3])
        go_median = np.median(go[3])
        assert py_median == pytest.approx(go_median, rel=0.15)

    def test_cri_upper_bounds_agree_within_tolerance(self):
        py = rr.risk_days_bs(**BS_KWARGS, use_go=False)
        go = rr.risk_days_bs(**BS_KWARGS, use_go=True)
        assert py[1][1] == pytest.approx(go[1][1], rel=0.20)


# ---------------------------------------------------------------------------
# mode_kde
# ---------------------------------------------------------------------------


class TestModeKde:
    """KDE-on-log-scale mode estimation for posterior distributions.

    Tests that the mode of the human and animal k posteriors matches
    the documented values.  These are pre-computed Bayesian posteriors
    stored in static/*.parquet that must not change; the expected modes
    serve as canary values to detect unintended breakage.
    """

    @staticmethod
    def _load_parquet(name: str) -> "np.ndarray":
        from pathlib import Path

        import pandas as pd

        static = Path(__file__).resolve().parent.parent / "static"
        return pd.read_parquet(static / name).iloc[:, 0].values

    def test_human_posterior_mode(self):
        """Human posterior k mode should be ≈ 0.000672 (document value 0.000673)."""
        k_human = self._load_parquet("k_param_human.parquet")
        mode = rr._kde_mode_log(k_human)
        assert mode == pytest.approx(0.000672, abs=1e-6)

    def test_animal_posterior_mode(self):
        """Animal posterior k mode should be ≈ 0.02086."""
        k_animal = self._load_parquet("k_param_animal.parquet")
        mode = rr._kde_mode_log(k_animal)
        assert mode == pytest.approx(0.020862, abs=1e-5)

    def test_mode_kde_public_wrapper(self):
        """The public mode_kde wrapper should match the private implementation."""
        import numpy as np

        rng = np.random.default_rng(seed=1)
        data = rng.lognormal(mean=-7.0, sigma=0.3, size=5000)
        assert rr.mode_kde(data, n_grid=5_000) == rr._kde_mode_log(data, n_grid=5_000)

    def test_all_zeros_or_negatives_raises(self):
        """Zero or negative values should raise ValueError (log undefined)."""
        import numpy as np

        with pytest.raises(ValueError, match="positive"):
            rr._kde_mode_log(np.array([0.0, 0.5, 1.0]))
        with pytest.raises(ValueError, match="positive"):
            rr._kde_mode_log(np.array([-0.1, 0.5, 1.0]))

    def test_result_is_positive_and_finite(self):
        """Mode estimate should always be a positive finite number."""
        import numpy as np

        rng = np.random.default_rng(seed=2)
        data = rng.lognormal(mean=-7.0, sigma=0.4, size=1000)
        mode = rr._kde_mode_log(data)
        assert np.isfinite(mode)
        assert mode > 0


# ---------------------------------------------------------------------------
# sample_invgamma
# ---------------------------------------------------------------------------


class TestSampleInvgamma:
    """Unit tests for sample_invgamma()."""

    def test_beta_parameterisation(self):
        """Direct (alpha, beta) produces samples with expected median."""
        samples = rr.sample_invgamma(50_000, alpha=2.0, beta=0.002019, seed=1)
        # InvGamma(2, 0.002019) theoretical median ≈ 0.001203
        assert np.median(samples) == pytest.approx(0.001203, rel=0.03)
        assert len(samples) == 50_000
        assert np.all(samples > 0)

    def test_mode_parameterisation(self):
        """(alpha, mode) auto-calculates beta = mode * (alpha + 1)."""
        beta_samples = rr.sample_invgamma(
            50_000, alpha=3.0, beta=0.002692, seed=1
        )
        mode_samples = rr.sample_invgamma(
            50_000, alpha=3.0, mode=0.000673, seed=1
        )
        # beta = 0.000673 * (3 + 1) = 0.002692 — same distribution
        assert np.array_equal(beta_samples, mode_samples)

    def test_mode_and_beta_mutually_exclusive(self):
        """Providing both beta and mode raises ValueError."""
        with pytest.raises(ValueError, match="not both"):
            rr.sample_invgamma(100, alpha=2.0, beta=0.002, mode=0.000673)

    def test_neither_beta_nor_mode_raises(self):
        """Providing neither beta nor mode raises ValueError."""
        with pytest.raises(ValueError, match="Exactly one"):
            rr.sample_invgamma(100, alpha=2.0)

    def test_reproducible_with_seed(self):
        """Same seed produces identical samples."""
        s1 = rr.sample_invgamma(500, alpha=2.0, mode=0.000673, seed=42)
        s2 = rr.sample_invgamma(500, alpha=2.0, mode=0.000673, seed=42)
        assert np.array_equal(s1, s2)

    def test_different_seeds_differ(self):
        """Different seeds produce different samples."""
        s1 = rr.sample_invgamma(500, alpha=2.0, mode=0.000673, seed=1)
        s2 = rr.sample_invgamma(500, alpha=2.0, mode=0.000673, seed=2)
        assert not np.array_equal(s1, s2)


# ---------------------------------------------------------------------------
# InvGamma IWP agreement
# ---------------------------------------------------------------------------


class TestInvgammaIwpAgreement:
    """The IWP point estimate from InvGamma(k_mode=human_mode) should
    closely match the IWP point estimate from the human posterior when
    both use mode as the k point estimate and mode as the IWP point
    estimate."""

    @staticmethod
    def _load_human_posterior() -> "np.ndarray":
        from pathlib import Path

        import pandas as pd

        static = Path(__file__).resolve().parent.parent / "static"
        return pd.read_parquet(static / "k_param_human.parquet").iloc[
            :, 0
        ].values

    @staticmethod
    def _common_params() -> dict:
        return dict(
            doubling_time=20.5 / 24,
            doubling_time_norm_sd=1.33 / 24,
            lod50=2.73,
            lod50_sd=0.193,
            lod95_lod50_ratio=12.33 / 2.73,
            volume_transfused=20,
            volume_transfused_range=(15, 30),
            pool_size=16,
            retests=1,
            n_bs=500,
            seed=42,
            threads=2,
            point_estimate="mode",
            use_go=False,
        )

    def test_iwp_mode_agrees_with_human_posterior(self):
        """IWP mode from InvGamma ≈ IWP mode from human posterior."""
        k_human = self._load_human_posterior()
        k_mode_human = rr.mode_kde(k_human, n_grid=5_000, cap=50_000)
        params = self._common_params()

        # Human posterior
        r_human = rr.risk_days_bs(
            k=k_mode_human,
            k_posterior_sample=k_human,
            **params,
        )
        iwp_human = r_human[0]

        # Inverse Gamma with same mode
        r_ig = rr.risk_days_bs(
            k=0.000673,
            k_invgamma_alpha=2.0,
            k_invgamma_mode=0.000673,
            **params,
        )
        iwp_ig = r_ig[0]

        # Both are positive risk-days estimates
        assert iwp_human > 0
        assert iwp_ig > 0

        # Point estimates should agree within 25% (n_bs=500 is noisy)
        assert iwp_human == pytest.approx(iwp_ig, rel=0.25)

    def test_iwp_mode_agrees_with_human_posterior_go(self):
        """Same as above, but using the Go implementation."""
        k_human = self._load_human_posterior()
        k_mode_human = rr.mode_kde(k_human, n_grid=5_000, cap=50_000)
        params = {**self._common_params(), "use_go": True}

        # Human posterior
        r_human = rr.risk_days_bs(
            k=k_mode_human,
            k_posterior_sample=k_human,
            **params,
        )
        iwp_human = r_human[0]

        # Inverse Gamma with same mode
        r_ig = rr.risk_days_bs(
            k=0.000673,
            k_invgamma_alpha=2.0,
            k_invgamma_mode=0.000673,
            **params,
        )
        iwp_ig = r_ig[0]

        assert iwp_human > 0
        assert iwp_ig > 0
        assert iwp_human == pytest.approx(iwp_ig, rel=0.25)


# ---------------------------------------------------------------------------
# residual_risk_rd
# ---------------------------------------------------------------------------


class TestResidualRiskRd:
    # Simulate a plausible iwp_bs (500 draws centred on ~4 days)
    _IWP_BS = list(_RNG.gamma(shape=2.0, scale=2.0, size=500))

    def test_basic_calculation(self):
        # incidence = 10/100 000 PY = 1e-4/year; iwp_pe = 4 days
        # rr_pe = 1e-4 * 4 / 365.25 * 1e6 ≈ 1.095 per million
        rr_pe, _, _ = rr.residual_risk_rd(
            iwp_pe=4.0,
            iwp_bs=self._IWP_BS,
            incidence=1e-4,
            incidence_norm_sd=1e-5,
            per=1e6,
            seed=42,
        )
        expected = 1e-4 * 4.0 / 365.25 * 1e6
        assert rr_pe == pytest.approx(expected, rel=1e-9)

    def test_one_in_x_inversion(self):
        # one_in_x=True should give 1 / (incidence * iwp_pe / 365.25)
        incidence = 1e-4
        iwp_pe = 4.0
        rr_pe, _, _ = rr.residual_risk_rd(
            iwp_pe=iwp_pe,
            iwp_bs=self._IWP_BS,
            incidence=incidence,
            incidence_norm_sd=1e-5,
            seed=42,
            one_in_x=True,
        )
        expected = 1.0 / (incidence * iwp_pe / 365.25)
        assert rr_pe == pytest.approx(expected, rel=1e-9)

    def test_cri_is_ordered(self):
        _, cri, _ = rr.residual_risk_rd(
            iwp_pe=4.0,
            iwp_bs=self._IWP_BS,
            incidence=1e-4,
            incidence_norm_sd=1e-5,
            seed=42,
        )
        assert cri[0] < cri[1]

    def test_higher_incidence_gives_higher_rr(self):
        common = dict(iwp_pe=4.0, iwp_bs=self._IWP_BS, incidence_norm_sd=1e-5, seed=42)
        rr_low, _, _ = rr.residual_risk_rd(incidence=1e-4, **common)
        rr_high, _, _ = rr.residual_risk_rd(incidence=5e-4, **common)
        assert rr_high > rr_low

    def test_zero_incidence_raises(self):
        with pytest.raises(ValueError):
            rr.residual_risk_rd(
                iwp_pe=4.0,
                iwp_bs=self._IWP_BS,
                incidence=0.0,
                incidence_norm_sd=1e-5,
            )

    def test_negative_incidence_raises(self):
        with pytest.raises(ValueError):
            rr.residual_risk_rd(
                iwp_pe=4.0,
                iwp_bs=self._IWP_BS,
                incidence=-1e-4,
                incidence_norm_sd=1e-5,
            )

    def test_zero_iwp_pe_raises(self):
        with pytest.raises(ValueError):
            rr.residual_risk_rd(
                iwp_pe=0.0,
                iwp_bs=self._IWP_BS,
                incidence=1e-4,
                incidence_norm_sd=1e-5,
            )


# ---------------------------------------------------------------------------
# risk_days_bs with InvGamma k distribution (Python and Go)
# ---------------------------------------------------------------------------

# Shared kwargs for InvGamma bootstrap tests
_INVGAMMA_BS_KWARGS = dict(
    k=0.000673,  # InvGamma(2, 0.002019) mode
    doubling_time=DEFAULTS["doubling_time"],
    doubling_time_norm_sd=DEFAULTS["doubling_time_norm_sd"],
    lod50=DEFAULTS["lod50"],
    lod50_sd=DEFAULTS["lod50_sd"],
    lod95_lod50_ratio=DEFAULTS["lod95_lod50_ratio"],
    volume_transfused=DEFAULTS["volume_transfused"],
    volume_transfused_range=DEFAULTS["volume_transfused_range"],
    pool_size=DEFAULTS["pool_size"],
    retests=DEFAULTS["retests"],
    k_invgamma_alpha=2.0,
    k_invgamma_beta=0.002019,
    n_bs=500,
    seed=42,
    threads=2,
)


class TestRiskDaysBsPythonInvGamma:
    """InvGamma k distribution via the Python backend."""

    def test_sanity(self):
        result = rr.risk_days_bs(**_INVGAMMA_BS_KWARGS, use_go=False)
        _assert_bs_result_sane(result, _INVGAMMA_BS_KWARGS["n_bs"])

    def test_mode_parameterisation_equivalent(self):
        """k_invgamma_mode should give same result as the equivalent k_invgamma_beta."""
        kwargs_mode = {
            k: v for k, v in _INVGAMMA_BS_KWARGS.items() if k != "k_invgamma_beta"
        }
        kwargs_mode["k_invgamma_mode"] = 0.000673  # beta = 0.000673 * 3 = 0.002019
        r_beta = rr.risk_days_bs(**_INVGAMMA_BS_KWARGS, use_go=False)
        r_mode = rr.risk_days_bs(**kwargs_mode, use_go=False)
        assert sorted(r_beta[3]) == pytest.approx(sorted(r_mode[3]), rel=1e-6)

    def test_reproducible(self):
        r1 = rr.risk_days_bs(**_INVGAMMA_BS_KWARGS, use_go=False)
        r2 = rr.risk_days_bs(**_INVGAMMA_BS_KWARGS, use_go=False)
        assert sorted(r1[3]) == sorted(r2[3])

    def test_point_estimate_matches_risk_days(self):
        result = rr.risk_days_bs(
            **_INVGAMMA_BS_KWARGS, use_go=False, point_estimate="primary parameters"
        )
        expected_pe = rr._risk_days(
            copies_per_virion=DEFAULTS["copies_per_virion"],
            C0=DEFAULTS["C0"],
            doubling_time=DEFAULTS["doubling_time"],
            volume_transfused=DEFAULTS["volume_transfused"],
            k=_INVGAMMA_BS_KWARGS["k"],
            pool_size=DEFAULTS["pool_size"],
            lod50=DEFAULTS["lod50"],
            lod95_lod50_ratio=DEFAULTS["lod95_lod50_ratio"],
            retests=DEFAULTS["retests"],
        )
        assert result[0] == pytest.approx(expected_pe, rel=1e-6)


class TestRiskDaysBsGoInvGamma:
    """InvGamma k distribution via the Go backend."""

    def test_sanity(self):
        result = rr.risk_days_bs(**_INVGAMMA_BS_KWARGS, use_go=True)
        _assert_bs_result_sane(result, _INVGAMMA_BS_KWARGS["n_bs"])

    def test_point_estimate_matches_risk_days(self):
        result = rr.risk_days_bs(
            **_INVGAMMA_BS_KWARGS, use_go=True, point_estimate="primary parameters"
        )
        expected_pe = rr._risk_days(
            copies_per_virion=DEFAULTS["copies_per_virion"],
            C0=DEFAULTS["C0"],
            doubling_time=DEFAULTS["doubling_time"],
            volume_transfused=DEFAULTS["volume_transfused"],
            k=_INVGAMMA_BS_KWARGS["k"],
            pool_size=DEFAULTS["pool_size"],
            lod50=DEFAULTS["lod50"],
            lod95_lod50_ratio=DEFAULTS["lod95_lod50_ratio"],
            retests=DEFAULTS["retests"],
        )
        assert result[0] == pytest.approx(expected_pe, rel=1e-6)

    def test_mode_parameterisation_accepted(self):
        """k_invgamma_mode should be accepted by the Go path."""
        kwargs_mode = {
            k: v for k, v in _INVGAMMA_BS_KWARGS.items() if k != "k_invgamma_beta"
        }
        kwargs_mode["k_invgamma_mode"] = 0.000673
        result = rr.risk_days_bs(**kwargs_mode, use_go=True)
        _assert_bs_result_sane(result, _INVGAMMA_BS_KWARGS["n_bs"])

    def test_simulation_medians_agree_with_python(self):
        """Python and Go InvGamma paths should agree within 15% on the median."""
        py = rr.risk_days_bs(**_INVGAMMA_BS_KWARGS, use_go=False)
        go = rr.risk_days_bs(**_INVGAMMA_BS_KWARGS, use_go=True)
        assert np.median(py[3]) == pytest.approx(np.median(go[3]), rel=0.15)

