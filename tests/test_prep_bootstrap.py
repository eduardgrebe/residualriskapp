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

"""
Integration tests for the PrEP bootstrap function risk_days_prep_bs().

These tests exercise the full bootstrap pipeline (parameter sampling →
integration → statistics) for each k-distribution path, result structure
validation, reproducibility, point-estimate methods, and PrEP-specific
parameter effects.

Bootstrap tests require ProcessPoolExecutor, which fails in the macOS
sandbox with PermissionError.  Run outside the sandbox:

    pytest tests/test_prep_bootstrap.py -v
"""

import unittest

import numpy as np
import polars as pl
import pytest

from residualrisk._go import find_go_binary
from residualrisk.prep import _risk_days_prep, risk_days_prep_bs


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

# Realistic PrEP bootstrap parameters (small n_bs for speed)
PREP_BS_KWARGS = dict(
    k=0.000673,
    doubling_time=0.8542,
    doubling_time_norm_sd=0.2813,
    lod50=2.73,
    lod50_sd=0.53,
    lod95_lod50_ratio=3.5,
    volume_transfused=200,
    volume_transfused_range=(100, 340),
    pool_size=16,
    retests=1,
    set_point=336,
    set_point_dist_uniform=(19.1, 2265),
    eclipse=7.0,
    eclipse_dist_uniform=(4.0, 10.0),
    a=0.7,
    b=0.6,
    offset=1,
    ser_min=28.7,
    ser_max=250,
    ser_alpha=50.49434,
    ser_beta=1.15062,
    n_bs=30,
    seed=42,
    threads=1,
    point_estimate="median",
)


# ---------------------------------------------------------------------------
# Result structure & statistics
# ---------------------------------------------------------------------------

@pytest.mark.multiprocessing
class TestPrepBsResultStructure(unittest.TestCase):
    """Verify the shape and ordering of risk_days_prep_bs() outputs."""

    def setUp(self):
        self.rd_pe, self.rd_cri, self.rd_range, self.rdests, self.sim_df = (
            risk_days_prep_bs(
                **PREP_BS_KWARGS,
                k_invgamma_alpha=2.0,
                k_invgamma_beta=0.002019,
            )
        )

    def test_pe_positive(self):
        self.assertGreater(self.rd_pe, 0)

    def test_returns_five_tuple(self):
        result = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
        )
        self.assertEqual(len(result), 5)

    def test_correct_number_of_simulations(self):
        self.assertEqual(len(self.rdests), PREP_BS_KWARGS["n_bs"])

    def test_all_simulations_non_negative(self):
        self.assertTrue(all(r >= 0 for r in self.rdests))

    def test_cri_ordered(self):
        self.assertLessEqual(self.rd_cri[0], self.rd_cri[1])

    def test_range_contains_cri(self):
        self.assertLessEqual(self.rd_range[0], self.rd_cri[0])
        self.assertGreaterEqual(self.rd_range[1], self.rd_cri[1])

    def test_pe_within_range(self):
        self.assertGreaterEqual(self.rd_pe, self.rd_range[0])
        self.assertLessEqual(self.rd_pe, self.rd_range[1])

    def test_sim_df_none_by_default(self):
        self.assertIsNone(self.sim_df)


# ---------------------------------------------------------------------------
# return_sim_df
# ---------------------------------------------------------------------------

@pytest.mark.multiprocessing
class TestPrepBsSimDf(unittest.TestCase):
    """Verify the simulation DataFrame when return_sim_df=True."""

    def setUp(self):
        self.rd_pe, self.rd_cri, self.rd_range, self.rdests, self.sim_df = (
            risk_days_prep_bs(
                **PREP_BS_KWARGS,
                k_invgamma_alpha=2.0,
                k_invgamma_beta=0.002019,
                return_sim_df=True,
            )
        )

    def test_sim_df_is_dataframe(self):
        self.assertIsInstance(self.sim_df, pl.DataFrame)

    def test_sim_df_row_count(self):
        self.assertEqual(len(self.sim_df), PREP_BS_KWARGS["n_bs"])

    def test_sim_df_has_iwp(self):
        self.assertIn("iwp", self.sim_df.columns)

    def test_sim_df_has_prep_columns(self):
        for col in ("set_point", "eclipse", "a", "b", "offset",
                     "ser_min", "ser_max", "ser_alpha", "ser_beta"):
            self.assertIn(col, self.sim_df.columns, f"Missing column: {col}")

    def test_sim_df_has_baseline_columns(self):
        for col in ("k", "doubling_time", "lod50", "volume_transfused"):
            self.assertIn(col, self.sim_df.columns, f"Missing column: {col}")

    def test_sim_df_iwp_matches_rdests(self):
        np.testing.assert_array_almost_equal(
            self.sim_df["iwp"].to_numpy(), self.rdests, decimal=10
        )

    def test_sim_df_set_points_sampled(self):
        """set_point values should be drawn from the uniform distribution."""
        sp = self.sim_df["set_point"].to_numpy()
        lo, hi = PREP_BS_KWARGS["set_point_dist_uniform"]
        self.assertTrue(np.all(sp >= lo))
        self.assertTrue(np.all(sp <= hi))
        # Should not all be identical (they're sampled uniformly)
        self.assertGreater(len(np.unique(sp)), 1)

    def test_sim_df_eclipses_sampled(self):
        """eclipse values should be drawn from the uniform distribution."""
        ec = self.sim_df["eclipse"].to_numpy()
        lo, hi = PREP_BS_KWARGS["eclipse_dist_uniform"]
        self.assertTrue(np.all(ec >= lo))
        self.assertTrue(np.all(ec <= hi))
        self.assertGreater(len(np.unique(ec)), 1)


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

@pytest.mark.multiprocessing
class TestPrepBsReproducibility(unittest.TestCase):
    """Same seed → same results; different seed → different results."""

    def test_same_seed_same_results(self):
        r1 = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
        )
        r2 = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
        )
        self.assertEqual(r1[0], r2[0])  # PE
        np.testing.assert_array_equal(r1[3], r2[3])  # rdests

    def test_different_seed_different_results(self):
        r1 = risk_days_prep_bs(
            **{**PREP_BS_KWARGS, "seed": 1},
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
        )
        r2 = risk_days_prep_bs(
            **{**PREP_BS_KWARGS, "seed": 2},
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
        )
        self.assertNotEqual(r1[0], r2[0])


# ---------------------------------------------------------------------------
# Point-estimate methods
# ---------------------------------------------------------------------------

@pytest.mark.multiprocessing
class TestPrepBsPointEstimates(unittest.TestCase):
    """Verify each point_estimate method produces sensible results."""

    def _run(self, method):
        return risk_days_prep_bs(
            **{**PREP_BS_KWARGS, "point_estimate": method},
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
        )

    def test_median(self):
        rd_pe, _, _, rdests, _ = self._run("median")
        self.assertAlmostEqual(rd_pe, float(np.median(rdests)), places=10)

    def test_mean(self):
        rd_pe, _, _, rdests, _ = self._run("mean")
        self.assertAlmostEqual(rd_pe, float(np.mean(rdests)), places=10)

    def test_primary_parameters(self):
        rd_pe, rd_cri, _, rdests, _ = self._run("primary parameters")
        self.assertGreater(rd_pe, 0)
        # Primary-params PE uses fixed (non-sampled) values, so it need not
        # lie within the CrI, but it should be in a plausible range.
        self.assertGreater(rd_pe, 0)


# ---------------------------------------------------------------------------
# k-distribution paths — full bootstrap
# ---------------------------------------------------------------------------

class TestPrepBsKDistPaths(unittest.TestCase):
    """Exercise every k-distribution path through the full bootstrap."""

    def _assert_valid(self, result, n_bs=PREP_BS_KWARGS["n_bs"]):
        rd_pe, rd_cri, rd_range, rdests, _ = result
        self.assertGreater(rd_pe, 0)
        self.assertEqual(len(rdests), n_bs)
        self.assertTrue(all(r >= 0 for r in rdests))
        self.assertLessEqual(rd_cri[0], rd_cri[1])

    @pytest.mark.multiprocessing
    def test_posterior_sample(self):
        posterior = np.random.default_rng(99).exponential(0.001, size=200)
        result = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_posterior_sample=posterior,
        )
        self._assert_valid(result)

    @pytest.mark.multiprocessing
    def test_invgamma_alpha_beta(self):
        result = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
        )
        self._assert_valid(result)

    @pytest.mark.multiprocessing
    def test_invgamma_alpha_mode(self):
        result = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_invgamma_alpha=2.0,
            k_invgamma_mode=0.000673,
        )
        self._assert_valid(result)

    @pytest.mark.multiprocessing
    def test_lnmix(self):
        result = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_lnmix_w=0.90,
            k_lnmix_mu1=-7.2403,
            k_lnmix_sigma1=0.3241,
            k_lnmix_mu2=-3.7423,
            k_lnmix_sigma2=0.5258,
        )
        self._assert_valid(result)

    @pytest.mark.multiprocessing
    def test_legacy_gamma(self):
        result = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_gamma_shape=2.0,
            k_gamma_scale=0.001,
        )
        self._assert_valid(result)

    def test_no_distribution_raises(self):
        with self.assertRaisesRegex(ValueError, "k-distribution"):
            risk_days_prep_bs(**PREP_BS_KWARGS)


# ---------------------------------------------------------------------------
# PrEP-specific parameter effects
# ---------------------------------------------------------------------------

@unittest.skipUnless(find_go_binary(), "Go binary not available")
class TestPrepBsParameterEffects(unittest.TestCase):
    """Verify that PrEP-specific parameters affect results as expected.

    Uses n_bs=200 with Go acceleration so directional effects reliably
    exceed Monte Carlo noise.
    """

    def _run(self, **overrides):
        kw = {**PREP_BS_KWARGS, "n_bs": 200, "use_go": True, **overrides}
        return risk_days_prep_bs(
            **kw,
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
        )

    def test_more_retests_higher_risk(self):
        """More retests → higher RDE.

        Each retest is an extra opportunity for a false-negative result
        to release an initially-reactive donation, so retests=1 should
        give higher risk than retests=0.

        Uses pool_size=1 so that individual-donation NAT sensitivity is
        high enough for the retest effect to be detectable.
        """
        rd_retests0, _, _, _, _ = self._run(retests=0, pool_size=1)
        rd_retests1, _, _, _, _ = self._run(retests=1, pool_size=1)
        self.assertGreater(rd_retests1, rd_retests0)

    def test_larger_pool_higher_risk(self):
        """Larger pool → higher RDE (diluted sample → harder to detect)."""
        rd_pool1, _, _, _, _ = self._run(pool_size=1)
        rd_pool16, _, _, _, _ = self._run(pool_size=16)
        self.assertGreater(rd_pool16, rd_pool1)

    def test_shorter_eclipse_higher_risk(self):
        """Shorter eclipse → higher RDE (viremia starts earlier)."""
        rd_e2, _, _, _, _ = self._run(
            eclipse=2.0, eclipse_dist_uniform=(1.0, 3.0)
        )
        rd_e14, _, _, _, _ = self._run(
            eclipse=14.0, eclipse_dist_uniform=(12.0, 16.0)
        )
        self.assertGreater(rd_e2, rd_e14)

    def test_wider_serology_window_higher_risk(self):
        """Later seroconversion (wider serology window) → higher RDE."""
        # Tight serology: seroconversion detects early (small ser_max, fast decay)
        rd_tight, _, _, _, _ = self._run(
            ser_min=10, ser_max=50, ser_alpha=5.0, ser_beta=3.0
        )
        # Wide serology: seroconversion detects late (large ser_max, slow decay)
        rd_wide, _, _, _, _ = self._run(
            ser_min=28.7, ser_max=500, ser_alpha=100.0, ser_beta=1.5
        )
        self.assertGreater(rd_wide, rd_tight)


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestPrepBsEdgeCases(unittest.TestCase):
    """Edge cases and error handling."""

    def test_n_bs_zero_raises(self):
        with self.assertRaisesRegex(ValueError, "n_bs must be greater than zero"):
            risk_days_prep_bs(
                **{**PREP_BS_KWARGS, "n_bs": 0},
                k_invgamma_alpha=2.0,
                k_invgamma_beta=0.002019,
            )

    def test_n_bs_negative_raises(self):
        with self.assertRaisesRegex(ValueError, "n_bs must be greater than zero"):
            risk_days_prep_bs(
                **{**PREP_BS_KWARGS, "n_bs": -5},
                k_invgamma_alpha=2.0,
                k_invgamma_beta=0.002019,
            )

    @pytest.mark.multiprocessing
    def test_n_bs_one(self):
        """Single bootstrap iteration should still work."""
        rd_pe, rd_cri, rd_range, rdests, _ = risk_days_prep_bs(
            **{**PREP_BS_KWARGS, "n_bs": 1},
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
        )
        self.assertEqual(len(rdests), 1)
        self.assertGreater(rd_pe, 0)


# ---------------------------------------------------------------------------
# Go acceleration (use_go=True)
# ---------------------------------------------------------------------------

@unittest.skipIf(find_go_binary() is None, "Go binary not available")
class TestPrepBsGoDispatch(unittest.TestCase):
    """Verify that use_go=True dispatches to Go and returns valid results."""

    def test_use_go_invgamma(self):
        rd_pe, rd_cri, rd_range, rdests, sim_df = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
            use_go=True,
        )
        self.assertGreater(rd_pe, 0)
        self.assertEqual(len(rdests), PREP_BS_KWARGS["n_bs"])
        self.assertLessEqual(rd_cri[0], rd_cri[1])

    def test_use_go_return_sim_df(self):
        rd_pe, rd_cri, rd_range, rdests, sim_df = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
            use_go=True,
            return_sim_df=True,
        )
        self.assertIsInstance(sim_df, pl.DataFrame)
        self.assertEqual(len(sim_df), PREP_BS_KWARGS["n_bs"])
        self.assertIn("set_point", sim_df.columns)
        self.assertIn("eclipse", sim_df.columns)

    def test_use_go_reproducible(self):
        r1 = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
            use_go=True,
        )
        r2 = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_invgamma_alpha=2.0,
            k_invgamma_beta=0.002019,
            use_go=True,
        )
        self.assertEqual(r1[0], r2[0])
        np.testing.assert_array_equal(r1[3], r2[3])

    def test_use_go_lnmix(self):
        rd_pe, _, _, rdests, _ = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_lnmix_w=0.90,
            k_lnmix_mu1=-7.2403,
            k_lnmix_sigma1=0.3241,
            k_lnmix_mu2=-3.7423,
            k_lnmix_sigma2=0.5258,
            use_go=True,
        )
        self.assertGreater(rd_pe, 0)
        self.assertEqual(len(rdests), PREP_BS_KWARGS["n_bs"])

    def test_use_go_posterior(self):
        posterior = np.random.default_rng(99).exponential(0.001, size=200)
        rd_pe, _, _, rdests, _ = risk_days_prep_bs(
            **PREP_BS_KWARGS,
            k_posterior_sample=posterior,
            use_go=True,
        )
        self.assertGreater(rd_pe, 0)
        self.assertEqual(len(rdests), PREP_BS_KWARGS["n_bs"])


# ---------------------------------------------------------------------------
# Integration method (gauss-legendre default vs quad)
# ---------------------------------------------------------------------------

# Standard single-call PrEP params for _risk_days_prep. Serology is supplied
# per-test so the two regimes (narrow vs wide active window) can be contrasted.
_PREP_SINGLE = dict(
    copies_per_virion=2,
    C0=0.00025,
    doubling_time=0.8542,
    set_point=336,
    eclipse=7.0,
    a=0.7,
    b=0.6,
    offset=1.0,
    volume_transfused=200.0,
    k=0.000673,
    pool_size=16,
    lod50=2.73,
    lod95_lod50_ratio=3.5,
    retests=1,
    z=1.6449,
)
# Narrow active window (~[8.7, 22.5] days): sharp Weibull serology cutoff. This
# is the regime where adaptive quad silently misses the peak and returns ~0.
_SER_NARROW = dict(ser_min=10, ser_max=500, ser_alpha=9.1, ser_beta=5.2)
# Wide active window (~[10, 169] days): the production serology defaults, where
# quad already integrates correctly.
_SER_PROD = dict(ser_min=28.7, ser_max=250, ser_alpha=50.49434, ser_beta=1.15062)

# Fine-grid Simpson reference values (0.01-day grid over [-100, 500]).
_TRUTH_NARROW = 1.00864
_TRUTH_PROD = 3.09187


class TestPrepIntegrationMethod(unittest.TestCase):
    """PrEP risk-days integration defaults to a fixed 1000-point Gauss-Legendre
    rule (matching the Go backend), with adaptive scipy quad selectable for
    reproducing prior analyses.

    Unlike the baseline integrand, the PrEP integrand has *compact support* (it
    is exactly zero before the eclipse phase and after the serology cutoff), so
    adaptive quad can silently miss a narrow active window and return ~0 — which
    Gauss-Legendre fixes.

    These call _risk_days_prep directly (a single deterministic integration), so
    they run without ProcessPoolExecutor and are safe inside the sandbox.
    """

    def test_default_is_gauss_legendre(self):
        self.assertEqual(
            _risk_days_prep(**_PREP_SINGLE, **_SER_PROD),
            _risk_days_prep(
                **_PREP_SINGLE, **_SER_PROD, integration_method="gauss-legendre"
            ),
        )

    def test_gl_recovers_narrow_window_where_quad_fails(self):
        # The key regression: at narrow serology GL recovers the true ~1.0086,
        # while quad collapses to ~0 (the silent compact-support failure).
        gl = _risk_days_prep(
            **_PREP_SINGLE, **_SER_NARROW, integration_method="gauss-legendre"
        )
        qd = _risk_days_prep(
            **_PREP_SINGLE, **_SER_NARROW, integration_method="quad"
        )
        assert gl == pytest.approx(_TRUTH_NARROW, rel=1e-3)
        self.assertGreater(gl, 1.0)
        # quad is catastrophically wrong here (orders of magnitude below truth)
        self.assertLess(qd, 1e-3)

    def test_gl_matches_quad_at_production(self):
        # At production serology the window is wide; GL and quad agree and both
        # match the Simpson reference — confirming the default switch does not
        # shift standard-scenario results.
        gl = _risk_days_prep(
            **_PREP_SINGLE, **_SER_PROD, integration_method="gauss-legendre"
        )
        qd = _risk_days_prep(
            **_PREP_SINGLE, **_SER_PROD, integration_method="quad"
        )
        assert gl == pytest.approx(qd, rel=1e-3)
        assert gl == pytest.approx(_TRUTH_PROD, rel=1e-3)

    def test_invalid_method_raises(self):
        with self.assertRaises(ValueError):
            _risk_days_prep(
                **_PREP_SINGLE, **_SER_PROD, integration_method="simpson"
            )

    def test_use_go_with_quad_raises(self):
        # Conflict guard: the Go backend only implements gauss-legendre. Raises
        # before any bootstrap runs, so this is sandbox-safe.
        with self.assertRaises(ValueError):
            risk_days_prep_bs(
                **PREP_BS_KWARGS,
                k_invgamma_alpha=2.0,
                k_invgamma_beta=0.002019,
                use_go=True,
                integration_method="quad",
            )


if __name__ == "__main__":
    unittest.main()
