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
Cross-validation tests: Python PrEP ↔ Go PrEP parity.

These tests require the Go binary (go/bin/riskdays_go) and
ProcessPoolExecutor (fails in macOS sandbox with PermissionError).
Run outside the sandbox: pytest tests/test_prep_go_parity.py -v
"""

import unittest

import numpy as np
import pytest

from residualrisk._go import find_go_binary, risk_days_prep_bs_go
from residualrisk.prep import risk_days_prep_bs


COMMON_KWARGS = dict(
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
    k_invgamma_alpha=2.0,
    k_invgamma_beta=0.002019,
    n_bs=500,
    seed=42,
    point_estimate="median",
    return_sim_df=True,
)


@unittest.skipIf(find_go_binary() is None, "Go binary not available")
class TestPrepGoParity(unittest.TestCase):
    """Go-side sanity, reproducibility, and dispatch checks.

    These exercise only the Go path (binary subprocess), so they are
    sandbox-safe. The actual numerical Python↔Go cross-validation lives in
    TestPrepPythonGoAgreement below.
    """

    def test_go_prep_returns_valid(self):
        """Go PrEP wrapper returns valid results."""
        rd_pe, rd_cri, rd_range, rdests, sim_df = risk_days_prep_bs_go(
            **COMMON_KWARGS
        )
        self.assertGreater(rd_pe, 0)
        self.assertEqual(len(rdests), 500)
        self.assertLessEqual(rd_cri[0], rd_pe)
        self.assertGreaterEqual(rd_cri[1], rd_pe)
        self.assertIsNotNone(sim_df)
        self.assertIn("set_point", sim_df.columns)
        self.assertIn("eclipse", sim_df.columns)

    def test_go_prep_reproducible(self):
        """Same seed → same results from Go."""
        r1 = risk_days_prep_bs_go(**COMMON_KWARGS)
        r2 = risk_days_prep_bs_go(**COMMON_KWARGS)
        self.assertEqual(r1[0], r2[0])  # PE
        np.testing.assert_array_equal(r1[3], r2[3])  # rdests

    def test_go_prep_different_k_dists(self):
        """Go PrEP works with different k distributions."""
        # Posterior sample
        rng = np.random.default_rng(99)
        kw = {**COMMON_KWARGS, "n_bs": 100}
        del kw["k_invgamma_alpha"]
        del kw["k_invgamma_beta"]
        kw["k_posterior_sample"] = rng.exponential(0.001, size=200)
        rd_pe, _, _, rdests, _ = risk_days_prep_bs_go(**kw)
        self.assertGreater(rd_pe, 0)
        self.assertEqual(len(rdests), 100)

    def test_use_go_dispatcher(self):
        """risk_days_prep_bs(use_go=True) dispatches to Go."""
        rd_pe, rd_cri, rd_range, rdests, sim_df = risk_days_prep_bs(
            **COMMON_KWARGS, use_go=True
        )
        self.assertGreater(rd_pe, 0)
        self.assertEqual(len(rdests), 500)
        self.assertIsNotNone(sim_df)

    def test_go_prep_lnmix(self):
        """Go PrEP works with lognormal mixture k distribution."""
        kw = {**COMMON_KWARGS, "n_bs": 100}
        del kw["k_invgamma_alpha"]
        del kw["k_invgamma_beta"]
        kw["k_lnmix_w"] = 0.90
        kw["k_lnmix_mu1"] = -7.2403
        kw["k_lnmix_sigma1"] = 0.3241
        kw["k_lnmix_mu2"] = -3.7423
        kw["k_lnmix_sigma2"] = 0.5258
        rd_pe, _, _, rdests, _ = risk_days_prep_bs_go(**kw)
        self.assertGreater(rd_pe, 0)
        self.assertEqual(len(rdests), 100)

    def test_go_prep_fixed_ab_default(self):
        """Without ranges, Go holds the sinusoidal a and b fixed at the scalars."""
        _, _, _, _, sim_df = risk_days_prep_bs_go(**{**COMMON_KWARGS, "n_bs": 100})
        self.assertEqual(sim_df["a"].n_unique(), 1)
        self.assertEqual(sim_df["b"].n_unique(), 1)
        self.assertAlmostEqual(sim_df["a"][0], COMMON_KWARGS["a"])
        self.assertAlmostEqual(sim_df["b"][0], COMMON_KWARGS["b"])

    def test_go_prep_varied_ab(self):
        """With ranges, Go samples a and b uniformly per iteration; sim_df
        carries the real per-iteration values, within range and a <= offset."""
        kw = {**COMMON_KWARGS, "n_bs": 300,
              "a_dist_uniform": (0.5, 0.9), "b_dist_uniform": (0.4, 0.8)}
        _, _, _, _, sim_df = risk_days_prep_bs_go(**kw)
        self.assertGreater(sim_df["a"].n_unique(), 1)
        self.assertGreater(sim_df["b"].n_unique(), 1)
        self.assertGreaterEqual(sim_df["a"].min(), 0.5)
        self.assertLessEqual(sim_df["a"].max(), 0.9)
        self.assertGreaterEqual(sim_df["b"].min(), 0.4)
        self.assertLessEqual(sim_df["b"].max(), 0.8)
        self.assertLessEqual(sim_df["a"].max(), COMMON_KWARGS["offset"])

    def test_go_prep_drug_effect_fixed_default(self):
        """Without a range, Go holds drug_effect fixed at 1.0 (no reduction)."""
        _, _, _, _, sim_df = risk_days_prep_bs_go(**{**COMMON_KWARGS, "n_bs": 100})
        self.assertEqual(sim_df["drug_effect"].n_unique(), 1)
        self.assertAlmostEqual(sim_df["drug_effect"][0], 1.0)

    def test_go_prep_drug_effect_varied(self):
        """With a range, Go samples drug_effect uniformly per iteration; sim_df
        carries the real per-iteration values, within range."""
        kw = {**COMMON_KWARGS, "n_bs": 300, "drug_effect_dist_uniform": (0.5, 0.9)}
        _, _, _, _, sim_df = risk_days_prep_bs_go(**kw)
        self.assertGreater(sim_df["drug_effect"].n_unique(), 1)
        self.assertGreaterEqual(sim_df["drug_effect"].min(), 0.5)
        self.assertLessEqual(sim_df["drug_effect"].max(), 0.9)

    def test_go_prep_drug_effect_linear(self):
        """drug_effect linearly scales the deterministic primary-parameters PE."""
        kw = {**COMMON_KWARGS, "n_bs": 100, "point_estimate": "primary parameters"}
        pe_full, _, _, _, _ = risk_days_prep_bs_go(**kw)
        pe_half, _, _, _, _ = risk_days_prep_bs_go(**{**kw, "drug_effect": 0.5})
        self.assertAlmostEqual(pe_half, 0.5 * pe_full, delta=abs(pe_full) * 1e-9)


@pytest.mark.multiprocessing
@unittest.skipIf(find_go_binary() is None, "Go binary not available")
class TestPrepPythonGoAgreement(unittest.TestCase):
    """Numerically cross-validate the Python and Go PrEP bootstrap.

    Modeled on the baseline TestPythonGoAgreement. The Python path uses
    ProcessPoolExecutor (hence @pytest.mark.multiprocessing — excluded from
    sandboxed/fast runs; run outside the sandbox). Python and Go use
    independent RNGs, so sampled-parameter draws differ and we compare
    distributional summaries within tolerance; the deterministic
    primary-parameters point estimate agrees far more tightly.

    Both sides run once in setUpClass with point_estimate="primary parameters",
    which returns the deterministic PE plus the full per-iteration sample.
    """

    # n_bs is bumped above COMMON_KWARGS's 500: the k prior is InvGamma(α=2)
    # (infinite variance), so the upper-CrI (97.5th pct) is too noisy to compare
    # at n_bs=500 (Python-vs-Go ~29%). At n_bs=2000 the agreement is stable
    # (median ~1%, CrI bounds ≤~10%), well within the tolerances below.
    N_BS = 2000

    @classmethod
    def setUpClass(cls):
        kwargs = {**COMMON_KWARGS, "point_estimate": "primary parameters", "n_bs": cls.N_BS}
        cls.py = risk_days_prep_bs(**kwargs, use_go=False)
        cls.go = risk_days_prep_bs(**kwargs, use_go=True)

    def test_primary_parameters_pe_agree(self):
        # Deterministic single integration on both sides: identical analytic
        # tcrit (Go FindTcrit / Python _find_tcrit) + a fixed 1000-point
        # Gauss-Legendre rule, so they agree to machine precision (residual
        # ~1e-13, limited only by numpy-vs-gonum GL node/weight roundoff). The
        # 1e-9 bound locks that in while still catching any real integration
        # discrepancy (e.g. a quad-style missed-window collapse).
        self.assertAlmostEqual(self.py[0], self.go[0], delta=abs(self.go[0]) * 1e-9)

    def test_medians_agree(self):
        # Independent RNGs → compare the bootstrap median within tolerance.
        py_median = float(np.median(self.py[3]))
        go_median = float(np.median(self.go[3]))
        self.assertAlmostEqual(py_median, go_median, delta=go_median * 0.20)

    def test_cri_bounds_agree(self):
        # 95% credible-interval bounds within a (wider) RNG-driven tolerance.
        for bound in (0, 1):
            self.assertAlmostEqual(
                self.py[1][bound], self.go[1][bound],
                delta=abs(self.go[1][bound]) * 0.25,
            )


@pytest.mark.multiprocessing
@unittest.skipIf(find_go_binary() is None, "Go binary not available")
class TestPrepPythonGoAgreementVariedAB(TestPrepPythonGoAgreement):
    """Same Python↔Go cross-validation as the parent, but with the sinusoidal
    amplitude (a) and frequency (b) *varied* uniformly — confirming the two
    backends sample and integrate the oscillation parameters equivalently.

    Inherits the parent's three checks. The primary-parameters PE test is
    unaffected (the PE uses the fixed scalar a/b, not the ranges), so it still
    agrees to ~1e-9. At n_bs=2000 the varied-a/b distributional agreement is
    median ~1.4% / CrI bounds ≲12%, well within the inherited 20%/25% bounds.
    """

    A_DIST = (0.5, 0.9)
    B_DIST = (0.4, 0.8)

    @classmethod
    def setUpClass(cls):
        kwargs = {
            **COMMON_KWARGS,
            "point_estimate": "primary parameters",
            "n_bs": cls.N_BS,
            "a_dist_uniform": cls.A_DIST,
            "b_dist_uniform": cls.B_DIST,
        }
        cls.py = risk_days_prep_bs(**kwargs, use_go=False)
        cls.go = risk_days_prep_bs(**kwargs, use_go=True)


@pytest.mark.multiprocessing
@unittest.skipIf(find_go_binary() is None, "Go binary not available")
class TestPrepPythonGoAgreementDrugEffect(TestPrepPythonGoAgreement):
    """Same Python<->Go cross-validation as the parent, but with the drug-effect
    transmissibility factor *varied* uniformly — confirming the two backends
    sample and apply it equivalently.

    Inherits the parent's three checks. The primary-parameters PE test is
    unaffected: the PE uses the scalar drug_effect (left at its 1.0 default
    here), not the range, so it still agrees to ~1e-9. The varied-drug-effect
    distributional agreement stays within the inherited 20%/25% bounds.
    """

    DRUG_EFFECT_DIST = (0.5, 1.0)

    @classmethod
    def setUpClass(cls):
        kwargs = {
            **COMMON_KWARGS,
            "point_estimate": "primary parameters",
            "n_bs": cls.N_BS,
            "drug_effect_dist_uniform": cls.DRUG_EFFECT_DIST,
        }
        cls.py = risk_days_prep_bs(**kwargs, use_go=False)
        cls.go = risk_days_prep_bs(**kwargs, use_go=True)


if __name__ == "__main__":
    unittest.main()
