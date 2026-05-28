# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

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

    @classmethod
    def setUpClass(cls):
        kwargs = {**COMMON_KWARGS, "point_estimate": "primary parameters"}
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


if __name__ == "__main__":
    unittest.main()
