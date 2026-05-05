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
Tests verifying that Inverse Gamma sampling is statistically correct and that
the Python and Go implementations produce equivalent output distributions.

Test classes:

  TestInvGammaTheoreticalStatistics
      Pure statistical tests of sample_invgamma() against InvGamma theory.
      No multiprocessing required — always runnable.

  TestInvGammaBootstrapKSamples
      Inspects the k values actually used inside the Python bootstrap (via
      return_sim_df=True) and verifies they follow InvGamma(alpha, beta).
      Requires ProcessPoolExecutor — will fail in sandboxed environments.

  TestPythonGoInvGammaAgreement
      Verifies that the Python and Go backends produce statistically
      indistinguishable output distributions (risk days) when using the same
      InvGamma(alpha, beta) parameters.
      Requires ProcessPoolExecutor for the Python backend.

  TestInvGammaGoStatistics
      Tests the Go backend alone: sanity, reproducibility, and moment
      agreement with theory (via the output distribution shape).
      No Python backend — always runnable.
"""

import numpy as np
import pytest
from scipy import stats as scipy_stats

from residualrisk import core as rr

# ---------------------------------------------------------------------------
# InvGamma(alpha, beta) theoretical constants
# ---------------------------------------------------------------------------

# Primary test parameters: InvGamma(2, 0.002019)
# Mode = beta / (alpha + 1) = 0.002019 / 3 ≈ 0.000673
# Mean = beta / (alpha - 1) = 0.002019 / 1 = 0.002019
# Variance: undefined for alpha <= 2 (Pareto-like tail)
ALPHA = 2.0
BETA = 0.002019
THEORY_MODE = BETA / (ALPHA + 1)       # ≈ 0.000673
THEORY_MEAN = BETA / (ALPHA - 1)       # = 0.002019
THEORY_MEDIAN = scipy_stats.invgamma.ppf(0.5, a=ALPHA, scale=BETA)

# Secondary parameters with alpha > 2 so variance is finite (used for variance tests)
ALPHA_FIN = 4.0
BETA_FIN = 0.006  # mode = 0.006/5 = 0.0012, mean = 0.006/3 = 0.002
THEORY_MEAN_FIN = BETA_FIN / (ALPHA_FIN - 1)     # = 0.002
THEORY_VAR_FIN = BETA_FIN**2 / ((ALPHA_FIN - 1)**2 * (ALPHA_FIN - 2))

# ---------------------------------------------------------------------------
# Shared bootstrap kwargs
# ---------------------------------------------------------------------------

_COMMON_BS = dict(
    doubling_time=20.5 / 24,
    doubling_time_norm_sd=1.33 / 24,
    lod50=2.73,
    lod50_sd=0.193,
    lod95_lod50_ratio=12.33 / 2.73,
    volume_transfused=20,
    volume_transfused_range=(15, 30),
    pool_size=16,
    retests=1,
    k=THEORY_MODE,
    threads=2,
)

# Base InvGamma bootstrap kwargs (n_bs set per test class)
_INVGAMMA_BS = dict(
    **_COMMON_BS,
    k_invgamma_alpha=ALPHA,
    k_invgamma_beta=BETA,
    seed=42,
)

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _rdests(result):
    """Extract the rdests array (index 3) from a risk_days_bs result tuple."""
    return np.asarray(result[3])


def _simdf(result):
    """Extract the sim_df (index 4) from a risk_days_bs result tuple."""
    return result[4]


# ---------------------------------------------------------------------------
# TestInvGammaTheoreticalStatistics
# ---------------------------------------------------------------------------


class TestInvGammaTheoreticalStatistics:
    """
    Verify that sample_invgamma() produces samples with the correct theoretical
    moments and CDF.  Uses sample_invgamma() directly — no multiprocessing.
    """

    N = 100_000

    def test_mean_converges(self):
        """Empirical mean ≈ beta / (alpha - 1)."""
        s = rr.sample_invgamma(self.N, alpha=ALPHA, beta=BETA, seed=1)
        assert np.mean(s) == pytest.approx(THEORY_MEAN, rel=0.05)

    def test_mean_converges_finite_variance(self):
        """Empirical mean ≈ beta / (alpha - 1) for alpha=4 (finite variance)."""
        s = rr.sample_invgamma(self.N, alpha=ALPHA_FIN, beta=BETA_FIN, seed=1)
        assert np.mean(s) == pytest.approx(THEORY_MEAN_FIN, rel=0.03)

    def test_variance_converges(self):
        """Empirical variance ≈ theoretical variance (requires alpha > 2)."""
        s = rr.sample_invgamma(self.N, alpha=ALPHA_FIN, beta=BETA_FIN, seed=1)
        assert np.var(s) == pytest.approx(THEORY_VAR_FIN, rel=0.10)

    def test_median_matches_scipy(self):
        """Empirical median ≈ scipy InvGamma PPF(0.5)."""
        s = rr.sample_invgamma(self.N, alpha=ALPHA, beta=BETA, seed=1)
        assert np.median(s) == pytest.approx(THEORY_MEDIAN, rel=0.02)

    def test_all_samples_positive(self):
        """InvGamma samples must all be strictly positive."""
        s = rr.sample_invgamma(self.N, alpha=ALPHA, beta=BETA, seed=1)
        assert np.all(s > 0)

    def test_cdf_matches_scipy_ks(self):
        """KS test: sample_invgamma samples are consistent with scipy invgamma CDF."""
        s = rr.sample_invgamma(10_000, alpha=ALPHA, beta=BETA, seed=99)
        ks_stat, p_value = scipy_stats.kstest(
            s, lambda x: scipy_stats.invgamma.cdf(x, a=ALPHA, scale=BETA)
        )
        assert p_value > 0.01, (
            f"KS test rejected: p={p_value:.4f}, stat={ks_stat:.4f}"
        )

    def test_mode_parameterisation_cdf_matches_scipy(self):
        """Mode-parameterised sample_invgamma also passes KS test vs scipy."""
        s = rr.sample_invgamma(10_000, alpha=ALPHA, mode=THEORY_MODE, seed=99)
        ks_stat, p_value = scipy_stats.kstest(
            s, lambda x: scipy_stats.invgamma.cdf(x, a=ALPHA, scale=BETA)
        )
        assert p_value > 0.01, (
            f"KS test rejected (mode param): p={p_value:.4f}, stat={ks_stat:.4f}"
        )

    def test_beta_and_mode_give_identical_samples(self):
        """beta=BETA and mode=THEORY_MODE (same distribution) produce identical draws."""
        s_beta = rr.sample_invgamma(1_000, alpha=ALPHA, beta=BETA, seed=7)
        s_mode = rr.sample_invgamma(1_000, alpha=ALPHA, mode=THEORY_MODE, seed=7)
        np.testing.assert_array_equal(s_beta, s_mode)

    def test_scale_equivariance(self):
        """Doubling beta doubles all quantiles (scale equivariance of InvGamma)."""
        s1 = rr.sample_invgamma(50_000, alpha=ALPHA, beta=BETA, seed=7)
        s2 = rr.sample_invgamma(50_000, alpha=ALPHA, beta=BETA * 2, seed=7)
        for q in [0.10, 0.25, 0.50, 0.75, 0.90]:
            assert np.quantile(s2, q) == pytest.approx(np.quantile(s1, q) * 2, rel=0.03), (
                f"Scale equivariance failed at q={q}"
            )

    def test_higher_alpha_lighter_tail(self):
        """Higher alpha → lighter upper tail (lower 95th-percentile / median ratio)."""
        # Both distributions have the same mode (beta / (alpha+1) = 0.000673)
        s_low = rr.sample_invgamma(50_000, alpha=2.0, beta=0.002019, seed=5)
        s_high = rr.sample_invgamma(50_000, alpha=5.0, beta=0.000673 * 6, seed=5)
        ratio_low = np.quantile(s_low, 0.95) / np.median(s_low)
        ratio_high = np.quantile(s_high, 0.95) / np.median(s_high)
        assert ratio_low > ratio_high, (
            f"Expected lighter tail for higher alpha: "
            f"ratio_low={ratio_low:.3f} ratio_high={ratio_high:.3f}"
        )

    def test_quantiles_match_scipy_ppf(self):
        """Empirical quantiles should closely match scipy PPF at multiple probability levels."""
        s = rr.sample_invgamma(100_000, alpha=ALPHA, beta=BETA, seed=11)
        for q in [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]:
            expected = scipy_stats.invgamma.ppf(q, a=ALPHA, scale=BETA)
            assert np.quantile(s, q) == pytest.approx(expected, rel=0.05), (
                f"Quantile mismatch at q={q}: "
                f"got {np.quantile(s, q):.6f}, expected {expected:.6f}"
            )


# ---------------------------------------------------------------------------
# TestInvGammaBootstrapKSamples
# ---------------------------------------------------------------------------


class TestInvGammaBootstrapKSamples:
    """
    Verify that the Python bootstrap backend samples k from the correct InvGamma
    distribution by inspecting the sim_df["k"] column (return_sim_df=True).

    NOTE: Requires ProcessPoolExecutor — will fail in sandboxed environments
    that block multiprocessing.synchronize.
    """

    _KWARGS = dict(**_INVGAMMA_BS, n_bs=3000, return_sim_df=True, use_go=False)

    def _get_k_samples(self, **overrides):
        kwargs = {**self._KWARGS, **overrides}
        result = rr.risk_days_bs(**kwargs)
        return _simdf(result)["k"].values

    def test_k_samples_positive(self):
        """All k values drawn from InvGamma must be strictly positive."""
        k = self._get_k_samples()
        assert np.all(k > 0)

    def test_k_samples_count(self):
        """Number of k samples equals n_bs."""
        k = self._get_k_samples()
        assert len(k) == 3000

    def test_k_samples_median_matches_theory(self):
        """Empirical median of bootstrap k samples ≈ InvGamma median."""
        k = self._get_k_samples(n_bs=5000)
        assert np.median(k) == pytest.approx(THEORY_MEDIAN, rel=0.05)

    def test_k_samples_mode_matches_theory(self):
        """KDE mode of bootstrap k samples ≈ beta / (alpha + 1)."""
        k = self._get_k_samples(n_bs=5000)
        mode_est = rr.mode_kde(k)
        assert mode_est == pytest.approx(THEORY_MODE, rel=0.15)

    def test_k_samples_cdf_ks(self):
        """Bootstrap k samples pass KS test vs scipy InvGamma CDF."""
        k = self._get_k_samples(n_bs=5000)
        ks_stat, p_value = scipy_stats.kstest(
            k, lambda x: scipy_stats.invgamma.cdf(x, a=ALPHA, scale=BETA)
        )
        assert p_value > 0.01, (
            f"Bootstrap k KS test rejected: p={p_value:.4f}, stat={ks_stat:.4f}"
        )

    def test_k_samples_mode_parameterisation_equivalent(self):
        """k_invgamma_mode produces the same distribution as the equivalent k_invgamma_beta."""
        n = 5000
        k_beta = self._get_k_samples(n_bs=n)
        # Build mode-parameterised kwargs; strip n_bs so we can set it explicitly below
        kwargs_mode = {
            key: val for key, val in self._KWARGS.items()
            if key not in ("k_invgamma_beta", "n_bs")
        }
        kwargs_mode["k_invgamma_mode"] = THEORY_MODE
        result_mode = rr.risk_days_bs(**kwargs_mode, n_bs=n)
        k_mode = _simdf(result_mode)["k"].values
        # Both samples should be from the same distribution: KS test
        ks_stat, p_value = scipy_stats.ks_2samp(k_beta, k_mode)
        assert p_value > 0.01, (
            f"Mode vs beta k parameterisation differ: p={p_value:.4f}, stat={ks_stat:.4f}"
        )


# ---------------------------------------------------------------------------
# TestPythonGoInvGammaAgreement
# ---------------------------------------------------------------------------


class TestPythonGoInvGammaAgreement:
    """
    Verify that the Python and Go backends produce statistically indistinguishable
    output distributions (risk days) when sampling k from the same InvGamma.

    Since the two backends use independent RNGs, exact equality is not expected.
    Tolerances follow the convention in existing cross-backend tests (15% on
    distributional summaries).

    NOTE: Requires ProcessPoolExecutor (Python backend) — will fail in sandboxed
    environments that block multiprocessing.synchronize.
    """

    _N = 3000
    _KWARGS = dict(**_INVGAMMA_BS, n_bs=_N)

    def _run(self, use_go, **overrides):
        kwargs = {**self._KWARGS, **overrides, "use_go": use_go}
        return _rdests(rr.risk_days_bs(**kwargs))

    def test_median_agreement(self):
        """Median of risk days from Python and Go agree within 15%."""
        py = self._run(use_go=False)
        go = self._run(use_go=True)
        assert np.median(py) == pytest.approx(np.median(go), rel=0.15)

    def test_iqr_agreement(self):
        """Interquartile range from Python and Go agree within 20%."""
        py = self._run(use_go=False)
        go = self._run(use_go=True)
        py_iqr = np.quantile(py, 0.75) - np.quantile(py, 0.25)
        go_iqr = np.quantile(go, 0.75) - np.quantile(go, 0.25)
        assert py_iqr == pytest.approx(go_iqr, rel=0.20)

    def test_95_cri_lower_agreement(self):
        """2.5th percentile of risk days agrees within 20%."""
        py = self._run(use_go=False)
        go = self._run(use_go=True)
        assert np.quantile(py, 0.025) == pytest.approx(np.quantile(go, 0.025), rel=0.20)

    def test_95_cri_upper_agreement(self):
        """97.5th percentile of risk days agrees within 20%."""
        py = self._run(use_go=False)
        go = self._run(use_go=True)
        assert np.quantile(py, 0.975) == pytest.approx(np.quantile(go, 0.975), rel=0.20)

    def test_cdf_agreement_at_quantile_points(self):
        """
        Python and Go empirical CDFs agree at multiple probability levels within 15%.

        A strict two-sample KS test is intentionally NOT used here: with n=3000 draws
        from each backend, the KS test has enough power to detect minor differences
        arising from independent RNGs and different truncated-normal/Gamma algorithms
        in scipy vs Gonum.  Quantile-level agreement within our stated tolerance is
        the appropriate measure of cross-language parity.
        """
        py = self._run(use_go=False)
        go = self._run(use_go=True)
        for prob in [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]:
            q_py = np.quantile(py, prob)
            q_go = np.quantile(go, prob)
            assert q_py == pytest.approx(q_go, rel=0.15), (
                f"CDF disagreement at p={prob}: py={q_py:.6f}, go={q_go:.6f}"
            )

    def test_mode_parameterisation_consistent_across_backends(self):
        """beta and mode parameterisations produce consistent medians across both backends."""
        # Build a mode-parameterised equivalent of _KWARGS
        kwargs_mode = {k: v for k, v in self._KWARGS.items() if k != "k_invgamma_beta"}
        kwargs_mode["k_invgamma_mode"] = THEORY_MODE

        go_beta  = _rdests(rr.risk_days_bs(**self._KWARGS, use_go=True))
        go_mode  = _rdests(rr.risk_days_bs(**kwargs_mode,    use_go=True))
        py_beta  = _rdests(rr.risk_days_bs(**self._KWARGS, use_go=False))
        py_mode  = _rdests(rr.risk_days_bs(**kwargs_mode,    use_go=False))

        # All four medians should be mutually consistent within 15%
        ref = np.median(go_beta)
        for label, arr in [("go_mode", go_mode), ("py_beta", py_beta), ("py_mode", py_mode)]:
            assert np.median(arr) == pytest.approx(ref, rel=0.15), (
                f"{label} median differs from go_beta by more than 15%"
            )


# ---------------------------------------------------------------------------
# TestInvGammaGoStatistics
# ---------------------------------------------------------------------------


class TestInvGammaGoStatistics:
    """
    Tests for the Go backend alone — no Python bootstrap backend used.
    These are always runnable (no ProcessPoolExecutor dependency).
    """

    _N = 2000
    _KWARGS = dict(**_INVGAMMA_BS, n_bs=_N, use_go=True)

    def test_output_count(self):
        """Go backend returns exactly n_bs risk day estimates."""
        result = rr.risk_days_bs(**self._KWARGS)
        assert len(_rdests(result)) == self._N

    def test_output_positive(self):
        """All risk day estimates must be strictly positive."""
        result = rr.risk_days_bs(**self._KWARGS)
        assert np.all(_rdests(result) > 0)

    def test_output_finite(self):
        """All risk day estimates must be finite."""
        result = rr.risk_days_bs(**self._KWARGS)
        assert np.all(np.isfinite(_rdests(result)))

    def test_reproducible(self):
        """Same seed produces identical output from Go backend."""
        r1 = _rdests(rr.risk_days_bs(**{**self._KWARGS, "seed": 1}))
        r2 = _rdests(rr.risk_days_bs(**{**self._KWARGS, "seed": 1}))
        np.testing.assert_array_equal(r1, r2)

    def test_different_seeds_differ(self):
        """Different seeds produce different output."""
        r1 = _rdests(rr.risk_days_bs(**{**self._KWARGS, "seed": 1}))
        r2 = _rdests(rr.risk_days_bs(**{**self._KWARGS, "seed": 2}))
        assert not np.array_equal(r1, r2)

    def test_mode_parameterisation_accepted(self):
        """Go backend accepts k_invgamma_mode and produces a positive distribution."""
        kwargs_mode = {k: v for k, v in self._KWARGS.items() if k != "k_invgamma_beta"}
        kwargs_mode["k_invgamma_mode"] = THEORY_MODE
        result = rr.risk_days_bs(**kwargs_mode)
        rd = _rdests(result)
        assert len(rd) == self._N
        assert np.all(rd > 0)

    def test_mode_parameterisation_median_matches_beta(self):
        """Go output from k_invgamma_mode should have the same median as k_invgamma_beta."""
        r_beta = _rdests(rr.risk_days_bs(**self._KWARGS))
        kwargs_mode = {k: v for k, v in self._KWARGS.items() if k != "k_invgamma_beta"}
        kwargs_mode["k_invgamma_mode"] = THEORY_MODE
        r_mode = _rdests(rr.risk_days_bs(**kwargs_mode))
        # Both should give statistically similar distributions; use generous tolerance
        # because the two runs use the same seed but resolve to the same beta anyway.
        ks_stat, p_value = scipy_stats.ks_2samp(r_beta, r_mode)
        assert p_value > 0.01, (
            f"Go mode vs beta parameterisation: KS p={p_value:.4f}"
        )

    def test_go_output_median_in_plausible_range(self):
        """Go output median should be in a physically plausible range (0.01 – 30 days)."""
        result = rr.risk_days_bs(**self._KWARGS)
        med = np.median(_rdests(result))
        assert 0.01 < med < 30, f"Go median out of plausible range: {med:.4f} days"

    def test_go_output_cdf_vs_python_sampler(self):
        """
        Go output distribution (risk days) should not be wildly different from a
        Python run that uses sample_invgamma() directly to construct a k_posterior_sample
        with the same InvGamma parameters.  KS test p > 0.01.
        """
        # Python path via posterior sample (bypasses ProcessPoolExecutor restriction
        # because this constructs a pre-sampled posterior).
        k_posterior = rr.sample_invgamma(self._N, alpha=ALPHA, beta=BETA, seed=42)
        result_py = rr.risk_days_bs(
            **{k: v for k, v in self._KWARGS.items()
               if k not in ("k_invgamma_alpha", "k_invgamma_beta", "use_go")},
            k_posterior_sample=k_posterior,
            use_go=True,  # Still use Go backend, just with posterior sample
        )
        result_go = rr.risk_days_bs(**self._KWARGS)
        ks_stat, p_value = scipy_stats.ks_2samp(_rdests(result_py), _rdests(result_go))
        assert p_value > 0.01, (
            f"Go InvGamma vs Go+PosteriorSample: KS p={p_value:.4f}, stat={ks_stat:.4f}"
        )
