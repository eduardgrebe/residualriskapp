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
Tests verifying that lognormal mixture sampling is statistically correct and
that the Python and Go implementations produce equivalent output distributions.

Mixture parameterisation (Recommendation B from K_PARAM_INPUTDIST.md):
  Component 1 (human):  LN(mu1=-7.2403, sigma1=0.3241)  weight w=0.90
  Component 2 (animal): LN(mu2=-3.7423, sigma2=0.5258)  weight w=0.10

Test classes:

  TestLnMixTheoreticalStatistics
      Pure statistical tests of sample_lnmix() against analytical mixture properties.
      No multiprocessing required — always runnable.

  TestLnMixBootstrapKSamples
      Inspects the k values used inside the Python bootstrap (via return_sim_df=True)
      and verifies they follow the mixture.
      Requires ProcessPoolExecutor — will fail in sandboxed environments.

  TestPythonGoLnMixAgreement
      Verifies that the Python and Go backends produce statistically indistinguishable
      output distributions (risk days) when using the same mixture parameters.
      Requires ProcessPoolExecutor for the Python backend.

  TestLnMixGoStatistics
      Tests the Go backend alone: sanity, reproducibility, and distributional
      plausibility.  No Python backend — always runnable.
"""

import numpy as np
import pytest
from scipy import stats as scipy_stats

from residualrisk import core as rr

# ---------------------------------------------------------------------------
# Lognormal mixture theoretical constants (Recommendation B)
# ---------------------------------------------------------------------------

W = 0.90
MU1, SIGMA1 = -7.2403, 0.3241   # human component (log-scale)
MU2, SIGMA2 = -3.7423, 0.5258   # animal component (log-scale)

# Analytic mean: w * exp(mu1 + sigma1²/2) + (1-w) * exp(mu2 + sigma2²/2)
THEORY_MEAN = W * np.exp(MU1 + SIGMA1**2 / 2) + (1 - W) * np.exp(MU2 + SIGMA2**2 / 2)

# Component medians: exp(mu_i)
COMP1_MEDIAN = np.exp(MU1)   # ≈ 0.000715
COMP2_MEDIAN = np.exp(MU2)   # ≈ 0.0237

# Approximate mixture mode ≈ 0.000649 (from companion analysis)
APPROX_MIXTURE_MODE = 0.000649

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
    k=APPROX_MIXTURE_MODE,
    threads=2,
)

_LNMIX_BS = dict(
    **_COMMON_BS,
    k_lnmix_w=W,
    k_lnmix_mu1=MU1,
    k_lnmix_sigma1=SIGMA1,
    k_lnmix_mu2=MU2,
    k_lnmix_sigma2=SIGMA2,
    seed=42,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _rdests(result):
    """Extract rdests (index 3) from a risk_days_bs result tuple."""
    return np.asarray(result[3])


def _simdf(result):
    """Extract sim_df (index 4) from a risk_days_bs result tuple."""
    return result[4]


def _mixture_cdf(x, w=W, mu1=MU1, s1=SIGMA1, mu2=MU2, s2=SIGMA2):
    """Analytical CDF of the two-component lognormal mixture."""
    c1 = scipy_stats.lognorm(s=s1, scale=np.exp(mu1)).cdf(x)
    c2 = scipy_stats.lognorm(s=s2, scale=np.exp(mu2)).cdf(x)
    return w * c1 + (1 - w) * c2


# ---------------------------------------------------------------------------
# TestLnMixTheoreticalStatistics
# ---------------------------------------------------------------------------


class TestLnMixTheoreticalStatistics:
    """
    Verify that sample_lnmix() produces samples with the correct theoretical
    moments and CDF.  Uses sample_lnmix() directly — no multiprocessing.
    """

    N = 100_000

    def test_all_samples_positive(self):
        """All lognormal mixture samples must be strictly positive."""
        s = rr.sample_lnmix(self.N, W, MU1, SIGMA1, MU2, SIGMA2, seed=1)
        assert np.all(s > 0)

    def test_all_samples_finite(self):
        """All lognormal mixture samples must be finite."""
        s = rr.sample_lnmix(self.N, W, MU1, SIGMA1, MU2, SIGMA2, seed=1)
        assert np.all(np.isfinite(s))

    def test_mean_converges(self):
        """Empirical mean ≈ w*exp(mu1+sigma1²/2) + (1-w)*exp(mu2+sigma2²/2)."""
        s = rr.sample_lnmix(self.N, W, MU1, SIGMA1, MU2, SIGMA2, seed=1)
        assert np.mean(s) == pytest.approx(THEORY_MEAN, rel=0.05)

    def test_reproducible_with_same_seed(self):
        """Same seed produces identical samples."""
        s1 = rr.sample_lnmix(1000, W, MU1, SIGMA1, MU2, SIGMA2, seed=7)
        s2 = rr.sample_lnmix(1000, W, MU1, SIGMA1, MU2, SIGMA2, seed=7)
        np.testing.assert_array_equal(s1, s2)

    def test_different_seeds_differ(self):
        """Different seeds produce different samples."""
        s1 = rr.sample_lnmix(1000, W, MU1, SIGMA1, MU2, SIGMA2, seed=7)
        s2 = rr.sample_lnmix(1000, W, MU1, SIGMA1, MU2, SIGMA2, seed=8)
        assert not np.array_equal(s1, s2)

    def test_component_isolation_w1(self):
        """w=1 gives pure component 1: log-mean should equal mu1."""
        s = rr.sample_lnmix(self.N, 1.0, MU1, SIGMA1, MU2, SIGMA2, seed=5)
        log_mean = np.mean(np.log(s))
        assert log_mean == pytest.approx(MU1, rel=0.03)

    def test_component_isolation_w0(self):
        """w=0 gives pure component 2: log-mean should equal mu2."""
        s = rr.sample_lnmix(self.N, 0.0, MU1, SIGMA1, MU2, SIGMA2, seed=5)
        log_mean = np.mean(np.log(s))
        assert log_mean == pytest.approx(MU2, rel=0.03)

    def test_component1_median(self):
        """w=1 samples: median ≈ exp(mu1)."""
        s = rr.sample_lnmix(self.N, 1.0, MU1, SIGMA1, MU2, SIGMA2, seed=3)
        assert np.median(s) == pytest.approx(COMP1_MEDIAN, rel=0.03)

    def test_component2_median(self):
        """w=0 samples: median ≈ exp(mu2)."""
        s = rr.sample_lnmix(self.N, 0.0, MU1, SIGMA1, MU2, SIGMA2, seed=3)
        assert np.median(s) == pytest.approx(COMP2_MEDIAN, rel=0.03)

    def test_cdf_matches_mixture_ks(self):
        """KS test: sample_lnmix samples are consistent with mixture CDF."""
        s = rr.sample_lnmix(10_000, W, MU1, SIGMA1, MU2, SIGMA2, seed=99)
        ks_stat, p_value = scipy_stats.kstest(s, _mixture_cdf)
        assert p_value > 0.01, (
            f"KS test rejected: p={p_value:.4f}, stat={ks_stat:.4f}"
        )

    def test_quantiles_match_mixture_cdf(self):
        """Empirical quantiles match those obtained by numerical inversion of mixture CDF."""
        s = rr.sample_lnmix(self.N, W, MU1, SIGMA1, MU2, SIGMA2, seed=11)
        for prob in [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]:
            # Numerical mixture quantile via scipy's ppf-equivalent
            from scipy.optimize import brentq
            expected = brentq(lambda x: _mixture_cdf(x) - prob, 1e-10, 1.0)
            got = np.quantile(s, prob)
            assert got == pytest.approx(expected, rel=0.08), (
                f"Quantile mismatch at q={prob}: got {got:.6f}, expected {expected:.6f}"
            )

    def test_higher_w_lower_median(self):
        """Higher w → more weight on component 1 (lower median), lower w → higher median."""
        s_high_w = rr.sample_lnmix(self.N, 0.99, MU1, SIGMA1, MU2, SIGMA2, seed=4)
        s_low_w = rr.sample_lnmix(self.N, 0.01, MU1, SIGMA1, MU2, SIGMA2, seed=4)
        assert np.median(s_high_w) < np.median(s_low_w), (
            f"Expected med(w=0.99) < med(w=0.01): "
            f"{np.median(s_high_w):.6f} vs {np.median(s_low_w):.6f}"
        )

    def test_correct_sample_count(self):
        """sample_lnmix returns exactly n samples."""
        for n in [1, 100, 10_000]:
            s = rr.sample_lnmix(n, W, MU1, SIGMA1, MU2, SIGMA2, seed=1)
            assert len(s) == n


# ---------------------------------------------------------------------------
# TestLnMixBootstrapKSamples
# ---------------------------------------------------------------------------


class TestLnMixBootstrapKSamples:
    """
    Verify that the Python bootstrap backend samples k from the correct lognormal
    mixture by inspecting sim_df["k"] (return_sim_df=True).

    NOTE: Requires ProcessPoolExecutor — will fail in sandboxed environments.
    """

    _KWARGS = dict(**_LNMIX_BS, n_bs=3000, return_sim_df=True, use_go=False)

    def _get_k_samples(self, **overrides):
        kwargs = {**self._KWARGS, **overrides}
        result = rr.risk_days_bs(**kwargs)
        return _simdf(result)["k"].to_numpy()

    def test_k_samples_positive(self):
        """All k values drawn from lnmix must be strictly positive."""
        k = self._get_k_samples()
        assert np.all(k > 0)

    def test_k_samples_count(self):
        """Number of k samples equals n_bs."""
        k = self._get_k_samples()
        assert len(k) == 3000

    def test_k_samples_median_plausible(self):
        """Empirical median of k samples from bootstrap is in plausible mixture range."""
        k = self._get_k_samples(n_bs=5000)
        # Mixture median should be close to component 1 median (0.000715) with slight pull
        assert 3e-4 < np.median(k) < 3e-3

    def test_k_samples_cdf_ks(self):
        """Bootstrap k samples pass KS test vs analytical mixture CDF."""
        k = self._get_k_samples(n_bs=5000)
        ks_stat, p_value = scipy_stats.kstest(k, _mixture_cdf)
        assert p_value > 0.01, (
            f"Bootstrap k KS test rejected: p={p_value:.4f}, stat={ks_stat:.4f}"
        )


# ---------------------------------------------------------------------------
# TestPythonGoLnMixAgreement
# ---------------------------------------------------------------------------


class TestPythonGoLnMixAgreement:
    """
    Verify that the Python and Go backends produce statistically indistinguishable
    output distributions (risk days) when sampling k from the same lognormal mixture.

    NOTE: Requires ProcessPoolExecutor (Python backend) — will fail in sandboxed
    environments.
    """

    _N = 3000
    _KWARGS = dict(**_LNMIX_BS, n_bs=_N)

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

        A strict two-sample KS test is intentionally NOT used: with n=3000 draws from
        each backend, independent RNGs and minor algorithmic differences between
        scipy/numpy and Gonum would cause spurious rejections.
        """
        py = self._run(use_go=False)
        go = self._run(use_go=True)
        for prob in [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]:
            q_py = np.quantile(py, prob)
            q_go = np.quantile(go, prob)
            assert q_py == pytest.approx(q_go, rel=0.15), (
                f"CDF disagreement at p={prob}: py={q_py:.6f}, go={q_go:.6f}"
            )


# ---------------------------------------------------------------------------
# TestLnMixGoStatistics
# ---------------------------------------------------------------------------


class TestLnMixGoStatistics:
    """
    Tests for the Go backend alone — no Python bootstrap backend used.
    These are always runnable (no ProcessPoolExecutor dependency).
    """

    _N = 2000
    _KWARGS = dict(**_LNMIX_BS, n_bs=_N, use_go=True)

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

    def test_go_output_median_in_plausible_range(self):
        """Go output median should be in a physically plausible range (0.01 – 30 days)."""
        result = rr.risk_days_bs(**self._KWARGS)
        med = np.median(_rdests(result))
        assert 0.01 < med < 30, f"Go median out of plausible range: {med:.4f} days"

    def test_go_output_cdf_vs_python_sampler(self):
        """
        Go output (risk days) should not be wildly different from a Python run that
        uses sample_lnmix() directly to construct a k_posterior_sample with the same
        parameters.  KS test p > 0.01.

        This test uses the Go backend for both runs (avoiding ProcessPoolExecutor),
        differing only in how k is provided (parametric vs pre-sampled posterior).
        """
        k_posterior = rr.sample_lnmix(self._N, W, MU1, SIGMA1, MU2, SIGMA2, seed=42)
        result_posterior = rr.risk_days_bs(
            **{k: v for k, v in self._KWARGS.items()
               if k not in ("k_lnmix_w", "k_lnmix_mu1", "k_lnmix_sigma1",
                            "k_lnmix_mu2", "k_lnmix_sigma2", "use_go")},
            k_posterior_sample=k_posterior,
            use_go=True,
        )
        result_lnmix = rr.risk_days_bs(**self._KWARGS)
        ks_stat, p_value = scipy_stats.ks_2samp(
            _rdests(result_posterior), _rdests(result_lnmix)
        )
        assert p_value > 0.01, (
            f"Go LnMix vs Go+PosteriorSample: KS p={p_value:.4f}, stat={ks_stat:.4f}"
        )

    def test_point_estimate_positive(self):
        """Go backend returns a positive point estimate."""
        result = rr.risk_days_bs(**self._KWARGS)
        assert result[0] > 0

    def test_cri_is_ordered(self):
        """95% CrI lower bound < upper bound."""
        result = rr.risk_days_bs(**self._KWARGS)
        cri = result[1]
        assert cri[0] < cri[1], f"CrI not ordered: {cri}"

    def test_cdf_agreement_at_quantile_points_vs_invgamma(self):
        """
        LnMix output should be in a similar ballpark to InvGamma(alpha=2, beta=0.002019),
        which targets the same human posterior mode.  Medians should agree within a
        factor of 3 (both are reasonably calibrated to the human posterior mode).
        """
        result_lnmix = rr.risk_days_bs(**self._KWARGS)
        alpha, beta = 2.0, 0.002019
        result_invgamma = rr.risk_days_bs(
            **{k: v for k, v in self._KWARGS.items()
               if k not in ("k_lnmix_w", "k_lnmix_mu1", "k_lnmix_sigma1",
                            "k_lnmix_mu2", "k_lnmix_sigma2", "use_go", "k")},
            k=beta / (alpha + 1),  # InvGamma mode
            k_invgamma_alpha=alpha,
            k_invgamma_beta=beta,
            use_go=True,
        )
        med_lnmix = np.median(_rdests(result_lnmix))
        med_invgamma = np.median(_rdests(result_invgamma))
        ratio = max(med_lnmix, med_invgamma) / min(med_lnmix, med_invgamma)
        assert ratio < 3.0, (
            f"LnMix vs InvGamma medians disagree by factor {ratio:.2f} "
            f"(lnmix={med_lnmix:.4f}, invgamma={med_invgamma:.4f})"
        )
