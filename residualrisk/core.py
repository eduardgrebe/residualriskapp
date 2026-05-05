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

import math
import statistics

import numpy as np
import pandas as pd
import scipy.stats as stats


def _concentration(C0, doubling_time, t):
    concentration = C0 * 2 ** (t / doubling_time)
    return concentration


def _prob_infectious_copies(n_copies, k):
    prob = 1.0 - math.exp(-k * n_copies)
    return max(0.0, min(1.0, prob))


def _prob_infectious(t, C0, doubling_time, volume_transfused, k, copies_per_virion=2):
    C = _concentration(C0, doubling_time, t)
    n_copies = C * copies_per_virion * volume_transfused
    prob = _prob_infectious_copies(n_copies, k)
    return prob


def _prob_infectious_copies_wc(n_copies):
    if n_copies < 2:
        return 0.0
    elif n_copies >= 2:
        return 1.0


def _prob_infectious_wc(t, C0, doubling_time, volume_transfused, copies_per_virion=2):
    C = _concentration(C0, doubling_time, t)
    n_copies = C * copies_per_virion * volume_transfused
    prob = _prob_infectious_copies_wc(n_copies)
    return prob


def _prob_pos_init(C, doubling_time, pool_size, lod50, lod95_lod50_ratio, z):
    if (not isinstance(pool_size, int)) or pool_size < 1:
        raise Exception("pool_size must be an integer of at least 1")

    # C is in copies copies_per_virion * C when C in virions
    X = z * (math.log10(((C) / (pool_size * lod50))) / math.log10(lod95_lod50_ratio))
    # print(X)
    from scipy.stats import norm

    prob = norm.cdf(X)
    return prob


def _prob_neg_retest(C, doubling_time, pool_size, lod50, lod95_lod50_ratio, retests, z):
    if (not isinstance(pool_size, int)) or pool_size < 1:
        raise Exception("pool_size must be an integer of at least 1")

    if (not isinstance(retests, int)) or retests < 0:
        raise Exception("retests must be a positive integer")
    elif retests == 0:
        return 0
    elif retests >= 1:
        # C is in copies copies_per_virion * C when C in virions
        X = z * (math.log10(((C) / lod50)) / math.log10(lod95_lod50_ratio))
        # print(X)
        from scipy.stats import norm

        prob = (1 - norm.cdf(X)) ** retests
        return prob


def _prob_nondetection(
    t,
    copies_per_virion,
    C0,
    doubling_time,
    pool_size,
    lod50,
    lod95_lod50_ratio,
    retests,
    z=1.6449,
):
    Cv = _concentration(C0, doubling_time, t)
    Cc = copies_per_virion * Cv
    p_pos_init = _prob_pos_init(
        Cc, doubling_time, pool_size, lod50, lod95_lod50_ratio, z
    )
    p_neg_retest = _prob_neg_retest(
        Cc, doubling_time, pool_size, lod50, lod95_lod50_ratio, retests, z
    )
    prob = 1 - p_pos_init * (1 - p_neg_retest)
    return prob


def _prob_infectious_nondetection(
    t,
    copies_per_virion,
    C0,
    doubling_time,
    volume_transfused,
    k,
    pool_size,
    lod50,
    lod95_lod50_ratio,
    retests,
    z=1.6449,
):
    product = _prob_infectious(
        t, C0, doubling_time, volume_transfused, k
    ) * _prob_nondetection(
        t,
        copies_per_virion,
        C0,
        doubling_time,
        pool_size,
        lod50,
        lod95_lod50_ratio,
        retests,
        z,
    )
    return product


def _prob_infectious_nondetection_wc(
    t,
    copies_per_virion,
    C0,
    doubling_time,
    volume_transfused,
    pool_size,
    lod50,
    lod95_lod50_ratio,
    retests,
    z=1.6449,
):
    product = _prob_infectious_wc(
        t, C0, doubling_time, volume_transfused
    ) * _prob_nondetection(
        t,
        copies_per_virion,
        C0,
        doubling_time,
        pool_size,
        lod50,
        lod95_lod50_ratio,
        retests,
        z,
    )
    return product


def _risk_days(
    copies_per_virion,
    C0,
    doubling_time,
    volume_transfused,
    k,
    pool_size,
    lod50,
    lod95_lod50_ratio,
    retests,
    z=1.6449,
    limits=(-100, 500),
):
    # Ideally we would integrate from -np.inf to np.inf, but that causes an
    # overflow error, so we choose safe limits instead
    from scipy.integrate import quad

    rd = quad(
        _prob_infectious_nondetection,
        limits[0],
        limits[1],
        args=(
            copies_per_virion,
            C0,
            doubling_time,
            volume_transfused,
            k,
            pool_size,
            lod50,
            lod95_lod50_ratio,
            retests,
            z,
        ),
    )[0]
    return rd


def get_cpu_core_count() -> int:
    """
    Return the number of logical CPU cores available on this machine.

    Uses `multiprocessing.cpu_count()` which works on Windows, macOS,
    Linux, and most other platforms.  If the call fails for any reason
    (e.g., in a restricted environment), it falls back to 1.
    """
    import multiprocessing
    import os

    try:
        return multiprocessing.cpu_count()
    except NotImplementedError:
        # Fallback: use os.cpu_count() if available
        return os.cpu_count() or 8


def mode_rounded(list, precision=6):
    return stats.mode(np.array(list).round(precision)).mode


def _kde_mode_log(data, n_grid=100_000, cap=50_000):
    """Estimate the mode of a positive, right-skewed distribution via
    KDE on the log scale.

    Applies Silverman's rule for bandwidth selection to log(k), then
    maps the density back to the original scale via the change-of-
    variables  f(k) = f_logk(log k) / k  and finds the maximum.

    This is the methodologically correct approach for a log-
    approximately-normal quantity such as the k posterior.

    Parameters
    ----------
    data : array-like
        Positive-valued posterior samples.
    n_grid : int
        Number of log-spaced grid points for density evaluation
        (default 100 000).
    cap : int or None
        Maximum number of samples to use.  If *data* exceeds *cap*,
        a random subset of size *cap* is drawn before fitting.
        Pass ``None`` to use all samples regardless of size.
        Default 50 000.

    Returns
    -------
    float
        Mode estimate on the original scale.
    """
    import warnings
    from scipy.stats import gaussian_kde

    data = np.asarray(data, dtype=float)
    if np.any(data <= 0):
        raise ValueError("All values must be positive for log-scale KDE.")

    # Cap to avoid O(n_data × n_grid) blow-up with large samples
    if cap is not None and len(data) > cap:
        rng = np.random.default_rng(42)
        data = rng.choice(data, size=cap, replace=False)

    log_data = np.log(data)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kde = gaussian_kde(log_data, bw_method="silverman")

    grid = np.logspace(np.log10(data.min()), np.log10(data.max()), n_grid)
    # Density on original scale: f(k) = f_logk(log k) / k
    fk = kde(np.log(grid)) / grid
    mode = grid[np.argmax(fk)]
    return float(mode)


def mode_kde(data, n_grid=100_000, cap=50_000):
    """Public wrapper for _kde_mode_log — estimate the mode of a
    positive posterior distribution via KDE on the log scale.

    See _kde_mode_log for full documentation.
    """
    return _kde_mode_log(data, n_grid=n_grid, cap=cap)


def sample_invgamma(n, alpha, beta=None, mode=None, seed=None):
    """Sample from an Inverse Gamma distribution.

    Supports two parameterisations:

    1. **alpha + beta** (direct)::

          sample_invgamma(n, alpha=2.0, beta=0.002019)

    2. **alpha + mode** (beta auto-calculated)::

          sample_invgamma(n, alpha=2.0, mode=0.000673)

       Beta is computed as ``mode * (alpha + 1)`` so that the resulting
       InvGamma(alpha, beta) has its mode at the specified value.

    Parameters
    ----------
    n : int
        Number of samples.
    alpha : float
        Shape parameter (must be > 0).
    beta : float, optional
        Scale parameter.  Maps to ``scipy.stats.invgamma(a=alpha, scale=beta)``.
        Exactly one of *beta* or *mode* must be provided.
    mode : float, optional
        Target mode.  Beta is computed as ``mode * (alpha + 1)``.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    np.ndarray
        Array of *n* samples from InvGamma(alpha, beta).

    Raises
    ------
    ValueError
        If neither or both of *beta* and *mode* are provided.
    """
    if beta is None and mode is None:
        raise ValueError(
            "Exactly one of 'beta' or 'mode' must be provided."
        )
    if beta is not None and mode is not None:
        raise ValueError(
            "Provide 'beta' or 'mode', not both."
        )
    if mode is not None:
        beta = mode * (alpha + 1)
    rng = np.random.default_rng(seed)
    from scipy.stats import invgamma
    return invgamma.rvs(alpha, scale=beta, size=n, random_state=rng)


def sample_lnmix(n, w, mu1, sigma1, mu2, sigma2, seed=None):
    """Sample from a two-component lognormal mixture distribution.

    Each sample is drawn from component 1 (LN(mu1, sigma1)) with probability w,
    or from component 2 (LN(mu2, sigma2)) with probability 1-w.

    Parameters follow the numpy/scipy lognormal convention:
    - ``mu`` is the mean of the underlying normal (log-scale mean)
    - ``sigma`` is the std of the underlying normal (log-scale std)

    This corresponds to ``scipy.stats.lognorm(s=sigma, scale=np.exp(mu))``.

    Parameters
    ----------
    n : int
        Number of samples.
    w : float
        Weight of component 1; must be in [0, 1].
    mu1 : float
        Log-scale mean of component 1 (e.g. -7.2403 for human posterior fit).
    sigma1 : float
        Log-scale std of component 1 (e.g. 0.3241 for human posterior fit).
    mu2 : float
        Log-scale mean of component 2 (e.g. -3.7423 for animal posterior fit).
    sigma2 : float
        Log-scale std of component 2 (e.g. 0.5258 for animal posterior fit).
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    np.ndarray
        Array of *n* positive samples from the mixture distribution.

    Raises
    ------
    ValueError
        If w is not in [0, 1].

    Examples
    --------
    Default 90/10 human-weighted mixture (Recommendation B)::

        samples = sample_lnmix(10000, w=0.90,
                                mu1=-7.2403, sigma1=0.3241,
                                mu2=-3.7423, sigma2=0.5258)
    """
    if not (0.0 <= w <= 1.0):
        raise ValueError(f"w must be in [0, 1], got {w}")
    rng = np.random.default_rng(seed)
    component = rng.random(n) < w
    comp1 = rng.lognormal(mean=mu1, sigma=sigma1, size=n)
    comp2 = rng.lognormal(mean=mu2, sigma=sigma2, size=n)
    return np.where(component, comp1, comp2)


def _risk_days_bs_python(
    k,
    doubling_time,
    doubling_time_norm_sd,
    lod50,
    lod50_sd,
    lod95_lod50_ratio,
    volume_transfused,
    volume_transfused_range,
    pool_size,
    retests,
    C0=0.00025,
    copies_per_virion=2,
    alpha=0.05,
    z=1.6449,
    k_posterior_sample=None,
    k_gamma_shape=None,
    k_gamma_scale=None,
    k_invgamma_alpha=None,
    k_invgamma_beta=None,
    k_invgamma_mode=None,
    k_lnmix_w=None,
    k_lnmix_mu1=None,
    k_lnmix_sigma1=None,
    k_lnmix_mu2=None,
    k_lnmix_sigma2=None,
    n_bs=10000,
    seed=126887,
    threads=get_cpu_core_count() - 1,
    point_estimate="primary parameters",
    mode_precision=2,
    progress=None,
    return_sim_df=False,
):
    if n_bs <= 0:
        raise ValueError("n_bs must be greater than zero to perform simulations.")

    np.random.seed(seed)
    if k_posterior_sample is not None:
        ks = np.random.choice(k_posterior_sample, size=n_bs, replace=True)
    elif (
        k_posterior_sample is None
        and k_gamma_shape is not None
        and k_gamma_scale is not None
    ):
        ks = np.random.gamma(k_gamma_shape, k_gamma_scale, n_bs)
    elif k_invgamma_alpha is not None:
        _beta = k_invgamma_beta
        if _beta is None:
            if k_invgamma_mode is not None:
                _beta = k_invgamma_mode * (k_invgamma_alpha + 1)
            else:
                raise ValueError(
                    "k_invgamma_alpha requires k_invgamma_beta or k_invgamma_mode"
                )
        # Uses the legacy numpy global state set above for reproducibility.
        ks = stats.invgamma.rvs(k_invgamma_alpha, scale=_beta, size=n_bs)
    elif k_lnmix_w is not None:
        if any(p is None for p in [k_lnmix_mu1, k_lnmix_sigma1, k_lnmix_mu2, k_lnmix_sigma2]):
            raise ValueError(
                "All lnmix parameters (k_lnmix_w, mu1, sigma1, mu2, sigma2) must be provided together."
            )
        ks = sample_lnmix(n_bs, k_lnmix_w, k_lnmix_mu1, k_lnmix_sigma1,
                           k_lnmix_mu2, k_lnmix_sigma2, seed=seed)
    else:
        raise ValueError(
            "k_posterior_sample and k_gamma parameters must not both be 'None'."
        )
    doubling_times = stats.truncnorm.rvs(
        0, np.inf, doubling_time, doubling_time_norm_sd, n_bs
    )
    lod50s = stats.truncnorm.rvs(0, np.inf, lod50, lod50_sd, n_bs)
    volumes_transfused = np.random.uniform(
        volume_transfused_range[0], volume_transfused_range[1], n_bs
    )
    print("Starting parallel risk days calculation on ", threads, " cores...")
    rdests = []
    args_list = [
        (
            copies_per_virion,
            C0,
            doubling_times[i],
            volumes_transfused[i],
            ks[i],
            pool_size,
            lod50s[i],
            lod95_lod50_ratio,
            retests,
            z,
            (-100, 500),
        )
        for i in range(n_bs)
    ]
    from concurrent.futures import ProcessPoolExecutor, as_completed

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(_risk_days, *args) for args in args_list]
        completed_count = 0
        for future in as_completed(futures):
            rdests.append(future.result())
            completed_count += 1
            # Update progress bar only when percentage changes (reduces warnings from multiprocessing)
            # Note: Streamlit warnings about missing ScriptRunContext are expected and harmless when using ProcessPoolExecutor
            if progress is not None:
                new_percent = int((completed_count / n_bs) * 100)
                if completed_count == 1 or new_percent > getattr(
                    progress, "_last_percent", 0
                ):
                    progress._last_percent = new_percent
                    progress_percentage = completed_count / n_bs
                    progress.progress(
                        progress_percentage,
                        text=f"Completed {completed_count}/{n_bs}...",
                    )
    rd_range = [np.min(rdests), np.max(rdests)]
    rd_cri = np.quantile(rdests, (alpha / 2, 1 - alpha / 2))
    if return_sim_df:
        sim_df = pd.DataFrame(
            args_list,
            columns=[
                "copies_per_virion",
                "C0",
                "doubling_time",
                "volume_transfused",
                "k",
                "pool_size",
                "lod50",
                "lod95_lod50_ratio",
                "retests",
                "z",
                "limits",
            ],
        )
        sim_df["lod95"] = (
            sim_df["lod50"] * sim_df["lod95_lod50_ratio"]
        )  # Convert ratio to actual lod95
        sim_df["iwp"] = rdests
        sim_df["random_seed"] = np.repeat(seed, n_bs)

    if point_estimate == "primary parameters":
        rd_pe = _risk_days(
            copies_per_virion,
            C0,
            doubling_time,
            volume_transfused,
            k,
            pool_size,
            lod50,
            lod95_lod50_ratio,
            retests,
        )
    elif point_estimate == "median":
        rd_pe = statistics.median(rdests)
    elif point_estimate == "mean":
        rd_pe = statistics.mean(rdests)
    elif point_estimate == "mode":
        rd_pe = _kde_mode_log(rdests, cap=None)
    else:
        rd_pe = None

    if return_sim_df:
        return (rd_pe, rd_cri, rd_range, rdests, sim_df)
    else:
        return (rd_pe, rd_cri, rd_range, rdests, None)


def risk_days_bs(
    k,
    doubling_time,
    doubling_time_norm_sd,
    lod50,
    lod50_sd,
    lod95_lod50_ratio,
    volume_transfused,
    volume_transfused_range,
    pool_size,
    retests,
    C0=0.00025,
    copies_per_virion=2,
    alpha=0.05,
    z=1.6449,
    k_posterior_sample=None,
    k_gamma_shape=None,
    k_gamma_scale=None,
    k_invgamma_alpha=None,
    k_invgamma_beta=None,
    k_invgamma_mode=None,
    k_lnmix_w=None,
    k_lnmix_mu1=None,
    k_lnmix_sigma1=None,
    k_lnmix_mu2=None,
    k_lnmix_sigma2=None,
    n_bs=10000,
    seed=126887,
    threads=get_cpu_core_count() - 1,
    point_estimate="primary parameters",
    mode_precision=2,
    progress=None,
    return_sim_df=False,
    use_go=False,
):
    """
    Risk days bootstrap calculation with optional Go acceleration.

    Parameters
    ----------
    use_go : bool, default=False
        If True, uses the high-performance Go implementation.
        If False (default), uses the Python implementation.
        If Go implementation fails, automatically falls back to Python.

    All other parameters are passed through to the underlying implementation.
    """
    if use_go:
        try:
            from ._go import risk_days_bs_go

            return risk_days_bs_go(
                k,
                doubling_time,
                doubling_time_norm_sd,
                lod50,
                lod50_sd,
                lod95_lod50_ratio,
                volume_transfused,
                volume_transfused_range,
                pool_size,
                retests,
                C0,
                copies_per_virion,
                alpha,
                z,
                k_posterior_sample,
                k_gamma_shape,
                k_gamma_scale,
                k_invgamma_alpha,
                k_invgamma_beta,
                k_invgamma_mode,
                k_lnmix_w,
                k_lnmix_mu1,
                k_lnmix_sigma1,
                k_lnmix_mu2,
                k_lnmix_sigma2,
                n_bs,
                seed,
                threads,
                point_estimate,
                mode_precision,
                progress,
                return_sim_df,
            )
        except Exception as e:
            print(f"Warning: Go implementation failed ({e}), falling back to Python")
            # Fall through to Python implementation

    # Use Python implementation
    return _risk_days_bs_python(
        k,
        doubling_time,
        doubling_time_norm_sd,
        lod50,
        lod50_sd,
        lod95_lod50_ratio,
        volume_transfused,
        volume_transfused_range,
        pool_size,
        retests,
        C0,
        copies_per_virion,
        alpha,
        z,
        k_posterior_sample,
        k_gamma_shape,
        k_gamma_scale,
        k_invgamma_alpha,
        k_invgamma_beta,
        k_invgamma_mode,
        k_lnmix_w,
        k_lnmix_mu1,
        k_lnmix_sigma1,
        k_lnmix_mu2,
        k_lnmix_sigma2,
        n_bs,
        seed,
        threads,
        point_estimate,
        mode_precision,
        progress,
        return_sim_df,
    )


def iwp_from_lookback_data(
    n_transmissions,
    intervals,
    negative_diagnostic_delay,
    positive_diagnostic_delay,
    alpha=0.05,
    n_bs=10000,
    seed=126887,
):
    """
    Estimate the infectious window period (IWP) from lookback investigation data.

    Models transfusion transmissions as a Poisson process.  Each seroconverting
    donor with a prior donation contributes 1/adjusted_IDI to the total exposure,
    where adjusted_IDI = IDI + negative_diagnostic_delay - positive_diagnostic_delay.

    Parameters
    ----------
    n_transmissions : int
        Number of confirmed transfusion transmissions from prior donations.
    intervals : array-like
        Raw inter-donation intervals (days) for each seroconverting donor.
    negative_diagnostic_delay : float
        Diagnostic delay (days) of the most sensitive test applied at the
        prior (negative) donation.
    positive_diagnostic_delay : float
        Diagnostic delay (days) of the least sensitive test that was positive
        at the seroconversion donation.
    alpha : float
        Significance level for the confidence interval (default 0.05 → 95% CI).
    n_bs : int
        Number of Gamma posterior samples to draw for use with
        residual_risk_rd() (default 10000).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    iwp_pe : float
        Point estimate of the IWP (days).  Zero when n_transmissions == 0.
    iwp_ci : tuple of float
        (lower, upper) frequentist CI derived from the chi-squared distribution.
    iwp_samples : np.ndarray
        Samples from Gamma(n_transmissions + 0.5, scale=1/total_exposure),
        suitable for passing as iwp_bs to residual_risk_rd().
    """
    from scipy import stats as scipy_stats

    adjusted = [
        x + negative_diagnostic_delay - positive_diagnostic_delay for x in intervals
    ]
    if any(adj <= 0 for adj in adjusted):
        raise ValueError(
            "All adjusted IDIs (IDI + negative_diagnostic_delay - "
            "positive_diagnostic_delay) must be positive. Check diagnostic "
            "delay parameters."
        )

    total_exposure = sum(1.0 / adj for adj in adjusted)

    iwp_pe = n_transmissions / total_exposure if n_transmissions > 0 else 0.0

    if n_transmissions > 0:
        iwp_ci_lb = (
            scipy_stats.chi2.ppf(alpha / 2, df=2 * n_transmissions) / 2 / total_exposure
        )
    else:
        iwp_ci_lb = 0.0
    iwp_ci_ub = (
        scipy_stats.chi2.ppf(1.0 - alpha / 2, df=2 * (n_transmissions + 1))
        / 2
        / total_exposure
    )

    # Gamma(n + 0.5, 1/T) is the posterior under Jeffreys prior on the
    # Poisson rate.  Valid for n_transmissions == 0 (shape = 0.5).
    np.random.seed(seed)
    iwp_samples = np.random.gamma(
        shape=n_transmissions + 0.5,
        scale=1.0 / total_exposure,
        size=n_bs,
    )

    return iwp_pe, (iwp_ci_lb, iwp_ci_ub), iwp_samples


def residual_risk_rd(
    iwp_pe,
    iwp_bs,
    incidence,
    incidence_norm_sd,
    per=1e6,
    seed=126887,
    alpha=0.05,
    one_in_x=False,
):
    # Validate inputs to prevent division by zero
    if incidence <= 0:
        raise ValueError(f"incidence must be positive, got {incidence}")
    if iwp_pe <= 0:
        raise ValueError(f"iwp_pe must be positive, got {iwp_pe}")

    if one_in_x:
        rr_pe = 1 / (incidence * iwp_pe / 365.25)
    else:
        rr_pe = incidence * iwp_pe / 365.25 * per
    n_bs = len(iwp_bs)
    np.random.seed(seed)
    incidence_draws = stats.truncnorm.rvs(0, np.inf, incidence, incidence_norm_sd, n_bs)
    rr = []
    for i in range(n_bs):
        # Skip iterations where the product would be zero or negative
        product = incidence_draws[i] * iwp_bs[i] / 365.25
        if product <= 0:
            continue
        if one_in_x:
            rr.append(1 / product)
        else:
            rr.append(product * per)
    rr_cri = np.quantile(rr, (alpha / 2, 1 - alpha / 2))
    rr_sd = np.std(rr)
    return (rr_pe, rr_cri, rr_sd)
