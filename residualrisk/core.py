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
from functools import lru_cache

import numpy as np
import polars as pl
import scipy.stats as stats

# Cap the viral-growth exponent to avoid float OverflowError when t is large
# relative to a small doubling_time. 2**700 (~5e210) is well within the float64
# range, yet large enough that the modelled viral load saturates every
# downstream detection probability — so the risk-days integrand is zero there
# regardless of the exact (capped) value. This only ever triggers far out in
# the integration tail; it never affects standard scenarios.
_MAX_GROWTH_EXP2 = 700.0


def _concentration(C0, doubling_time, t):
    exponent = t / doubling_time
    if exponent > _MAX_GROWTH_EXP2:
        exponent = _MAX_GROWTH_EXP2
    return C0 * 2.0**exponent


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


@lru_cache(maxsize=None)
def _gauss_legendre_rule(n):
    """Cached n-point Gauss-Legendre nodes and weights on [-1, 1].

    Matches gonum's quad.Fixed (used by the Go backend), so the Python and Go
    "gauss-legendre" results agree to machine precision.
    """
    return np.polynomial.legendre.leggauss(n)


def _integrate_gauss_legendre(func, a, b, args, n=1000):
    """Fixed n-point Gauss-Legendre quadrature of func over [a, b]."""
    nodes, weights = _gauss_legendre_rule(n)
    x = 0.5 * (b - a) * nodes + 0.5 * (a + b)
    y = np.fromiter((func(xi, *args) for xi in x), dtype=float, count=x.size)
    return 0.5 * (b - a) * float(np.dot(weights, y))


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
    integration_method="gauss-legendre",
):
    # Ideally we would integrate from -np.inf to np.inf, but that causes an
    # overflow error, so we choose safe limits instead.
    args = (
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
    )
    if integration_method == "gauss-legendre":
        # Fixed 1000-point Gauss-Legendre, matching the Go backend. The default:
        # robust to narrow / compact-support integrands where adaptive quad can
        # silently miss the active window.
        rd = _integrate_gauss_legendre(
            _prob_infectious_nondetection, limits[0], limits[1], args
        )
    elif integration_method == "quad":
        # Adaptive Gauss-Kronrod (scipy). Retained so prior analyses computed
        # with quad can be reproduced exactly via the Python API.
        from scipy.integrate import quad

        rd = quad(
            _prob_infectious_nondetection, limits[0], limits[1], args=args
        )[0]
    else:
        raise ValueError(
            "integration_method must be 'gauss-legendre' or 'quad', "
            f"got {integration_method!r}"
        )
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


def _kde_mode_log(data, n_grid=5_000, cap=50_000):
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


def mode_kde(data, n_grid=5_000, cap=50_000):
    """Public wrapper for _kde_mode_log — estimate the mode of a
    positive posterior distribution via KDE on the log scale.

    See _kde_mode_log for full documentation.
    """
    return _kde_mode_log(data, n_grid=n_grid, cap=cap)


def _sample_k(
    n_bs,
    seed,
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
):
    """Sample *n_bs* values of *k* from the specified distribution.

    Dispatches in priority order:
    1. ``k_posterior_sample`` — resample with replacement from an array
    2. ``k_gamma_shape`` / ``k_gamma_scale`` — legacy Gamma (deprecated)
    3. ``k_invgamma_alpha`` + (``k_invgamma_beta`` or ``k_invgamma_mode``) — Inverse Gamma
    4. ``k_lnmix_*`` — two-component lognormal mixture

    Assumes ``np.random.seed(seed)`` has already been called by the caller
    (legacy global-state RNG contract used throughout the bootstrap functions).
    """
    if k_posterior_sample is not None:
        return np.random.choice(k_posterior_sample, size=n_bs, replace=True)
    elif k_gamma_shape is not None and k_gamma_scale is not None:
        return np.random.gamma(k_gamma_shape, k_gamma_scale, n_bs)
    elif k_invgamma_alpha is not None:
        _beta = k_invgamma_beta
        if _beta is None:
            if k_invgamma_mode is not None:
                _beta = k_invgamma_mode * (k_invgamma_alpha + 1)
            else:
                raise ValueError(
                    "k_invgamma_alpha requires k_invgamma_beta or k_invgamma_mode"
                )
        return stats.invgamma.rvs(k_invgamma_alpha, scale=_beta, size=n_bs)
    elif k_lnmix_w is not None:
        if any(p is None for p in [k_lnmix_mu1, k_lnmix_sigma1, k_lnmix_mu2, k_lnmix_sigma2]):
            raise ValueError(
                "All lnmix parameters (k_lnmix_w, mu1, sigma1, mu2, sigma2) must be provided together."
            )
        return sample_lnmix(n_bs, k_lnmix_w, k_lnmix_mu1, k_lnmix_sigma1,
                            k_lnmix_mu2, k_lnmix_sigma2, seed=seed)
    else:
        raise ValueError(
            "At least one k-distribution must be specified: k_posterior_sample, "
            "k_gamma_shape/scale, k_invgamma_alpha, or k_lnmix_w."
        )


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


def _sample_positive_normal(mean, sd, n):
    """Sample ``n`` values from Normal(``mean``, ``sd``) truncated to positive
    values (> 0).

    Matches the Go backend's ``GenerateTruncatedNormal`` (truncation at 0).

    NOTE the ``scipy.stats.truncnorm`` convention: its ``a``/``b`` bounds are in
    *standard deviations from* ``loc``, so the lower truncation point 0 maps to
    ``a = (0 - mean) / sd`` — **not** ``a = 0``, which truncates at the mean.
    Passing ``a = 0`` (the long-standing bug this helper replaces) discards the
    entire lower half of the distribution and inflates the sampled mean by
    ``≈ 0.8 * sd``.

    Uses the legacy NumPy global RNG state (seeded by the caller) for
    reproducibility, consuming ``n`` draws regardless of the truncation point.
    """
    if sd <= 0:
        return np.full(n, float(mean))
    return stats.truncnorm.rvs(-mean / sd, np.inf, mean, sd, n)


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
    integration_method="gauss-legendre",
):
    if n_bs <= 0:
        raise ValueError("n_bs must be greater than zero to perform simulations.")

    np.random.seed(seed)
    ks = _sample_k(
        n_bs, seed,
        k_posterior_sample=k_posterior_sample,
        k_gamma_shape=k_gamma_shape,
        k_gamma_scale=k_gamma_scale,
        k_invgamma_alpha=k_invgamma_alpha,
        k_invgamma_beta=k_invgamma_beta,
        k_invgamma_mode=k_invgamma_mode,
        k_lnmix_w=k_lnmix_w,
        k_lnmix_mu1=k_lnmix_mu1,
        k_lnmix_sigma1=k_lnmix_sigma1,
        k_lnmix_mu2=k_lnmix_mu2,
        k_lnmix_sigma2=k_lnmix_sigma2,
    )
    doubling_times = _sample_positive_normal(doubling_time, doubling_time_norm_sd, n_bs)
    lod50s = _sample_positive_normal(lod50, lod50_sd, n_bs)
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
    from functools import partial

    _rd = partial(_risk_days, integration_method=integration_method)
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(_rd, *args) for args in args_list]
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
        _col_names = [
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
        ]
        sim_df = pl.DataFrame(
            {name: list(col) for name, col in zip(_col_names, zip(*args_list))}
        ).with_columns(
            (pl.col("lod50") * pl.col("lod95_lod50_ratio")).alias("lod95"),  # Convert ratio to actual lod95
            pl.Series("iwp", rdests),
            pl.lit(seed).alias("random_seed"),
        )

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
            integration_method=integration_method,
        )
    elif point_estimate == "median":
        rd_pe = statistics.median(rdests)
    elif point_estimate == "mean":
        rd_pe = statistics.mean(rdests)
    elif point_estimate == "mode":
        rd_pe = _kde_mode_log(rdests, n_grid=1_000_000, cap=None) # Accurate but impractically slow
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
    integration_method="gauss-legendre",
):
    """
    Risk days bootstrap calculation with optional Go acceleration.

    Parameters
    ----------
    use_go : bool, default=False
        If True, uses the high-performance Go implementation.
        If False (default), uses the Python implementation.
        If Go implementation fails, automatically falls back to Python.
    integration_method : {"gauss-legendre", "quad"}, default="gauss-legendre"
        Numerical integration scheme for the risk-days integrand.
        "gauss-legendre" uses a fixed 1000-point Gauss-Legendre rule matching
        the Go backend (robust default). "quad" uses scipy's adaptive
        Gauss-Kronrod quadrature and is provided for reproducing prior
        analyses computed with quad; it is only available on the Python path
        (use_go=False), since the Go backend always uses Gauss-Legendre.

    All other parameters are passed through to the underlying implementation.
    """
    if use_go and integration_method != "gauss-legendre":
        raise ValueError(
            "Go acceleration only implements 'gauss-legendre' integration; "
            "set use_go=False to use integration_method='quad'."
        )
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
        integration_method,
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
    incidence_draws = _sample_positive_normal(incidence, incidence_norm_sd, n_bs)
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


def total_residual_risk_rd(
    components,
    per=1e6,
    seed=126887,
    alpha=0.05,
    return_samps=False,
):
    """Additive total residual risk across populations, with a *joint* credible
    interval.

    ``components`` is a non-empty sequence of
    ``(iwp_pe, iwp_bs, incidence, incidence_norm_sd)`` tuples — one per population
    (e.g. baseline, oral-PrEP, injectable-PrEP). ``iwp_bs`` is that population's
    bootstrap sample of risk-day equivalents (IWP).

    The per-iteration residual-risk *probability* for component ``c`` is

        p_c[b] = incidence_draw_c[b] * iwp_bs_c[b] / 365.25

    and the total per iteration is ``T[b] = sum_c p_c[b]``. The credible interval
    is the empirical ``(alpha/2, 1 - alpha/2)`` quantiles of ``T`` — i.e. we sum
    the components *within each iteration* and then take quantiles of the summed
    distribution.

    Two modelling assumptions are baked in, both deliberate:

    1. **Shared parameters are aligned across components.** The ``iwp_bs_c``
       arrays must already share the same per-iteration draws of the parameters
       common to every component (infectivity ``k``, viral doubling time, LOD,
       transfused volume). This holds when they are produced by the Go backend
       with a common seed, because those parameters are drawn *before* the
       baseline/PrEP branch, so the PrEP-specific draws do not perturb them.
       Pairing the arrays by iteration index then preserves the positive
       correlation those shared parameters induce between components — which is
       what makes the summed-quantile interval a valid joint CrI rather than an
       independence approximation. (The pure-Python bootstrap draws in a
       different order and is **not** aligned; see the UI caveat.)

    2. **Incidence is independent across populations.** Each component's
       incidence is drawn with its own seed (``seed + i``), reflecting the
       assumption that the populations' incidence-rate uncertainties are
       independent.

    Returns ``(rr_pe, rr_cri, onein_pe, onein_cri)``; with ``return_samps=True``
    the per-iteration total-probability sample is appended. ``rr_pe``/``rr_cri``
    are per ``per`` transfusions; ``onein_pe``/``onein_cri`` are in the
    "1 in N transfusions" representation, with ``onein_cri`` following the same
    ``(alpha/2, 1 - alpha/2)`` quantile ordering as
    ``residual_risk_rd(one_in_x=True)``.
    """
    components = list(components)
    if not components:
        raise ValueError("components must be a non-empty sequence")

    n = len(components[0][1])
    if n == 0:
        raise ValueError("iwp_bs arrays must be non-empty")
    if any(len(c[1]) != n for c in components):
        raise ValueError(
            "all components must have the same number of bootstrap iterations; "
            "iwp_bs arrays must be equal length and per-iteration aligned"
        )

    total_prob_pe = 0.0
    total_prob_samp = np.zeros(n, dtype=float)
    for i, (iwp_pe, iwp_bs, incidence, incidence_norm_sd) in enumerate(components):
        if incidence <= 0:
            raise ValueError(f"incidence must be positive, got {incidence}")
        if iwp_pe <= 0:
            raise ValueError(f"iwp_pe must be positive, got {iwp_pe}")
        total_prob_pe += incidence * iwp_pe / 365.25
        # Independent incidence draws per population (seed + i). Uses the legacy
        # global RNG via _sample_positive_normal, exactly like residual_risk_rd,
        # so a single-component call reproduces residual_risk_rd bit-for-bit.
        np.random.seed(seed + i)
        incidence_draws = _sample_positive_normal(incidence, incidence_norm_sd, n)
        total_prob_samp += incidence_draws * np.asarray(iwp_bs, dtype=float) / 365.25

    if total_prob_pe <= 0 or np.any(total_prob_samp <= 0):
        raise ValueError("total residual-risk probability must be positive")

    rr_pe = total_prob_pe * per
    rr_cri = tuple(np.quantile(total_prob_samp * per, (alpha / 2, 1 - alpha / 2)))
    onein_pe = 1.0 / total_prob_pe
    onein_cri = tuple(np.quantile(1.0 / total_prob_samp, (alpha / 2, 1 - alpha / 2)))

    if return_samps:
        return (rr_pe, rr_cri, onein_pe, onein_cri, total_prob_samp)
    return (rr_pe, rr_cri, onein_pe, onein_cri)
