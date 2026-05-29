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

"""PrEP breakthrough infection residual risk model."""

import math
import statistics
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

import numpy as np
import polars as pl
from scipy.integrate import quad

from .core import (
    _integrate_gauss_legendre,
    _kde_mode_log,
    _prob_infectious_copies,
    _prob_neg_retest,
    _prob_pos_init,
    _sample_k,
    _sample_positive_normal,
    get_cpu_core_count,
)


def _sin_varied(t, a, b, offset):
    """Sinusoidal oscillation factor for the post-breakthrough viral-load plateau.

    Returns ``offset + a * sin(b * t)``. In the PrEP breakthrough model the
    plateau viral load is ``set_point * _sin_varied(t - tcrit, a, b, offset)``
    (see :func:`_vl_postbt`), so this factor scales the set-point: the plateau
    oscillates between ``(offset - a)`` and ``(offset + a)`` times ``set_point``.

    Parameters
    ----------
    t : float or array
        Time **since the plateau onset** ``tcrit`` (days); callers pass
        ``t - tcrit``.
    a : float
        Oscillation amplitude, as a fraction of ``set_point``.
    b : float
        Angular frequency (radians/day); the oscillation period is ``2*pi / b``.
    offset : float
        Central level of the oscillation, as a multiple of ``set_point``
        (``offset = 1`` centres the plateau on ``set_point``).

    Notes
    -----
    If ``a > offset`` the factor goes negative on its downswings, which would
    imply a negative viral load; :func:`_vl_postbt` clamps the result to ``0``.
    """
    return offset + a * np.sin(b * t)


def _find_tcrit(eclipse, C0, doubling_time, set_point):
    """Time at which exponential growth first reaches the set-point.

    Closed-form solution of ``C0 * 2**((t - eclipse) / doubling_time) = set_point``:

        tcrit = eclipse + doubling_time * log2(set_point / C0)

    Mirrors the Go backend's ``FindTcrit`` (``go/riskdays/prep.go``) exactly, so
    the two implementations agree to machine precision. Replaces the former
    grid-search ``_vl_postbt_vec`` — this is exact (not rounded to a 0.1-day
    grid), O(1), and immune to the empty-``argmin`` crash that occurred when
    ``tcrit`` exceeded the search grid.
    """
    return eclipse + doubling_time * math.log2(set_point / C0)


def _vl_postbt(t, eclipse, C0, doubling_time, set_point, a, b, offset, tcrit):
    if t < eclipse:
        vl = 0.0
    elif t <= tcrit:
        # Exponential growth phase. The exponent is bounded here (at t=tcrit it
        # equals log2(set_point/C0)), so no overflow.
        vl = C0 * 2 ** ((t - eclipse) / doubling_time)
    else:
        # Oscillating plateau. The exponential is deliberately NOT evaluated in
        # this branch: for large t with a small (now-sampleable) doubling_time it
        # would overflow, and the value would only be discarded — cf. the
        # growth-exponent cap in core._concentration.
        vl = set_point * _sin_varied(t=t - tcrit, a=a, b=b, offset=offset)
    # Modelled viral load can dip below zero when the sinusoidal set-point
    # oscillation amplitude exceeds its offset (a > offset); clamp to a
    # physical floor of zero.
    return max(0.0, vl)


def _prob_infectious_prep(
    t,
    eclipse,
    C0,
    doubling_time,
    set_point,
    a,
    b,
    offset,
    volume_transfused,
    k,
    copies_per_virion=2.0,
):
    tcrit = _find_tcrit(eclipse, C0, doubling_time, set_point)
    C = _vl_postbt(
        t=t,
        eclipse=eclipse,
        C0=C0,
        doubling_time=doubling_time,
        set_point=set_point,
        a=a,
        b=b,
        offset=offset,
        tcrit=tcrit,
    )
    n_copies = C * copies_per_virion * volume_transfused
    prob = _prob_infectious_copies(n_copies, k)
    return prob


def _prob_nondetection_serology_prep(t, min, max, alpha, beta):
    if t < min:
        p = 1.0
    elif t > max:
        p = 0.0
    else:
        p = math.exp(-((t - min) / alpha) ** beta)
    return p


def _prob_nondetection_prep(
    t,
    copies_per_virion,
    C0,
    doubling_time,
    eclipse,
    set_point,
    a,
    b,
    offset,
    pool_size,
    lod50,
    lod95_lod50_ratio,
    retests,
    z=1.6449,
    seroconversion_delay_median=45,
):
    tcrit = _find_tcrit(eclipse, C0, doubling_time, set_point)
    Cv = _vl_postbt(
        t=t,
        eclipse=eclipse,
        C0=C0,
        doubling_time=doubling_time,
        set_point=set_point,
        a=a,
        b=b,
        offset=offset,
        tcrit=tcrit,
    )
    Cc = copies_per_virion * Cv
    if Cc == 0.0:
        return 1.0
    elif Cc > 0.0:
        p_pos_init = _prob_pos_init(Cc, doubling_time, pool_size, lod50, lod95_lod50_ratio, z)
        p_neg_retest = _prob_neg_retest(
            Cc, doubling_time, pool_size, lod50, lod95_lod50_ratio, retests, z
        )
        prob = 1 - p_pos_init * (1 - p_neg_retest)
        return prob


def _prob_infectious_nondetection_prep(
    t,
    eclipse,
    C0,
    doubling_time,
    set_point,
    a,
    b,
    offset,
    volume_transfused,
    k,
    copies_per_virion,
    pool_size,
    lod50,
    lod95_lod50_ratio,
    retests,
    ser_min,
    ser_max,
    ser_alpha,
    ser_beta,
    z=1.6449,
):
    product = (
        _prob_infectious_prep(
            t=t,
            eclipse=eclipse,
            C0=C0,
            doubling_time=doubling_time,
            set_point=set_point,
            a=a,
            b=b,
            offset=offset,
            volume_transfused=volume_transfused,
            k=k,
        )
        * _prob_nondetection_prep(
            t=t,
            copies_per_virion=copies_per_virion,
            C0=C0,
            doubling_time=doubling_time,
            eclipse=eclipse,
            set_point=set_point,
            a=a,
            b=b,
            offset=offset,
            pool_size=pool_size,
            lod50=lod50,
            lod95_lod50_ratio=lod95_lod50_ratio,
            retests=retests,
            z=z,
        )
        * _prob_nondetection_serology_prep(
            t,
            min=ser_min,
            max=ser_max,
            alpha=ser_alpha,
            beta=ser_beta,
        )
    )
    return product


def _risk_days_prep(
    copies_per_virion,
    C0,
    doubling_time,
    set_point,
    eclipse,
    a,
    b,
    offset,
    volume_transfused,
    k,
    pool_size,
    lod50,
    lod95_lod50_ratio,
    retests,
    ser_min,
    ser_max,
    ser_alpha,
    ser_beta,
    z,
    limits=(-100, 500),
    integration_method="gauss-legendre",
):
    # Ideally we would integrate from -np.inf to np.inf, but that causes an
    # overflow error, so we choose safe limits instead.
    args = (
        eclipse, C0, doubling_time, set_point, a, b, offset,
        volume_transfused, k, copies_per_virion, pool_size,
        lod50, lod95_lod50_ratio, retests, ser_min, ser_max,
        ser_alpha, ser_beta,
    )
    if integration_method == "gauss-legendre":
        # Fixed 1000-point Gauss-Legendre, matching the Go backend. The default:
        # the PrEP integrand has compact support (it is exactly zero before the
        # eclipse phase and after the serology cutoff), so adaptive quad can
        # silently miss the active window and return ~0 — the GL rule samples
        # the whole interval and cannot.
        rd = _integrate_gauss_legendre(
            _prob_infectious_nondetection_prep, limits[0], limits[1], args
        )
    elif integration_method == "quad":
        # Adaptive Gauss-Kronrod (scipy). Retained so prior analyses computed
        # with quad can be reproduced exactly via the Python API. limit=500
        # matches the historical PrEP call.
        rd = quad(
            _prob_infectious_nondetection_prep,
            limits[0],
            limits[1],
            limit=500,
            args=args,
        )[0]
    else:
        raise ValueError(
            "integration_method must be 'gauss-legendre' or 'quad', "
            f"got {integration_method!r}"
        )
    return rd


def risk_days_prep_bs(
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
    set_point=336,
    set_point_dist_uniform=(19.1, 2265),
    eclipse=7.0,
    eclipse_dist_uniform=(4.0, 10.0),
    a=0.7,
    b=0.6,
    offset=1,
    a_dist_uniform=None,
    b_dist_uniform=None,
    ser_min=28.7,
    ser_max=250,
    ser_alpha=50.49434,
    ser_beta=1.15062,
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
    """Bootstrap the PrEP-breakthrough risk-day-equivalents (RDE) distribution.

    Each iteration samples ``k``, ``doubling_time`` / ``lod50`` (positive
    truncated normal), and ``set_point`` / ``eclipse`` / ``volume_transfused``
    (uniform over their ``*_dist_uniform`` / ``*_range``). The sinusoidal
    set-point oscillation parameters ``a`` and ``b`` are held fixed at their
    scalar values unless ``a_dist_uniform`` / ``b_dist_uniform`` are given, in
    which case they are sampled ``Uniform(lo, hi)`` per iteration (``None`` =
    fixed, the default — back-compatible). ``offset`` is always fixed. ``a``
    (and the upper bound of ``a_dist_uniform``) must be ``<= offset``, or the
    plateau viral load would go negative.

    Returns ``(rd_pe, rd_cri, rd_range, rdests, sim_df)``; ``sim_df`` is ``None``
    unless ``return_sim_df=True``.
    """
    if use_go and integration_method != "gauss-legendre":
        raise ValueError(
            "Go acceleration only implements 'gauss-legendre' integration; "
            "set use_go=False to use integration_method='quad'."
        )
    if use_go:
        try:
            from ._go import risk_days_prep_bs_go

            return risk_days_prep_bs_go(
                k=k,
                doubling_time=doubling_time,
                doubling_time_norm_sd=doubling_time_norm_sd,
                lod50=lod50,
                lod50_sd=lod50_sd,
                lod95_lod50_ratio=lod95_lod50_ratio,
                volume_transfused=volume_transfused,
                volume_transfused_range=volume_transfused_range,
                pool_size=pool_size,
                retests=retests,
                set_point=set_point,
                set_point_dist_uniform=set_point_dist_uniform,
                eclipse=eclipse,
                eclipse_dist_uniform=eclipse_dist_uniform,
                a=a,
                b=b,
                offset=offset,
                a_dist_uniform=a_dist_uniform,
                b_dist_uniform=b_dist_uniform,
                ser_min=ser_min,
                ser_max=ser_max,
                ser_alpha=ser_alpha,
                ser_beta=ser_beta,
                C0=C0,
                copies_per_virion=copies_per_virion,
                alpha=alpha,
                z=z,
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
                n_bs=n_bs,
                seed=seed,
                threads=threads,
                point_estimate=point_estimate,
                mode_precision=mode_precision,
                progress=progress,
                return_sim_df=return_sim_df,
            )
        except Exception:
            pass  # fall back to Python

    if n_bs <= 0:
        raise ValueError("n_bs must be greater than zero to perform simulations.")
    # The sinusoidal amplitude must not exceed the offset, or the plateau viral
    # load would go negative (it would otherwise be clamped to 0 in _vl_postbt).
    if a > offset:
        raise ValueError(f"a ({a}) must be <= offset ({offset}).")
    if a_dist_uniform is not None and a_dist_uniform[1] > offset:
        raise ValueError(
            f"a_dist_uniform upper bound ({a_dist_uniform[1]}) must be <= offset ({offset})."
        )

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
    set_points = np.random.uniform(set_point_dist_uniform[0], set_point_dist_uniform[1], n_bs)
    eclipses = np.random.uniform(eclipse_dist_uniform[0], eclipse_dist_uniform[1], n_bs)
    lod50s = _sample_positive_normal(lod50, lod50_sd, n_bs)
    volumes_transfused = np.random.uniform(
        volume_transfused_range[0], volume_transfused_range[1], n_bs
    )
    # Sinusoidal oscillation parameters: fixed at the scalar value unless a
    # uniform range is given. np.full draws no RNG, so reproducibility of the
    # other parameters is unchanged when these are not varied. offset is never
    # varied; a <= offset is enforced above.
    a_s = (
        np.random.uniform(a_dist_uniform[0], a_dist_uniform[1], n_bs)
        if a_dist_uniform is not None else np.full(n_bs, a)
    )
    b_s = (
        np.random.uniform(b_dist_uniform[0], b_dist_uniform[1], n_bs)
        if b_dist_uniform is not None else np.full(n_bs, b)
    )

    args_list = [
        (
            copies_per_virion, C0, doubling_times[i], set_points[i], eclipses[i],
            a_s[i], b_s[i], offset, volumes_transfused[i], ks[i], pool_size, lod50s[i],
            lod95_lod50_ratio, retests, ser_min, ser_max, ser_alpha, ser_beta, z,
            (-100, 500),
        )
        for i in range(n_bs)
    ]

    rdests = []
    _rd = partial(_risk_days_prep, integration_method=integration_method)
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
                if completed_count == 1 or new_percent > getattr(progress, "_last_percent", 0):
                    progress._last_percent = new_percent
                    progress.progress(completed_count / n_bs, text=f"Progress: {new_percent}%")

    rd_range = [np.min(rdests), np.max(rdests)]
    rd_cri = np.quantile(rdests, (alpha / 2, 1 - alpha / 2))

    if return_sim_df:
        _col_names = [
            "copies_per_virion",
            "C0",
            "doubling_time",
            "set_point",
            "eclipse",
            "a",
            "b",
            "offset",
            "volume_transfused",
            "k",
            "pool_size",
            "lod50",
            "lod95_lod50_ratio",
            "retests",
            "ser_min",
            "ser_max",
            "ser_alpha",
            "ser_beta",
            "z",
            "limits",
        ]
        sim_df = pl.DataFrame(
            {name: list(col) for name, col in zip(_col_names, zip(*args_list))}
        ).with_columns(
            (pl.col("lod50") * pl.col("lod95_lod50_ratio")).alias("lod95"),
            pl.Series("iwp", rdests),
            pl.lit(seed).alias("random_seed"),
        )

    if point_estimate == "primary parameters":
        rd_pe = _risk_days_prep(
            copies_per_virion, C0, doubling_time, set_point, eclipse,
            a, b, offset, volume_transfused, k, pool_size, lod50,
            lod95_lod50_ratio, retests, ser_min, ser_max, ser_alpha, ser_beta, z,
            (-100, 500),
            integration_method=integration_method,
        )
    elif point_estimate == "median":
        rd_pe = statistics.median(rdests)
    elif point_estimate == "mean":
        rd_pe = statistics.mean(rdests)
    elif point_estimate == "mode":
        # KDE-log mode — matches Go's KDEModeLog and the baseline Python path
        # (_risk_days_bs_python). Accurate but slow. mode_precision is unused
        # here (kept for API / Go-bridge compatibility), as in the baseline.
        rd_pe = _kde_mode_log(rdests, n_grid=1_000_000, cap=None)
    else:
        rd_pe = None

    if return_sim_df:
        return (rd_pe, rd_cri, rd_range, rdests, sim_df)
    else:
        return (rd_pe, rd_cri, rd_range, rdests, None)
