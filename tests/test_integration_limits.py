# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025-2026  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the exposed integration domain (``limits``).

The domain used to be duplicated — a Python default plus four hardcoded
``-100``/``500`` literals in the Go ``riskdays.go`` — so a caller could not change
it and the two sides could silently drift. ``limits`` is now a public parameter on
``risk_days_bs`` / ``risk_days_prep_bs``, threaded through the JSON bridge to Go
and validated identically on both paths.

The parity tests compare Python's *single integration* against Go's
``"primary parameters"`` point estimate — the same integral — so they need no
``ProcessPoolExecutor`` and stay sandbox-safe. The one test that must exercise the
Python *bootstrap* carries the ``multiprocessing`` marker.
"""

import pytest

import residualrisk as rr
from residualrisk import core as C
from residualrisk import prep as P

_BASE = dict(
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
    k_invgamma_alpha=2.0,
    k_invgamma_beta=0.002019,
    n_bs=200,
    seed=42,
)

# Inverted, zero-width, and non-finite domains. `not (lo < hi)` catches NaN because
# every NaN comparison is False.
_BAD_LIMITS = [
    (500, -100),
    (5, 5),
    (float("nan"), 500),
    (float("-inf"), 500),
    (-100, float("inf")),
]


def _py_baseline(limits):
    """The baseline integral at the point values — no multiprocessing."""
    return C._risk_days(
        2, 0.00025, 0.8542, 200, 0.000673, 16, 2.73, 3.5, 1, z=1.6449, limits=limits
    )


def _py_prep(limits):
    """The PrEP integral at the point values (PrEP defaults) — no multiprocessing."""
    return P._risk_days_prep(
        2, 0.00025, 0.8542, 336, 7.0, 0.7, 0.6, 1,
        200, 0.000673, 16, 2.73, 3.5, 1,
        28.7, 250, 50.49434, 1.15062, 1.6449, 1.0, limits,
    )


def _go_pe(limits, prep=False, **over):
    """Go's "primary parameters" PE — the same single integral, via the binary."""
    fn = rr.risk_days_prep_bs if prep else rr.risk_days_bs
    pe, *_ = fn(
        **{**_BASE, **over},
        use_go=True,
        point_estimate="primary parameters",
        limits=limits,
    )
    return pe


class TestLimitsValidation:
    """Both backends must reject the same degenerate domains, before dispatch."""

    @pytest.mark.parametrize("bad", _BAD_LIMITS)
    @pytest.mark.parametrize("use_go", [False, True])
    def test_baseline_rejects(self, bad, use_go):
        with pytest.raises(ValueError, match="limits"):
            rr.risk_days_bs(**_BASE, limits=bad, use_go=use_go)

    @pytest.mark.parametrize("bad", _BAD_LIMITS)
    @pytest.mark.parametrize("use_go", [False, True])
    def test_prep_rejects(self, bad, use_go):
        with pytest.raises(ValueError, match="limits"):
            rr.risk_days_prep_bs(**_BASE, limits=bad, use_go=use_go)


class TestLimitsReachGo:
    """A custom domain must actually reach the Go backend.

    Go previously hardcoded (-100, 500) and silently ignored the caller's domain, so
    it would return the *default* answer while Python honoured the request. These
    lock that divergence out: Go must reproduce Python's value for each domain.
    """

    @pytest.mark.parametrize("limits", [(-100, 500), (-50, 300), (0, 10)])
    def test_baseline_python_go_agree(self, limits):
        assert _go_pe(limits) == pytest.approx(_py_baseline(limits), rel=1e-9)

    @pytest.mark.parametrize("limits", [(-100, 500), (0, 60), (0, 30)])
    def test_prep_python_go_agree(self, limits):
        assert _go_pe(limits, prep=True) == pytest.approx(_py_prep(limits), rel=1e-9)

    def test_truncating_domain_changes_the_answer(self):
        """Guards the parity tests against passing vacuously.

        The default (-100, 500) is deliberately generous — the integrand is ~0 well
        inside it — so a domain must actually bite for the parity checks to mean
        anything. (0, 10) cuts the mass that sits at negative t: ~-91% in practice.
        """
        full = _py_baseline((-100, 500))
        cut = _py_baseline((0, 10))
        assert cut < 0.5 * full
        assert _go_pe((0, 10)) == pytest.approx(cut, rel=1e-9)

    def test_omitting_limits_uses_the_default(self):
        """Callers that never pass `limits` must be unaffected."""
        assert _go_pe((-100, 500)) == pytest.approx(
            rr.risk_days_bs(**_BASE, use_go=True, point_estimate="primary parameters")[0],
            rel=1e-12,
        )


class TestPointEstimateHonoursZ:
    """Regression: the baseline "primary parameters" PE omitted ``z`` (and ``limits``)
    when calling ``_risk_days``, so it silently used that function's own defaults
    (``z=1.6449``) while the bootstrap used the caller's ``z``. A custom ``z``
    therefore produced a PE inconsistent with its own credible interval — and with
    Go, which does honour ``z``.
    """

    def test_custom_z_actually_moves_the_answer(self):
        """Otherwise the parity test below would be vacuous."""
        assert abs(_py_baseline_z(2.5) - _py_baseline_z(1.6449)) > 1e-6

    @pytest.mark.multiprocessing
    def test_python_pe_honours_custom_z_and_matches_go(self):
        pe_py, *_ = rr.risk_days_bs(
            **_BASE, z=2.5, use_go=False, point_estimate="primary parameters", threads=2
        )
        pe_go, *_ = rr.risk_days_bs(
            **_BASE, z=2.5, use_go=True, point_estimate="primary parameters"
        )
        assert pe_py == pytest.approx(pe_go, rel=1e-9)


def _py_baseline_z(z):
    return C._risk_days(
        2, 0.00025, 0.8542, 200, 0.000673, 16, 2.73, 3.5, 1, z=z, limits=(-100, 500)
    )
