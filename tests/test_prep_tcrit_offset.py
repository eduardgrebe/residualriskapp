# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025-2026  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""``tcrit`` targets the plateau's central level, so the PrEP viral-load trajectory
is continuous at the growth→plateau crossover for **every** offset.

``_find_tcrit`` used to solve for the time growth reaches the bare set-point, while
the plateau branch begins at ``offset * set_point`` (the sinusoid starts at zero
phase). The two agreed only at ``offset == 1``; at any other offset the modelled
viral load jumped by exactly a factor of ``offset`` at ``tcrit`` — an instantaneous,
physically impossible change. The offset was UI-reachable over 0.05–2.0 at the time.

Retargeting ``tcrit`` to ``offset * set_point / copies_per_virion`` restores
continuity everywhere and is a **bit-for-bit no-op at offset = 1**, so no prior
result moves — `test_offset_one_is_bit_for_bit_unchanged` pins exactly that.

The parity tests compare Python's *single integration* against Go's
``"primary parameters"`` point estimate — the same integral — so they need no
``ProcessPoolExecutor`` and stay sandbox-safe.
"""

import json
import math
import subprocess

import pytest

import residualrisk as rr
from residualrisk import prep as P
from residualrisk._go import find_go_binary

# PrEP defaults, matching tests/test_integration_limits.py::_py_prep.
_ECLIPSE, _C0, _DT, _SP, _CPV = 7.0, 0.00025, 0.8542, 336, 2
_A, _B = 0.7, 0.6

_OFFSETS = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0]

_BASE = dict(
    k=0.000673,
    doubling_time=_DT,
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
    n_bs=50,
    seed=42,
)


def _tcrit(offset, set_point=_SP):
    return P._find_tcrit(_ECLIPSE, _C0, _DT, set_point, _CPV, offset=offset)


def _vl(t, offset, tcrit, set_point=_SP, a=_A):
    return P._vl_postbt(t, _ECLIPSE, _C0, _DT, set_point, a, _B, offset, tcrit, _CPV)


def _py_rde(offset, set_point=_SP, a=_A):
    """The PrEP integral at the point values — single integration, no pool."""
    return P._risk_days_prep(
        _CPV, _C0, _DT, set_point, _ECLIPSE, a, _B, offset,
        200, 0.000673, 16, 2.73, 3.5, 1,
        28.7, 250, 50.49434, 1.15062, 1.6449,
    )


def _amplitude_for(offset):
    """a must be <= offset, or the plateau viral load would go negative."""
    return min(_A, offset)


def _go_rde(offset, **over):
    """Go's "primary parameters" PE — the same single integral, via the binary."""
    pe, *_ = rr.risk_days_prep_bs(
        **{**_BASE, **over},
        offset=offset,
        a=_amplitude_for(offset),
        b=_B,
        use_go=True,
        point_estimate="primary parameters",
    )
    return pe


class TestContinuityAtTcrit:
    """The defect: a jump of exactly `offset` in viral load at the crossover."""

    @pytest.mark.parametrize("offset", _OFFSETS)
    def test_trajectory_is_continuous(self, offset):
        tc = _tcrit(offset)
        eps = 1e-9
        left = _vl(tc - eps, offset, tc)   # end of exponential growth
        right = _vl(tc + eps, offset, tc)  # start of the oscillating plateau
        assert left == pytest.approx(right, rel=1e-6), (
            f"viral load jumps by {right / left:.2f}x at tcrit for offset={offset}"
        )

    @pytest.mark.parametrize("offset", _OFFSETS)
    def test_growth_ends_at_the_plateau_centre(self, offset):
        """The analytic statement of the same thing: growth must reach
        offset * set_point / copies_per_virion, which is where the plateau begins."""
        tc = _tcrit(offset)
        assert _vl(tc, offset, tc) == pytest.approx(offset * _SP / _CPV, rel=1e-9)

    @pytest.mark.parametrize("offset", _OFFSETS)
    def test_tcrit_matches_closed_form(self, offset):
        expected = _ECLIPSE + _DT * math.log2((offset * _SP / _CPV) / _C0)
        assert _tcrit(offset) == pytest.approx(expected, rel=1e-12)

    def test_higher_offset_delays_tcrit(self):
        """Sanity on direction: a higher plateau takes longer to grow to."""
        tcrits = [_tcrit(o) for o in _OFFSETS]
        assert tcrits == sorted(tcrits)


class TestNoRegressionAtTheDefault:
    """offset = 1 must be untouched — every shipped and published result used it."""

    def test_offset_one_is_bit_for_bit_unchanged(self):
        """The retarget multiplies the target level by offset, so at offset = 1 the
        expression is literally the old one. Pinned exactly, not approximately."""
        old_formula = _ECLIPSE + _DT * math.log2((_SP / _CPV) / _C0)
        assert _tcrit(1.0) == old_formula

    def test_default_offset_is_one(self):
        """_find_tcrit's default must stay 1.0, so callers that never heard of the
        offset keep the old behaviour."""
        assert P._find_tcrit(_ECLIPSE, _C0, _DT, _SP, _CPV) == _tcrit(1.0)


class TestOffsetIsRedundantWithSetPoint:
    """Why the UI control was removed: after the retarget, the offset is *exactly* a
    set-point multiplier. (set_point, offset, a) and (set_point*offset, 1, a/offset)
    are the same model — same tcrit, same plateau, same RDE. The set point is the
    parameter to vary: it is in clinical units and has its own bootstrap range."""

    @pytest.mark.parametrize("offset", [0.5, 0.8, 1.2, 1.5, 2.0])
    def test_reparameterisation_identity_tcrit(self, offset):
        assert _tcrit(offset) == pytest.approx(
            _tcrit(1.0, set_point=_SP * offset), rel=1e-12
        )

    @pytest.mark.parametrize("offset", [0.5, 0.8, 1.2, 1.5, 2.0])
    def test_reparameterisation_identity_rde(self, offset):
        a_scaled = _A / offset
        if a_scaled > 1.0:
            pytest.skip("a/offset > 1 would violate the a <= offset constraint")
        assert _py_rde(offset) == pytest.approx(
            _py_rde(1.0, set_point=_SP * offset, a=a_scaled), rel=1e-9
        )


class TestPythonGoParity:
    """Both engines must retarget identically — otherwise the backends disagree at
    any offset != 1, which is exactly where the bug lived."""

    @pytest.mark.parametrize("offset", _OFFSETS)
    def test_rde_matches(self, offset):
        if find_go_binary() is None:
            pytest.skip("no Go binary available")
        expected = _py_rde(offset, a=_amplitude_for(offset))
        assert _go_rde(offset) == pytest.approx(expected, rel=1e-6)


class TestOffsetValidation:
    """offset <= 0 makes the target level non-positive and tcrit -inf."""

    @pytest.mark.parametrize("bad", [0.0, -0.5])
    @pytest.mark.parametrize("use_go", [False, True])
    def test_rejected(self, bad, use_go):
        if use_go and find_go_binary() is None:
            pytest.skip("no Go binary available")
        with pytest.raises(ValueError, match="offset must be positive"):
            rr.risk_days_prep_bs(**_BASE, offset=bad, a=0.0, use_go=use_go)

    def test_go_binary_rejects_nonpositive_offset(self):
        """The Go Validate() backstop for a raw JSON payload. Note PrEP scalars are
        deliberately not defaulted from a 0 sentinel (see models.go SetDefaults), so
        an *omitted* offset lands here as 0 and is correctly rejected rather than
        silently producing a zero plateau."""
        binary = find_go_binary()
        if binary is None:
            pytest.skip("no Go binary available")
        payload = {
            "k": 0.000673, "doubling_time": _DT, "doubling_time_norm_sd": 0.2813,
            "lod50": 2.73, "lod50_sd": 0.53, "lod95_lod50_ratio": 3.5,
            "volume_transfused": 200, "volume_transfused_min": 100,
            "volume_transfused_max": 340, "pool_size": 16, "retests": 1,
            "k_invgamma_alpha": 2.0, "k_invgamma_beta": 0.002019,
            "n_bs": 50, "seed": 42, "threads": 1,
            "prep_mode": True, "set_point": _SP, "eclipse": _ECLIPSE,
            "a": 0.0, "b": _B, "drug_effect": 1.0,
            "ser_min": 28.7, "ser_max": 250, "ser_alpha": 50.49434, "ser_beta": 1.15062,
            # "offset" omitted -> 0 -> must be rejected
        }
        result = subprocess.run(
            [binary], input=json.dumps(payload), capture_output=True, text=True, timeout=60
        )
        assert result.returncode != 0
        assert "offset must be positive" in (result.stderr + result.stdout).lower()
