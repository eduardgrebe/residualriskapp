# Residual HIV Transfusion Transmission Risk Estimator
# Copyright (C) 2025-2026 Vitalant and Eduard Grebe Consulting
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

"""The four *k* input modes are mutually exclusive, not a priority cascade.

Both engines used to dispatch on the first non-``None`` k-parameter in a fixed order
(posterior > gamma > invgamma > lnmix), so specifying two silently ran the
higher-priority one. The failure was quiet and scientific: a library caller reusing a
config dict that still carried ``k_posterior_sample`` while setting
``k_invgamma_alpha``/``beta`` got the *posterior*, and their "InvGamma sensitivity
analysis" was silently the distribution they were trying to vary away from. A partial
spec was equally silent — ``k_gamma_shape`` without ``k_gamma_scale`` fell through to
whichever *other* mode happened to be populated.

Validation now happens pre-dispatch in Python (so both backends reject identically and
in-process, not inside a worker) and in the Go ``Validate()`` (so a hand-written JSON
payload is rejected too). AGENTS.md already described these modes as mutually
exclusive; the code now matches the documentation.

These tests all raise *before* any bootstrap runs, so none needs a process pool — the
suite stays sandbox-safe.
"""

import json
import subprocess

import numpy as np
import pytest

import residualrisk as rr
from residualrisk import core as C
from residualrisk._go import find_go_binary

_POSTERIOR = np.array([0.0005, 0.001, 0.0015])

# One fully-specified instance of each mode, keyed by the name the validator reports.
_MODES = {
    "k_posterior_sample": dict(k_posterior_sample=_POSTERIOR),
    "k_gamma_shape/k_gamma_scale": dict(k_gamma_shape=2.0, k_gamma_scale=0.001),
    "k_invgamma_*": dict(k_invgamma_alpha=2.0, k_invgamma_beta=0.002019),
    "k_lnmix_*": dict(
        k_lnmix_w=0.90,
        k_lnmix_mu1=-7.2403,
        k_lnmix_sigma1=0.3241,
        k_lnmix_mu2=-3.7423,
        k_lnmix_sigma2=0.5258,
    ),
}

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
    n_bs=50,
    seed=42,
)


class TestValidatorDirectly:
    """core._validate_k_inputs — the single source of the rule."""

    @pytest.mark.parametrize("mode", list(_MODES))
    def test_each_mode_alone_is_valid(self, mode):
        C._validate_k_inputs(**_MODES[mode])  # must not raise

    @pytest.mark.parametrize(
        "a, b",
        [
            ("k_posterior_sample", "k_invgamma_*"),
            ("k_posterior_sample", "k_lnmix_*"),
            ("k_posterior_sample", "k_gamma_shape/k_gamma_scale"),
            ("k_gamma_shape/k_gamma_scale", "k_invgamma_*"),
            ("k_gamma_shape/k_gamma_scale", "k_lnmix_*"),
            ("k_invgamma_*", "k_lnmix_*"),
        ],
    )
    def test_every_pair_raises(self, a, b):
        """All six pairs — not just the posterior-wins case the old tests pinned."""
        with pytest.raises(ValueError, match="Exactly one k-distribution"):
            C._validate_k_inputs(**_MODES[a], **_MODES[b])

    def test_all_four_raises(self):
        merged = {k: v for m in _MODES.values() for k, v in m.items()}
        with pytest.raises(ValueError, match="got 4"):
            C._validate_k_inputs(**merged)

    def test_no_mode_raises(self):
        with pytest.raises(ValueError, match="A k-distribution must be specified"):
            C._validate_k_inputs()

    def test_error_names_the_offending_modes(self):
        """The message must say *which* modes collided — otherwise the caller has to
        go hunting through a config dict for the stale key."""
        with pytest.raises(ValueError) as exc:
            C._validate_k_inputs(
                **_MODES["k_posterior_sample"], **_MODES["k_invgamma_*"]
            )
        assert "k_posterior_sample" in str(exc.value)
        assert "k_invgamma_*" in str(exc.value)

    # Partial specs: a half-populated mode still *counts* as specified, so it collides
    # rather than silently falling through to the next branch of the cascade.
    def test_partial_gamma_alone_raises(self):
        with pytest.raises(ValueError, match="must be given together"):
            C._validate_k_inputs(k_gamma_shape=2.0)

    def test_partial_gamma_plus_another_mode_raises(self):
        with pytest.raises(ValueError, match="Exactly one k-distribution"):
            C._validate_k_inputs(k_gamma_scale=0.001, **_MODES["k_invgamma_*"])

    def test_partial_lnmix_raises(self):
        with pytest.raises(ValueError, match="must be "):
            C._validate_k_inputs(k_lnmix_w=0.9, k_lnmix_mu1=-7.24)

    def test_invgamma_beta_without_alpha_raises(self):
        with pytest.raises(ValueError, match="require k_invgamma_alpha"):
            C._validate_k_inputs(k_invgamma_beta=0.002019)

    def test_invgamma_alpha_alone_raises(self):
        with pytest.raises(ValueError, match="exactly one of k_invgamma_beta"):
            C._validate_k_inputs(k_invgamma_alpha=2.0)

    def test_invgamma_beta_and_mode_together_raise(self):
        """Aligns _sample_k with the public sample_invgamma(), which already raised on
        both — the same library disagreed with itself about the same argument pair."""
        with pytest.raises(ValueError, match="exactly one of k_invgamma_beta"):
            C._validate_k_inputs(
                k_invgamma_alpha=2.0, k_invgamma_beta=0.002019, k_invgamma_mode=0.000673
            )

    def test_invgamma_mode_parameterisation_is_valid(self):
        C._validate_k_inputs(k_invgamma_alpha=2.0, k_invgamma_mode=0.000673)


class TestBaselineBootstrapEntryPoint:
    """risk_days_bs rejects pre-dispatch, on both backends."""

    @pytest.mark.parametrize("use_go", [False, True])
    def test_two_modes_raise(self, use_go):
        if use_go and find_go_binary() is None:
            pytest.skip("no Go binary available")
        with pytest.raises(ValueError, match="Exactly one k-distribution"):
            rr.risk_days_bs(
                **_BASE,
                use_go=use_go,
                **_MODES["k_posterior_sample"],
                **_MODES["k_invgamma_*"],
            )

    @pytest.mark.parametrize("use_go", [False, True])
    def test_no_mode_raises(self, use_go):
        if use_go and find_go_binary() is None:
            pytest.skip("no Go binary available")
        with pytest.raises(ValueError, match="A k-distribution must be specified"):
            rr.risk_days_bs(**_BASE, use_go=use_go)

    def test_raises_before_dispatch(self, monkeypatch):
        """The guard must fire before the engine is entered at all — otherwise the Go
        path would only fail in the subprocess and the Python path inside a worker.

        core imports risk_days_bs_go lazily (inside the function), so the patch has to
        land on the _go module, not on core.
        """
        import residualrisk._go as _go_mod

        monkeypatch.setattr(
            _go_mod,
            "risk_days_bs_go",
            lambda *a, **kw: pytest.fail("dispatched despite invalid k"),
        )
        with pytest.raises(ValueError, match="Exactly one k-distribution"):
            rr.risk_days_bs(
                **_BASE,
                use_go=True,
                **_MODES["k_posterior_sample"],
                **_MODES["k_lnmix_*"],
            )


class TestGoBinaryRejectsPayload:
    """The Go Validate() is the backstop for a hand-written JSON payload, which never
    passes through the Python validator at all."""

    def _payload(self, **k_params):
        # The Go CLI schema takes volume_transfused_min/_max, not the Python API's
        # volume_transfused_range tuple (see README's JSON example).
        return {
            "k": 0.000673,
            "doubling_time": 0.8542,
            "doubling_time_norm_sd": 0.2813,
            "lod50": 2.73,
            "lod50_sd": 0.53,
            "lod95_lod50_ratio": 3.5,
            "volume_transfused": 200,
            "volume_transfused_min": 100,
            "volume_transfused_max": 340,
            "pool_size": 16,
            "retests": 1,
            "n_bs": 50,
            "seed": 42,
            "threads": 1,
            **k_params,
        }

    def _run(self, payload):
        binary = find_go_binary()
        if binary is None:
            pytest.skip("no Go binary available")
        return subprocess.run(
            [binary],
            input=json.dumps(payload),
            capture_output=True,
            text=True,
            timeout=60,
        )

    def test_two_modes_rejected(self):
        result = self._run(
            self._payload(
                k_posterior_sample=[0.0005, 0.001, 0.0015],
                k_invgamma_alpha=2.0,
                k_invgamma_beta=0.002019,
            )
        )
        assert result.returncode != 0
        assert "exactly one k distribution" in (result.stderr + result.stdout).lower()

    def test_single_mode_still_accepted(self):
        """Guard against the new check rejecting a *valid* payload."""
        result = self._run(
            self._payload(k_invgamma_alpha=2.0, k_invgamma_beta=0.002019)
        )
        assert result.returncode == 0, result.stderr
        assert json.loads(result.stdout)["point_estimate"] > 0
