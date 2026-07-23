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

"""Golden regression for the bootstrap **headline** output (point estimate + CrI).

The pre-existing goldens (``test_golden_value_default_params``, ``3.7207``) pin
``_risk_days`` — the *single integration*. Nothing pinned the **bootstrap's** headline
PE and credible interval, so a change to the parameter sampling, the aggregation, or
the RNG plumbing could slip through every other test.

Why there are two literals rather than one: the backends draw from **different RNG
streams** (numpy vs Go), so at a fixed seed they produce different *realisations* of
the same distributions. There is no shared value to pin — each backend gets its own
golden. A third test then asserts the two agree within Monte-Carlo error, which is
what actually guards against one backend's *distribution* drifting from the other's
(a per-backend golden alone cannot see that: it would happily pin the drifted value
once regenerated).

``threads`` is pinned only so the goldens reproduce by construction. The result is in
fact thread-count independent — the bootstrap pre-draws every parameter from the seed
and fills results by index (the ``core.py:640`` indexed-fill fix), so neither worker
count nor completion order can influence it (verified on the Go path across
threads=1,2,4,8).
"""

import pytest

import residualrisk as rr

_KW = dict(
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
    n_bs=2000,
    seed=42,
    threads=4,
    point_estimate="median",
)

# Pinned from the current implementation at the config above. A *legitimate* change to
# either RNG (a numpy upgrade, a change to Go's sampling) invalidates these — they
# should then be regenerated deliberately, never reflexively to make a red test green.
_GO_PE = 3.5129519456102267
_GO_CRI = (1.0941703700582464, 7.525835605798333)

_PY_PE = 3.523607492842414
_PY_CRI = (1.0595265225746342, 7.669053212149041)

# Tight enough to catch any real change (the computation is deterministic), loose
# enough to absorb last-ULP floating-point noise.
_REL = 1e-6


class TestBootstrapGolden:
    def test_go_backend_golden(self):
        pe, cri, _rng, rdests, _ = rr.risk_days_bs(**_KW, use_go=True)
        assert len(rdests) == _KW["n_bs"]
        assert pe == pytest.approx(_GO_PE, rel=_REL)
        assert cri[0] == pytest.approx(_GO_CRI[0], rel=_REL)
        assert cri[1] == pytest.approx(_GO_CRI[1], rel=_REL)

    @pytest.mark.multiprocessing
    def test_python_backend_golden(self):
        pe, cri, _rng, rdests, _ = rr.risk_days_bs(**_KW, use_go=False)
        assert len(rdests) == _KW["n_bs"]
        assert pe == pytest.approx(_PY_PE, rel=_REL)
        assert cri[0] == pytest.approx(_PY_CRI[0], rel=_REL)
        assert cri[1] == pytest.approx(_PY_CRI[1], rel=_REL)

    @pytest.mark.multiprocessing
    def test_backends_agree_within_monte_carlo_error(self):
        """The semantic guard the per-backend goldens cannot provide.

        Each golden only pins its own backend, so if one backend's *distribution*
        drifted (a mis-sampled parameter, a wrong integration domain) its golden could
        simply be regenerated and the divergence would go unnoticed. At n_bs=2000 the
        two independent realisations agree to well under a percent on the PE and a few
        percent on the interval bounds; a genuine divergence blows straight past this.
        """
        py_pe, py_cri, *_ = rr.risk_days_bs(**_KW, use_go=False)
        go_pe, go_cri, *_ = rr.risk_days_bs(**_KW, use_go=True)
        assert py_pe == pytest.approx(go_pe, rel=0.05)
        assert py_cri[0] == pytest.approx(go_cri[0], rel=0.10)
        assert py_cri[1] == pytest.approx(go_cri[1], rel=0.10)
