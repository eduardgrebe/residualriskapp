# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025-2026  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""UI-level checks for the Streamlit estimator page via ``streamlit.testing``.

These render the page in-process. The bootstrap only runs on a button click,
which these tests do not trigger, so they never reach ``ProcessPoolExecutor``
and are safe inside the sandbox (no ``@pytest.mark.multiprocessing`` marker).
"""

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from residualrisk import find_go_binary

_ESTIMATOR = str(Path(__file__).resolve().parent.parent / "estimator.py")


def _run_with_prep():
    """Render the estimator with both PrEP scenarios enabled so the oPrEP/iPrEP
    serology parameter widgets are present."""
    at = AppTest.from_file(_ESTIMATOR, default_timeout=120).run()
    assert not at.exception, at.exception
    for cb in at.checkbox:
        if "PrEP breakthrough" in (cb.label or ""):
            cb.set_value(True)
    at.run()
    assert not at.exception, at.exception
    return at


class TestSerologyWeibullLabels:
    """The serology Weibull term is ``S(t) = exp(-((t - min) / alpha) ** beta)``,
    so ``alpha`` is the Weibull **scale** and ``beta`` is the **shape** (and the
    defaults — alpha ~= 50-91 days, beta ~= 1-3 — confirm it). The UI labels and
    help text must match; they were previously swapped.
    """

    def test_alpha_labelled_scale_and_beta_shape(self):
        at = _run_with_prep()
        labels = [n.label for n in at.number_input]
        # One oPrEP + one iPrEP widget for each parameter.
        assert labels.count("Seroconversion Weibull scale (α)") == 2, labels
        assert labels.count("Seroconversion Weibull shape (β)") == 2, labels
        # The old swapped labels must be gone.
        assert "Seroconversion Weibull shape (α)" not in labels
        assert "Seroconversion Weibull scale (β)" not in labels


class TestOffsetNotExposed:
    """The Offset control is deliberately gone. After tcrit was retargeted to the
    plateau's central level, the offset is *exactly* a set-point multiplier —
    (set_point, o, a) is the same model as (set_point*o, 1, a/o) — so it added no
    expressive power while being a first-order lever on the answer (sweeping it
    0.5–2.0 moved the RDE ~5x). The set point is the parameter to vary: clinical
    units, its own bootstrap range, and it drives tcrit correctly.

    Removing it also makes two former UI bugs structurally unreachable: offset=0 built
    the amplitude-range slider with min==max (StreamlitAPIException on render), and
    a > offset raised ValueError from the library on Run.
    """

    def test_no_offset_input_rendered(self):
        at = _run_with_prep()
        assert not [n for n in at.number_input if n.label == "Offset"], (
            "the Offset control is back — see this class's docstring before restoring it"
        )

    def test_amplitude_capped_at_one(self):
        """With the offset pinned at 1, a > offset is unreachable: the amplitude
        widget must not let a exceed it."""
        at = _run_with_prep()
        amps = [n for n in at.number_input if n.label == "Amplitude (a)"]
        assert amps, "Amplitude (a) input not found"
        assert all(a.max == 1.0 for a in amps), [a.max for a in amps]

    def test_vary_sin_renders_without_crash(self):
        at = _run_with_prep()
        for cb in at.checkbox:
            if "Vary sinusoidal" in (cb.label or ""):
                cb.set_value(True)
        at.run()
        assert not at.exception, at.exception
        # the amplitude-range slider was built (min < max) rather than crashing
        sliders = [s for s in at.slider if "Amplitude (a) range" in (s.label or "")]
        assert sliders
        assert all(s.max == 1.0 for s in sliders), [s.max for s in sliders]


class TestInvalidPrepInputsCaught:
    """A library ValueError on Run must surface as st.sidebar.error, not a full-page
    traceback — the mechanistic branch catches it, mirroring the lookback branch
    (regression for estimator.py:816).

    This used to be provoked with a > offset. The Offset control is gone (see
    TestOffsetNotExposed) and the amplitude is capped at 1, so that input is no longer
    reachable; an inverted seroconversion window (min > max) is, and exercises the same
    catch — the point of the test is the error *path*, not the specific parameter.
    """

    def test_invalid_prep_params_show_error_not_traceback(self):
        if find_go_binary() is None:
            pytest.skip("needs the Go binary so the baseline runs without ProcessPoolExecutor")
        at = _run_with_prep()
        for s in at.select_slider:
            if "number of simulations" in (s.label or ""):
                s.set_value(1000)  # smallest count → fast baseline before the PrEP error
        for n in at.number_input:
            # ser_max <= ser_min → the library raises on Run.
            if n.label == "Time to seroconversion min (days)":
                n.set_value(200)
            elif n.label == "Time to seroconversion max (days)":
                n.set_value(100)
        at.run()
        for b in at.button:
            if (b.label or "") in ("Run simulations", "Calculate RDEs"):
                b.click()
        at.run()
        assert not at.exception, at.exception  # no full-page traceback
        errors = [e.value for e in at.sidebar.error]
        assert any("Invalid parameters" in m for m in errors), errors


class TestPointEstimateDefault:
    """The reported RDE point estimate defaults to the bootstrap mode (always inside
    the CrI), not 'primary parameters' — whose plug-in value can land in the far
    right tail of a skewed RDE distribution (even above the upper CrI). Picking
    'primary parameters' explicitly surfaces a tail-caveat warning.
    """

    def _pe_selectbox(self, at):
        for s in at.selectbox:
            if "point estimate of RDEs" in (s.label or ""):
                return s
        raise AssertionError("point-estimate selectbox not found")

    def test_default_method_is_mode(self):
        at = AppTest.from_file(_ESTIMATOR, default_timeout=120).run()
        assert not at.exception, at.exception
        assert self._pe_selectbox(at).value == "mode"

    def test_primary_parameters_shows_tail_warning(self):
        at = AppTest.from_file(_ESTIMATOR, default_timeout=120).run()
        assert not at.exception, at.exception
        self._pe_selectbox(at).set_value("primary parameters")
        at.run()
        assert not at.exception, at.exception
        assert any("far right tail" in (w.value or "") for w in at.warning), [
            w.value for w in at.warning
        ]


class TestPrereleaseBanner:
    """The header banner (``estimator.py`` ``_prerelease_notice``) warns on
    non-stable builds: ``.dev`` is checked first (so it wins over a coexisting
    ``bN``), then ``aN``/``bN``/``rcN`` → alpha/beta/release candidate; a clean
    ``X.Y.Z`` release shows no banner."""

    @pytest.mark.parametrize(
        "version, needle",
        [
            ("1.1.0b3.dev1", "unstable test build"),
            ("1.1.0a7", "alpha release"),
            ("1.1.0b3", "beta release"),
            ("1.1.0rc1", "release candidate"),
        ],
    )
    def test_prerelease_version_shows_banner(self, monkeypatch, version, needle):
        import residualrisk

        monkeypatch.setattr(residualrisk, "__version__", version)
        at = AppTest.from_file(_ESTIMATOR, default_timeout=120).run()
        assert not at.exception, at.exception
        assert any(needle in (w.value or "") for w in at.warning), [
            w.value for w in at.warning
        ]

    def test_stable_version_shows_no_banner(self, monkeypatch):
        import residualrisk

        monkeypatch.setattr(residualrisk, "__version__", "1.1.0")
        at = AppTest.from_file(_ESTIMATOR, default_timeout=120).run()
        assert not at.exception, at.exception
        assert not any(
            "Use at your own risk" in (w.value or "") for w in at.warning
        ), [w.value for w in at.warning]
