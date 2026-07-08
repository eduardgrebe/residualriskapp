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


class TestOffsetFloorPreventsSliderCrash:
    """Offset=0 with 'Vary sinusoidal oscillation parameters' enabled built the
    amplitude-range slider with min==max (== the offset), raising
    StreamlitAPIException on render. Offset is now floored at 0.05 so that state
    is unreachable (regression for estimator.py:855).
    """

    def test_offset_input_floored_above_zero(self):
        at = _run_with_prep()
        offsets = [n for n in at.number_input if n.label == "Offset"]
        assert offsets, "Offset input not found"
        assert all(o.min and o.min > 0.0 for o in offsets), [o.min for o in offsets]

    def test_vary_sin_at_minimum_offset_renders_without_crash(self):
        at = _run_with_prep()
        for cb in at.checkbox:
            if "Vary sinusoidal" in (cb.label or ""):
                cb.set_value(True)
        for n in at.number_input:
            if n.label == "Offset":
                n.set_value(n.min)  # smallest allowed offset
        at.run()
        assert not at.exception, at.exception
        # the amplitude-range slider was built (min < max) rather than crashing
        assert any("Amplitude (a) range" in (s.label or "") for s in at.slider)


class TestInvalidPrepInputsCaught:
    """a > offset makes the library raise ValueError on Run. The mechanistic
    branch now catches it and shows st.sidebar.error instead of a full-page
    traceback — mirroring the lookback branch (regression for estimator.py:816).
    """

    def test_amplitude_over_offset_shows_error_not_traceback(self):
        if find_go_binary() is None:
            pytest.skip("needs the Go binary so the baseline runs without ProcessPoolExecutor")
        at = _run_with_prep()
        for s in at.select_slider:
            if "number of simulations" in (s.label or ""):
                s.set_value(1000)  # smallest count → fast baseline before the PrEP error
        for n in at.number_input:
            if n.label == "Offset":
                n.set_value(0.5)
            elif n.label == "Amplitude (a)":
                n.set_value(0.9)  # a (0.9) > offset (0.5) → invalid
        at.run()
        for b in at.button:
            if (b.label or "") in ("Run simulations", "Calculate RDEs"):
                b.click()
        at.run()
        assert not at.exception, at.exception  # no full-page traceback
        errors = [e.value for e in at.sidebar.error]
        assert any("Invalid parameters" in m for m in errors), errors
