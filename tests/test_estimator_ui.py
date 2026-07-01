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

from streamlit.testing.v1 import AppTest

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
