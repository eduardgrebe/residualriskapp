# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025-2026  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Resolution of the wheel-bundled Go binary in ``residualrisk._go``: the
platform → (GOOS, GOARCH) mapping (including Windows' uppercase machine strings,
which a naive mapping misses) and the version-skew guard."""

import pytest

import residualrisk._go as _go


@pytest.mark.parametrize(
    "system, machine, expected",
    [
        ("Linux", "x86_64", ("linux", "amd64")),
        ("Linux", "aarch64", ("linux", "arm64")),
        ("Darwin", "arm64", ("darwin", "arm64")),
        ("Darwin", "x86_64", ("darwin", "amd64")),
        ("Windows", "AMD64", ("windows", "amd64")),  # Windows reports uppercase
        ("Windows", "ARM64", ("windows", "arm64")),  # Windows on ARM
        ("SunOS", "sparc", (None, None)),  # unsupported → no bundle
    ],
)
def test_current_go_platform(monkeypatch, system, machine, expected):
    monkeypatch.setattr("platform.system", lambda: system)
    monkeypatch.setattr("platform.machine", lambda: machine)
    assert _go._current_go_platform() == expected


def test_binary_version_matches_real():
    """A resolvable binary runs and its --version equals the library version."""
    binary = _go.find_go_binary()
    if binary is None:
        pytest.skip("no Go binary available")
    assert _go._binary_version_matches(binary) is True


def test_binary_version_matches_rejects_missing():
    """A missing/non-runnable path is rejected (→ pure-Python fallback)."""
    assert _go._binary_version_matches("/nonexistent/riskdays_go") is False
