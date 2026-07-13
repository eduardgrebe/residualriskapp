# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025-2026  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Resolution of the Go binary in ``residualrisk._go``: the platform →
(GOOS, GOARCH) mapping (including Windows' uppercase machine strings, which a naive
mapping misses), the version-skew guard on the bundled binary, and the runnability
smoke test applied to every candidate."""

import logging
import stat

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


# ---------------------------------------------------------------------------
# Runnability smoke test (_binary_runs / find_go_binary)
#
# A binary can EXIST yet be unable to run: built for the wrong architecture (an
# amd64 build on arm64), corrupt, non-executable, missing a shared library.
# find_go_binary() previously returned such a file on mere existence, so the app
# reported Go acceleration as available (estimator.py's used_go gate is
# `find_go_binary() is not None`) and only the first dispatch discovered otherwise —
# silently degrading to the 10-50x slower Python engine. Every candidate is now
# smoke-tested with `--version`.
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _clear_binary_cache():
    """_binary_runs caches per path for the life of the process; isolate tests."""
    _go._binary_ok.clear()
    yield
    _go._binary_ok.clear()


def _fake_binary(path, script, executable):
    path.write_text(script)
    if executable:
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)


def test_binary_runs_true_for_real_binary():
    binary = _go.find_go_binary()
    if binary is None:
        pytest.skip("no Go binary available")
    assert _go._binary_runs(binary) is True


def test_binary_runs_false_for_missing_path():
    assert _go._binary_runs("/nonexistent/riskdays_go") is False


def test_binary_runs_false_for_non_executable(tmp_path):
    """Present but not executable — e.g. a botched copy or a lost permission bit."""
    p = _fake_binary(tmp_path / "riskdays_go", "not a real binary", executable=False)
    assert _go._binary_runs(p) is False


def test_binary_runs_false_when_it_exits_nonzero(tmp_path):
    """Stands in for the wrong-architecture / missing-shared-library case: the file is
    executable, but actually running it fails."""
    p = _fake_binary(tmp_path / "riskdays_go", "#!/bin/sh\nexit 1\n", executable=True)
    assert _go._binary_runs(p) is False


def test_binary_runs_is_cached(tmp_path, monkeypatch):
    """One subprocess per path per process: find_go_binary() is called at load, at
    dispatch, and for the UI backend gate, so an uncached smoke test would be costly."""
    p = _fake_binary(tmp_path / "riskdays_go", "#!/bin/sh\nexit 0\n", executable=True)
    calls = []
    real_run = _go.subprocess.run

    def counting_run(*args, **kwargs):
        calls.append(args)
        return real_run(*args, **kwargs)

    monkeypatch.setattr(_go.subprocess, "run", counting_run)
    assert _go._binary_runs(p) is True
    assert _go._binary_runs(p) is True
    assert len(calls) == 1, "smoke test should run once per path, then be cached"


def test_find_go_binary_skips_broken_override(tmp_path, monkeypatch, caplog):
    """The regression this guards: an existing-but-unrunnable binary must NOT be
    handed back. It is skipped (a working candidate wins instead), and because the
    override was set *explicitly*, we say so rather than silently substituting."""
    broken = _fake_binary(
        tmp_path / "riskdays_go", "#!/bin/sh\nexit 1\n", executable=True
    )
    monkeypatch.setenv("RESIDUALRISK_GO_BINARY", broken)
    with caplog.at_level(logging.WARNING, logger="residualrisk._go"):
        found = _go.find_go_binary()
    assert found != broken
    assert "does not run" in caplog.text


def test_find_go_binary_returns_none_when_nothing_runs(tmp_path, monkeypatch):
    """With no runnable binary anywhere, degrade *honestly* — return None so the app
    reports Go as unavailable, rather than a path that will fail at first dispatch."""
    broken = _fake_binary(
        tmp_path / "riskdays_go", "#!/bin/sh\nexit 1\n", executable=True
    )
    monkeypatch.setenv("RESIDUALRISK_GO_BINARY", broken)
    monkeypatch.setattr(_go, "_bundled_go_binary", lambda: None)
    # Only the broken path exists, and nothing executes successfully.
    monkeypatch.setattr(_go.Path, "exists", lambda self: str(self) == broken)
    monkeypatch.setattr(
        _go.subprocess,
        "run",
        lambda *a, **kw: type("R", (), {"returncode": 1, "stdout": ""})(),
    )
    assert _go.find_go_binary() is None
