#!/usr/bin/env python3
# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025-2026  Vitalant and Eduard Grebe Consulting
"""Smoke-test the wheel-bundled Go binary for the *current* platform.

Extracts the platform's binary out of a built wheel (dist/*.whl), then checks it
(1) runs and reports the expected version, and (2) completes a small bootstrap and
emits a result. Used by the GitHub Actions cross-OS job to validate the
Linux-cross-compiled macOS/Windows binaries actually execute on real macOS/Windows
runners (which Codeberg's hosted Linux runners cannot do).

Usage:  python scripts/smoke_bundled_binary.py [path/to/wheel]
"""

import glob
import json
import platform
import subprocess
import sys
import zipfile
from pathlib import Path

wheel = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob("dist/*.whl"))[-1]
# residualrisk-<version>-py3-none-any.whl
expected_version = Path(wheel).name.split("-")[1]

goos = {"linux": "linux", "darwin": "darwin", "windows": "windows"}[
    platform.system().lower()
]
goarch = {
    "x86_64": "amd64", "amd64": "amd64", "arm64": "arm64", "aarch64": "arm64",
}[platform.machine().lower()]
suffix = ".exe" if goos == "windows" else ""
member = f"residualrisk/_bin/riskdays_go-{goos}-{goarch}{suffix}"

dest = Path("smoke_extracted")
with zipfile.ZipFile(wheel) as z:
    z.extract(member, dest)
binary = dest / member
if suffix == "":
    binary.chmod(0o755)

# 1. Runs and matches the wheel version.
ver = subprocess.run([str(binary), "--version"], capture_output=True, text=True)
assert ver.returncode == 0, f"--version failed: {ver.stderr}"
assert ver.stdout.strip() == expected_version, (
    f"binary version {ver.stdout.strip()!r} != wheel version {expected_version!r}"
)

# 2. Completes a small computation (last stdout line is the result JSON).
payload = json.dumps({
    "k": 0.000673, "doubling_time": 0.8542, "doubling_time_norm_sd": 0.2813,
    "lod50": 2.73, "lod50_sd": 0.53, "lod95_lod50_ratio": 3.5,
    "volume_transfused": 200, "volume_transfused_min": 100,
    "volume_transfused_max": 340, "pool_size": 16, "retests": 1,
    "k_invgamma_alpha": 2.0, "k_invgamma_beta": 0.002019, "n_bs": 1000,
})
run = subprocess.run([str(binary)], input=payload, capture_output=True, text=True, timeout=120)
assert run.returncode == 0, f"computation failed: {run.stderr}"
results = [
    json.loads(ln)
    for ln in run.stdout.strip().splitlines()
    if ln.strip() and '"progress"' not in ln
]
assert results, "computation emitted no non-progress result line"

print(f"OK  {goos}/{goarch}  version {expected_version}  (bundled binary runs and computes)")
