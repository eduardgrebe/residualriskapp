#!/usr/bin/env python3
# This script written by Deepseek V4.1 Flash
# Convenience script only, not a core part of the software
# Copyright (C) 2026 Eduard Grebe Consulting
# Under supervision of Eduard Grebe <eduard@grebe.consulting>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version. See <https://www.gnu.org/licenses/>.

"""Verify that every file in an installed virtual environment still matches the
hash recorded in its own ``.dist-info/RECORD``.

Why this exists
---------------
``uv`` (like pip) tracks which *versions* are installed, not the *content* of the
installed files. If an upgrade leaves an old file behind under a name that is
unchanged between the two versions, the environment stays silently broken while
``uv sync`` reports success and ``uv lock --check`` passes --
``scripts/audit_deps.sh`` checks the *lockfile*, never the installed tree.

This is not hypothetical. In August 2026 a pyarrow 24 -> 25 upgrade left the
24.x extension modules (``lib.cpython-314-darwin.so`` and friends) installed
under the 25.0.1 dist-info, so ``import pyarrow`` died with::

    Library not loaded: @rpath/libarrow_python.2400.dylib

Every ``uv sync`` since -- including ``uv sync --upgrade`` -- reported success,
and the breakage surfaced a month later as an unrelated-looking Streamlit UI
test failure (Streamlit hashes polars frames via ``to_arrow()``, which imports
pyarrow). The same event left stale files in nine other packages (pandas, idna,
packaging, starlette, traitlets, ...) and, because uv's *cache* carried the same
contamination, ``uv sync --reinstall`` would not have repaired it -- only
``uv cache clean`` followed by a rebuild did.

Hash-checking the tree takes ~2 s and turns that class of silent corruption into
a loud failure.

Usage
-----
    python scripts/verify_venv.py                # ./.venv beside the repo root
    python scripts/verify_venv.py --venv PATH    # any other environment

Exit status: 0 = clean, 1 = problems found, 2 = no environment to check.
"""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import sys
from dataclasses import dataclass, field
from pathlib import Path

# Finder metadata is not package content. A few wheels ship a `.DS_Store`
# (observed in pydeck 0.9.1), and macOS rewrites it whenever the directory is
# browsed, which would otherwise show up as a permanent, harmless mismatch.
IGNORED_NAMES = frozenset({".DS_Store", ".localized"})

_CHUNK = 1024 * 1024


def find_site_packages(venv: Path) -> Path | None:
    """Return the venv's site-packages directory, or ``None`` if there is none."""
    for pattern in ("lib/python*/site-packages", "Lib/site-packages"):
        for candidate in sorted(venv.glob(pattern)):
            if candidate.is_dir():
                return candidate
    return None


def hash_file(path: Path) -> str:
    """Return the urlsafe-base64 sha256 of *path*, in RECORD's encoding."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_CHUNK), b""):
            digest.update(chunk)
    return base64.urlsafe_b64encode(digest.digest()).rstrip(b"=").decode()


def recorded_hash(value: str) -> str | None:
    """Extract the digest from a RECORD hash cell (``sha256=<b64>``).

    Returns ``None`` for entries RECORD deliberately leaves unhashed (``.pyc``
    files, ``RECORD`` itself) and for algorithms we do not implement.
    """
    algorithm, _, encoded = value.partition("=")
    if algorithm == "sha256" and encoded:
        return encoded
    return None


@dataclass
class PackageReport:
    """Problems found in one installed distribution."""

    package: str
    mismatched: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)

    @property
    def total(self) -> int:
        return len(self.mismatched) + len(self.missing)


def verify(site_packages: Path) -> tuple[list[PackageReport], int]:
    """Hash every RECORD-listed file under *site_packages*.

    Returns the reports for packages with problems (empty if clean) and the
    number of files whose hash was actually checked.
    """
    root = site_packages.resolve()
    reports: list[PackageReport] = []
    checked = 0

    for dist_info in sorted(site_packages.glob("*.dist-info")):
        record = dist_info / "RECORD"
        if not record.is_file():
            continue
        report = PackageReport(package=dist_info.name[: -len(".dist-info")])
        with record.open(newline="", encoding="utf-8") as handle:
            for row in csv.reader(handle):
                if len(row) < 2:
                    continue
                relative, digest = row[0], recorded_hash(row[1])
                if digest is None or Path(relative).name in IGNORED_NAMES:
                    continue
                target = (site_packages / relative).resolve()
                # RECORD paths are relative to site-packages; never follow one
                # that escapes it.
                if not target.is_relative_to(root):
                    continue
                if not target.is_file():
                    report.missing.append(relative)
                    continue
                checked += 1
                if hash_file(target) != digest:
                    report.mismatched.append(relative)
        if report.total:
            reports.append(report)

    return reports, checked


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Verify an installed venv against each package's RECORD hashes. "
            "Detects environments left stale by an upgrade (see module docstring)."
        )
    )
    parser.add_argument(
        "--venv",
        type=Path,
        default=Path(__file__).resolve().parents[1] / ".venv",
        help="virtual environment to check (default: ./.venv beside the repo root)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="print only the summary (and any problems), not progress",
    )
    parser.add_argument(
        "--max-per-package",
        type=int,
        default=5,
        help="how many offending paths to list per package (default: 5)",
    )
    args = parser.parse_args(argv)

    venv: Path = args.venv
    if not venv.is_dir():
        print(f"ERROR: no virtual environment at {venv}", file=sys.stderr)
        return 2
    site_packages = find_site_packages(venv)
    if site_packages is None:
        print(f"ERROR: no site-packages directory under {venv}", file=sys.stderr)
        return 2

    if not args.quiet:
        print(f"Verifying {site_packages}")
        print("against the hashes in each package's RECORD ...")

    reports, checked = verify(site_packages)

    if not reports:
        print(f"OK: {checked} files verified — all match their RECORD.")
        return 0

    problems = sum(report.total for report in reports)
    print(
        f"FAIL: {problems} problem file(s) in {len(reports)} package(s) "
        f"({checked} files verified)."
    )
    for report in reports:
        print(
            f"  {report.package}: {len(report.mismatched)} with stale/wrong "
            f"content, {len(report.missing)} missing"
        )
        for relative in report.mismatched[: args.max_per_package]:
            print(f"      WRONG    {relative}")
        for relative in report.missing[: args.max_per_package]:
            print(f"      MISSING  {relative}")
        hidden = report.total - min(report.total, 2 * args.max_per_package)
        if hidden > 0:
            print(f"      ... and {hidden} more")

    print(
        "\nThe environment does not match uv.lock's packages. uv will not notice\n"
        "by itself (it tracks versions, not file content), and a contaminated\n"
        "cache reproduces the same files, so clean it before rebuilding:\n"
        "\n"
        "    uv cache clean && rm -rf .venv && uv sync\n"
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
