"""Hatchling build hook — cross-compile the Go accelerator for each target platform
and bundle the binaries into the wheel (``residualrisk/_bin/``), so
``pip install residualrisk`` ships a ready-to-run accelerator with no Go toolchain
and no separate build step on the user's machine.

This produces a *fat* wheel: all targets in one ``py3-none-any`` wheel, and
``residualrisk._go.find_go_binary`` selects (and version-checks) the right one at
runtime. For per-platform wheels instead, build a single target and set
``build_data['tag']`` + ``build_data['pure_python'] = False`` (see the note at the
end of ``initialize``).

Cross-compilation is pure-Go (``CGO_ENABLED=0``), so one machine builds every
target and the Linux binaries are static (portable across distros; manylinux is
trivial). The Go toolchain is needed at *build* time only. If ``go`` is absent the
hook skips bundling and produces a pure-Python wheel; ``find_go_binary`` then falls
back to a ``$RESIDUALRISK_GO_BINARY`` / system binary or the pure-Python engine.
"""

import os
import shutil
import subprocess
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface

# (GOOS, GOARCH, exe-suffix). Static builds → no libc dependency → portable.
_TARGETS = [
    ("linux", "amd64", ""),
    ("linux", "arm64", ""),
    ("darwin", "amd64", ""),
    ("darwin", "arm64", ""),
    ("windows", "amd64", ".exe"),
    ("windows", "arm64", ".exe"),
]


class GoBinaryBuildHook(BuildHookInterface):
    PLUGIN_NAME = "custom"

    def initialize(self, version, build_data):
        if shutil.which("go") is None:
            self.app.display_warning(
                "Go toolchain not found — building a pure-Python wheel with no "
                "bundled accelerator (find_go_binary falls back at runtime)."
            )
            return

        root = Path(self.root)
        go_dir = root / "go"
        bin_dir = root / "residualrisk" / "_bin"
        bin_dir.mkdir(parents=True, exist_ok=True)

        for goos, goarch, suffix in _TARGETS:
            name = f"riskdays_go-{goos}-{goarch}{suffix}"
            out = bin_dir / name
            env = {**os.environ, "GOOS": goos, "GOARCH": goarch, "CGO_ENABLED": "0"}
            subprocess.run(
                ["go", "build", "-trimpath", "-ldflags", "-s -w", "-o", str(out), "."],
                cwd=str(go_dir),
                env=env,
                check=True,
            )
            out.chmod(0o755)
            # Force the binary into the wheel even though residualrisk/_bin/ is
            # git-ignored (hatchling's default file selection follows VCS).
            build_data["force_include"][str(out)] = f"residualrisk/_bin/{name}"

        # Fat wheel: keep the default py3-none-any tag. For per-platform wheels,
        # build a single (GOOS, GOARCH) above and set:
        #     build_data["pure_python"] = False
        #     build_data["tag"] = "py3-none-manylinux2014_x86_64"   # etc.
