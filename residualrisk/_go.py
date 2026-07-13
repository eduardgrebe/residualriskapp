# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025  Vitalant and Eduard Grebe Consulting
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

"""
Python wrapper for calling the Go implementation of risk_days_bs()
"""

import json
import logging
import queue
import struct
import subprocess
import threading
from pathlib import Path

import numpy as np
import polars as pl


# Cache of bundled-binary validity checks (abs path -> bool), so the --version
# handshake runs at most once per process.
logger = logging.getLogger(__name__)

_bundled_validity: dict = {}

# path -> bool: does this binary actually execute? (see _binary_runs). Cached so the
# smoke test costs one subprocess per path per process, not one per call — the app
# calls find_go_binary() at load, at dispatch, and for the UI backend gate.
_binary_ok: dict = {}


def _current_go_platform():
    """(GOOS, GOARCH) for the running interpreter, or (None, None) if we don't ship a
    binary for it. platform.machine() varies by OS and case (x86_64 / AMD64 / arm64 /
    ARM64 / aarch64), so normalise via lower()."""
    import platform

    goos = {"linux": "linux", "darwin": "darwin", "windows": "windows"}.get(
        platform.system().lower()
    )
    goarch = {
        "x86_64": "amd64", "amd64": "amd64",
        "arm64": "arm64", "aarch64": "arm64",
    }.get(platform.machine().lower())
    return goos, goarch


def _bundled_go_binary():
    """Path to the platform's binary bundled in the wheel (residualrisk/_bin/), or
    None. It is used only if it actually runs and its version matches the library —
    a stale/mismatched/incompatible bundle falls through to the pure-Python path
    rather than silently producing wrong results (result cached per path)."""
    import os

    goos, goarch = _current_go_platform()
    if not (goos and goarch):
        return None
    suffix = ".exe" if goos == "windows" else ""
    path = Path(__file__).parent / "_bin" / f"riskdays_go-{goos}-{goarch}{suffix}"
    if not path.exists():
        return None
    key = str(path)
    if key not in _bundled_validity:
        if not os.access(path, os.X_OK):
            try:
                path.chmod(0o755)
            except OSError:
                pass
        _bundled_validity[key] = _binary_version_matches(key)
    return key if _bundled_validity[key] else None


def _binary_version_matches(path):
    """True iff ``path --version`` runs and equals the library ``__version__``."""
    try:
        import residualrisk

        result = subprocess.run(
            [path, "--version"], capture_output=True, text=True, timeout=10
        )
        return (
            result.returncode == 0
            and result.stdout.strip() == residualrisk.__version__
        )
    except Exception:
        return False


def _binary_runs(path):
    """True iff ``path --version`` actually executes (exit 0). Cached per path.

    Existence is not enough. A binary can be present but unrunnable: built for the
    wrong architecture (an amd64 build on arm64), truncated/corrupt, non-executable,
    or missing a shared library. Such a file is still "found", so
    ``find_go_binary()`` hands it back, ``estimator.py``'s ``used_go`` gate reports
    Go acceleration as available, and only the *first dispatch* discovers otherwise —
    silently degrading to the 10-50x slower Python engine.

    Smoke-testing each candidate once means an unusable binary is skipped in favour
    of the next one (or None), so the app degrades *honestly*: it says Go is
    unavailable, warns, and uses Python.

    Note this deliberately does NOT require a version match — unlike the bundled
    binary (see ``_bundled_go_binary``), where a version mismatch signals a stale
    wheel and must be rejected. A developer's ``go/bin/`` build may legitimately lag
    the library version; refusing to run it would be hostile.
    """
    if path not in _binary_ok:
        try:
            result = subprocess.run(
                [path, "--version"], capture_output=True, text=True, timeout=10
            )
            _binary_ok[path] = result.returncode == 0
        except Exception:
            _binary_ok[path] = False
    return _binary_ok[path]


def find_go_binary():
    """
    Find the riskdays_go Go binary.

    Searches in order:
    1. $RESIDUALRISK_GO_BINARY env var (explicit override)
    2. go/bin/riskdays_go (repo root, relative to this package)
    3. residualrisk/_bin/riskdays_go-<goos>-<goarch> (bundled in the wheel)
    4. ~/.local/bin/riskdays_go (sudo-free user-install location)
    5. riskdays_go in PATH

    Returns:
        Path to binary or None if not found
    """
    import os

    # Every candidate is smoke-tested with `--version` (_binary_runs, cached): a
    # binary that exists but cannot execute — wrong architecture, corrupt,
    # non-executable, missing shared library — must NOT be returned, or the app
    # reports Go acceleration as available and only discovers otherwise on the first
    # dispatch, silently degrading to the far slower Python engine. Skipping it lets
    # the next candidate win, and if none run we return None and degrade honestly.
    env_override = os.environ.get("RESIDUALRISK_GO_BINARY")
    if env_override and Path(env_override).exists():
        if _binary_runs(env_override):
            return env_override
        # An *explicit* override that cannot run is worth saying out loud: silently
        # falling through to some other binary would be baffling for the operator who
        # deliberately pointed us at this one.
        logger.warning(
            "RESIDUALRISK_GO_BINARY=%s exists but does not run (`--version` failed) — "
            "ignoring it and continuing the search.",
            env_override,
        )

    # Package lives at <repo>/residualrisk/_go.py; binary at <repo>/go/bin/
    repo_root = Path(__file__).parent.parent
    relative_binary = repo_root / "go" / "bin" / "riskdays_go"
    if relative_binary.exists() and _binary_runs(str(relative_binary)):
        return str(relative_binary)

    # Bundled binary shipped inside the wheel (residualrisk/_bin/), selected by
    # platform and validated against the library version. See hatch_build.py.
    # Stricter than the others: it must run AND match __version__, since a mismatch
    # means a stale wheel rather than a dev build that merely lags.
    bundled = _bundled_go_binary()
    if bundled:
        return bundled

    # User-local install — the recommended manual-install location (no sudo).
    user_binary = Path.home() / ".local" / "bin" / "riskdays_go"
    if user_binary.exists() and _binary_runs(str(user_binary)):
        return str(user_binary)

    # Try PATH
    try:
        result = subprocess.run(
            ["which", "riskdays_go"], capture_output=True, text=True, check=False
        )
        if result.returncode == 0:
            found = result.stdout.strip()
            if found and _binary_runs(found):
                return found
    except Exception:
        pass

    return None


def mode_kde_go(
    data: "np.ndarray",
    n_grid: int = 100_000,
    cap: int | None = None,
    threads: int = 0,
) -> float:
    """
    Estimate the mode of a positive right-skewed distribution via KDE on the
    log scale, using the Go binary for speed (typically 30× faster than the
    pure-Python implementation for large posteriors).

    Parameters
    ----------
    data:
        1-D array of positive values (e.g. a k-parameter posterior sample).
    n_grid:
        Number of log-spaced grid points for the KDE, passed to the Go binary.
        Default 100 000 — the app's grid. The Go KDE is FFT-based, so this is fast,
        and the mode matches a 1 000 000-grid result to ~4 sig figs. The pure-Python
        fallback (``_kde_mode_log`` / ``mode_kde``) defaults to a coarser 5 000 for
        speed, since it is O(n_data × n_grid); the mode differs by only ~0.1%.
    cap:
        Maximum number of samples to use.  Data is pre-subsampled in Python
        (seed=42) before serialisation to keep the JSON payload small.
        ``None`` (default) means no cap — all samples are used.
    threads:
        Parallel goroutines; 0 → all CPU cores.

    Returns
    -------
    float
        Estimated mode, or ``None`` if the Go binary is unavailable.
    """
    go_bin = find_go_binary()
    if go_bin is None:
        return None

    arr = np.asarray(data, dtype=float)
    if cap is not None and len(arr) > cap:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(arr), size=cap, replace=False)
        arr = arr[idx]

    payload = json.dumps({"data": arr.tolist(), "n_grid": n_grid, "cap": 0, "threads": threads})
    try:
        result = subprocess.run(
            [go_bin, "--kde-mode"],
            input=payload,
            capture_output=True,
            text=True,
            check=True,
        )
        return json.loads(result.stdout)["mode"]
    except Exception:
        return None


def mode_hsm_go(data: "np.ndarray") -> float:
    """
    Estimate the mode using the Half-Sample Mode (HSM) algorithm via the Go
    binary.

    HSM iteratively finds the shortest interval containing half the data,
    recursing until 1–3 points remain.  It is bandwidth-free and outlier-robust.

    Parameters
    ----------
    data:
        1-D array of positive values.

    Returns
    -------
    float
        Estimated mode, or ``None`` if the Go binary is unavailable.
    """
    go_bin = find_go_binary()
    if go_bin is None:
        return None

    arr = np.asarray(data, dtype=float)
    payload = json.dumps({"data": arr.tolist()})
    try:
        result = subprocess.run(
            [go_bin, "--hsm-mode"],
            input=payload,
            capture_output=True,
            text=True,
            check=True,
        )
        return json.loads(result.stdout)["mode"]
    except Exception:
        return None


def risk_days_bs_go(
    k,
    doubling_time,
    doubling_time_norm_sd,
    lod50,
    lod50_sd,
    lod95_lod50_ratio,
    volume_transfused,
    volume_transfused_range,
    pool_size,
    retests,
    C0=0.00025,
    copies_per_virion=2,
    alpha=0.05,
    z=1.6449,
    k_posterior_sample=None,
    k_gamma_shape=None,
    k_gamma_scale=None,
    k_invgamma_alpha=None,
    k_invgamma_beta=None,
    k_invgamma_mode=None,
    k_lnmix_w=None,
    k_lnmix_mu1=None,
    k_lnmix_sigma1=None,
    k_lnmix_mu2=None,
    k_lnmix_sigma2=None,
    n_bs=10000,
    seed=126887,
    threads=None,
    point_estimate="primary parameters",
    mode_precision=2,
    progress=None,
    return_sim_df=False,
    # Keep `limits` LAST: both call sites in core.py pass everything up to
    # return_sim_df positionally, so inserting it earlier would silently shift them.
    limits=(-100, 500),
):
    """
    Go implementation of risk_days_bs() via CLI subprocess.

    This function has the same signature as the Python version in residualrisk.py
    and returns results in the same format.

    Parameters match the Python implementation exactly.

    Returns:
        Tuple: (rd_pe, rd_cri, rd_range, rdests, sim_df) or (rd_pe, rd_cri, rd_range, rdests, None)

    Raises:
        RuntimeError: If Go binary not found or execution fails
    """
    # Find Go binary
    binary_path = find_go_binary()
    if not binary_path:
        raise RuntimeError(
            "Go binary 'riskdays_go' not found. "
            "Please build it with 'cd go && make build' "
            "or install it with 'cd go && sudo make install'"
        )

    # Determine number of threads
    if threads is None:
        import multiprocessing

        threads = max(1, multiprocessing.cpu_count() - 1)

    # Build input JSON
    input_data = {
        "k": k,
        "doubling_time": doubling_time,
        "doubling_time_norm_sd": doubling_time_norm_sd,
        "lod50": lod50,
        "lod50_sd": lod50_sd,
        "lod95_lod50_ratio": lod95_lod50_ratio,
        "volume_transfused": volume_transfused,
        "volume_transfused_min": volume_transfused_range[0],
        "volume_transfused_max": volume_transfused_range[1],
        "pool_size": pool_size,
        "retests": retests,
        "c0": C0,
        "copies_per_virion": copies_per_virion,
        "alpha": alpha,
        "z": z,
        "n_bs": n_bs,
        "seed": seed,
        "threads": threads,
        "point_estimate": point_estimate,
        "mode_precision": mode_precision,
        "limits": [limits[0], limits[1]],
    }

    # Add k distribution parameters
    if k_posterior_sample is not None:
        input_data["k_posterior_sample"] = (
            k_posterior_sample.tolist()
            if hasattr(k_posterior_sample, "tolist")
            else list(k_posterior_sample)
        )
    elif k_gamma_shape is not None and k_gamma_scale is not None:
        input_data["k_gamma_shape"] = k_gamma_shape
        input_data["k_gamma_scale"] = k_gamma_scale
    elif k_invgamma_alpha is not None:
        # Resolve mode → beta before sending to Go
        _beta = k_invgamma_beta
        if _beta is None:
            if k_invgamma_mode is not None:
                _beta = k_invgamma_mode * (k_invgamma_alpha + 1)
            else:
                raise ValueError(
                    "k_invgamma_alpha requires k_invgamma_beta or k_invgamma_mode"
                )
        input_data["k_invgamma_alpha"] = k_invgamma_alpha
        input_data["k_invgamma_beta"] = _beta
    elif k_lnmix_w is not None:
        if any(p is None for p in [k_lnmix_mu1, k_lnmix_sigma1, k_lnmix_mu2, k_lnmix_sigma2]):
            raise ValueError(
                "All lnmix parameters (k_lnmix_w, mu1, sigma1, mu2, sigma2) must be provided together."
            )
        input_data["k_lnmix_w"] = k_lnmix_w
        input_data["k_lnmix_mu1"] = k_lnmix_mu1
        input_data["k_lnmix_sigma1"] = k_lnmix_sigma1
        input_data["k_lnmix_mu2"] = k_lnmix_mu2
        input_data["k_lnmix_sigma2"] = k_lnmix_sigma2
    else:
        raise ValueError(
            "Either k_posterior_sample, k_gamma parameters, k_invgamma parameters, or k_lnmix parameters must be provided"
        )

    # If sim_df requested, tell Go to return per-iteration params via binary format
    if return_sim_df:
        input_data["return_params"] = True

    # Run Go binary
    try:
        # Binary mode on stdout when return_params=True to read raw float64 arrays.
        # Text mode (current behaviour) otherwise.
        use_binary = return_sim_df
        process = subprocess.Popen(
            [binary_path],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=False,   # always bytes — we decode stderr lines manually
            bufsize=0,    # unbuffered for binary; progress lines are small
        )

        # Send input and close stdin
        input_json = json.dumps(input_data)
        process.stdin.write(input_json.encode())
        process.stdin.close()

        # Read both stdout and stderr in background threads to avoid deadlock
        # (With 1M+ simulations, both pipes can fill up and block the Go process)
        error_msg = "Unknown error"
        stderr_queue = queue.Queue()
        stdout_data = []

        def read_stderr():
            """Background thread to read stderr for progress updates (bytes → UTF-8 lines)"""
            try:
                for line_bytes in iter(process.stderr.readline, b""):
                    if not line_bytes:
                        break
                    stderr_queue.put(line_bytes.decode("utf-8", errors="replace"))
            except Exception as e:
                stderr_queue.put(f"STDERR_ERROR: {e}")
            finally:
                process.stderr.close()

        def read_stdout():
            """Background thread to read stdout as raw bytes"""
            try:
                stdout_data.append(process.stdout.read())
            except Exception as e:
                stdout_data.append(f"STDOUT_ERROR: {e}".encode())
            finally:
                process.stdout.close()

        # Start both reader threads
        stderr_thread = threading.Thread(target=read_stderr, daemon=True)
        stdout_thread = threading.Thread(target=read_stdout, daemon=True)
        stderr_thread.start()
        stdout_thread.start()

        # Process stderr messages from queue in real-time
        while stderr_thread.is_alive() or not stderr_queue.empty():
            try:
                line = stderr_queue.get(timeout=0.1)
                if line.startswith("STDERR_ERROR:"):
                    break
                try:
                    msg = json.loads(line.strip())
                    if msg.get("type") == "progress" and progress is not None:
                        # Update Streamlit progress bar in real-time
                        percent_value = int(msg["percent"] * 100)
                        progress.progress(
                            msg["percent"], text=f"Progress: {percent_value}%"
                        )
                    elif msg.get("type") == "error":
                        error_msg = msg.get("message", error_msg)
                except json.JSONDecodeError:
                    # Ignore non-JSON lines
                    pass
            except queue.Empty:
                continue

        # Wait for stdout thread to finish and process to complete
        stdout_thread.join()
        process.wait()

        # Check for errors
        if process.returncode != 0:
            raise RuntimeError(f"Go binary failed: {error_msg}")

        # Get stdout bytes
        raw = stdout_data[0] if stdout_data else b""
        if isinstance(raw, str) or (isinstance(raw, bytes) and raw.startswith(b"STDOUT_ERROR:")):
            raise RuntimeError(f"Failed to read stdout: {raw}")

        if use_binary:
            # Parse binary wire format:
            #   [8 bytes]  uint64 LE — length of JSON header
            #   [N bytes]  JSON header
            #   [rest]     column-major float64 LE arrays (iwp, k, doubling_time, lod50, volume_transfused)
            header_len = struct.unpack_from("<Q", raw, 0)[0]
            header = json.loads(raw[8:8 + header_len])
            n_cols = len(header["columns"])
            n_bs_actual = header["n_bs"]
            arrays = np.frombuffer(raw[8 + header_len:], dtype="<f8").reshape(n_cols, n_bs_actual)
            col_idx = {name: i for i, name in enumerate(header["columns"])}

            rd_pe = header["point_estimate"]
            rd_cri = tuple(header["credible_interval"])
            rd_range = tuple(header["range"])
            rdests = arrays[col_idx["iwp"]].tolist()

            # Build sim_df from the REAL per-iteration parameters returned by Go.
            sim_df = pl.DataFrame({
                "k":                 arrays[col_idx["k"]],
                "doubling_time":     arrays[col_idx["doubling_time"]],
                "lod50":             arrays[col_idx["lod50"]],
                "volume_transfused": arrays[col_idx["volume_transfused"]],
                "iwp":               arrays[col_idx["iwp"]],
            }).with_columns(
                pl.lit(copies_per_virion).alias("copies_per_virion"),
                pl.lit(C0).alias("C0"),
                pl.lit(pool_size).alias("pool_size"),
                pl.lit(lod95_lod50_ratio).alias("lod95_lod50_ratio"),
                (pl.col("lod50") * lod95_lod50_ratio).alias("lod95"),
                pl.lit(retests).alias("retests"),
                pl.lit(z).alias("z"),
                # The integration domain. Every other parameter is recorded, and
                # `limits` is caller-settable, so a Go-path frame must record it too
                # — the Python sim_df has always carried it. Pass the raw values (no
                # float cast) so the dtype follows the input types, as Python's does.
                pl.lit([limits[0], limits[1]]).alias("limits"),
                pl.lit(seed).alias("random_seed"),
            )
            return (rd_pe, rd_cri, rd_range, rdests, sim_df)
        else:
            # Standard JSON output path (return_sim_df=False)
            output = json.loads(raw.decode("utf-8"))
            rd_pe = output["point_estimate"]
            rd_cri = tuple(output["credible_interval"])
            rd_range = tuple(output["range"])
            rdests = output["simulations"]
            return (rd_pe, rd_cri, rd_range, rdests, None)

    except subprocess.SubprocessError as e:
        raise RuntimeError(f"Failed to run Go binary: {e}")
    except json.JSONDecodeError as e:
        raise RuntimeError(f"Failed to parse Go output: {e}. Output was: {raw!r:.200}")
    except Exception as e:
        raise RuntimeError(f"Unexpected error calling Go implementation: {e}")


def risk_days_prep_bs_go(
    k,
    doubling_time,
    doubling_time_norm_sd,
    lod50,
    lod50_sd,
    lod95_lod50_ratio,
    volume_transfused,
    volume_transfused_range,
    pool_size,
    retests,
    set_point=336,
    set_point_dist_uniform=(19.1, 2265),
    eclipse=7.0,
    eclipse_dist_uniform=(4.0, 10.0),
    a=0.7,
    b=0.6,
    offset=1,
    a_dist_uniform=None,
    b_dist_uniform=None,
    drug_effect=1.0,
    drug_effect_dist_uniform=None,
    ser_min=28.7,
    ser_max=250,
    ser_alpha=50.49434,
    ser_beta=1.15062,
    C0=0.00025,
    copies_per_virion=2,
    alpha=0.05,
    z=1.6449,
    k_posterior_sample=None,
    k_gamma_shape=None,
    k_gamma_scale=None,
    k_invgamma_alpha=None,
    k_invgamma_beta=None,
    k_invgamma_mode=None,
    k_lnmix_w=None,
    k_lnmix_mu1=None,
    k_lnmix_sigma1=None,
    k_lnmix_mu2=None,
    k_lnmix_sigma2=None,
    n_bs=10000,
    seed=126887,
    threads=None,
    point_estimate="primary parameters",
    mode_precision=2,
    progress=None,
    return_sim_df=False,
    # Keep `limits` LAST: both call sites in core.py pass everything up to
    # return_sim_df positionally, so inserting it earlier would silently shift them.
    limits=(-100, 500),
):
    """
    Go implementation of risk_days_prep_bs() via CLI subprocess.

    Signature matches the Python version in residualrisk/prep.py.
    Sets ``prep_mode: true`` in the Go binary input so the PrEP 3-phase viral
    dynamics and serology non-detection are used.

    Returns:
        Tuple: (rd_pe, rd_cri, rd_range, rdests, sim_df)
        sim_df is None unless return_sim_df=True.
    """
    binary_path = find_go_binary()
    if not binary_path:
        raise RuntimeError(
            "Go binary 'riskdays_go' not found. "
            "Please build it with 'cd go && make build' "
            "or install it with 'cd go && sudo make install'"
        )

    if threads is None:
        import multiprocessing
        threads = max(1, multiprocessing.cpu_count() - 1)

    # Build input JSON — baseline fields + prep_mode + PrEP-specific fields
    input_data = {
        "k": k,
        "doubling_time": doubling_time,
        "doubling_time_norm_sd": doubling_time_norm_sd,
        "lod50": lod50,
        "lod50_sd": lod50_sd,
        "lod95_lod50_ratio": lod95_lod50_ratio,
        "volume_transfused": volume_transfused,
        "volume_transfused_min": volume_transfused_range[0],
        "volume_transfused_max": volume_transfused_range[1],
        "pool_size": pool_size,
        "retests": retests,
        "c0": C0,
        "copies_per_virion": copies_per_virion,
        "alpha": alpha,
        "z": z,
        "n_bs": n_bs,
        "seed": seed,
        "threads": threads,
        "point_estimate": point_estimate,
        "mode_precision": mode_precision,
        "limits": [limits[0], limits[1]],
        # PrEP-specific
        "prep_mode": True,
        "set_point": set_point,
        "set_point_dist_uniform": list(set_point_dist_uniform),
        "eclipse": eclipse,
        "eclipse_dist_uniform": list(eclipse_dist_uniform),
        "a": a,
        "b": b,
        "offset": offset,
        "drug_effect": drug_effect,
        "ser_min": ser_min,
        "ser_max": ser_max,
        "ser_alpha": ser_alpha,
        "ser_beta": ser_beta,
    }
    # Optional uniform ranges for the sinusoidal params; omitted → Go treats
    # them as fixed at a / b (its [0,0] sentinel).
    if a_dist_uniform is not None:
        input_data["a_dist_uniform"] = list(a_dist_uniform)
    if b_dist_uniform is not None:
        input_data["b_dist_uniform"] = list(b_dist_uniform)
    # Omitted → Go treats drug_effect as fixed at the scalar (its [0,0] sentinel).
    if drug_effect_dist_uniform is not None:
        input_data["drug_effect_dist_uniform"] = list(drug_effect_dist_uniform)

    # Add k distribution parameters (same logic as baseline wrapper)
    if k_posterior_sample is not None:
        input_data["k_posterior_sample"] = (
            k_posterior_sample.tolist()
            if hasattr(k_posterior_sample, "tolist")
            else list(k_posterior_sample)
        )
    elif k_gamma_shape is not None and k_gamma_scale is not None:
        input_data["k_gamma_shape"] = k_gamma_shape
        input_data["k_gamma_scale"] = k_gamma_scale
    elif k_invgamma_alpha is not None:
        _beta = k_invgamma_beta
        if _beta is None:
            if k_invgamma_mode is not None:
                _beta = k_invgamma_mode * (k_invgamma_alpha + 1)
            else:
                raise ValueError(
                    "k_invgamma_alpha requires k_invgamma_beta or k_invgamma_mode"
                )
        input_data["k_invgamma_alpha"] = k_invgamma_alpha
        input_data["k_invgamma_beta"] = _beta
    elif k_lnmix_w is not None:
        if any(p is None for p in [k_lnmix_mu1, k_lnmix_sigma1, k_lnmix_mu2, k_lnmix_sigma2]):
            raise ValueError(
                "All lnmix parameters (k_lnmix_w, mu1, sigma1, mu2, sigma2) must be provided together."
            )
        input_data["k_lnmix_w"] = k_lnmix_w
        input_data["k_lnmix_mu1"] = k_lnmix_mu1
        input_data["k_lnmix_sigma1"] = k_lnmix_sigma1
        input_data["k_lnmix_mu2"] = k_lnmix_mu2
        input_data["k_lnmix_sigma2"] = k_lnmix_sigma2
    else:
        raise ValueError(
            "Either k_posterior_sample, k_gamma parameters, k_invgamma parameters, or k_lnmix parameters must be provided"
        )

    if return_sim_df:
        input_data["return_params"] = True

    try:
        use_binary = return_sim_df
        process = subprocess.Popen(
            [binary_path],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=False,
            bufsize=0,
        )

        input_json = json.dumps(input_data)
        process.stdin.write(input_json.encode())
        process.stdin.close()

        error_msg = "Unknown error"
        stderr_queue = queue.Queue()
        stdout_data = []

        def read_stderr():
            try:
                for line_bytes in iter(process.stderr.readline, b""):
                    if not line_bytes:
                        break
                    stderr_queue.put(line_bytes.decode("utf-8", errors="replace"))
            except Exception as e:
                stderr_queue.put(f"STDERR_ERROR: {e}")
            finally:
                process.stderr.close()

        def read_stdout():
            try:
                stdout_data.append(process.stdout.read())
            except Exception as e:
                stdout_data.append(f"STDOUT_ERROR: {e}".encode())
            finally:
                process.stdout.close()

        stderr_thread = threading.Thread(target=read_stderr, daemon=True)
        stdout_thread = threading.Thread(target=read_stdout, daemon=True)
        stderr_thread.start()
        stdout_thread.start()

        while stderr_thread.is_alive() or not stderr_queue.empty():
            try:
                line = stderr_queue.get(timeout=0.1)
                if line.startswith("STDERR_ERROR:"):
                    break
                try:
                    msg = json.loads(line.strip())
                    if msg.get("type") == "progress" and progress is not None:
                        percent_value = int(msg["percent"] * 100)
                        progress.progress(
                            msg["percent"], text=f"Progress: {percent_value}%"
                        )
                    elif msg.get("type") == "error":
                        error_msg = msg.get("message", error_msg)
                except json.JSONDecodeError:
                    pass
            except queue.Empty:
                continue

        stdout_thread.join()
        process.wait()

        if process.returncode != 0:
            raise RuntimeError(f"Go binary failed: {error_msg}")

        raw = stdout_data[0] if stdout_data else b""
        if isinstance(raw, str) or (isinstance(raw, bytes) and raw.startswith(b"STDOUT_ERROR:")):
            raise RuntimeError(f"Failed to read stdout: {raw}")

        if use_binary:
            header_len = struct.unpack_from("<Q", raw, 0)[0]
            header = json.loads(raw[8:8 + header_len])
            n_cols = len(header["columns"])
            n_bs_actual = header["n_bs"]
            arrays = np.frombuffer(raw[8 + header_len:], dtype="<f8").reshape(n_cols, n_bs_actual)
            col_idx = {name: i for i, name in enumerate(header["columns"])}

            rd_pe = header["point_estimate"]
            rd_cri = tuple(header["credible_interval"])
            rd_range = tuple(header["range"])
            rdests = arrays[col_idx["iwp"]].tolist()

            sim_df = pl.DataFrame({
                "k":                 arrays[col_idx["k"]],
                "doubling_time":     arrays[col_idx["doubling_time"]],
                "set_point":         arrays[col_idx["set_point"]],
                "eclipse":           arrays[col_idx["eclipse"]],
                "a":                 arrays[col_idx["a"]],
                "b":                 arrays[col_idx["b"]],
                "drug_effect":       arrays[col_idx["drug_effect"]],
                "lod50":             arrays[col_idx["lod50"]],
                "volume_transfused": arrays[col_idx["volume_transfused"]],
                "iwp":               arrays[col_idx["iwp"]],
            }).with_columns(
                pl.lit(copies_per_virion).alias("copies_per_virion"),
                pl.lit(C0).alias("C0"),
                pl.lit(pool_size).alias("pool_size"),
                pl.lit(lod95_lod50_ratio).alias("lod95_lod50_ratio"),
                (pl.col("lod50") * lod95_lod50_ratio).alias("lod95"),
                pl.lit(retests).alias("retests"),
                pl.lit(z).alias("z"),
                pl.lit(offset).alias("offset"),
                pl.lit(ser_min).alias("ser_min"),
                pl.lit(ser_max).alias("ser_max"),
                pl.lit(ser_alpha).alias("ser_alpha"),
                pl.lit(ser_beta).alias("ser_beta"),
                # The integration domain. Every other parameter is recorded, and
                # `limits` is caller-settable, so a Go-path frame must record it too
                # — the Python sim_df has always carried it. Pass the raw values (no
                # float cast) so the dtype follows the input types, as Python's does.
                pl.lit([limits[0], limits[1]]).alias("limits"),
                pl.lit(seed).alias("random_seed"),
            )
            return (rd_pe, rd_cri, rd_range, rdests, sim_df)
        else:
            output = json.loads(raw.decode("utf-8"))
            rd_pe = output["point_estimate"]
            rd_cri = tuple(output["credible_interval"])
            rd_range = tuple(output["range"])
            rdests = output["simulations"]
            return (rd_pe, rd_cri, rd_range, rdests, None)

    except subprocess.SubprocessError as e:
        raise RuntimeError(f"Failed to run Go binary: {e}")
    except json.JSONDecodeError as e:
        raise RuntimeError(f"Failed to parse Go output: {e}. Output was: {raw!r:.200}")
    except Exception as e:
        raise RuntimeError(f"Unexpected error calling Go implementation: {e}")
