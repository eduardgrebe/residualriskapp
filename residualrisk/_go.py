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
import queue
import subprocess
import threading
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats


def find_go_binary():
    """
    Find the riskdays_go Go binary.

    Searches in order:
    1. $RESIDUALRISK_GO_BINARY env var (explicit override)
    2. go/bin/riskdays_go (repo root, relative to this package)
    3. /usr/local/bin/riskdays_go
    4. riskdays_go in PATH

    Returns:
        Path to binary or None if not found
    """
    import os

    env_override = os.environ.get("RESIDUALRISK_GO_BINARY")
    if env_override and Path(env_override).exists():
        return env_override

    # Package lives at <repo>/residualrisk/_go.py; binary at <repo>/go/bin/
    repo_root = Path(__file__).parent.parent
    relative_binary = repo_root / "go" / "bin" / "riskdays_go"
    if relative_binary.exists():
        return str(relative_binary)

    # Try system installation
    system_binary = Path("/usr/local/bin/riskdays_go")
    if system_binary.exists():
        return str(system_binary)

    # Try PATH
    try:
        result = subprocess.run(
            ["which", "riskdays_go"], capture_output=True, text=True, check=False
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass

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

    # Run Go binary
    try:
        process = subprocess.Popen(
            [binary_path],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,  # Line buffering for real-time output
        )

        # Send input and close stdin
        input_json = json.dumps(input_data)
        process.stdin.write(input_json)
        process.stdin.close()

        # Read both stdout and stderr in background threads to avoid deadlock
        # (With 1M+ simulations, both pipes can fill up and block the Go process)
        error_msg = "Unknown error"
        stderr_queue = queue.Queue()
        stdout_data = []

        def read_stderr():
            """Background thread to read stderr for progress updates"""
            try:
                for line in iter(process.stderr.readline, ""):
                    if not line:
                        break
                    stderr_queue.put(line)
            except Exception as e:
                stderr_queue.put(f"STDERR_ERROR: {e}")
            finally:
                process.stderr.close()

        def read_stdout():
            """Background thread to read stdout (potentially large JSON)"""
            try:
                stdout_data.append(process.stdout.read())
            except Exception as e:
                stdout_data.append(f"STDOUT_ERROR: {e}")
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

        # Get stdout data
        stdout = stdout_data[0] if stdout_data else ""
        if stdout.startswith("STDOUT_ERROR:"):
            raise RuntimeError(f"Failed to read stdout: {stdout}")

        # Parse output
        output = json.loads(stdout)

        # Extract results
        rd_pe = output["point_estimate"]
        rd_cri = tuple(output["credible_interval"])
        rd_range = tuple(output["range"])
        rdests = output["simulations"]

        # Return in same format as Python version
        if return_sim_df:
            # Construct sim_df matching Python implementation
            # IMPORTANT: Regenerate the same random samples using the same seed.
            # Note: Due to differences between Go's and Python's RNG implementations,
            # the regenerated parameters here may not exactly match what the Go
            # binary used internally, even with the same seed. The iwp results are
            # correct (from Go), but the associated parameter samples are approximations
            # based on Python's RNG with the same seed and distributions.
            np.random.seed(seed)

            # Generate k samples — mirrors the distribution used by Go.
            # Note: Python and Go use independent RNGs; values won't match
            # exactly even with the same seed. See docstring for details.
            if k_posterior_sample is not None:
                ks = np.random.choice(k_posterior_sample, size=n_bs, replace=True)
            elif k_gamma_shape is not None and k_gamma_scale is not None:
                ks = np.random.gamma(k_gamma_shape, k_gamma_scale, n_bs)
            elif k_invgamma_alpha is not None:
                ks = scipy_stats.invgamma.rvs(k_invgamma_alpha, scale=_beta, size=n_bs)
            elif k_lnmix_w is not None:
                from .core import sample_lnmix
                ks = sample_lnmix(n_bs, k_lnmix_w, k_lnmix_mu1, k_lnmix_sigma1,
                                   k_lnmix_mu2, k_lnmix_sigma2, seed=seed)

            # Generate doubling time samples (truncated normal)
            doubling_times = scipy_stats.truncnorm.rvs(
                0, np.inf, doubling_time, doubling_time_norm_sd, n_bs
            )

            # Generate lod50 samples (truncated normal)
            lod50s = scipy_stats.truncnorm.rvs(0, np.inf, lod50, lod50_sd, n_bs)

            # Generate volume transfused samples (uniform)
            volumes_transfused = np.random.uniform(
                volume_transfused_range[0], volume_transfused_range[1], n_bs
            )

            # Build args_list matching Python format
            args_list = [
                (
                    copies_per_virion,
                    C0,
                    doubling_times[i],
                    volumes_transfused[i],
                    ks[i],
                    pool_size,
                    lod50s[i],
                    lod95_lod50_ratio,
                    retests,
                    z,
                    (-100, 500),
                )
                for i in range(n_bs)
            ]

            # Create DataFrame
            sim_df = pd.DataFrame(
                args_list,
                columns=[
                    "copies_per_virion",
                    "C0",
                    "doubling_time",
                    "volume_transfused",
                    "k",
                    "pool_size",
                    "lod50",
                    "lod95_lod50_ratio",
                    "retests",
                    "z",
                    "limits",
                ],
            )

            # Convert lod95 from ratio to actual value
            sim_df["lod95"] = sim_df["lod50"] * sim_df["lod95_lod50_ratio"]

            # Add iwp results
            sim_df["iwp"] = rdests

            # Add random seed
            sim_df["random_seed"] = np.repeat(seed, n_bs)

            return (rd_pe, rd_cri, rd_range, rdests, sim_df)
        else:
            return (rd_pe, rd_cri, rd_range, rdests, None)

    except subprocess.SubprocessError as e:
        raise RuntimeError(f"Failed to run Go binary: {e}")
    except json.JSONDecodeError as e:
        raise RuntimeError(f"Failed to parse Go output: {e}. Output was: {stdout}")
    except Exception as e:
        raise RuntimeError(f"Unexpected error calling Go implementation: {e}")
