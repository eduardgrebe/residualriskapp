"""Quick check: compare IWP point estimates from human posterior vs InvGamma.

Run from the repo root:
    uv run python scripts/check_invgamma_iwp.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

from residualrisk import core as rr

# ── Load human posterior ─────────────────────────────────────────────────────
static = Path("static")
k_human = pd.read_parquet(static / "k_param_human.parquet").iloc[:, 0].values
k_mode_human = rr.mode_kde(k_human, n_grid=5_000, cap=50_000)
print(f"Human posterior: n={len(k_human)}, mode_kde={k_mode_human:.8f}")

# ── Load artificial posterior ─────────────────────────────────────────────────────
static = Path("static")
k_expdecay = pd.read_parquet(static / "k_param_expdecay.parquet").iloc[:, 0].values
k_mode_exp = rr.mode_kde(k_expdecay, n_grid=5_000, cap=100_000)
print(f"Artificial posterior: n={len(k_expdecay)}, mode_kde={k_mode_exp:.8f}")


# ── Common parameters (match app defaults) ───────────────────────────────────
params = dict(
    doubling_time=20.5 / 24,
    doubling_time_norm_sd=1.33 / 24,
    lod50=2.73,
    lod50_sd=0.193,
    lod95_lod50_ratio=12.33 / 2.73,
    volume_transfused=20,
    volume_transfused_range=(15, 30),
    pool_size=16,
    retests=1,
    n_bs=1000000,
    seed=42,
    threads=15,
    point_estimate="mode",
    use_go=True,
)

# ── Human posterior ──────────────────────────────────────────────────────────
r_human = rr.risk_days_bs(
    k=k_mode_human,
    k_posterior_sample=k_human,
    **params,
)
iwp_human, cri_human, rng_human, rdests_human, _ = r_human
print(f"\nHuman posterior:")
print(f"  IWP mode:   {iwp_human:.4f}")
print(f"  IWP 95% CrI: ({cri_human[0]:.4f}, {cri_human[1]:.4f})")
print(f"  IWP range:   ({rng_human[0]:.4f}, {rng_human[1]:.4f})")

# ── Inverse Gamma (a=2, mode = human posterior mode) ─────────────────────────
r_ig = rr.risk_days_bs(
    k=0.000673,
    k_invgamma_a=2.0,
    k_invgamma_mode=0.000673,
    **params,
)
iwp_ig, cri_ig, rng_ig, rdests_ig, _ = r_ig
print(f"\nInvGamma(a=2, mode=0.000673):")
print(f"  IWP mode:   {iwp_ig:.4f}")
print(f"  IWP 95% CrI: ({cri_ig[0]:.4f}, {cri_ig[1]:.4f})")
print(f"  IWP range:   ({rng_ig[0]:.4f}, {rng_ig[1]:.4f})")

# ── Comparison ───────────────────────────────────────────────────────────────
print(
    f"\nDifference: {abs(iwp_human - iwp_ig):.4f} "
    f"({abs(iwp_human - iwp_ig) / iwp_human * 100:.1f}%)"
)
