# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Residual HIV transfusion transmission risk estimation."""

from .core import (
    get_cpu_core_count,
    iwp_from_lookback_data,
    mode_rounded,
    residual_risk_rd,
    risk_days_bs,
)
from ._go import find_go_binary

__version__ = "0.1.6"

__all__ = [
    "risk_days_bs",
    "iwp_from_lookback_data",
    "residual_risk_rd",
    "get_cpu_core_count",
    "mode_rounded",
    "find_go_binary",
]
