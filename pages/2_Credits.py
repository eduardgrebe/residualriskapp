# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from pathlib import Path

import streamlit as st

st.set_page_config(page_title="Credits — Residual HIV-TT Risk Estimator", layout="wide")

_STATIC_DIR = Path(__file__).parent.parent / "static"

st.markdown((Path(__file__).parent.parent / "docs" / "credits.md").read_text())

# Vitalant Research Institute logo — larger, centred at the bottom of the page.
# Streamlit elements must live here in the page module; placed inside the markdown
# content file (credits.md) they render as literal text rather than running.
st.divider()
_, _logo_col, _ = st.columns([1, 2, 1])
_logo_col.image(str(_STATIC_DIR / "vri_logo_white_ri.png"), width="stretch")
