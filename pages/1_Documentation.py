# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import re
from pathlib import Path

import streamlit as st

st.set_page_config(page_title="Documentation — Residual HIV-TT Risk Estimator")

DOCS = Path(__file__).parent.parent / "docs"
text = (DOCS / "theory.md").read_text()

# Render theory.md. Text (including $…$ / $$…$$ LaTeX) is rendered with
# st.markdown; figure references are split out and rendered with st.image,
# because st.markdown does not load relative local image paths. re.split with a
# capturing group yields [text, path, text, path, …] — odd indices are paths.
for i, part in enumerate(re.split(r"!\[[^\]]*\]\(([^)]+)\)", text)):
    if i % 2 == 1:  # captured image path, e.g. "figures/fig1_dose_response.png"
        img = DOCS / part
        if img.exists():
            st.image(str(img), width="stretch")
        else:
            st.warning(f"Figure not found: {part}")
    elif part.strip():
        st.markdown(part)
