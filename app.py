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

# Multipage entry point / router.
#   Run with: streamlit run app.py
# The estimator UI lives in estimator.py; the Documentation and Credits pages
# live in pages/. Using st.navigation (rather than the automatic pages/ folder
# nav) lets us set explicit page labels — in particular "Estimator" for the
# main page instead of the entry-script filename.
from pathlib import Path

import streamlit as st

import residualrisk as rr

APP_VERSION = "1.1.0a8"

_STATIC_DIR = Path(__file__).parent / "static"

# Shared page chrome and branding. Must be the first Streamlit command in the
# entry script; individual pages may additively override the title.
st.set_page_config(
    page_title="Residual HIV-TT Risk Estimator",
    page_icon=str(_STATIC_DIR / "transfusion_bag_purple.png"),
)

# Widen the centered main column from Streamlit's default (~730px) to ~1100px
# (≈50% wider) — enough for the results table to fit on one line, without the
# stretched look of layout="wide". Injected from the router so it applies to
# every page; the centered layout keeps the column horizontally centred.
st.markdown(
    "<style>.block-container, [data-testid='stMainBlockContainer']"
    " { max-width: 1100px !important; }</style>",
    unsafe_allow_html=True,
)

navigation = st.navigation([
    st.Page("estimator.py", title="Estimator", default=True),
    st.Page("pages/1_Documentation.py", title="Documentation"),
    st.Page("pages/2_Credits.py", title="Credits"),
])
navigation.run()

# Shared sidebar footer — rendered on every page, after the active page's own
# sidebar content, so it stays pinned at the bottom of the sidebar. The logo sits
# above the centred app/library version caption.
#
# Only the Estimator page has its own sidebar controls above the logo, so it gets
# a divider to separate them. The Documentation/Credits pages have no sidebar
# controls — the navigation menu's own separator already sits directly above the
# logo — so a second divider there would just double up.
if navigation.title == "Estimator":
    st.sidebar.divider()
st.sidebar.image(str(_STATIC_DIR / "vri_logo_white.png"), width="stretch")
st.sidebar.divider()
st.sidebar.caption(
    f"App v{APP_VERSION} · Library v{rr.__version__}",
    text_alignment="center",
)
