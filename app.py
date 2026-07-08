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
import base64
from pathlib import Path

import streamlit as st

import residualrisk as rr

APP_VERSION = "1.1.0a9"

_STATIC_DIR = Path(__file__).parent / "static"


@st.cache_data
def _img_data_uri(path: str) -> str:
    """Base64 ``data:`` URI for a PNG, so the logo images can be embedded directly
    in the theme-switching component's HTML below. Cached: read/encoded once."""
    return "data:image/png;base64," + base64.b64encode(Path(path).read_bytes()).decode()


# Self-contained HTML for the theme-aware sidebar logo, rendered via `st.iframe`.
# Why an iframe: the logo must follow Streamlit's *active* theme — including a
# choice made in its Settings menu — which Streamlit applies via emotion with no
# CSS hook and which `prefers-color-scheme` cannot see (that only reflects the OS
# theme). Detecting it needs JavaScript, and st.markdown/st.html strip <script>;
# an iframe runs JS. Streamlit's iframe sandbox is `allow-same-origin
# allow-scripts`, so this script can
# read the parent app's theme: it takes light/dark from the app's `color-scheme`
# (falling back to background luminance) and swaps the logo, re-checking on a
# MutationObserver plus a short poll, so it switches live on ANY theme change (menu
# or OS) with no reload. The link opens a new tab via `window.parent.open` (the
# sandbox has no allow-popups, so a plain target="_blank" inside the frame is
# blocked); VRI-stats attribution rides in the `utm_source` query param. Logos are
# inlined as data: URIs (the iframe HTML is not sanitised, so no static serving).
_LOGO_COMPONENT_HTML = """<!doctype html><html><head><meta charset="utf-8"><style>
  html,body{margin:0;padding:0;background:transparent}
  a{display:block;cursor:pointer;text-decoration:none}
  img{width:100%;display:block;border:0}
</style></head><body>
<a href="__URL__" target="_blank" rel="noopener"
   onclick="try{window.parent.open('__URL__','_blank','noopener');return false;}catch(e){return true;}">
  <img id="vri" src="__COLOUR__" alt="Vitalant Research Institute">
</a>
<script>
(function(){
  var img=document.getElementById("vri"), COLOUR="__COLOUR__", WHITE="__WHITE__", cur=null;
  function isDark(){
    try{
      var pd=window.parent.document, app=pd.querySelector('[data-testid="stApp"]')||pd.body;
      var cs=(getComputedStyle(app).colorScheme||"");
      if(/dark/.test(cs) && !/light/.test(cs)) return true;
      if(/light/.test(cs) && !/dark/.test(cs)) return false;
      var el=app, bg="";
      while(el){ bg=getComputedStyle(el).backgroundColor;
        if(bg && bg!=="transparent" && bg!=="rgba(0, 0, 0, 0)") break; el=el.parentElement; }
      var m=bg.match(/[0-9.]+/g);
      if(m && m.length>=3) return (0.299*(+m[0])+0.587*(+m[1])+0.114*(+m[2])) < 128;
    }catch(e){}
    return !!(window.matchMedia && matchMedia("(prefers-color-scheme: dark)").matches);
  }
  function apply(){ var d=isDark(); if(d!==cur){ cur=d; img.src=d?WHITE:COLOUR; } }
  function fit(){ try{ if(window.frameElement) window.frameElement.style.height=Math.ceil(img.getBoundingClientRect().height)+"px"; }catch(e){} }
  img.addEventListener("load", fit);
  apply(); fit();
  try{
    var pd=window.parent.document, app=pd.querySelector('[data-testid="stApp"]');
    new MutationObserver(apply).observe(pd.head, {childList:true, subtree:true, attributes:true});
    if(app) new MutationObserver(apply).observe(app, {attributes:true, attributeFilter:["style","class"]});
  }catch(e){}
  if(window.matchMedia) matchMedia("(prefers-color-scheme: dark)").addEventListener("change", apply);
  window.addEventListener("resize", fit);
  setInterval(apply, 400);
  setTimeout(fit, 200); setTimeout(fit, 600);
})();
</script>
</body></html>"""


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
# Theme-aware VRI logo (see _LOGO_COMPONENT_HTML above) — an iframe whose JS reads
# the parent app's active theme and swaps the logo live, so it tracks a Settings-menu
# theme change as well as an OS one. `utm_source` gives VRI's stats attribution (the
# Referer header is unavailable: Streamlit force-strips it on new-tab links, and an
# iframe navigation wouldn't carry the app's URL anyway).
_VRI_URL = "https://research.vitalant.org/?utm_source=residualrisk.org"
_logo_html = (
    _LOGO_COMPONENT_HTML
    .replace("__COLOUR__", _img_data_uri(str(_STATIC_DIR / "vri_logo_web.png")))
    .replace("__WHITE__", _img_data_uri(str(_STATIC_DIR / "vri_logo_white_web.png")))
    .replace("__URL__", _VRI_URL)
)
st.sidebar.iframe(_logo_html, height=120)
st.sidebar.markdown(f"A project of [Vitalant Research Institute]({_VRI_URL}).")
st.sidebar.divider()
st.sidebar.caption(
    f"App v{APP_VERSION} · Library v{rr.__version__}",
    text_alignment="center",
)
