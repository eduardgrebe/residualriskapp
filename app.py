# Residual HIV Transfusion Transmission Risk Estimator
# Copyright (C) 2025-2026 Vitalant and Eduard Grebe Consulting
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

APP_VERSION = "1.1.0b6.dev0"

_STATIC_DIR = Path(__file__).parent / "static"


@st.cache_data
def _img_data_uri(path: str) -> str:
    """Base64 ``data:`` URI for a PNG, so the logo images can be embedded directly
    in the theme-switching component's HTML below. Cached: read/encoded once."""
    return "data:image/png;base64," + base64.b64encode(Path(path).read_bytes()).decode()


# Self-contained HTML for the sidebar footer — the theme-aware VRI logo plus the
# "A project of…" attribution and the app/library version line — all rendered in a
# single `st.iframe(height="content")`.
#
# Why an iframe: the logo must follow Streamlit's *active* theme (including a
# Settings-menu choice), which Streamlit applies via emotion with no CSS hook and
# which `prefers-color-scheme` cannot see (OS only). Detecting it needs JavaScript,
# and st.markdown/st.html strip <script>; an iframe runs JS. Its sandbox is
# `allow-same-origin allow-scripts`, so the script reads the parent theme (from the
# app's `color-scheme`, falling back to background luminance), swaps the logo, and
# paints the frame body with the parent sidebar's background + inherits its
# font/colour — so the frame blends in (also hiding the opaque canvas Chromium
# renders behind a transparent iframe) and the text reads like native sidebar text.
# It re-checks on a MutationObserver plus a short poll, switching live on ANY theme
# change (menu or OS), no reload.
#
# Why the attribution + version live INSIDE the frame: st.iframe(height="content")
# floors the frame height at ~150px, so a short logo (narrow sidebar) left a large
# gap above a *separate* attribution rendered below the frame. Bundling the text in
# makes the content taller than the floor and lets CSS control the logo→text spacing
# tightly at any width.
#
# Links open a new tab via `window.parent.open` (the sandbox has no allow-popups, so
# a plain target="_blank" inside the frame is blocked); VRI-stats attribution rides
# in the `utm_source` query param. Assets inlined as data: URIs (the iframe HTML is
# not sanitised, so no static serving).
_LOGO_COMPONENT_HTML = """<!doctype html><html><head><meta charset="utf-8"><style>
  html,body{margin:0;padding:0;background:transparent}
  #wrap{display:flex;flex-direction:column;align-items:center;text-align:center;
    font-family:sans-serif;-webkit-font-smoothing:antialiased}
  a.logo{display:block;width:100%;cursor:pointer;text-decoration:none}
  img{width:100%;display:block;border:0}
  .attr{margin-top:.8rem;font-size:.875em;line-height:1.45}
  .attr a{color:var(--lnk);text-decoration:none}
  .attr a:hover{text-decoration:underline}
  hr{width:100%;border:0;border-top:1px solid currentColor;opacity:.2;margin:.8rem 0}
  .ver{font-size:.875em;opacity:.55;line-height:1.4}
</style></head><body>
<div id="wrap">
  <a class="logo" href="__URL__" target="_blank" rel="noopener"
     onclick="try{window.parent.open('__URL__','_blank','noopener');return false;}catch(e){return true;}">
    <img id="vri" src="__COLOUR__" alt="Vitalant Research Institute">
  </a>
  <div class="attr">A project of <a href="__URL__" target="_blank" rel="noopener"
     onclick="try{window.parent.open('__URL__','_blank','noopener');return false;}catch(e){return true;}">Vitalant Research Institute</a>.</div>
  <hr>
  <div class="ver">App v__APP_VER__ · Library v__LIB_VER__</div>
</div>
<script>
(function(){
  var img=document.getElementById("vri"), wrap=document.getElementById("wrap"),
      COLOUR="__COLOUR__", WHITE="__WHITE__", cur=null;
  function firstOpaqueBg(el){
    for(; el; el=el.parentElement){
      var bg=getComputedStyle(el).backgroundColor;
      if(bg && bg!=="transparent" && bg!=="rgba(0, 0, 0, 0)") return bg;
    }
    return "";
  }
  function parentEl(sel){ return window.parent.document.querySelector(sel); }
  function isDark(){
    try{
      var app=parentEl('[data-testid="stApp"]')||window.parent.document.body;
      var cs=(getComputedStyle(app).colorScheme||"");
      if(/dark/.test(cs) && !/light/.test(cs)) return true;
      if(/light/.test(cs) && !/dark/.test(cs)) return false;
      var m=firstOpaqueBg(app).match(/[0-9.]+/g);
      if(m && m.length>=3) return (0.299*(+m[0])+0.587*(+m[1])+0.114*(+m[2])) < 128;
    }catch(e){}
    return !!(window.matchMedia && matchMedia("(prefers-color-scheme: dark)").matches);
  }
  function apply(){
    var d=isDark();
    if(d!==cur){ cur=d; img.src=d?WHITE:COLOUR; }
    // Match the sidebar: blend the frame background (also covers the opaque canvas
    // Chromium paints on dark) and inherit its font + text colour so the attribution
    // and version read like native sidebar text. Firefox is already transparent.
    try{
      var sb=parentEl('[data-testid="stSidebar"]')||parentEl('[data-testid="stApp"]')||window.parent.document.body;
      var cs=getComputedStyle(sb);
      document.body.style.background=firstOpaqueBg(sb)||"transparent";
      wrap.style.color=cs.color;
      wrap.style.fontFamily=cs.fontFamily;
      wrap.style.fontSize=cs.fontSize;
      wrap.style.setProperty("--lnk", d?"#6cb4ff":"#0068c9");
    }catch(e){}
  }
  apply();
  try{
    var pd=window.parent.document, app=parentEl('[data-testid="stApp"]');
    new MutationObserver(apply).observe(pd.head, {childList:true, subtree:true, attributes:true});
    if(app) new MutationObserver(apply).observe(app, {attributes:true, attributeFilter:["style","class"]});
  }catch(e){}
  if(window.matchMedia) matchMedia("(prefers-color-scheme: dark)").addEventListener("change", apply);
  setInterval(apply, 400);
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
# Theme-aware VRI logo + "A project of…" attribution + app/library version, all
# rendered in one iframe (see _LOGO_COMPONENT_HTML above): its JS reads the app's
# active theme and swaps the logo live (tracking a Settings-menu change as well as an
# OS one). Bundling the text into the frame keeps the logo→text spacing tight at any
# sidebar width — st.iframe floors the frame height, which otherwise orphaned a
# separate caption below a short logo. `utm_source` gives VRI's stats attribution
# (the Referer header is unavailable: Streamlit force-strips it on new-tab links, and
# an iframe navigation wouldn't carry the app's URL anyway).
_VRI_URL = "https://research.vitalant.org/?utm_source=residualrisk.org"
_logo_html = (
    _LOGO_COMPONENT_HTML
    .replace("__COLOUR__", _img_data_uri(str(_STATIC_DIR / "vri_logo_web.png")))
    .replace("__WHITE__", _img_data_uri(str(_STATIC_DIR / "vri_logo_white_web.png")))
    .replace("__APP_VER__", APP_VERSION)
    .replace("__LIB_VER__", rr.__version__)
    .replace("__URL__", _VRI_URL)
)
st.sidebar.iframe(_logo_html, height="content")
