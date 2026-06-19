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

st.set_page_config(page_title="Documentation — Residual HIV-TT Risk Estimator", layout="wide")

DOCS = Path(__file__).parent.parent / "docs"
text = (DOCS / "theory.md").read_text()

# Pattern for a markdown image reference: ![alt](path).
_IMG = re.compile(r"!\[[^\]]*\]\(([^)]+)\)")
# Pattern for a top-level section heading ("## 3. ...") at line start. Two
# hashes only — "### " subsections stay within their parent section.
_SECTION = re.compile(r"(?m)^## (.+)$")


def render_markdown_with_figures(md: str) -> None:
    """Render markdown text (including ``$…$`` / ``$$…$$`` KaTeX), splitting out
    figure references to ``st.image`` because ``st.markdown`` does not load
    relative local image paths. ``re.split`` with a capturing group yields
    ``[text, path, text, path, …]`` — odd indices are captured paths."""
    for i, part in enumerate(_IMG.split(md)):
        if i % 2 == 1:  # captured image path, e.g. "figures/fig1_dose_response.png"
            img = DOCS / part
            if img.exists():
                st.image(str(img), width="stretch")
            else:
                st.warning(f"Figure not found: {part}")
        elif part.strip():
            st.markdown(part)


# Peel a short trailing footer (the note after the document's final horizontal
# rule) so it renders below the accordion instead of being buried in the
# collapsed References expander. Guarded so it only fires on a genuine short
# trailing note (no section heading, modest length).
body_text, sep, footer = text.rpartition("\n---\n")
if not sep or "## " in footer or len(footer) > 800:
    body_text, footer = text, ""

# Split the body into a preamble (title + intro, before the first "## " section)
# and the numbered sections. The preamble stays visible; each section is rendered
# inside a collapsible expander titled by its heading, so the page opens as a
# scannable table of contents. theory.md itself remains a single document (it
# still renders as one page on GitHub) — the accordion is purely a presentation
# layer here.
matches = list(_SECTION.finditer(body_text))
preamble = body_text[: matches[0].start()] if matches else body_text
render_markdown_with_figures(preamble)

for idx, m in enumerate(matches):
    title = m.group(1).strip()
    body_start = m.end()
    body_end = matches[idx + 1].start() if idx + 1 < len(matches) else len(body_text)
    body = body_text[body_start:body_end]
    with st.expander(title, expanded=(idx == 0)):  # §1 open by default
        render_markdown_with_figures(body)

if footer.strip():
    st.markdown("---")
    render_markdown_with_figures(footer)
