# -*- coding: utf-8 -*-
"""Inserts the 4 new SOC-normalised plots (+ commentary) into
FINDINGS_2026-09-20.pdf as brand-new pages, right at the end of Finding 1
(before the "Finding 2" heading). Every existing page is left byte-for-byte
untouched -- this only splices in new pages, it does not re-render the
document from markdown.
"""
import os

import fitz  # pymupdf
from reportlab.lib.pagesizes import LETTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from PIL import Image as PILImage

HERE = os.path.dirname(os.path.abspath(__file__))
ORIG_PDF = os.path.join(HERE, "FINDINGS_2026-09-20.pdf")
INSERT_PDF = os.path.join(HERE, "_soc_insert.pdf")
FINAL_PDF = os.path.join(HERE, "FINDINGS_2026-09-20.pdf")

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="Body2", parent=styles["Normal"], fontSize=10.3,
                           leading=15, spaceAfter=8, alignment=TA_JUSTIFY,
                           fontName="Helvetica"))
styles.add(ParagraphStyle(name="Bullet2", parent=styles["Body2"], leftIndent=14,
                           bulletIndent=2, spaceAfter=10))
styles.add(ParagraphStyle(name="Cap2", parent=styles["Body2"], fontSize=9,
                           leading=12, alignment=TA_JUSTIFY, spaceBefore=4,
                           spaceAfter=14, textColor="#222222"))


def fig(fname, max_width=6.3 * inch, max_height=4.3 * inch):
    path = os.path.join(HERE, fname)
    with PILImage.open(path) as im:
        w, h = im.size
    ar = w / h
    width, height = max_width, max_width / ar
    if height > max_height:
        height, width = max_height, max_height * ar
    return Image(path, width=width, height=height)


story = [
    Paragraph(
        "<b>SOC (depth-of-discharge)-normalised versions, added 2026-09-23:</b> "
        "the plots above compare model and real curves at matched <i>EFC</i>, "
        "which requires picking a single “nearest” real RPT for each "
        "model cycle — inherently approximate, since the model's own RPT "
        "cadence drifts from the real RPT schedule post-knee (a given EFC on "
        "the model side rarely lands on the exact real RPT's EFC, and the two "
        "curves don't necessarily reach the same discharge endpoint either). "
        "Renormalising both curves' discharge capacity onto a common [0, 1] "
        "depth-of-discharge (DoD) axis before comparing removes this "
        "EFC/endpoint-matching sensitivity entirely and isolates the thing "
        "actually being assessed: whether the model's voltage <i>shape</i> "
        "across the full discharge tracks the real curve's shape, independent "
        "of exactly how much capacity either curve reached. This makes "
        "DoD-normalised V(Q) and dV/dQ the more honest way to track fit "
        "quality across RPT1/2/4/5, since even a perfectly-shaped model curve "
        "would show a spurious gap under raw-EFC/Q matching if its "
        "capacity-fade rate is even slightly off from the real RPT's, whereas "
        "DoD-normalisation is blind to that mismatch and grades the shape "
        "alone.",
        styles["Body2"],
    ),
    Spacer(1, 6),
    fig("old_baseline_soc_4rpt.png"),
    Paragraph(
        "<font name=\"Courier\">old_baseline_soc_4rpt.png</font> — old "
        "baseline (pre LAM/LLI-split retune), DoD-normalised discharge "
        "voltage for RPT1, RPT2, RPT4, RPT5 each against its real "
        "counterpart. RPT1/RPT2 (pre-knee) already overlap closely; the "
        "same large post-knee sag visible in the un-normalised plots above "
        "persists on RPT4/RPT5 — so the mid-DoD (~0.2-0.8) gap is not an "
        "artefact of raw-Q/EFC mismatch, it is a genuine discharge-shape "
        "error.",
        styles["Cap2"],
    ),
    fig("new_baseline_ocp_soc_4rpt.png"),
    Paragraph(
        "<font name=\"Courier\">new_baseline_ocp_soc_4rpt.png</font> — "
        "new baseline with the Si OCP aging-deformation fix enabled "
        "(throughput-driven, RPT4 excluded from the fitted anchor set — "
        "see the OCP-aging-deformation work in the CELL064 overpotential "
        "investigation). RPT4/RPT5's mid-DoD gap shrinks substantially "
        "relative to the old baseline, confirming the OCP-deformation fix "
        "genuinely improves discharge-curve <i>shape</i>, not just an "
        "EFC-alignment coincidence.",
        styles["Cap2"],
    ),
    fig("old_baseline_soc_4rpt_dvdq.png"),
    fig("new_baseline_ocp_soc_4rpt_dvdq.png"),
    Paragraph(
        "<font name=\"Courier\">old_baseline_soc_4rpt_dvdq.png</font> / "
        "<font name=\"Courier\">new_baseline_ocp_soc_4rpt_dvdq.png</font> "
        "— the same comparison for DoD-normalised dV/dQ, the more "
        "curvature-sensitive view: the old baseline's RPT4/RPT5 dV/dQ "
        "deviates from the real curve's shape across most of the discharge, "
        "while the new (OCP-deformation-on) baseline's dV/dQ tracks the "
        "real curve's shape much more closely, particularly through the "
        "mid-discharge region where the raw-Q plots above showed the "
        "largest visual overpotential gap.",
        styles["Cap2"],
    ),
]

doc = SimpleDocTemplate(
    INSERT_PDF, pagesize=LETTER,
    leftMargin=0.85 * inch, rightMargin=0.85 * inch,
    topMargin=0.85 * inch, bottomMargin=0.85 * inch,
)
doc.build(story)

base = fitz.open(ORIG_PDF)
insert = fitz.open(INSERT_PDF)
n_insert_pages = len(insert)

# Locate the "Finding 2" heading page (0-indexed) so the new pages land
# right after the end of Finding 1 and before Finding 2 starts.
insert_at = None
for i, page in enumerate(base):
    if "Finding 2: post-knee LAM" in page.get_text():
        insert_at = i
        break
if insert_at is None:
    raise RuntimeError("Could not locate 'Finding 2' heading page")

base.insert_pdf(insert, start_at=insert_at)
tmp_out = os.path.join(HERE, "_findings_spliced.pdf")
base.save(tmp_out)
base.close()
insert.close()

os.replace(tmp_out, FINAL_PDF)
os.remove(INSERT_PDF)
print(f"Inserted {n_insert_pages} pages before original page index {insert_at}; "
      f"saved {FINAL_PDF}")
