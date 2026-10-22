"""Figure style for the paper (Journal of The Electrochemical Society).

JES prints figures at one-column width (3.25 in, 8.2 cm) wherever possible,
with lettering at least 1.6 mm high and 300 dpi, and asks for a lowercase
letter in the body of each part of a multi-part figure. Every paper figure is
therefore drawn at its final printed size (SINGLE or DOUBLE width), so the
font sizes below are the printed sizes: 8 pt text, nothing below 7 pt (7 pt
Arial capitals are 1.8 mm high). Panels carry letters instead of titles; the
caption describes them.
"""
import os
import string

import matplotlib as mpl

SINGLE_IN = 3.25  # one column (8.2 cm)
DOUBLE_IN = 6.75  # two columns (17.1 cm)
MIN_FONT_PT = 7.0
LABEL_PT = 9.0  # panel letters

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 7.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
    "legend.title_fontsize": 7,
    "figure.titlesize": 8,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "axes.linewidth": 0.6,
    "axes.titlepad": 3,
    "axes.labelpad": 2,
    "lines.linewidth": 1.1,
    "lines.markersize": 3.5,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.major.pad": 2,
    "ytick.major.pad": 2,
    "grid.linewidth": 0.4,
    "grid.alpha": 0.35,
    "axes.grid": True,
    "legend.frameon": True,
    "legend.framealpha": 0.85,
    "legend.edgecolor": "0.8",
    "legend.fancybox": False,
    "legend.borderpad": 0.3,
    "legend.labelspacing": 0.25,
    "legend.handlelength": 1.6,
    "legend.handletextpad": 0.4,
    "legend.borderaxespad": 0.4,
    "legend.columnspacing": 0.9,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def use():
    """Apply the paper style (call before creating figures)."""
    mpl.rcParams.update(RC)


def letter(i):
    return f"{string.ascii_lowercase[i]})"


def panel_label(ax, i, extra=None):
    """Bold 'a)', 'b)', ... above the top-left corner of the axes (out of the
    data), optionally with a short plain label above the right end (e.g. an RPT).
    Set as axes titles so the layout engine leaves room for them."""
    ax.set_title(letter(i), loc="left", fontsize=LABEL_PT, fontweight="bold")
    if extra:
        ax.set_title(extra, loc="right", fontsize=7.5)


def small_fonts(fig, min_pt=MIN_FONT_PT):
    """Texts in a figure printed below min_pt (should be empty)."""
    out = []
    for t in fig.findobj(mpl.text.Text):
        if t.get_visible() and t.get_text().strip() and t.get_fontsize() < min_pt - 1e-6:
            out.append((round(t.get_fontsize(), 2), t.get_text()[:40]))
    return out


def save(fig, out_dir, name, formats=("png", "pdf")):
    """Save at the drawn (printed) size: <out_dir>/<name>.png (300 dpi) and .pdf.
    Returns the paths. Fails if any text is below the minimum font size."""
    bad = small_fonts(fig)
    if bad:
        raise ValueError(f"{name}: text below {MIN_FONT_PT} pt: {bad[:5]}")
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for fmt in formats:
        p = os.path.join(out_dir, f"{name}.{fmt}")
        fig.savefig(p)
        paths.append(p)
    return paths
