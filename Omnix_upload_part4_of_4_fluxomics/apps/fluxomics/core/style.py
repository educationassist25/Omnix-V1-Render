"""Shared visual style: colors, matplotlib defaults and figure export."""
from __future__ import annotations

import io

import matplotlib
import matplotlib.colors as mcolors

# categorical slots (validated colorblind-safe order) — groups take them in this fixed order
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID, MUTED = "#000000", "#000000", "#e7e6ec", "#b9b8c0"   # all text solid black (v2.4.2)
M0_GRAY = "#d9d8d3"
# diverging (Z-scores, differences): muted blue - neutral - muted red, centered at 0
DIVERGING = mcolors.LinearSegmentedColormap.from_list(
    "flux_div", ["#2f5d9c", "#7fa6d3", "#f3f2ef", "#dc8b70", "#a8402c"])
# sequential (incorporation %): light lavender -> deep indigo
SEQUENTIAL = mcolors.LinearSegmentedColormap.from_list("flux_seq", ["#f1effa", "#b9b0e6", "#6f5fc7", "#2c2068"])
EXPORT_FORMATS = ["PNG", "TIFF", "PDF", "SVG", "JPEG"]
MIME = {"PNG": "image/png", "TIFF": "image/tiff", "PDF": "application/pdf", "SVG": "image/svg+xml",
        "JPEG": "image/jpeg"}


def group_colors(groups) -> dict:
    return {g: CATEGORICAL[i % len(CATEGORICAL)] for i, g in enumerate(groups)}


def iso_colors(n: int) -> list:
    """M+0 gray, M+1..M+n along the sequential ramp."""
    if n <= 0:
        return [M0_GRAY]
    return [M0_GRAY] + [mcolors.to_hex(SEQUENTIAL(0.25 + 0.75 * k / max(1, n - 1))) for k in range(n)]


def apply_mpl_style():
    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5, "axes.edgecolor": "#9b9aa3", "axes.linewidth": 0.6,
        "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "axes.titleweight": "bold",
        "text.color": INK, "axes.titlecolor": INK, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
        "legend.labelcolor": INK,
        "axes.titlesize": 10, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
        "ps.fonttype": 42, "svg.fonttype": "none", "figure.dpi": 110, "savefig.facecolor": "white",
    })


def apply_plotly_style():
    """Solid black text in every Plotly figure (the app's figures use the plotly_white template)."""
    try:
        import plotly.io as pio
        for name in ("plotly_white", "plotly"):
            t = pio.templates[name]
            t.layout.font.color = INK
            t.layout.title.font.color = INK
            for ax in ("xaxis", "yaxis"):
                a = t.layout[ax]
                a.tickfont.color = INK
                a.title.font.color = INK
            t.layout.legend.font.color = INK
            t.layout.coloraxis.colorbar.tickfont.color = INK
            t.layout.coloraxis.colorbar.title.font.color = INK
            t.layout.annotationdefaults.font.color = INK
    except Exception:
        pass


def fig_bytes(fig, fmt: str = "PNG", dpi: int = 300) -> bytes:
    """Export a matplotlib figure. Vector formats (PDF/SVG) keep text editable."""
    fmt = fmt.upper()
    buf = io.BytesIO()
    kw = dict(format=fmt.lower() if fmt != "JPEG" else "jpeg", dpi=dpi, bbox_inches="tight", facecolor="white")
    if fmt == "JPEG":
        kw["pil_kwargs"] = {"quality": 95}
    if fmt == "TIFF":
        # journals and Word / Prism / Illustrator expect a flat RGB TIFF (no alpha channel) with the DPI tag set
        from PIL import Image
        kw["format"] = "png"
        fig.savefig(buf, **kw)
        im = Image.open(io.BytesIO(buf.getvalue()))
        if im.mode in ("RGBA", "LA", "P"):
            bg = Image.new("RGB", im.size, "white")
            im = im.convert("RGBA")
            bg.paste(im, mask=im.split()[-1])
            im = bg
        out = io.BytesIO()
        im.convert("RGB").save(out, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
        return out.getvalue()
    fig.savefig(buf, **kw)
    return buf.getvalue()


apply_mpl_style()
apply_plotly_style()
