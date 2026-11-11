"""
Figure files that look exactly like the chart on screen.

The interactive charts in Integrix are Plotly figures. Their downloads are made from the same figure by Plotly's
own image engine (Kaleido, which draws the chart with the same plotly.js the page uses), in the light screen style
with black text, at the chosen size and resolution:

* PNG, SVG, PDF and JPEG come straight from Kaleido (SVG and PDF are vector files with real text).
* TIFF (flat RGB, LZW, DPI tag set) is made from the full-resolution PNG, so it looks the same.
* EPS is converted from the SVG (svglib + ReportLab): a small vector file with the same layout and colours; its
  text uses the standard PostScript font Helvetica.
* Raster files are capped at 60 million pixels (server memory), as in every Omnix platform.

Kaleido needs a headless Chrome. On Render it is downloaded once at start-up (render_start.sh); if it is not
available, `available()` is False and the app offers the publication-style figure instead.
"""
from __future__ import annotations

import io
import zipfile

FORMATS = {"PNG": "png", "TIFF": "tiff", "JPEG": "jpg", "SVG": "svg", "PDF": "pdf", "EPS": "eps"}
MIME = {"png": "image/png", "tiff": "image/tiff", "jpg": "image/jpeg", "svg": "image/svg+xml",
        "pdf": "application/pdf", "eps": "application/postscript", "zip": "application/zip"}
MAX_PIXELS = 60e6
SCREEN_DPI = 96.0            # 1 inch = 96 CSS pixels, as in the browser

_STATE = {"ok": None, "error": "", "checked": 0.0}
RECHECK_S = 60.0             # Chrome may still be downloading right after a deploy: look again a minute later


def available() -> bool:
    """True when Kaleido and its Chrome can draw a figure. A success is kept for the server process; a failure
    is checked again after RECHECK_S seconds."""
    import time
    if _STATE["ok"] or (_STATE["ok"] is False and time.time() - _STATE["checked"] < RECHECK_S):
        return bool(_STATE["ok"])
    try:
        import plotly.graph_objects as go
        go.Figure(go.Scatter(x=[0, 1], y=[0, 1])).to_image(format="png", width=60, height=40, scale=1)
        _STATE["ok"] = True
    except Exception as e:                          # no kaleido, no Chrome yet, or Chrome cannot start
        _STATE["ok"], _STATE["error"] = False, f"{type(e).__name__}: {e}"[:300]
    _STATE["checked"] = time.time()
    return bool(_STATE["ok"])


def error() -> str:
    return _STATE["error"]


def screen_figure(fig, theme_figure=None, blacken=None, transparent: bool = False):
    """The figure as the light-mode screen shows it: same theme, pure-black text, white (or clear) background."""
    import plotly.graph_objects as go
    out = theme_figure(fig, "light") if theme_figure else go.Figure(fig)
    if blacken:
        blacken(out)
    bg = "rgba(0,0,0,0)" if transparent else "white"
    out.update_layout(paper_bgcolor=bg, plot_bgcolor=bg)
    return vector_traces(out)


def vector_traces(fig):
    """WebGL layers (Scattergl, used on screen for charts with many dots) drawn as ordinary SVG layers.

    The server's headless Chrome has no graphics card, and at print resolution (several thousand pixels) WebGL
    draws nothing: the file would show the axes and legend but no dots. The SVG version looks the same, is drawn
    at any resolution, and stays vector in SVG, PDF and EPS."""
    import plotly.graph_objects as go
    if not any(t.type == "scattergl" for t in fig.data):
        return fig
    data = []
    for t in fig.data:
        if t.type == "scattergl":
            j = t.to_plotly_json()
            j.pop("type", None)
            data.append(go.Scatter(j))
        else:
            data.append(t)
    return go.Figure(data=data, layout=fig.layout)


def effective_dpi(w_in: float, h_in: float, dpi: int) -> int:
    if w_in * h_in * dpi * dpi <= MAX_PIXELS:
        return int(dpi)
    return max(72, int((MAX_PIXELS / (w_in * h_in)) ** 0.5))


def figure_bytes(fig, ext: str, dpi: int, w_in: float, h_in: float, transparent: bool = False) -> bytes:
    """One file of `fig` (already styled with screen_figure) at w_in x h_in inches."""
    w_px, h_px = max(80, int(round(w_in * SCREEN_DPI))), max(60, int(round(h_in * SCREEN_DPI)))
    if ext in ("svg", "pdf"):
        return fig.to_image(format=ext, width=w_px, height=h_px, scale=1)
    if ext == "eps":
        svg = fig.to_image(format="svg", width=w_px, height=h_px, scale=1)
        return svg_to_eps(svg, w_in, h_in)
    eff = effective_dpi(w_in, h_in, dpi)
    png = fig.to_image(format="png", width=w_px, height=h_px, scale=eff / SCREEN_DPI)
    from PIL import Image
    im = Image.open(io.BytesIO(png))
    im.load()
    out = io.BytesIO()
    if ext == "png":
        im.save(out, format="PNG", dpi=(eff, eff))
        return out.getvalue()
    if im.mode in ("RGBA", "LA", "P"):                # flatten on white: TIFF / JPEG / EPS have no transparency
        im = im.convert("RGBA")
        bg = Image.new("RGB", im.size, "white")
        bg.paste(im, mask=im.split()[-1])
        im = bg
    im = im.convert("RGB")
    if ext == "tiff":
        im.save(out, format="TIFF", compression="tiff_lzw", dpi=(eff, eff))
    elif ext == "jpg":
        im.save(out, format="JPEG", quality=95, dpi=(eff, eff), subsampling=0)
    else:
        raise ValueError(ext)
    return out.getvalue()


def svg_to_eps(svg: bytes, w_in: float, h_in: float) -> bytes:
    """Vector EPS of a Plotly SVG, w_in x h_in inches."""
    import tempfile
    from reportlab.graphics import renderPS
    from svglib.svglib import svg2rlg
    import re

    def _render(src: bytes) -> bytes:
        d = svg2rlg(io.BytesIO(src))
        sx, sy = w_in * 72.0 / d.width, h_in * 72.0 / d.height       # Plotly's SVG is in CSS px; EPS in points
        d.scale(sx, sy)
        d.width, d.height = w_in * 72.0, h_in * 72.0
        with tempfile.NamedTemporaryFile(suffix=".eps") as tmp:
            renderPS.drawToFile(d, tmp.name)
            return open(tmp.name, "rb").read()

    svg = flatten_gradients(svg)
    try:
        return _render(svg)
    except Exception:
        # ReportLab's PostScript writer cannot clip some shapes (e.g. the cell image of a heatmap). Plotly only
        # clips to the plot area, where those shapes already fit, so draw them without the clip.
        return _render(re.sub(rb'\s?clip-path="[^"]*"', b"", svg))


def flatten_gradients(svg: bytes, steps: int = 64) -> bytes:
    """Replace each gradient-filled rectangle (Plotly colour bars) by `steps` solid bands of the same colours.
    PostScript output cannot shade gradients; at print size 64 bands look the same and stay vector."""
    import re
    text = svg.decode("utf-8")
    grads = {}
    for m in re.finditer(r'<linearGradient([^>]*)>(.*?)</linearGradient>', text, re.S):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        stops = []
        for st_ in re.finditer(r'<stop([^>]*)/?>', m.group(2)):
            a = dict(re.findall(r'([\w-]+)="([^"]*)"', st_.group(1)))
            off = float(a.get("offset", "0").rstrip("%")) / (100.0 if "%" in a.get("offset", "") else 1.0)
            rgb = [float(v) for v in re.findall(r"[\d.]+", a.get("stop-color", "rgb(0,0,0)"))[:3]]
            stops.append((off, rgb))
        if attrs.get("id") and stops:
            grads[attrs["id"]] = (attrs, sorted(stops))

    def colour(stops, t):
        for (o0, c0), (o1, c1) in zip(stops, stops[1:]):
            if o0 <= t <= o1:
                k = 0 if o1 == o0 else (t - o0) / (o1 - o0)
                return [c0[i] + k * (c1[i] - c0[i]) for i in range(3)]
        return stops[0][1] if t < stops[0][0] else stops[-1][1]

    def repl(m):
        tag = m.group(0)
        gid = m.group(1)
        if gid not in grads:
            return tag
        attrs, stops = grads[gid]
        g = {k: float(v) for k, v in re.findall(r'\s(x|y|width|height)="([-\d.]+)"', tag)}
        if not {"x", "y", "width", "height"} <= set(g):
            return tag
        x1, y1, x2, y2 = (float(attrs.get(k, d)) for k, d in (("x1", 0), ("y1", 0), ("x2", 1), ("y2", 0)))
        vertical = abs(x2 - x1) < 1e-9
        out = []
        for i in range(steps):
            t0, t1 = i / steps, (i + 1) / steps
            tm = (t0 + t1) / 2
            # position along the gradient vector -> offset
            off = (y1 + (y2 - y1) * tm) if vertical else (x1 + (x2 - x1) * tm)
            off = (off - (y1 if vertical else x1)) / ((y2 - y1) if vertical else (x2 - x1))
            r, gg, b = colour(stops, tm)
            if vertical:   # bbox units: y runs top (0) to bottom (1); the gradient goes from y1 to y2
                ya = y1 + (y2 - y1) * t0
                yb = y1 + (y2 - y1) * t1
                top, bot = min(ya, yb), max(ya, yb)
                xx, yy, ww, hh = g["x"], g["y"] + top * g["height"], g["width"], (bot - top) * g["height"] + 0.3
            else:
                xa_ = x1 + (x2 - x1) * t0
                xb_ = x1 + (x2 - x1) * t1
                lft, rgt = min(xa_, xb_), max(xa_, xb_)
                xx, yy, ww, hh = g["x"] + lft * g["width"], g["y"], (rgt - lft) * g["width"] + 0.3, g["height"]
            out.append(f'<rect x="{xx:.3f}" y="{yy:.3f}" width="{ww:.3f}" height="{hh:.3f}" '
                       f'style="fill: rgb({r:.0f}, {gg:.0f}, {b:.0f}); stroke: none;"/>')
        return "".join(out)

    text = re.sub(r"""<rect[^>]*fill:\s*url\(['"]?#([^'")]+)['"]?\)[^>]*/>""", repl, text)
    return text.encode("utf-8")


def all_formats_zip(fig, name: str, dpi: int, w_in: float, h_in: float, transparent: bool = False) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for ext in FORMATS.values():
            z.writestr(f"{name}.{ext}", figure_bytes(fig, ext, dpi, w_in, h_in, transparent))
    return buf.getvalue()
