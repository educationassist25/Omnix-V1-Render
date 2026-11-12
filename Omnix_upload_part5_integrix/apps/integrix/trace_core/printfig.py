"""
Print typesetting of the Integrix charts: the same chart as on screen, laid out as a journal figure.

The screen layout is made for a browser window (wide canvas, tall rows, generous margins, 13 px text). Copied
unchanged into a file, that becomes a large canvas whose text is far too small once the figure is reduced to a
journal column, with empty space around the chart. `typeset()` keeps the chart itself (data, marks, colours,
order) and sets it for print:

* Width: a journal column. Single column 3.5 in (89 mm), 1.5 column 5.5 in (140 mm), double column 7.2 in (183 mm).
* Type: Arial (Liberation Sans on Linux servers) at real point sizes for the final printed size: body and axis
  titles 8 pt, tick labels and legend 7 pt, nothing below 6 pt. Black text, 0.75 pt axis lines.
* Height: computed from the content (rows, legend lines, label lengths), with tight margins, so the file holds
  the chart and nothing else. No title inside the image: a journal's figure legend carries it.

1 in = 96 px when the file is drawn (Plotly's unit), so 1 pt = 4/3 px.
"""
from __future__ import annotations

import math
import re

import numpy as np
import plotly.graph_objects as go

FAMILY = "Arial, 'Liberation Sans', Helvetica, sans-serif"
PX_PER_IN = 96.0
WIDTHS = {"Single column (3.5 in · 89 mm)": 3.5, "1.5 column (5.5 in · 140 mm)": 5.5,
          "Double column (7.2 in · 183 mm)": 7.2}
# default width for each chart (by the name the app gives it)
DEFAULT_WIDTH = {"layer_concordance": 7.2, "reaction_classification": 3.5, "reaction_patterns": 7.2,
                 "pathway_activity": 5.5, "mapping_coverage": 3.5, "canonical_pathways": 5.5, "concordance": 5.5,
                 "factor_variance": 5.5, "factor_scores": 3.5, "signature_roc": 3.5, "signature_permutation": 3.5,
                 "signature_panel": 3.5}
# widths a chart can be set at without its labels colliding (dense charts need the wider columns)
ALLOWED = {"reaction_patterns": (5.5, 7.2), "concordance": (5.5, 7.2), "layer_concordance": (7.2,),
           "factor_variance": (5.5, 7.2),
           "pathway_activity": (5.5, 7.2), "canonical_pathways": (3.5, 5.5, 7.2)}
MIN_PT = 6.0


def widths_for(name: str) -> dict:
    """The width choices offered for a chart, label -> inches."""
    ok = ALLOWED.get(name, (3.5, 5.5, 7.2))
    return {k: v for k, v in WIDTHS.items() if v in ok}


def px(pt: float) -> float:
    return pt * 4.0 / 3.0


_FONT_PROP = None


def _font():
    """The font used to measure text: Arial or its metric twin Liberation Sans when installed (the server's
    Chrome uses the same), otherwise DejaVu Sans, which is ~8 % wider, so the layout errs on the safe side."""
    global _FONT_PROP
    if _FONT_PROP is None:
        from matplotlib import font_manager as fm
        have = {f.name for f in fm.fontManager.ttflist}
        fam = next((n for n in ("Arial", "Liberation Sans", "Helvetica", "Arimo") if n in have), "DejaVu Sans")
        _FONT_PROP = fm.FontProperties(family=fam)
    return _FONT_PROP


@__import__("functools").lru_cache(maxsize=20000)
def _measure(s: str, pt: float) -> float:
    try:
        from matplotlib.textpath import TextPath
        tp = TextPath((0, 0), s, size=pt, prop=_font())
        return (tp.get_extents().width + 0.25 * pt) / 72.0             # + side bearings
    except Exception:
        return 0.56 * pt / 72.0 * len(s)


def text_w_in(s: str, pt: float) -> float:
    """Width (inches) of the longest line of `s` (Plotly <br> lines) at `pt` points, from the font outlines."""
    lines = re.split(r"<br\s*/?>", str(s))
    return max(_measure(re.sub(r"<[^>]+>", "", ln).strip(), round(float(pt), 2)) for ln in lines) * 1.03


def _split_word(w: str, max_in: float, pt: float) -> list:
    """A word longer than a line: break after - _ / . or, failing that, at the last character that fits."""
    out, cur = [], ""
    for ch in w:
        if cur and text_w_in(cur + ch, pt) > max_in:
            cut = max((cur.rfind(c) for c in "-_/."), default=-1)
            if cut > len(cur) // 3:
                out.append(cur[:cut + 1])
                cur = cur[cut + 1:] + ch
            else:
                out.append(cur)
                cur = ch
        else:
            cur += ch
    out.append(cur)
    return out


def wrap(s: str, max_in: float, pt: float, max_lines: int = 3) -> str:
    """Wrap `s` into lines (Plotly <br>) that each fit in max_in inches at `pt`. Long words are broken; if the
    text needs more than max_lines lines, the last line ends with an ellipsis. Nothing is ever wider than
    max_in, so neighbouring labels cannot run into each other."""
    words = []
    for w in re.sub(r"<br\s*/?>", " ", str(s)).split():
        words += _split_word(w, max_in, pt) if text_w_in(w, pt) > max_in else [w]
    lines, cur = [], ""
    for w in words:
        cand = f"{cur} {w}" if cur else w
        if cur and text_w_in(cand, pt) > max_in:
            lines.append(cur)
            cur = w
        else:
            cur = cand
    if cur:
        lines.append(cur)
    if len(lines) > max_lines:
        last = lines[max_lines - 1]
        while last and text_w_in(last + "…", pt) > max_in:
            last = last[:-1].rstrip()
        lines = lines[:max_lines - 1] + [last + "…"]
    return "<br>".join(lines)


def _n_lines(s) -> int:
    return str(s).count("<br>") + 1


def _legend_rows(names, W_in: float, pt: float, symbol_in: float = 0.45) -> int:
    """Text lines a horizontal legend needs at width W_in (Plotly packs entries left to right; an entry whose
    name has several lines makes its row that tall)."""
    lines, used, row_max = 0, 0.0, 1
    for n in names:
        w = text_w_in(n, pt) + symbol_in + 0.12
        if used and used + w > W_in:
            lines += row_max
            used, row_max = 0.0, 1
        used += w
        row_max = max(row_max, _n_lines(n))
    return lines + row_max


def _base(fig, W_in: float, pt: float):
    """Fonts, lines, background and margins shared by every chart."""
    f = go.Figure(fig)
    body, small = px(pt), px(max(MIN_PT, pt - 1))
    f.update_layout(title_text="", template="simple_white", autosize=False, width=int(round(W_in * PX_PER_IN)),
                    font=dict(family=FAMILY, size=body, color="#000000"),
                    paper_bgcolor="white", plot_bgcolor="white",
                    margin=dict(l=4, r=6, t=6, b=4, pad=0),
                    legend=dict(font=dict(size=small), title_font=dict(size=small), bgcolor="rgba(0,0,0,0)",
                                itemwidth=30, tracegroupgap=2, itemsizing="constant"))
    for sh in f.layout.shapes or ():          # the template draws shapes at 30 % opacity: reference lines vanish
        if sh.opacity is None:
            sh.opacity = 1
    axis = dict(tickfont=dict(size=small, color="#000000"), title_font=dict(size=body, color="#000000"),
                linewidth=1, linecolor="#000000", tickwidth=1, tickcolor="#000000", ticklen=3,
                title_standoff=4, automargin=True)
    f.update_xaxes(**axis)
    f.update_yaxes(**axis)
    for a in f.layout.annotations or ():
        o = a.font.size if a.font and a.font.size else 13
        a.font.size = max(px(MIN_PT), min(px(pt * 2.2), body * o / 13.0))
        a.font.family = FAMILY
    for t in f.data:
        tf = getattr(t, "textfont", None)
        if tf is not None and tf.size:
            t.textfont.size = max(px(MIN_PT), small * tf.size / 11.0)
        cb = getattr(t, "colorbar", None)
        if cb is not None and t.type in ("heatmap", "contour"):
            t.colorbar.update(thickness=10, outlinewidth=0.5, tickfont=dict(size=small),
                              title=dict(font=dict(size=small), side="right"), len=0.9)
    return f


def _arr(v):
    """A Plotly array value as numbers (handles Plotly's binary-encoded arrays {"dtype", "bdata"})."""
    if isinstance(v, dict) and "bdata" in v:
        import base64
        a = np.frombuffer(base64.b64decode(v["bdata"]), dtype=np.dtype(v["dtype"])).astype(float)
        if v.get("shape"):
            a = a.reshape([int(x) for x in str(v["shape"]).split(",")])
        return a
    return np.asarray(v, float)


def _scale_markers(f, k: float):
    for t in f.data:
        if t.type in ("scatter", "scattergl") and t.marker is not None and t.marker.size is not None:
            s = _arr(t.marker.size)
            t.marker.size = (s * k).tolist() if s.ndim else float(s) * k
            if t.marker.line is not None and t.marker.line.width:
                t.marker.line.width = max(0.4, t.marker.line.width * k)


def _legend_below(f, names, W_in, pt, gap_in=0.32):
    """Horizontal legend under the plot; returns the inches it takes."""
    small = max(MIN_PT, pt - 1)
    rows = _legend_rows(names, W_in - 0.1, small)
    h = rows * (small * 1.55 / 72.0) + 0.04
    f.update_layout(legend=dict(orientation="h", x=0.01, xanchor="left", xref="container", y=0, yref="container",
                                yanchor="bottom", entrywidth=0))
    return h, gap_in


def typeset(fig, name: str, W_in: float | None = None, pt: float = 8.0, H_in: float | None = None):
    """Return (figure, width_in, height_in) for print. H_in overrides the computed height."""
    W_in = float(W_in or DEFAULT_WIDTH.get(name, 5.5))
    pt = max(MIN_PT + 1, float(pt))
    small = max(MIN_PT, pt - 1)
    f = _base(fig, W_in, pt)
    line_in = pt * 1.25 / 72.0
    ok = ALLOWED.get(name)
    if ok and W_in < min(ok):                       # too narrow for this chart's labels
        W_in = min(ok)
        f.update_layout(width=int(round(W_in * PX_PER_IN)))
    for t in f.data:                                # legend entries never wider than the figure
        if t.name and t.showlegend is not False and text_w_in(t.name, small) > W_in - 0.7:
            t.name = wrap(t.name, W_in - 0.7, small, 3)
    fn = _LAYOUTS.get(name, _generic)
    h = fn(f, W_in, pt, small, line_in)
    H = float(H_in or h)
    _wrap_axis_titles(f, W_in, H, pt)
    f.update_layout(height=int(round(H * PX_PER_IN)))
    return f, W_in, H


def _wrap_axis_titles(f, W, H, pt):
    """Axis titles longer than their axis are set on two lines instead of being cut at the figure edge."""
    for ax_name in [k for k in f.layout if re.match(r"[xy]axis\d*$", k)]:
        ax = f.layout[ax_name]
        t = ax.title.text if ax.title else None
        if not t or "<br>" in t:
            continue
        room = (0.55 * W) if ax_name.startswith("x") else 0.5 * H      # about the plot's length
        if text_w_in(t, pt) > room:
            ax.title.text = wrap(t, room, pt, 2)


# ------------------------------------------------------------------------------------------------ per chart
def _generic(f, W, pt, small, line_in):
    names = [t.name for t in f.data if t.showlegend is not False and t.name]
    h = W * 0.72
    if f.layout.showlegend is not False and len(names) > 1:
        lh, gap = _legend_below(f, names, W, pt)
        f.update_layout(margin_b=int((lh + gap) * PX_PER_IN))
        h += lh
    return h


def _xy_scatter(f, W, pt, small, line_in, aspect=0.78):
    """Scatter maps (reaction patterns, concordance): square-ish plot, legend underneath."""
    _scale_markers(f, min(1.0, max(0.55, W / 11.5)))
    for a in f.layout.annotations or ():
        if a.xref != "paper":                       # gene / enzyme labels: at most 7.5 pt, so they stay clear
            a.font.size = min(a.font.size or px(small), px(min(small, 7.5)))
        if a.showarrow:
            a.arrowwidth = 0.75
            a.standoff = 2
            a.bgcolor = "rgba(255,255,255,0.85)"
    names = [t.name for t in f.data if t.name and t.showlegend is not False]
    lh, gap = _legend_below(f, names, W, pt)
    plot_h = (W - 0.75) * aspect
    ytitle = str(f.layout.yaxis.title.text or "")
    left_in = 0.42 + _n_lines(wrap(ytitle, 0.5 * (plot_h + 0.45 + lh + gap), pt, 2)) * pt * 1.25 / 72.0
    f.update_layout(margin=dict(l=int(left_in * PX_PER_IN), r=8, t=6, b=int((lh + gap + 0.45) * PX_PER_IN)))
    f.update_xaxes(title_standoff=2)
    f.update_yaxes(automargin=False)
    _relabel(f, W - left_in - 8 / PX_PER_IN, plot_h, min(small, 7.5))
    return plot_h + 0.45 + lh + gap


def _xy_points(f):
    xs, ys = [], []
    for t in f.data:
        if t.type in ("scatter", "scattergl") and t.x is not None and t.y is not None and t.mode != "lines":
            x, y = _arr(t.x), _arr(t.y)
            if x.ndim and len(x) == len(y):
                xs.append(x)
                ys.append(y)
    if not xs:
        return np.zeros((0, 2))
    P = np.c_[np.concatenate(xs), np.concatenate(ys)]
    return P[np.isfinite(P).all(1)]


def _relabel(f, plot_w: float, plot_h: float, lab_pt: float):
    """Place the gene / enzyme labels again for the print canvas (plot_w x plot_h inches): no label touches
    another, a quadrant caption or the frame; a label that cannot be placed cleanly is left out."""
    from .layout import place_labels
    xr, yr = list(f.layout.xaxis.range or []), list(f.layout.yaxis.range or [])
    if len(xr) != 2 or len(yr) != 2:
        return
    tx = lambda v: (np.asarray(v, float) - xr[0]) / (xr[1] - xr[0]) * plot_w
    ty = lambda v: (np.asarray(v, float) - yr[0]) / (yr[1] - yr[0]) * plot_h
    data_ann = [a for a in f.layout.annotations if a.xref in (None, "x") and a.yref in (None, "y")]
    keep = [a for a in f.layout.annotations if not (a.xref in (None, "x") and a.yref in (None, "y"))]
    if not data_ann:
        return
    anchors = np.array([[a.x, a.y] for a in data_ann], float)
    texts = [a.text for a in data_ann]
    A = np.c_[tx(anchors[:, 0]), ty(anchors[:, 1])]
    hl = lab_pt * 1.25 / 72.0 + 0.03
    sizes = np.array([[text_w_in(t_, lab_pt) + 0.05, hl] for t_ in texts])
    P = _xy_points(f)
    O = np.c_[tx(P[:, 0]), ty(P[:, 1])] if len(P) else np.zeros((0, 2))
    R = np.full(len(O), 0.03)
    # quadrant captions and other paper-anchored text: their boxes are occupied
    occ = []
    for a in keep:
        if a.xref == "paper" and a.yref == "paper" and a.text:
            w_, h_ = text_w_in(a.text, (a.font.size or 10) * 0.75) + 0.06, _n_lines(a.text) * hl + 0.04
            x0 = a.x * plot_w - (w_ if a.xanchor == "right" else (w_ / 2 if a.xanchor == "center" else 0))
            y0 = a.y * plot_h - (h_ if a.yanchor == "top" else (h_ / 2 if a.yanchor == "middle" else 0))
            gx, gy = np.meshgrid(np.arange(x0, x0 + w_ + 0.03, 0.03), np.arange(y0, y0 + h_ + 0.03, 0.03))
            occ.append(np.c_[gx.ravel(), gy.ravel()])
    if occ:
        O = np.r_[O, np.concatenate(occ)]
        R = np.r_[R, np.full(len(O) - len(R), 0.02)]
    C, placed, lead = place_labels(A, sizes, O, R, np.zeros((0, 2)), plot_w, plot_h, hide_unplaced=True,
                                   margin=0.02, anchor_r=np.full(len(A), 0.04))
    f.layout.annotations = keep
    for i, txt in enumerate(texts):
        if not placed[i]:
            continue
        cx = xr[0] + C[i, 0] / plot_w * (xr[1] - xr[0])
        cy = yr[0] + C[i, 1] / plot_h * (yr[1] - yr[0])
        common = dict(text=txt, font=dict(size=px(lab_pt), family=FAMILY, color="#000000"),
                      bgcolor="rgba(255,255,255,0.85)", borderpad=0)
        if lead[i]:
            f.add_annotation(x=anchors[i, 0], y=anchors[i, 1], ax=cx, ay=cy, axref="x", ayref="y", showarrow=True,
                             arrowhead=0, arrowwidth=0.75, arrowcolor="#7F8C8D", standoff=2, **common)
        else:
            f.add_annotation(x=cx, y=cy, showarrow=False, **common)


def _pattern(f, W, pt, small, line_in):
    return _xy_scatter(f, W, pt, small, line_in, aspect=0.70)


def _concordance(f, W, pt, small, line_in):
    return _xy_scatter(f, W, pt, small, line_in, aspect=0.85)


def _canonical(f, W, pt, small, line_in):
    """Horizontal bars, one row per pathway: labels wrapped to 40 % of the width, rows sized to their text."""
    bar = f.data[0]
    lab_in = 0.40 * W
    labs = [wrap(re.sub("<br>", " ", y), lab_in, small, 3) for y in bar.y]
    f.data[0].y = labs
    for t in f.data[1:]:
        if t.y is not None and len(t.y) == len(labs):
            t.y = labs
            t.marker.size = 4
            t.line.width = 1
    pitch = small * 1.25 / 72.0
    # categories are evenly spaced, so every row gets the height of the tallest label plus a gap
    rows_h = len(labs) * (max(_n_lines(l) for l in labs) * pitch + 0.9 * pitch)
    f.update_layout(bargap=0.28)
    f.update_yaxes(range=[-0.6, len(labs) - 0.4], autorange=False)
    f.update_xaxes(title_standoff=2)
    if "xaxis2" in f.layout:
        f.layout.xaxis2.update(tickfont=dict(size=px(small), color="#E67E22"),
                               title=dict(font=dict(size=px(pt), color="#E67E22"), standoff=2), nticks=5)
    return rows_h + 2 * (pt * 1.25 / 72 + small * 1.25 / 72 + 0.08) + 0.12


def _panel(f, W, pt, small, line_in):
    bar = f.data[0]
    labs = [wrap(y, 0.45 * W, small, 2) for y in bar.y]
    f.data[0].y = labs
    pitch = small * 1.25 / 72.0
    rows_h = len(labs) * (max(_n_lines(l) for l in labs) * pitch + 0.8 * pitch)
    f.update_layout(bargap=0.3)
    return rows_h + (pt * 1.25 + small * 1.25) / 72 + 0.2


def _dotplot(f, W, pt, small, line_in):
    """Pathway x layer dot grid (1:1 cells). Cell size follows the label line height; width holds labels, the
    grid, the class tiles and the colour key; height = rows + rotated column labels + key.

    The key is drawn here, across the full figure width (a Plotly legend would wrap to the narrow grid), and its
    size circles are drawn at exactly the diameters the grid uses for q = 0.05, 0.001 and <= 1e-6."""
    ya = f.layout.yaxis
    xa = f.layout.xaxis
    x0, x1 = xa.range
    units = x1 - x0
    lab_in = min(0.42 * W, 2.3)
    labs = [wrap(re.sub("<br>", " ", t), lab_in, small, 2) for t in ya.ticktext]
    lab_w = max(text_w_in(l, small) for l in sum((s_.split("<br>") for s_ in labs), [])) + 0.08
    cell = min(0.30, max(0.17, (W - lab_w - 0.25) / units))
    two = any("<br>" in l for l in labs)
    cell = max(cell, (2 if two else 1) * small * 1.2 / 72.0 + 0.04)
    grid_w = units * cell
    # angled headers at 45 degrees only when neighbours stay a full line apart (pitch x sin 45 >= line height)
    ang = 45 if cell * math.sin(math.radians(45)) >= small * 1.3 / 72.0 else 90
    col_top = max(text_w_in(t, small) for t in xa.ticktext) * math.sin(math.radians(ang)) + 0.12
    plot_h = (ya.range[1] - ya.range[0]) * cell

    # ---- key entries: class tiles present in the chart, then the three size circles
    classes = [t for t in f.data if t.name and t.marker is not None and t.marker.symbol == "square"]
    entries = [("sq", t.name, t.marker.color) for t in classes]
    for q, lab in ((0.05, "dot size: q = 0.05"), (1e-3, "q = 0.001"), (1e-6, "q \u2264 1e-6")):
        rad = 0.13 + 0.29 * min(1.0, -math.log10(q) / 6)        # same rule as the grid
        entries.append(("ci", lab, 2 * rad * cell))
    sym_in, gap_in, sep_in = 0.16, 0.06, 0.22
    row_h = max(small * 1.35 / 72.0, 2 * 0.42 * cell) + 0.04
    rows, cur, used = [], [], 0.0
    sym = lambda e: max(sym_in, e[2]) if e[0] == "ci" else sym_in
    for e in entries:
        w = sym(e) + gap_in + text_w_in(e[1], small)
        if cur and used + w > W - 0.1:
            rows.append(cur)
            cur, used = [], 0.0
        cur.append((e, w))
        used += w + sep_in
    rows.append(cur)
    key_h = len(rows) * row_h + 0.18                           # one blank line between grid and key

    H = col_top + plot_h + key_h
    left = int(lab_w * PX_PER_IN)
    right = max(4, int((W - lab_w - grid_w) * PX_PER_IN) - 2)
    top, bottom = int(col_top * PX_PER_IN), int(key_h * PX_PER_IN)
    f.update_yaxes(ticktext=labs, tickfont=dict(size=px(small)), automargin=False)
    f.update_xaxes(tickfont=dict(size=px(small)), automargin=False, constrain="range", tickangle=-ang)
    f.update_layout(margin=dict(l=left, r=right, t=top, b=bottom), showlegend=False)
    for a in f.layout.annotations or ():
        a.font.size = px(small)

    # inches from the figure's bottom-left corner -> paper coordinates of the plot area
    pw_in = W - (left + right) / PX_PER_IN
    ph_in = H - (top + bottom) / PX_PER_IN
    X = lambda xin: (xin - left / PX_PER_IN) / pw_in
    Y = lambda yin: (yin - bottom / PX_PER_IN) / ph_in
    y = key_h - 0.18 - row_h / 2
    for row in rows:
        x = 0.05
        for (kind, text, val), w in row:
            sw = sym((kind, text, val))
            cx = x + sw / 2
            if kind == "sq":
                hs = 0.055
                f.add_shape(type="rect", xref="paper", yref="paper", x0=X(cx - hs), x1=X(cx + hs), y0=Y(y - hs),
                            y1=Y(y + hs), fillcolor=val, line=dict(color="#666", width=0.6), opacity=1)
            else:
                r = val / 2
                f.add_shape(type="circle", xref="paper", yref="paper", x0=X(cx - r), x1=X(cx + r), y0=Y(y - r),
                            y1=Y(y + r), fillcolor="#BBBBBB", line=dict(color="#555", width=0.7), opacity=1)
            f.add_annotation(xref="paper", yref="paper", x=X(x + sw + gap_in), y=Y(y), text=text,
                             showarrow=False, xanchor="left", yanchor="middle", font=dict(size=px(small)))
            x += w + sep_in
        y -= row_h
    return H


def _layer_bars(f, W, pt, small, line_in):
    """Grouped bars per pathway: full pathway names, wrapped horizontally under each group."""
    n = len(f.data[0].x)
    per = (W - 0.6) / max(n, 1)
    full = [cd[0] for cd in f.data[0].customdata]
    labs = [wrap(p, per - 0.16, small, 6) for p in full]
    for t in f.data:
        t.x = labs
        t.marker.cornerradius = 2
    f.update_xaxes(tickangle=0, tickfont=dict(size=px(small)))
    f.update_layout(legend=dict(orientation="h", x=0.5, xanchor="center", y=1.0, yanchor="bottom", yref="paper"),
                    margin=dict(l=4, r=4, t=int(0.28 * PX_PER_IN), b=4), bargap=0.25)
    lines = max(_n_lines(l) for l in labs)
    return 2.0 + lines * small * 1.25 / 72 + 0.45


def _donut(f, W, pt, small, line_in):
    names = list(f.data[0].labels)
    lh, gap = _legend_below(f, names, W, pt, gap_in=0.08)
    f.update_layout(margin=dict(l=4, r=4, t=4, b=int((lh + gap) * PX_PER_IN)))
    f.data[0].marker.line.width = 1
    for a in f.layout.annotations or ():
        a.text = re.sub(r"font-size:\s*\d+px", f"font-size:{px(small):.0f}px", a.text)
        a.font.size = px(pt * 1.5)
    return min(W, 3.0) * 0.82 + lh + gap


def _coverage(f, W, pt, small, line_in):
    lh, gap = _legend_below(f, [t.name for t in f.data], W, pt, gap_in=0.3)
    f.update_layout(margin_b=int((lh + gap) * PX_PER_IN))
    f.update_yaxes(title_text="Features")
    return 2.1 + lh + gap


def _heatmap(f, W, pt, small, line_in):
    """Cell values drawn as annotations, black on light cells and white on dark ones (readable in print)."""
    from plotly.colors import sample_colorscale, get_colorscale
    hm = f.data[0]
    z = _arr(hm.z)
    nrow = z.shape[0]
    zmin = hm.zmin if hm.zmin is not None else np.nanmin(z)
    zmax = hm.zmax if hm.zmax is not None else np.nanmax(z)
    cs = [list(c) for c in (hm.colorscale if hm.colorscale is not None else get_colorscale("Teal"))]
    cs = ([[p_, c_ if c_.startswith("rgb") else
                                     "rgb(%d, %d, %d)" % tuple(int(c_.lstrip("#")[k:k + 2], 16) for k in (0, 2, 4))]
                                    for p_, c_ in cs])
    hm.texttemplate = None
    hm.text = None
    xs = list(hm.x) if hm.x is not None else list(range(z.shape[1]))
    ys = list(hm.y) if hm.y is not None else list(range(nrow))
    for i, yv in enumerate(ys):
        for j, xv in enumerate(xs):
            v = z[i, j]
            if not np.isfinite(v):
                continue
            c = sample_colorscale(cs, [float(np.clip((v - zmin) / max(zmax - zmin, 1e-9), 0, 1))])[0]
            rgb = [float(x) for x in re.findall(r"[\d.]+", c)[:3]]
            lum = (0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]) / 255
            f.add_annotation(x=xv, y=yv, text=f"{v:.1f}", showarrow=False,
                             font=dict(size=px(small), color="#000000" if lum > 0.5 else "#FFFFFF"))
    return nrow * 0.30 + 0.9


def _roc(f, W, pt, small, line_in):
    f.update_layout(legend=dict(x=0.98, xanchor="right", y=0.04, yanchor="bottom", bgcolor="rgba(255,255,255,0.8)"))
    for t in f.data:
        if t.line is not None and t.line.width and t.line.width > 1.5:
            t.line.width = 1.5
    return W * 0.92


def _perm(f, W, pt, small, line_in):
    for a in f.layout.annotations or ():
        a.font.size = px(small)
    for s_ in f.layout.shapes or ():
        s_.opacity = 1
        if s_.line is not None and s_.line.width and s_.line.width > 1.5:
            s_.line.width = 1.5
    f.update_yaxes(title_text="Permutations")
    return W * 0.70


def _scores(f, W, pt, small, line_in):
    _scale_markers(f, 0.7)
    names = [t.name for t in f.data if t.name]
    lh, gap = _legend_below(f, names, W, pt, gap_in=0.42)
    f.update_layout(margin_b=int((lh + gap) * PX_PER_IN))
    return W * 0.82 + lh + gap


_LAYOUTS = {"reaction_patterns": _pattern, "concordance": _concordance, "canonical_pathways": _canonical,
            "signature_panel": _panel, "pathway_activity": _dotplot, "layer_concordance": _layer_bars,
            "reaction_classification": _donut, "mapping_coverage": _coverage, "factor_variance": _heatmap,
            "signature_roc": _roc, "signature_permutation": _perm, "factor_scores": _scores}
