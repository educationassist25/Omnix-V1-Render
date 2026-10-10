"""
core.pathmap — pathway maps from the hand-drawn layouts in core.pathway_templates (v2.5).

Every detected metabolite of the layout becomes a CARD with two heatmaps (unchanged from v2.4.2):

    ┌──────────────────────────────────────────────┐
    │ Citrate                                      │
    │      Incorporation (%)        Z-score         │
    │      ▬▬▬▬▬▬  ▬▬▬▬▬▬           ▬▬▬▬▬▬  ▬▬▬▬▬▬   │  group colour strip
    │      ● ● ●   ● ● ●            ● ● ●   ● ● ●    │  time-point dots
    │ M+0  ■ ■ ■   ■ ■ ■            ■ ■ ■   ■ ■ ■    │  rows = isotopologues
    └──────────────────────────────────────────────┘

* Incorporation: mean isotopologue fraction (%) per group × time point.
* Z-score: row-wise Z-score of each isotopologue channel across all analysed samples from the active intensity
  dataset (log2 or log2 + IQR — the matrix the Heatmap step uses), averaged per group × time point.

Display: an SVG vector viewer (pan / zoom; text scales with the map) and the same figure for export.

Layout: the template coordinates set the drawing (panel, arrows, cofactors, links to other pathways, sub-panels,
the TCA ring). The x and y spacing are expanded until no card, label, link or cofactor label overlaps; arrows,
corner points, arcs and backgrounds follow the expanded coordinates, so the drawing keeps its design. Undetected
metabolites keep their name (italic). Detected members of the matching annotation pathways that are not on the
drawing are listed in a block underneath. All geometry is in cell units (1 unit = one heatmap row).
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .pathway_templates import T as TEMPLATES, Template
from .style import SEQUENTIAL, DIVERGING, group_colors

# ---- card geometry (cell units; v2.5.1: name on top of its heatmaps, tighter spacing, larger text)
# v2.5.5: smaller heatmap cells (0.72 × 0.75 units instead of 1 × 1), same text sizes
PAD = 0.4
FONT = 1.2            # metabolite name (bold)
LH = 1.55             # name row
CHAR = 0.6
PH = 1.0
PFONT = 0.74          # "Incorporation (%)" / "Z-score"
GS = 0.4
TD = 0.7
RL = 2.3
RFONT = 0.70          # M+0 … M+n
RH = 0.75             # height of one isotopologue row
CW = 0.72             # width of one sample column (values printed in the cells: CW_VAL)
CW_VAL = 2.7
DOT = 0.26            # time-point dot radius
GGAP = 0.22
PGAP = 0.7
VFONT = 0.52
# ---- drawing
MARGIN = 1.9          # free space between boxes (arrows run through it)
ROW_GAP = 3.0         # vertical gap between rows of the drawing (room for an arrow + its label)
NDFONT = 1.15         # undetected metabolite label
STFONT = 1.05         # link to another pathway (stub)
CFONT = 0.88          # cofactor / enzyme label on an arrow
TFONT = 2.0           # panel title
BFONT = 1.25          # sub-panel title
CIRC_FONT = 2.8       # text in the TCA circle
HEAD_L, HEAD_W = 0.85, 0.4
PANEL_PAD = 1.4
TITLE_BAND = 3.6

PANEL_TITLES = {"iso": "Incorporation (%)", "z": "Z-score"}
SEQ_CMAPS = {"Omnix purple": None, "Viridis": "viridis", "Magma": "magma", "Inferno": "inferno", "Cividis": "cividis",
             "Blues": "Blues", "Reds": "Reds", "Greens": "Greens", "Purples": "Purples", "Yellow–orange–red": "YlOrRd",
             "Yellow–green–blue": "YlGnBu", "Greys": "Greys"}
DIV_CMAPS = {"Omnix blue–red": None, "Red–blue (RdBu)": "RdBu_r", "Cool–warm": "coolwarm",
             "Purple–green (PRGn)": "PRGn_r", "Pink–green (PiYG)": "PiYG_r", "Brown–teal (BrBG)": "BrBG_r",
             "Red–yellow–blue": "RdYlBu_r", "Spectral": "Spectral_r", "Seismic": "seismic"}


# =============================================================================================================
# colours
# =============================================================================================================
def get_cmap(name: str, kind: str):
    import matplotlib.pyplot as plt
    n = (SEQ_CMAPS if kind == "iso" else DIV_CMAPS).get(name)
    if n is None:
        return SEQUENTIAL if kind == "iso" else DIVERGING
    return plt.get_cmap(n)


def _hex(cmap, v, vmin, vmax):
    import matplotlib.colors as mc
    if v is None or not np.isfinite(v):
        return "#f0eff3"
    t = (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5
    return mc.to_hex(cmap(float(min(max(t, 0.0), 1.0))))


def _ink_on(hex_color: str) -> str:
    import matplotlib.colors as mc
    r, g, b = mc.to_rgb(hex_color)
    return "#000000" if 0.299 * r + 0.587 * g + 0.114 * b > 0.55 else "#ffffff"


def time_colors(times) -> dict:
    """Time point -> dot colour (light → dark green with time; same ramp as the heatmap Time annotation)."""
    import matplotlib.colors as mc
    import matplotlib.pyplot as plt
    ts = sorted(dict.fromkeys(times))
    return {h: mc.to_hex(plt.get_cmap("Greens")(0.35 + 0.6 * k / max(1, len(ts) - 1))) for k, h in enumerate(ts)}


def _w(text: str, font: float) -> float:
    return len(str(text)) * font * CHAR


# =============================================================================================================
# data
# =============================================================================================================
def mid_table(F: pd.DataFrame, D: pd.DataFrame) -> pd.Series:
    """Mean isotopologue fraction per compound, group, hours, isotopologue."""
    from .design import analysed
    X = F.merge(analysed(D)[["Sample", "Group", "Hours"]], on="Sample")
    return X.groupby(["Compound", "Group", "Hours", "i"]).Fraction.mean()


def iso_blocks(iso: pd.Series | None, cols: list, comps=None) -> dict:
    """compound -> DataFrame (rows i, columns (group, hours)) of mean fraction (%)."""
    if iso is None or not len(iso):
        return {}
    U = iso.unstack(["Group", "Hours"])
    U = U.reindex(columns=pd.MultiIndex.from_tuples(cols, names=["Group", "Hours"]))
    out = {}
    for c, d in U.groupby(level="Compound"):
        if comps is not None and c not in comps:
            continue
        out[c] = d.droplevel("Compound").sort_index()
    return out


def zscore_blocks(X: pd.DataFrame, design: pd.DataFrame, cols: list, comps=None) -> dict:
    """Row-wise Z of each channel ('Compound | M+i') over all analysed samples, mean per group × time point."""
    if X is None or X.empty:
        return {}
    Z = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1, ddof=1), axis=0).replace([np.inf, -np.inf], np.nan)
    info = design.set_index("Sample").reindex(Z.columns)
    key = pd.MultiIndex.from_arrays([info.Group, info.Hours], names=["Group", "Hours"])
    G = Z.T.groupby(key).mean().T
    G = G.reindex(columns=pd.MultiIndex.from_tuples(cols, names=["Group", "Hours"]))
    parts = G.index.to_series().str.rsplit(" | ", n=1, expand=True)
    G.index = pd.MultiIndex.from_arrays([parts[0], parts[1].str.replace("M+", "", regex=False).astype(int)],
                                        names=["Compound", "i"])
    out = {}
    for c, d in G.groupby(level="Compound"):
        if comps is not None and c not in comps:
            continue
        out[c] = d.droplevel("Compound").sort_index()
    return out


_GREEK = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "⁺": "", "₂": "2", "₃": "3", "₄": "4"}


def _nk(s: str) -> str:
    s = str(s or "")
    for a, b in _GREEK.items():
        s = s.replace(a, b)
    return re.sub(r"[^a-z0-9]", "", s.lower())


def bind(tpl: Template, A: pd.DataFrame, detected) -> dict:
    """template node id -> data compound (KEGG first, then name). Generic class nodes need a user-confirmed ID."""
    detected = set(detected)
    ok = A[A.Confidence.isin(["High", "Medium", "User-confirmed"]) & A.Compound.isin(detected)]
    by_kegg, by_kegg_user = {}, {}
    for r in ok.itertuples():
        k = str(getattr(r, "KEGG", "") or "").strip()
        if k:
            by_kegg.setdefault(k, r.Compound)
            if r.Confidence == "User-confirmed":
                by_kegg_user.setdefault(k, r.Compound)
    by_name = {}
    for r in ok.itertuples():
        for n in (r.Compound, getattr(r, "Label", ""), getattr(r, "StandardName", "")):
            if n:
                by_name.setdefault(_nk(n), r.Compound)
    for c in detected:
        by_name.setdefault(_nk(c), c)
    out = {}
    for nid, n in tpl.nodes.items():
        if n["kind"] != "met":
            continue
        keggs = str(n["kegg"]).split("|")
        src = by_kegg_user if n.get("verify") else by_kegg
        c = next((src[k] for k in keggs if k in src), None)
        if c is None and not n.get("verify"):
            c = by_name.get(_nk(n["label"]))
        if c:
            out[nid] = c
    return out


def extras(tpl: Template, A: pd.DataFrame, detected, binding: dict) -> list:
    """Detected members of the template's annotation pathways that are not on the drawing."""
    if not tpl.annot:
        return []
    detected = set(detected)
    bound = set(binding.values())
    out = []
    for r in A.itertuples():
        paths = {x.strip() for x in str(r.Pathways).split(";") if x.strip()}
        if paths & set(tpl.annot) and r.Compound in detected and r.Compound not in bound:
            out.append(r.Compound)
    return list(dict.fromkeys(out))


# =============================================================================================================
# layout
# =============================================================================================================
@dataclass
class PathMap:
    tpl: Template
    nodes: pd.DataFrame          # id, kind ('card' | 'nd' | 'stub'), label, comp, cx, cy, w, h
    edges: list                  # dict pts, color, dashed, bi
    cofactors: list              # (x, y, text, ha)
    circles: list                # (cx, cy, rx, ry, [(x, y, text)])
    boxes: list                  # (x0, y0, x1, y1, title)
    notes: list                  # (x, y, text, rotation)
    panel: tuple                 # x0, y0, x1, y1
    extent: tuple                # x0, x1, y0, y1
    cols: list
    groups: list
    panels: list
    cw: float
    extra_title: tuple | None = None


def panel_width(cols: list, cw: float) -> float:
    groups = list(dict.fromkeys(g for g, _ in cols))
    return len(cols) * cw + GGAP * max(len(groups) - 1, 0)


def panel_offsets(cols, panels, cw) -> list:
    """x offset of each panel; a panel title wider than its heatmap pushes the next panel right (no overlap)."""
    pw_ = panel_width(cols, cw)
    out, x = [], 0.0
    for p in panels:
        out.append(x)
        x += max(pw_ + PGAP, _w(PANEL_TITLES[p], PFONT) + 0.45)
    return out


def _panel_span(cols, panels, cw):
    pw_ = panel_width(cols, cw)
    off = panel_offsets(cols, panels, cw)
    return off[-1] + max(pw_, _w(PANEL_TITLES[panels[-1]], PFONT))


def card_size(label: str, n_rows: int, cols: list, panels: list, cw: float) -> tuple:
    """Card = name on top, then [M+i labels | incorporation heatmap | Z-score heatmap]."""
    heat_w = RL + _panel_span(cols, panels, cw)
    return 2 * PAD + max(_w(label, FONT), heat_w), 2 * PAD + LH + PH + GS + TD + n_rows * RH


def _box_size(kind: str, label: str, nrows: int, cols, panels, cw):
    if kind == "card":
        return card_size(label, nrows, cols, panels, cw)
    if kind == "stub":
        return _w(label, STFONT) + 1.6, 1.9
    return _w(label, NDFONT) + 1.2, 2.0


def _clip_from(box, tx, ty, extra=0.0):
    """Point where the line from the centre of box (cx, cy, w, h) towards (tx, ty) leaves the box (+extra)."""
    cx, cy, w, h = box
    dx, dy = tx - cx, ty - cy
    if abs(dx) < 1e-12 and abs(dy) < 1e-12:
        return cx, cy
    t = min((w / 2 + extra) / abs(dx) if abs(dx) > 1e-12 else np.inf,
            (h / 2 + extra) / abs(dy) if abs(dy) > 1e-12 else np.inf)
    return cx + dx * t, cy + dy * t


def _inside(pt, box, extra=0.0):
    cx, cy, w, h = box
    return abs(pt[0] - cx) < w / 2 + extra and abs(pt[1] - cy) < h / 2 + extra


def layout(tpl: Template, binding: dict, rows: dict, cols: list, panels: list, labels: dict | None = None,
           extra: list | None = None, values_in_cells: bool = False, max_iter: int = 250) -> PathMap:
    """rows: compound -> isotopologue indices shown; labels: compound -> display name."""
    labels = labels or {}
    extra = [c for c in (extra or []) if c in rows]
    cw = CW_VAL if values_in_cells else CW
    groups = list(dict.fromkeys(g for g, _ in cols))
    # ---------- boxes
    recs = []
    for nid, n in tpl.nodes.items():
        comp = binding.get(nid)
        if n["kind"] == "stub":
            kind, lab = "stub", n["label"]
        elif comp in rows:
            kind, lab = "card", n["label"]
        else:
            kind, lab = "nd", n["label"]
        w, h = _box_size(kind, lab, len(rows.get(comp, [])), cols, panels, cw)
        recs.append(dict(id=nid, kind=kind, label=lab, comp=comp if kind == "card" else None,
                         bx=n["x"], by=-n["y"], w=w, h=h))
    for k, (cx, cy, r, lines) in enumerate(tpl.circles):       # text in a ring centre is kept free
        if lines:
            recs.append(dict(id=f"__circ{k}", kind="ctext", label="", comp=None, bx=cx, by=-cy,
                             w=max(_w(t, CIRC_FONT) for t in lines) + 1.0, h=len(lines) * CIRC_FONT * 1.25 + 0.8))
    N = pd.DataFrame(recs)
    bx, by = N.bx.to_numpy(float), N.by.to_numpy(float)
    W, H = N.w.to_numpy(float), N.h.to_numpy(float)
    idx = {r: i for i, r in enumerate(N.id)}
    # labelled straight edges need room for their label between the two boxes
    need = []
    for e in tpl.edges:                       # every straight arrow keeps a visible shaft (+ room for its label)
        if not e["via"] and "arc" not in e:
            i, j = idx[e["a"]], idx[e["b"]]
            heads = HEAD_L * (2 if e["bi"] else 1)
            gap = max(_w(e["label"], CFONT) + 1.2 + heads if e.get("label") else 0.0, 1.6 + heads)
            need.append((i, j, gap))
    # ---------- expand: grow x / y spacing until nothing overlaps
    sx = sy = 0.02
    for _ in range(max_iter):
        cx, cy = bx * sx, by * sy
        dx = np.abs(cx[:, None] - cx[None, :]); dy = np.abs(cy[:, None] - cy[None, :])
        nx = (W[:, None] + W[None, :]) / 2 + MARGIN; ny = (H[:, None] + H[None, :]) / 2 + MARGIN
        ov = np.triu((dx < nx - 1e-9) & (dy < ny - 1e-9), 1)
        grow_x = grow_y = 1.0
        if ov.any():
            i, j = np.where(ov)
            bdx, bdy = np.abs(bx[i] - bx[j]), np.abs(by[i] - by[j])
            fx = np.where(bdx > 1e-9, nx[i, j] / np.maximum(bdx * sx, 1e-12), np.inf)
            fy = np.where(bdy > 1e-9, ny[i, j] / np.maximum(bdy * sy, 1e-12), np.inf)
            use_x = fx <= fy
            if use_x.any():
                grow_x = max(grow_x, min(float(fx[use_x].max()), 3.0))
            if (~use_x).any():
                grow_y = max(grow_y, min(float(fy[~use_x].max()), 3.0))
        for i, j, lw_ in need:
            ddx, ddy = abs(bx[i] - bx[j]), abs(by[i] - by[j])
            if ddx >= ddy and ddx > 0:          # horizontal-ish: label sits above the arrow
                req = (W[i] + W[j]) / 2 + lw_
                if ddx * sx < req:
                    grow_x = max(grow_x, req / (ddx * sx))
            elif ddy > 0:                        # vertical-ish: label beside the arrow
                req = (H[i] + H[j]) / 2 + min(lw_, 3.0)
                if ddy * sy < req:
                    grow_y = max(grow_y, req / (ddy * sy))
        if grow_x <= 1.0 + 1e-9 and grow_y <= 1.0 + 1e-9:
            break
        sx *= grow_x * 1.0001 if grow_x > 1 else 1.0
        sy *= grow_y * 1.0001 if grow_y > 1 else 1.0
    if tpl.circles or any("arc" in e for e in tpl.edges):   # a ring stays round enough to read as a cycle
        sy = max(sy, sx / 1.3)
        sx = max(sx, sy / 1.3)
    cx, cy = bx * sx, by * sy
    # vertical compaction (layouts without a ring): every row of the drawing gets exactly the height of its tallest
    # box plus the gap, so rows of small labels do not inherit the spacing needed by rows of cards
    fy = lambda yy: -np.asarray(yy, float) * sy
    if not any("arc" in e for e in tpl.edges) and not tpl.circles:
        ty = -by                                           # template y (down)
        order = np.argsort(ty)
        lev = np.zeros(len(ty), int); k, last = 0, None
        for i_ in order:
            if last is not None and ty[i_] - last > 12:
                k += 1
            lev[i_] = k; last = ty[i_]
        n_lev = lev.max() + 1
        lev_y = np.array([ty[lev == q].mean() for q in range(n_lev)])
        lev_h = np.array([H[lev == q].max() for q in range(n_lev)])
        newy = np.zeros(n_lev)
        for q in range(1, n_lev):
            step = lev_h[q - 1] / 2 + lev_h[q] / 2 + ROW_GAP
            newy[q] = newy[q - 1] - step
        # linear extrapolation beyond the first / last row (corner points above or below the drawing)
        k0 = (newy[1] - newy[0]) / (lev_y[1] - lev_y[0]) if n_lev > 1 else -sy
        k1 = (newy[-1] - newy[-2]) / (lev_y[-1] - lev_y[-2]) if n_lev > 1 else -sy
        ly = np.r_[lev_y[0] - 1e4, lev_y, lev_y[-1] + 1e4]
        ny_ = np.r_[newy[0] - k0 * 1e4, newy, newy[-1] + k1 * 1e4]
        fy = lambda yy: np.interp(np.asarray(yy, float), ly, ny_)
        cy = fy(ty)
    N = N.assign(cx=cx, cy=cy)
    tf = lambda x, y: (x * sx, float(fy(y)))
    box = {r.id: (r.cx, r.cy, r.w, r.h) for r in N.itertuples()}
    real = N[N.kind != "ctext"]
    obstacles = [box[r.id] for r in real.itertuples()]

    # ---------- edges
    edges, cof = [], []
    for e in tpl.edges:
        a, b = box[e["a"]], box[e["b"]]
        col = e.get("color") or tpl.arrow
        if "arc" in e:
            acx, acy, r, a1, a2 = e["arc"]
            step = -0.5 if a2 < a1 else 0.5
            angs = np.arange(a1, a2 + step / 2, step)
            pts = [tf(acx + r * math.cos(math.radians(t)), acy - r * math.sin(math.radians(t))) for t in angs]
            k0 = 0
            while k0 < len(pts) and _inside(pts[k0], a, 0.3):
                k0 += 1
            k1 = len(pts) - 1
            while k1 >= 0 and _inside(pts[k1], b, 0.45):
                k1 -= 1
            pts = pts[k0:k1 + 1]
            if len(pts) < 3:
                continue
            edges.append(dict(pts=pts, color=col, dashed=False, bi=e["bi"]))
            continue
        tn_a, tn_b = tpl.nodes[e["a"]], tpl.nodes[e["b"]]
        vias = []
        for vx, vy in e["via"]:
            px_, py_ = tf(vx, vy)
            for tn, bb in ((tn_a, a), (tn_b, b)):        # keep orthogonal routing aligned with the nodes
                if abs(vx - tn["x"]) < 30:
                    px_ = bb[0]
                if abs(vy - tn["y"]) < 30:
                    py_ = bb[1]
            vias.append([px_, py_, vx, vy])
        for _ in range(4):                              # push corner points out of cards
            moved = False
            for v in vias:
                for ob in obstacles:
                    if _inside(v[:2], ob, 1.0):
                        ox = ob[2] / 2 + 1.0 - abs(v[0] - ob[0]); oy = ob[3] / 2 + 1.0 - abs(v[1] - ob[1])
                        if ox <= oy:
                            d = ox * (1 if v[0] >= ob[0] else -1)
                            for u in vias:
                                if abs(u[2] - v[2]) < 1:
                                    u[0] += d
                        else:
                            d = oy * (1 if v[1] >= ob[1] else -1)
                            for u in vias:
                                if abs(u[3] - v[3]) < 1:
                                    u[1] += d
                        moved = True
            if not moved:
                break
        mid = [(v[0], v[1]) for v in vias]
        nxt = mid[0] if mid else b[:2]
        prv = mid[-1] if mid else a[:2]
        p0 = _clip_from(a, *nxt, extra=0.18)
        pn = _clip_from(b, *prv, extra=0.28)
        pts = [p0] + mid + [pn]
        edges.append(dict(pts=pts, color=col, dashed=e["style"] == "dashed", bi=e["bi"], via=e["via"]))
        if e.get("label"):
            (x1, y1), (x2, y2) = pts[0], pts[1]
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2
            side = e.get("side", 1)
            L_ = math.hypot(x2 - x1, y2 - y1) or 1.0
            unx, uny = -(y2 - y1) / L_, (x2 - x1) / L_          # unit normal of the arrow
            if abs(uny) > 0.92:                                   # horizontal-ish: above (side 1) / below
                cof.append((mx, my + (0.62 if side > 0 else -0.62), e["label"], "center"))
            else:                                                 # vertical or diagonal: beside the arrow
                if unx < 0:
                    unx, uny = -unx, -uny
                sgn = 1 if side > 0 else -1
                cof.append((mx + sgn * unx * 0.6, my + sgn * uny * 0.6, e["label"], "left" if sgn > 0 else "right"))
    # ---------- circles, sub-panels, notes
    circles = []
    for k, (ccx, ccy, r, lines) in enumerate(tpl.circles):
        X0, Y0 = tf(ccx, ccy)
        texts = [(X0, Y0 + (len(lines) - 1) * CIRC_FONT * 0.62 - q * CIRC_FONT * 1.25, t) for q, t in enumerate(lines)]
        circles.append((X0, Y0, r * sx, r * sy, texts))
    boxes = []
    for x0, y0, bw, bh, title in tpl.boxes:
        mem = [nid for nid, n in tpl.nodes.items() if x0 <= n["x"] <= x0 + bw and y0 <= n["y"] <= y0 + bh]
        if not mem:
            continue
        L = min(box[m][0] - box[m][2] / 2 for m in mem) - 1.2
        Rr = max(box[m][0] + box[m][2] / 2 for m in mem) + 1.2
        B = min(box[m][1] - box[m][3] / 2 for m in mem) - 1.2
        Tt = max(box[m][1] + box[m][3] / 2 for m in mem) + 1.2 + BFONT * 1.8
        Rr = max(Rr, L + _w(title, BFONT) + 1.6)
        boxes.append((L, B, Rr, Tt, title))
    notes = []
    for x, y, text, rot in tpl.notes:
        X0, Y0 = tf(x, y)
        for e in edges:                          # follow a nearby (moved) routing line
            for vx_vy, pt in zip(e.get("via", []), e["pts"][1:-1]):
                if abs(vx_vy[0] - x) < 30:
                    X0 = pt[0] + 0.9
        notes.append((X0, Y0, text, rot))
    # ---------- extra block underneath
    xs = list(real.cx - real.w / 2) + list(real.cx + real.w / 2)
    ys = list(real.cy - real.h / 2) + list(real.cy + real.h / 2)
    for e in edges:
        xs += [p_[0] for p_ in e["pts"]]; ys += [p_[1] for p_ in e["pts"]]
    for c in circles:
        xs += [c[0] - c[2], c[0] + c[2]]; ys += [c[1] - c[3], c[1] + c[3]]
    for b in boxes:
        xs += [b[0], b[2]]; ys += [b[1], b[3]]
    for x, y, t_, ha in cof:
        w_ = _w(t_, CFONT)
        x0_ = x - w_ / 2 if ha == "center" else (x if ha == "left" else x - w_)
        xs += [x0_, x0_ + w_]; ys += [y - 0.5, y + 0.5]
    for x, y, t_, rot in notes:
        L_ = _w(t_, CFONT)
        xs += [x - 0.5, x + 0.5] if rot else [x, x + L_]
        ys += [y - L_ / 2, y + L_ / 2] if rot else [y - 0.5, y + 0.5]
    x_left, x_right, y_low = min(xs), max(xs), min(ys)
    extra_title = None
    if extra:
        X = []
        for comp in extra:
            lab = labels.get(comp, comp)
            w, h = card_size(lab, len(rows.get(comp, [])), cols, panels, cw)
            X.append(dict(id=f"x:{comp}", kind="card", label=lab, comp=comp, w=w, h=h))
        X = pd.DataFrame(X)
        cell_w = float(X.w.max()) + MARGIN
        span = max(x_right - x_left, cell_w * min(len(X), 3))
        ncol = max(1, min(len(X), int((span + MARGIN) // cell_w)))
        top = y_low - 3.0
        extra_title = (x_left, top, "Also detected in this pathway group (not on the drawing)")
        top -= 2.0
        yy = top
        for k0 in range(0, len(X), ncol):
            chunk = X.iloc[k0:k0 + ncol]
            for k, r in enumerate(chunk.itertuples()):
                X.loc[r.Index, "cx"] = x_left + k * cell_w + r.w / 2
                X.loc[r.Index, "cy"] = yy - r.h / 2
            yy -= float(chunk.h.max()) + MARGIN
        N = pd.concat([N, X], ignore_index=True)
        xs += list(X.cx - X.w / 2) + list(X.cx + X.w / 2) + [x_left + _w(extra_title[2], 0.95)]
        ys += list(X.cy - X.h / 2) + [top + 2.0]
    nodes = N[N.kind != "ctext"][["id", "kind", "label", "comp", "cx", "cy", "w", "h"]].reset_index(drop=True)
    x0, x1 = min(xs) - PANEL_PAD, max(xs) + PANEL_PAD
    y0, y1 = min(ys) - PANEL_PAD, max(ys) + TITLE_BAND
    x1 = max(x1, x0 + _w(tpl.title, TFONT) + 2 * PANEL_PAD)
    return PathMap(tpl, nodes, edges, cof, circles, boxes, notes, (x0, y0, x1, y1),
                   (x0 - 0.6, x1 + 0.6, y0 - 0.6, y1 + 0.6), list(cols), groups, list(panels), cw, extra_title)


# =============================================================================================================
# card geometry (shared by both renderers)
# =============================================================================================================
def card_parts(r, blocks: dict, rows: list, pm: PathMap):
    x0, y1 = r.cx - r.w / 2, r.cy + r.h / 2
    out = dict(label=[(x0 + PAD, y1 - PAD - LH / 2, r.label)], titles=[], strips=[], tdots=[], rowlabels=[],
               cells=[])
    if r.kind != "card":
        return out
    hx0 = x0 + PAD + RL
    offs = panel_offsets(pm.cols, pm.panels, pm.cw)
    t_top = y1 - PAD - LH                         # heatmaps directly under the name
    s_bot = t_top - PH - GS
    htop = s_bot - TD
    for ri, i in enumerate(rows):
        out["rowlabels"].append((f"M+{int(i)}", hx0 - 0.25, htop - ri * RH - RH / 2))
    for k, panel in enumerate(pm.panels):
        px = hx0 + offs[k]
        out["titles"].append((PANEL_TITLES[panel], px, t_top - PH / 2))
        B = blocks.get(panel, {}).get(r.comp)
        x = px
        for g in pm.groups:
            cs = [c for c in pm.cols if c[0] == g]
            if not cs:
                continue
            out["strips"].append((g, x + 0.04, s_bot + 0.1, x + len(cs) * pm.cw - 0.04, s_bot + GS - 0.06))
            for c in cs:
                out["tdots"].append((x + pm.cw / 2, s_bot - TD / 2, c))
                for ri, i in enumerate(rows):
                    v = np.nan
                    if B is not None and i in B.index and c in B.columns:
                        v = float(B.loc[i, c])
                    yt = htop - ri * RH
                    out["cells"].append((x + 0.04, yt - RH + 0.04, x + pm.cw - 0.04, yt - 0.04, v, int(i), c, panel))
                x += pm.cw
            x += GGAP
    return out


def _arrowheads(pts, bi):
    def head(tip, prev):
        dx, dy = tip[0] - prev[0], tip[1] - prev[1]
        L = math.hypot(dx, dy) or 1.0
        ux, uy = dx / L, dy / L
        bx_, by_ = tip[0] - ux * HEAD_L, tip[1] - uy * HEAD_L
        return [tip, (bx_ - uy * HEAD_W, by_ + ux * HEAD_W), (bx_ + uy * HEAD_W, by_ - ux * HEAD_W)]
    out = [head(pts[-1], pts[-2])]
    if bi:
        out.append(head(pts[0], pts[1]))
    return out


def _shorten(pts, bi):
    """Stop the line at the arrowhead base so heads stay crisp."""
    pts = [tuple(p) for p in pts]

    def cut(tip, prev):
        dx, dy = tip[0] - prev[0], tip[1] - prev[1]
        L = math.hypot(dx, dy) or 1.0
        k = min(HEAD_L * 0.7, L * 0.5) / L
        return (tip[0] - dx * k, tip[1] - dy * k)
    out = list(pts)
    out[-1] = cut(pts[-1], pts[-2])
    if bi:
        out[0] = cut(pts[0], pts[1])
    return out


def column_key(cols: list) -> str:
    groups = list(dict.fromkeys(g for g, _ in cols))
    return " | ".join(f"{g}: " + " · ".join(f"{h:g} h" for gg, h in cols if gg == g) for g in groups)


def _fmt(v, panel):
    if not np.isfinite(v):
        return "n.d."
    return f"{v:.1f}" if panel == "iso" else f"{v:+.2f}"


def _rounded(x0, y0, x1, y1, r):
    """Rounded rectangle as a list of points (for Plotly paths)."""
    r = min(r, (x1 - x0) / 2, (y1 - y0) / 2)
    pts = []
    for (cx_, cy_, a0) in ((x1 - r, y1 - r, 0), (x0 + r, y1 - r, 90), (x0 + r, y0 + r, 180), (x1 - r, y0 + r, 270)):
        for t in np.linspace(a0, a0 + 90, 7):
            pts.append((cx_ + r * math.cos(math.radians(t)), cy_ + r * math.sin(math.radians(t))))
    return pts


def _ellipse(cx, cy, rx, ry, n=120):
    return [(cx + rx * math.cos(t), cy + ry * math.sin(t)) for t in np.linspace(0, 2 * math.pi, n)]


def _lighten(hex_color, f):
    import matplotlib.colors as mc
    r, g, b = mc.to_rgb(hex_color)
    return mc.to_hex((r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f))


# =============================================================================================================
# matplotlib (publication)
# =============================================================================================================
def render_mpl(pm: PathMap, blocks: dict, rows: dict, *, iso_vmax=50.0, z_lim=2.0, iso_cmap="Omnix purple",
               z_cmap="Omnix blue–red", subtitle="", width=None, pt_per_unit=6.5, values_in_cells=False,
               z_label="Z-score", tips: dict | None = None):
    """
    Width in inches (fonts scale with the map, so the layout never overlaps); default 6.5 pt per heatmap row.
    tips: if a dict is given, every heatmap cell is drawn as its own element with an id and its hover text is stored
    in tips[id] (used by the interactive vector viewer).
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Polygon, Rectangle
    tpl = pm.tpl
    x0, x1, y0, y1 = pm.extent
    n_cb = int("iso" in pm.panels) + int("z" in pm.panels)
    leg_in = 0.45 + 1.15 * n_cb
    ptu = max((float(width) - leg_in) * 72 / (x1 - x0), 1.0) if width else pt_per_unit
    map_w, map_h = (x1 - x0) * ptu / 72, (y1 - y0) * ptu / 72
    sub_lines = []
    if subtitle:                                   # wrap the subtitle so it never runs past the figure
        max_chars = max(int((map_w + leg_in - 0.3) * 72 / (8 * 0.55)), 30)
        cur = ""
        for part in subtitle.split(" · "):
            nxt = part if not cur else cur + " · " + part
            if len(nxt) > max_chars and cur:
                sub_lines.append(cur); cur = part
            else:
                cur = nxt
        sub_lines.append(cur)
    sub_h = 0.17 * len(sub_lines) + (0.05 if sub_lines else 0)
    head_in = 0.3 + sub_h + 0.3
    W_in = map_w + leg_in
    cb_h = min(1.9, max(1.1, 0.3 * map_h))
    H_in = max(map_h, cb_h + 0.4) + head_in
    fig = plt.figure(figsize=(W_in, H_in))
    ax = fig.add_axes([0, 0, map_w / W_in, map_h / H_in])
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal"); ax.axis("off")
    gc = group_colors(pm.groups)
    tcol = time_colors([h for _, h in pm.cols])
    icm, zcm = get_cmap(iso_cmap, "iso"), get_cmap(z_cmap, "z")
    fs = lambda units: units * ptu * 0.95
    lw = max(0.6, 0.12 * ptu)
    px0, py0, px1, py1 = pm.panel
    ax.add_patch(FancyBboxPatch((px0, py0), px1 - px0, py1 - py0, boxstyle="round,pad=0,rounding_size=1.6",
                                fc=tpl.fill, ec=tpl.border, lw=1.0, zorder=0))
    ax.text(px0 + PANEL_PAD, py1 - 0.7, tpl.title, fontsize=fs(TFONT), fontweight="bold", color="black",
            ha="left", va="top", zorder=6)
    for (L, B, Rr, Tt, title) in pm.boxes:
        ax.add_patch(FancyBboxPatch((L, B), Rr - L, Tt - B, boxstyle="round,pad=0,rounding_size=1.1",
                                    fc=_lighten(tpl.border, 0.55), ec=tpl.border, lw=0.7, alpha=0.9, zorder=0.5))
        ax.text(L + 0.8, Tt - 0.4, title, fontsize=fs(BFONT), fontweight="bold", color="black", ha="left",
                va="top", zorder=6)
    for (ccx, ccy, rx, ry, texts) in pm.circles:
        ax.add_patch(Ellipse((ccx, ccy), 2 * rx, 2 * ry, fc=_lighten(tpl.border, 0.45), ec=tpl.border, lw=1.0,
                             alpha=0.75, zorder=0.6))
        for (tx, ty, tt) in texts:
            ax.text(tx, ty, tt, fontsize=fs(CIRC_FONT), fontweight="bold", color="black", ha="center",
                    va="center", zorder=0.7)
    for e in pm.edges:
        P_ = np.array(_shorten(e["pts"], e["bi"]))
        ax.add_line(Line2D(P_[:, 0], P_[:, 1], color=e["color"], lw=lw, ls=(0, (4, 3)) if e["dashed"] else "-",
                           solid_capstyle="round", zorder=1))
        for tri in _arrowheads(e["pts"], e["bi"]):
            ax.add_patch(Polygon(tri, closed=True, fc=e["color"], ec="none", zorder=1.1))
    for x, y, t_, ha in pm.cofactors:
        ax.text(x, y, t_, fontsize=fs(CFONT), style="italic", color="black", ha=ha, va="center", zorder=1.2)
    for x, y, t_, rot in pm.notes:
        ax.text(x, y, t_, fontsize=fs(CFONT), style="italic", color="black", ha="center", va="center",
                rotation=rot, zorder=1.2)
    if pm.extra_title:
        ex, ey, et = pm.extra_title
        ax.text(ex, ey, et, fontsize=fs(0.95), fontweight="bold", color="black", ha="left", va="top", zorder=2)
    rects, colors = [], []
    for r in pm.nodes.itertuples():
        if r.kind == "nd":
            ax.add_patch(FancyBboxPatch((r.cx - r.w / 2, r.cy - r.h / 2), r.w, r.h,
                                        boxstyle=f"round,pad=0,rounding_size={r.h / 2}", fc="white", alpha=0.55,
                                        ec="#8d8a99", lw=0.5, ls=":", zorder=2))
            ax.text(r.cx, r.cy, r.label, fontsize=fs(NDFONT), style="italic", color="black", ha="center",
                    va="center", zorder=4)
            continue
        if r.kind == "stub":
            ax.add_patch(FancyBboxPatch((r.cx - r.w / 2, r.cy - r.h / 2), r.w, r.h,
                                        boxstyle=f"round,pad=0,rounding_size={r.h / 2}", fc="white", alpha=0.75,
                                        ec=tpl.arrow, lw=0.6, ls=(0, (3, 2)), zorder=2))
            ax.text(r.cx, r.cy, r.label, fontsize=fs(STFONT), style="italic", color="black", ha="center",
                    va="center", zorder=4)
            continue
        ax.add_patch(FancyBboxPatch((r.cx - r.w / 2, r.cy - r.h / 2), r.w, r.h,
                                    boxstyle="round,pad=0,rounding_size=0.45", fc="white", ec=tpl.arrow, lw=0.7,
                                    zorder=2))
        rr = rows.get(r.comp, [])
        P = card_parts(r, blocks, rr, pm)
        for lx, ly, lt in P["label"]:
            t_obj = ax.text(lx, ly, lt, fontsize=fs(FONT), va="center", ha="left", color="black", fontweight="bold",
                            zorder=4)
            if tips is not None:
                gid = f"lab{len(tips)}"
                t_obj.set_gid(gid)
                tips[gid] = f"{r.label} — {r.comp}"
        for t_, x, y in P["titles"]:
            ax.text(x, y, t_, fontsize=fs(PFONT), color="black", va="center", ha="left", zorder=4)
        for x_, y_, (g_, h_) in P["tdots"]:
            ax.add_patch(Circle((x_, y_), DOT, fc=tcol[h_], ec="black", lw=0.2, zorder=4))
        for g, a, b, c, d in P["strips"]:
            ax.add_patch(Rectangle((a, b), c - a, d - b, fc=gc.get(g, "#888"), lw=0, zorder=4))
        for t_, x, y in P["rowlabels"]:
            ax.text(x, y, t_, fontsize=fs(RFONT), color="black", va="center", ha="right", zorder=4)
        for a, b, c, d, v, i, col, panel in P["cells"]:
            hx = _hex(icm, v, 0, iso_vmax) if panel == "iso" else _hex(zcm, v, -z_lim, z_lim)
            if tips is None:
                rects.append(Rectangle((a, b), c - a, d - b))
                colors.append(hx)
            else:
                gid = f"cell{len(tips)}"
                ax.add_patch(Rectangle((a, b), c - a, d - b, fc=hx, ec="none", zorder=3.5, gid=gid))
                val = ("n.d." if not np.isfinite(v) else f"Incorporation {v:.2f} %" if panel == "iso"
                       else f"{z_label} {v:+.2f}")
                tips[gid] = f"{r.label} ({r.comp}) · M+{i} · {col[0]} · {col[1]:g} h — {val}"
            if values_in_cells:
                ax.text((a + c) / 2, (b + d) / 2, _fmt(v, panel), fontsize=fs(VFONT), ha="center", va="center",
                        color=_ink_on(hx), zorder=5)
    if rects:
        ax.add_collection(PatchCollection(rects, facecolors=colors, edgecolors="none", zorder=3.5))
    fw = lambda inch: inch / W_in
    fh = lambda inch: inch / H_in
    k = 0
    for panel, cmap, lo, hi, lab in (("iso", icm, 0, iso_vmax, "Incorporation (%)"),
                                     ("z", zcm, -z_lim, z_lim, "Z-score")):
        if panel not in pm.panels:
            continue
        cax = fig.add_axes([fw(map_w + 0.35 + 1.15 * k), fh(0.15), fw(0.14), fh(cb_h)])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(lo, hi), cmap=cmap), cax=cax)
        cb.ax.tick_params(labelsize=7, colors="black", length=2)
        cb.outline.set_linewidth(0.5)
        cax.set_title(lab.replace(" (", "\n("), fontsize=7.5, color="black", loc="left", pad=4)
        k += 1
    if subtitle:
        fig.text(fw(0.1), 1 - fh(0.08), "\n".join(sub_lines), fontsize=8, va="top", color="black", linespacing=1.3)
    hs = [Rectangle((0, 0), 1, 1, color=gc.get(g, "#888")) for g in pm.groups] + \
         [Line2D([], [], marker="o", ls="", ms=6, mfc=c, mec="black", mew=0.3) for c in tcol.values()]
    ls_ = list(pm.groups) + [f"{h:g} h" for h in tcol]
    fig.legend(hs, ls_, loc="upper left", bbox_to_anchor=(fw(0.08), 1 - fh(0.1 + sub_h)),
               ncol=len(ls_), frameon=False, fontsize=7.5, handlelength=1.0, columnspacing=1.2, labelcolor="black",
               borderaxespad=0)
    return fig


# =============================================================================================================
# interactive vector viewer (text scales with zoom; hover = SVG <title>)
# =============================================================================================================
def render_svg(pm: PathMap, blocks: dict, rows: dict, **kw) -> str:
    """The publication figure as SVG (text kept as text), with a hover title on every heatmap cell and name."""
    import io
    import matplotlib.pyplot as plt
    from xml.sax.saxutils import escape
    tips = {}
    fig = render_mpl(pm, blocks, rows, tips=tips, **kw)
    buf = io.StringIO()
    with plt.rc_context({"svg.fonttype": "none", "svg.hashsalt": "omnix"}):
        fig.savefig(buf, format="svg", facecolor="white")
    plt.close(fig)
    svg = buf.getvalue()
    svg = svg[svg.index("<svg"):]
    def add_title(m):
        gid = m.group(1)
        return f'{m.group(0)}<title>{escape(tips[gid])}</title>' if gid in tips else m.group(0)
    svg = re.sub(r'<g id="((?:cell|lab)\d+)">', add_title, svg)
    return svg


VIEWER = """<div style="font-family:Arial,Helvetica,sans-serif">
<div id="pmv-bar" style="display:flex;gap:6px;align-items:center;margin:0 0 6px 2px">
<button data-a="in" title="Zoom in">+</button><button data-a="out" title="Zoom out">−</button>
<button data-a="fit" title="Fit to width">Fit</button><button data-a="one" title="Actual size">100 %</button>
<span id="pmv-z" style="font-size:12px;color:#000;margin-left:4px"></span>
<span style="font-size:12px;color:#000;margin-left:12px">Drag to move · wheel to zoom · hover a cell for its value</span></div>
<div id="pmv" style="position:relative;border:1px solid #d9d7e3;border-radius:8px;background:#fff;
height:__H__px;overflow:hidden;cursor:grab">
<div id="pmv-stage" style="transform-origin:0 0;position:absolute;left:0;top:0">__SVG__</div></div></div>
<style>#pmv-bar button{font:600 13px Arial;color:#000;background:#fff;border:1px solid #b9b2e3;border-radius:5px;
padding:3px 10px;cursor:pointer}#pmv-bar button:hover{background:#efedfa}</style>
<script>(function(){const w=document.getElementById('pmv'),st=document.getElementById('pmv-stage'),
zl=document.getElementById('pmv-z'),svg=st.querySelector('svg');let s=1,x=0,y=0,drag=null;
const W=svg.width.baseVal.value||svg.getBoundingClientRect().width;
function ap(){st.style.transform=`translate(${x}px,${y}px) scale(${s})`;zl.textContent=Math.round(s*100)+' %';}
function fit(){s=Math.min(1.6,(w.clientWidth-16)/W);x=8;y=8;ap();}
function zoom(f,cx,cy){const ns=Math.max(0.1,Math.min(12,s*f));x=cx-(cx-x)*ns/s;y=cy-(cy-y)*ns/s;s=ns;ap();}
w.addEventListener('wheel',e=>{e.preventDefault();const r=w.getBoundingClientRect();
zoom(e.deltaY<0?1.15:1/1.15,e.clientX-r.left,e.clientY-r.top);},{passive:false});
w.addEventListener('mousedown',e=>{drag=[e.clientX-x,e.clientY-y];w.style.cursor='grabbing';});
window.addEventListener('mousemove',e=>{if(!drag)return;x=e.clientX-drag[0];y=e.clientY-drag[1];ap();});
window.addEventListener('mouseup',()=>{drag=null;w.style.cursor='grab';});
document.getElementById('pmv-bar').addEventListener('click',e=>{const a=e.target.dataset.a;if(!a)return;
const cx=w.clientWidth/2,cy=w.clientHeight/2;if(a==='in')zoom(1.25,cx,cy);if(a==='out')zoom(0.8,cx,cy);
if(a==='fit')fit();if(a==='one'){s=1;x=8;y=8;ap();}});fit();})();</script>"""


def viewer_html(svg: str, height: int = 900) -> str:
    """Pan (drag) / zoom (wheel, buttons) viewer: the vector map scales as a whole, so text always fits."""
    return VIEWER.replace("__H__", str(int(height))).replace("__SVG__", svg)


def values_table(blocks: dict, rows: dict, cols: list, labels: dict | None = None) -> pd.DataFrame:
    labels = labels or {}
    recs = []
    for comp, rr in rows.items():
        for i in rr:
            rec = {"Metabolite": labels.get(comp, comp), "Compound": comp, "Isotopologue": f"M+{int(i)}"}
            for panel, name in (("iso", "Incorporation (%)"), ("z", "Z-score")):
                B = blocks.get(panel, {}).get(comp)
                for g, h in cols:
                    v = B.loc[i, (g, h)] if B is not None and i in B.index and (g, h) in B.columns else np.nan
                    rec[f"{g} · {h:g} h — {name}"] = v
            recs.append(rec)
    return pd.DataFrame(recs)
