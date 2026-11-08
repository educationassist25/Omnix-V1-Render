"""
rewiring_figures.py - Publication figures for the Correlation Rewiring Map.

The interactive view (rewiring_assets/rewiring_view.html) is for exploring; these functions
draw the same network and pair correlation plot with matplotlib so they can be downloaded in
any common format (PNG, JPEG, TIFF, SVG, PDF, EPS) at print resolution.

network_groups() is the single definition of which pathways (and therefore which metabolites)
the network shows: the top-N pathways of the pathway-level test after a p-value or FDR cutoff.
The interactive view implements the same rule in JavaScript (tested to agree).
"""

from __future__ import annotations

import io
import math
import textwrap

import numpy as np
import pandas as pd

FORMATS = {                     # label -> (matplotlib format, mime, extension)
    "PNG": ("png", "image/png", "png"),
    "JPEG": ("jpeg", "image/jpeg", "jpg"),
    "TIFF": ("tiff", "image/tiff", "tif"),
    "SVG": ("svg", "image/svg+xml", "svg"),
    "PDF": ("pdf", "application/pdf", "pdf"),
    "EPS": ("eps", "application/postscript", "eps"),
}
TOP_CHOICES = {"Top 5": 5, "Top 10": 10, "Top 15": 15, "Top 20": 20, "All pathways": 0}

# light-theme colours of the interactive view (solid colours: EPS has no transparency)
EDGE = {"lost": "#4a3aa7", "gained": "#1baf7a", "flipped": "#e87ba4", "intact": "#c9cfcc"}
DASH = {"lost": (0, (5, 3)), "gained": "solid", "flipped": (0, (1.5, 2)), "intact": "solid"}
GROUP = ("#2a78d6", "#eb6834")
INK, INK2, LINE = "#141a1b", "#4b5557", "#b9c1be"
DIV = ("#2a78d6", "#f0efec", "#e34948")


def display_label(name: str) -> str:
    return f"Relative Abundance ({name})"


# ---------------------------------------------------------------------------
# which pathways / metabolites the network shows
# ---------------------------------------------------------------------------
def network_groups(result, top: int = 10, metric: str = "p", cutoff: float = 0.05):
    """Returns (group_of, order, selected):
      group_of : ring group per feature, or None when the feature is hidden
      order    : ring groups in drawing order
      selected : DataFrame of the selected pathway tests (empty when showing all)
    With top = 0 (all pathways) or no pathway tests, every feature is shown under its ring
    pathway. Otherwise the pathways passing `cutoff` on `metric` ("p" or "q") are ranked by
    that metric (ties: larger statistic first) and the first `top` are kept; each metabolite is
    drawn once, under the best-ranked selected pathway that contains it."""
    P = result.prep
    pt = result.pathway_tests if result.pathway_tests is not None else pd.DataFrame()
    m = len(P.features)
    if not top or pt.empty:
        order = list(dict.fromkeys(P.pathways))
        return list(P.pathways), order, pd.DataFrame()
    col = "q_value" if metric == "q" else "p_value"
    sel = pt[pt[col] <= cutoff + 1e-12].copy()
    sel = sel.assign(_neg=-sel["statistic"]).sort_values([col, "_neg"], kind="stable").head(int(top))
    group_of = [None] * m
    order = []
    for r in sel.itertuples():
        placed = False
        for k in r.feature_idx:
            if group_of[k] is None:
                group_of[k] = r.pathway
                placed = True
        if placed:
            order.append(r.pathway)
    return group_of, order, sel.drop(columns=["_neg"])


def _layout(group_of, order, radius=1.0):
    """Ring positions (angle, x, y) in y-up coordinates, starting at the top and running clockwise,
    plus the arcs (group, start angle, end angle) and the angular step between neighbours."""
    pos = {}
    groups = [g for g in order if any(x == g for x in group_of)]
    members = {g: [i for i, x in enumerate(group_of) if x == g] for g in groups}
    n = sum(len(v) for v in members.values())
    if n == 0:
        return pos, [], 0.0
    gap = min(0.12, 0.9 / len(groups)) if len(groups) > 1 else 0.0
    step = (2 * math.pi - gap * len(groups)) / n
    a = math.pi / 2          # start at the top, clockwise
    arcs = []
    for g in groups:
        start = a
        for i in members[g]:
            ang = a - step / 2
            pos[i] = (ang, radius * math.cos(ang), radius * math.sin(ang))
            a -= step
        arcs.append((g, start, a))
        a -= gap
    return pos, arcs, step


def _div_colour(v):
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("div", list(DIV))
    return cmap((max(-1.5, min(1.5, v)) + 1.5) / 3.0)


# ---------------------------------------------------------------------------
# collision-free label placement (same rules as the interactive view)
# ---------------------------------------------------------------------------
# All geometry is in points (1 data unit = 1 pt in the saved figure), so measured text sizes
# and positions agree exactly.
NODE_FS, PW_FS, STAT_FS = 6.5, 7.5, 6.5      # font sizes (pt)
NODE_R = 3.6                                 # node marker radius (pt)
LAB_OFF, LAB_LINE = 7.5, 9.0                 # name starts LAB_OFF from the node centre; neighbours ≥ LAB_LINE apart
PW_LH = 9.5                                  # pathway label line height (pt)


def _text_size(text, size, weight="normal", family=None):
    """(width, height) of a single line of text in points."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath
    fp = FontProperties(size=size, weight=weight, family=family)
    ext = TextPath((0, 0), text, prop=fp).get_extents()
    return float(ext.x1) + 1.5, size * 1.2


def _seg_hits(seg, boxes, pad):
    (x1, y1), (x2, y2) = seg
    n = max(2, int(math.ceil(math.hypot(x2 - x1, y2 - y1) / 1.5)))
    for k in range(n + 1):
        x, y = x1 + (x2 - x1) * k / n, y1 + (y2 - y1) * k / n
        for b in boxes:
            if b[0] - pad < x < b[1] + pad and b[2] - pad < y < b[3] + pad:
                return True
    return False


def place_pathway_labels(items, r_arc, gap_x=6.0, gap_y=3.5, step=1.5):
    """items: dicts with 'a' (angle, y-up), 'w', 'h' (points). Adds 'box' = (x0, x1, y0, y1),
    'right' and 'lead' (leader segment or None). Each box starts just outside the arcs on the side
    facing away from the centre (so it cannot reach the ring); labels nearest the horizontal axis are
    placed first and a label that would touch an earlier one moves further up / down (away from the
    ring), with a leader line to its arc."""
    r0 = r_arc + 8.0
    placed, leaders = [], []

    def hit(b):
        return any(b[0] < p[1] + gap_x and b[1] > p[0] - gap_x and b[2] < p[3] + gap_y and b[3] > p[2] - gap_y
                   for p in placed)

    order = sorted(range(len(items)), key=lambda k: (abs(math.sin(items[k]["a"])), k))
    for k in order:
        o = items[k]
        c, s = math.cos(o["a"]), math.sin(o["a"])
        right = c >= 0
        ax_, ay_ = r0 * c, r0 * s
        centred = abs(s) < 0.2
        up = s >= 0 if not centred else False
        direction = 1 if s >= 0 else -1      # y-up: top half moves up

        def mk(y0):
            return (ax_ if right else ax_ - o["w"], ax_ + o["w"] if right else ax_, y0, y0 + o["h"])

        y00 = ay_ - o["h"] / 2 if centred else (ay_ if up else ay_ - o["h"])
        b, lead = mk(y00), None
        for _ in range(4000):
            moved = abs(b[2] - y00) > 0.5
            lead = (((r_arc + 2.5) * c, (r_arc + 2.5) * s),
                    ((b[0] - 2) if right else (b[1] + 2), (b[2] - 1.5) if direction > 0 else (b[3] + 1.5))) if moved else None
            if not hit(b) and not (lead and _seg_hits(lead, placed, 2)) and not any(_seg_hits(l, [b], 2) for l in leaders):
                break
            b = mk(b[2] + direction * step)
        o.update(box=b, right=right, lead=lead)
        placed.append(b)
        if lead:
            leaders.append(lead)
    return items


# ---------------------------------------------------------------------------
# network
# ---------------------------------------------------------------------------
def network_geometry(result, top: int = 10, metric: str = "p", cutoff: float = 0.05, show_labels: bool | None = None):
    """Everything the network figure draws, in points: node positions, ring radius, name extents,
    arcs and pathway label boxes. Separate from drawing so the no-overlap rules can be tested."""
    P = result.prep
    group_of, order, sel = network_groups(result, top, metric, cutoff)
    pos, arcs, step = _layout(group_of, order)
    shown = set(pos)
    labels = show_labels if show_labels is not None else len(shown) <= 200
    rn = max(150.0, (LAB_LINE / step - LAB_OFF) if step else 150.0)
    pos = {i: (a, rn * x, rn * y) for i, (a, x, y) in pos.items()}
    names = {i: _text_size(P.features[i], NODE_FS) for i in pos} if labels else {}
    max_w = max((w for w, _ in names.values()), default=0.0)
    r_arc = rn + (LAB_OFF + max_w + 8.0 if labels else 14.0)
    items = []
    for g, a0, a1 in arcs:
        lines = textwrap.wrap(g, 24, break_on_hyphens=False) or [g]
        stat = None
        if not sel.empty:
            row = sel[sel["pathway"] == g]
            if len(row):
                col = "q_value" if metric == "q" else "p_value"
                stat = f"{'FDR q-value' if metric == 'q' else 'p-value'} = {row.iloc[0][col]:.2g}"
        w = max([_text_size(ln, PW_FS, "bold")[0] for ln in lines] + ([_text_size(stat, STAT_FS)[0]] if stat else []))
        items.append({"pw": g, "a": (a0 + a1) / 2, "a0": a0, "a1": a1, "lines": lines, "stat": stat,
                      "w": w, "h": (len(lines) + (1 if stat else 0)) * PW_LH + 1.5})
    place_pathway_labels(items, r_arc)
    return {"group_of": group_of, "order": order, "sel": sel, "pos": pos, "rn": rn, "r_arc": r_arc,
            "labels": labels, "names": names, "items": items}


def name_polygon(geom, i):
    """Corners of metabolite i's (radial) name box, in points."""
    ang = geom["pos"][i][0]
    w, h = geom["names"][i]
    r1, r2 = geom["rn"] + LAB_OFF, geom["rn"] + LAB_OFF + w
    ux, uy = math.cos(ang), math.sin(ang)
    nx, ny = -uy, ux
    return [(r * ux + t * nx * h / 2, r * uy + t * ny * h / 2) for r, t in ((r1, -1), (r2, -1), (r2, 1), (r1, 1))]


def network_figure(result, top: int = 10, metric: str = "p", cutoff: float = 0.05, *, show_intact: bool = False,
                   show_reactions: bool = True, show_labels: bool | None = None, title: str | None = None,
                   size: float | None = None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Arc, Circle, PathPatch
    from matplotlib.path import Path

    P = result.prep
    geo = network_geometry(result, top, metric, cutoff, show_labels)
    pos, rn, r_arc, sel, order = geo["pos"], geo["rn"], geo["r_arc"], geo["sel"], geo["order"]
    shown = set(pos)

    # extent of everything drawn, then a figure where 1 data unit = 1 pt
    xs, ys = [-r_arc - 3, r_arc + 3], [-r_arc - 3, r_arc + 3]
    for it in geo["items"]:
        x0, x1, y0, y1 = it["box"]
        xs += [x0, x1]
        ys += [y0, y1]
    pad = 10.0
    x_lo, x_hi, y_lo, y_hi = min(xs) - pad, max(xs) + pad, min(ys) - pad, max(ys) + pad
    wpt, hpt = x_hi - x_lo, y_hi - y_lo
    fig = plt.figure(figsize=(wpt / 72.0, hpt / 72.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_aspect("equal")
    ax.axis("off")

    def curve(i, j):
        _, x1, y1 = pos[i]
        _, x2, y2 = pos[j]
        cx, cy = (x1 + x2) * 0.18, (y1 + y2) * 0.18
        return Path([(x1, y1), (cx, cy), (x2, y2)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])

    pairs = result.pairs
    if show_reactions and "rxn_steps" in pairs:
        for r in pairs[(pairs["rxn_steps"] > 0)].itertuples():
            if r.i in shown and r.j in shown:
                ax.add_patch(PathPatch(curve(r.i, r.j), fill=False, lw=0.6, color="#9aa4a6",
                                       ls="solid" if r.rxn_steps == 1 else (0, (1, 2)), zorder=1))
    order_cls = (["intact"] if show_intact else []) + ["lost", "gained", "flipped"]
    for c in order_cls:
        sub = pairs[pairs["class"] == c]
        for r in sub.itertuples():
            if r.i in shown and r.j in shown:
                lw = 0.6 if c == "intact" else 0.8 + min(1.4, abs(r.delta_r)) * 1.3
                ax.add_patch(PathPatch(curve(r.i, r.j), fill=False, lw=lw, color=EDGE[c], ls=DASH[c], zorder=2,
                                       capstyle="round"))
    lfc = result.features["log2fc"].values
    for i, (ang, x, y) in pos.items():
        ax.add_patch(Circle((x, y), NODE_R, facecolor=_div_colour(lfc[i]), edgecolor="white", lw=0.8, zorder=3))
        if geo["labels"]:
            rot = math.degrees(ang)
            right = math.cos(ang) >= 0
            ax.text((rn + LAB_OFF) * math.cos(ang), (rn + LAB_OFF) * math.sin(ang), P.features[i], fontsize=NODE_FS,
                    color=INK2, rotation=rot if right else rot + 180, rotation_mode="anchor",
                    ha="left" if right else "right", va="center_baseline", zorder=4)
    for it in geo["items"]:
        a0, a1 = it["a0"], it["a1"]
        ax.add_patch(Arc((0, 0), 2 * r_arc, 2 * r_arc, theta1=math.degrees(a1 + 0.004), theta2=math.degrees(a0 - 0.004),
                         lw=2.2, color=LINE, zorder=1))
        if it["lead"]:
            (lx1, ly1), (lx2, ly2) = it["lead"]
            ax.plot([lx1, lx2], [ly1, ly2], color=LINE, lw=0.6, zorder=1, solid_capstyle="butt")
        x0, x1, y0, y1 = it["box"]
        tx = x0 if it["right"] else x1
        rows = [(ln, PW_FS, "bold", INK) for ln in it["lines"]] + ([(it["stat"], STAT_FS, "normal", INK2)] if it["stat"] else [])
        for k, (ln, fs, fw, col) in enumerate(rows):
            ax.text(tx, y1 - 0.75 - k * PW_LH, ln, fontsize=fs, fontweight=fw, color=col,
                    ha="left" if it["right"] else "right", va="top", zorder=4)
    if not pos:
        ax.text(0, 0, "No pathway passes the current cutoff", ha="center", va="center", fontsize=9, color=INK)

    # legend, colour bar and caption sit below the network (outside it), title above
    ga, gb = P.group_labels
    handles = [Line2D([0], [0], color=EDGE[c], lw=2, ls=DASH[c]) for c in ("lost", "gained", "flipped")]
    names = ["Coupling lost", "Coupling gained", "Sign flipped"]
    if show_intact:
        handles.append(Line2D([0], [0], color=EDGE["intact"], lw=1.2))
        names.append("Intact")
    if show_reactions:
        handles += [Line2D([0], [0], color="#9aa4a6", lw=0.8), Line2D([0], [0], color="#9aa4a6", lw=0.8, ls=(0, (1, 2)))]
        names += ["Reaction (direct)", "Reaction (2 steps)"]
    ax.legend(handles, names, loc="upper left", bbox_to_anchor=(0.0, 0.0), frameon=False, fontsize=7.5, ncol=3,
              borderaxespad=0.3)
    import matplotlib.cm as cm
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    sm = cm.ScalarMappable(norm=Normalize(-1.5, 1.5), cmap=LinearSegmentedColormap.from_list("div", list(DIV)))
    cw = min(0.3, 150.0 / wpt)
    cax = fig.add_axes([1 - cw - 0.01, -34.0 / hpt, cw, 7.0 / hpt])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label(f"log2 fold change ({gb} / {ga})", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    sel_txt = ("all pathways" if not top or sel.empty else
               f"top {len(sel)} {result.prep.pathway_library} pathways by {'FDR q-value' if metric == 'q' else 'p-value'}"
               f" ≤ {cutoff:g}")
    ax.set_title(title or f"Correlation rewiring network: {ga} vs {gb} ({sel_txt})", fontsize=11, color=INK, pad=6)
    absorbed = [g for g in sel["pathway"]] if not sel.empty else []
    absorbed = [g for g in absorbed if g not in order]
    if absorbed:
        ax.text(0.0, -70.0 / hpt, textwrap.fill(
            f"Also selected; all their measured metabolites are drawn under a higher-ranked pathway: {'; '.join(absorbed)}",
            max(60, int(wpt / 3.6))), transform=ax.transAxes, fontsize=6.5, color=INK2, ha="left", va="top")
    return fig


# ---------------------------------------------------------------------------
# pair correlation plot
# ---------------------------------------------------------------------------
def pair_figure(result, i: int, j: int, *, size=(6.6, 4.4)):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    P = result.prep
    xi, xj = P.X[i], P.X[j]
    fig, ax = plt.subplots(figsize=size)
    ga, gb = P.group_labels
    stats = []
    for idx, col, lab in ((P.idx_a, GROUP[0], ga), (P.idx_b, GROUP[1], gb)):
        x, y = xi[idx], xj[idx]
        r = float(np.corrcoef(x, y)[0, 1]) if len(x) > 2 and x.std() > 0 and y.std() > 0 else float("nan")
        ax.scatter(x, y, s=26, color=col, edgecolors="white", linewidths=0.6, zorder=3, label=f"{lab} (r = {r:.2f})")
        if len(x) > 1 and x.std() > 0:
            b, a = np.polyfit(x, y, 1)
            xs = np.array([x.min(), x.max()])
            ax.plot(xs, a + b * xs, color=col, lw=2, zorder=2)
        stats.append(r)
    pr = result.pairs
    row = pr[(pr["i"] == min(i, j)) & (pr["j"] == max(i, j))]
    sub = ""
    if len(row):
        r0 = row.iloc[0]
        cls = r0["class"] if isinstance(r0["class"], str) else "not rewired"
        sub = f"Δr = {r0['delta_r']:+.2f} · FDR q-value = {r0['q_value']:.2g} · {cls}"
        if r0.get("rxn_steps", 0):
            sub += " · " + ("direct reaction" if r0["rxn_steps"] == 1 else f"2 steps via {r0['rxn_via']}")
    # full names, wrapped when long, so axis titles never run past the plot
    ax.set_xlabel(textwrap.fill(display_label(P.features[i]), 55, break_on_hyphens=False), fontsize=9, color=INK)
    ax.set_ylabel(textwrap.fill(display_label(P.features[j]), 45, break_on_hyphens=False), fontsize=9, color=INK)
    ax.set_title(textwrap.fill(f"{P.features[i]} ↔ {P.features[j]}", 70, break_on_hyphens=False), fontsize=10.5,
                 color=INK, loc="left", pad=16 if sub else 6)
    if sub:
        ax.text(0, 1.015, sub, transform=ax.transAxes, fontsize=7.5, color=INK2, va="bottom")
    ax.grid(color="#e3e7e5", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=8)
    # legend and note beside the plot, never on top of the points
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    ax.text(1.03, 0.72, textwrap.fill(f"Relative abundance on the analysed {P.scale_label} scale", 26),
            transform=ax.transAxes, fontsize=6.5, color=INK2, ha="left", va="top")
    fig.tight_layout()
    return fig


def figure_bytes(fig, fmt_label: str, dpi: int = 300) -> bytes:
    import matplotlib.pyplot as plt
    fmt = FORMATS[fmt_label][0]
    buf = io.BytesIO()
    kw = {"format": fmt, "bbox_inches": "tight", "facecolor": "white"}
    if fmt in ("png", "jpeg", "tiff"):
        kw["dpi"] = dpi
    if fmt == "jpeg":
        kw["pil_kwargs"] = {"quality": 95}
    if fmt == "tiff":
        kw["pil_kwargs"] = {"compression": "tiff_lzw"}
    import matplotlib
    # editable text in SVG (real <text>, not outlined glyphs) and TrueType fonts in PDF/EPS,
    # so labels stay editable in Illustrator / Inkscape
    with matplotlib.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42}):
        fig.savefig(buf, **kw)
    plt.close(fig)
    return buf.getvalue()
