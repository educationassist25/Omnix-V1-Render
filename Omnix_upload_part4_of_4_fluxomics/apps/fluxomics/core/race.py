"""
core.race — label arrival, all groups and time points side by side (static).

One figure, one panel per group × time point (columns ordered group by group, earliest time first), shared rows.
Each row is one isotopologue (M+1 and up) of one metabolite, grouped by pathway (bar colour). Bar = group-mean
isotopologue fraction (%) at that time, error bar = SEM or 95 % CI. Static figure (v2.3: no animation).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .style import INK, INK2, GRID, group_colors

FONT = dict(family="Inter, Segoe UI, Helvetica, Arial, sans-serif", size=12, color=INK)
# pathway colours deliberately avoid the group colours (blue / orange) used in the panel headers
PATH_COLORS = ["#1baf7a", "#8a5cc7", "#d4a300", "#e87ba4", "#2c8c8c", "#9c5a2e", "#7a7a86", "#b5a900",
               "#c94c4c", "#4a3aa7", "#5fb7d4", "#6b8e23"]


def race_table(SK: pd.DataFrame, members: dict, groups: list, times: list, label_map: dict | None = None,
               threshold: float = 2.0, max_per_pathway: int = 12, include_m0: bool = False) -> pd.DataFrame:
    """
    SK: summary per isotopologue channel — Compound ('Name | M+i'), Group, Hours, mean, sem, ci95.
    Returns long table: Key, Label, Pathway, Group, Hours, mean, sem, ci95 for the rows to show.
    """
    comp_pw = {}
    for p, cs in members.items():
        for c in cs:
            comp_pw.setdefault(c, p)
    X = SK[SK.Group.isin(groups) & SK.Hours.isin(times)].copy()
    X["Metabolite"] = X.Compound.str.rsplit(" | ", n=1).str[0]
    X["Iso"] = X.Compound.str.rsplit(" | ", n=1).str[1]
    if not include_m0:
        X = X[X.Iso != "M+0"]
    X = X[X.Metabolite.isin(comp_pw)]
    X["Pathway"] = X.Metabolite.map(comp_pw)
    peak = X.groupby("Compound")["mean"].max()
    keep = peak[peak >= threshold].index
    X = X[X.Compound.isin(keep)]
    # top N isotopologues per pathway (by their highest group-mean fraction)
    top = (X.groupby(["Pathway", "Compound"])["mean"].max().reset_index()
           .sort_values(["Pathway", "mean"], ascending=[True, False]).groupby("Pathway").head(max_per_pathway))
    X = X[X.Compound.isin(top.Compound)]
    lab = label_map or {}
    X["Label"] = [f"{lab.get(m, m)} {i}" for m, i in zip(X.Metabolite, X.Iso)]
    dup = X.drop_duplicates("Compound").Label.duplicated(keep=False)
    if dup.any():
        X["Label"] = [f"{m} {i}" for m, i in zip(X.Metabolite, X.Iso)]
    return X.rename(columns={"Compound": "Key"})[["Key", "Label", "Pathway", "Metabolite", "Iso", "Group", "Hours",
                                                  "mean", "sem", "ci95"]]


def _order(R: pd.DataFrame, times: list):
    last = R[R.Hours == max(times)].groupby(["Pathway", "Label"])["mean"].max().reset_index()
    allrows = R.drop_duplicates("Label")[["Pathway", "Label"]].merge(last, how="left").fillna({"mean": 0})
    allrows = allrows.sort_values(["Pathway", "mean"], ascending=[True, False])
    return list(allrows.Label), list(allrows.Pathway)


def _nice_ticks(vmax: float) -> tuple[float, list]:
    """Axis upper limit and 2–3 round ticks that stay clear of the panel edge (no clash with the next panel)."""
    vmax = max(float(vmax), 1e-6)
    raw = vmax / 2.5
    mag = 10 ** np.floor(np.log10(raw))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= raw), default=10 * mag)
    ticks = list(np.arange(0, vmax * 1.0001, step))
    upper = max(vmax, ticks[-1]) * 1.18            # room after the last tick label
    return upper, [float(t) for t in ticks]


def _fmt(t: float) -> str:
    return f"{t:g}"


def race_grid_plotly(R: pd.DataFrame, groups: list, times: list, unc: str = "SEM", xlab: str = "Fraction (%)",
                     short_group: dict | None = None):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    labels, rpw = _order(R, times)
    pws = list(dict.fromkeys(rpw))
    pcol = {p: PATH_COLORS[i % len(PATH_COLORS)] for i, p in enumerate(pws)}
    cols = [(g, h) for g in groups for h in times]
    ec = "sem" if unc == "SEM" else "ci95"
    vmax = float(np.nanmax(R["mean"] + R[ec].fillna(0))) if len(R) else 1.0
    upper, ticks = _nice_ticks(vmax)
    sg = short_group or {}
    nT = len(times)
    # wider gap between groups than between time points
    widths = []
    for gi in range(len(groups)):
        widths += [1.0] * nT
    fig = make_subplots(rows=1, cols=len(cols), shared_yaxes=True, horizontal_spacing=0.012,
                        column_widths=widths)
    rowcol = [pcol[p] for p in rpw]
    for k, (g, h) in enumerate(cols):
        d = R[(R.Group == g) & np.isclose(R.Hours, h)].set_index("Label").reindex(labels)
        fig.add_trace(go.Bar(x=d["mean"], y=labels, orientation="h", marker=dict(color=rowcol, line=dict(width=0)),
                             error_x=dict(type="data", array=d[ec].fillna(0), thickness=0.9, width=1.5,
                                          color="#4a4952"),
                             customdata=np.c_[rpw, d["mean"].round(2).astype(str), d[ec].round(2).astype(str)],
                             hovertemplate="%{y}<br>" + f"{g} · {h:g} h" + "<br>%{customdata[1]} % ± "
                                           "%{customdata[2]} (" + unc + ")<br>%{customdata[0]}<extra></extra>",
                             showlegend=False), row=1, col=k + 1)
    for k in range(len(cols)):
        fig.update_xaxes(range=[0, upper], tickvals=ticks, ticktext=[_fmt(t) for t in ticks], row=1, col=k + 1,
                         gridcolor=GRID, tickfont=dict(size=9), showline=True, linecolor="#c9c8d0")
    fig.update_yaxes(autorange="reversed", tickfont=dict(size=10), row=1, col=1, automargin=True)
    gc = group_colors(groups)
    H = max(480, 18 * len(labels) + 250)
    px = 1.0 / H                                        # one pixel in paper units (vertical)
    for k, (g, h) in enumerate(cols):                   # time label above each panel
        dom = fig.layout[f"xaxis{k + 1 if k else ''}"].domain
        fig.add_annotation(x=(dom[0] + dom[1]) / 2, xref="paper", y=1 + 6 * px, yref="paper", yanchor="bottom",
                           showarrow=False, text=f"{h:g} h", font=dict(size=11, color=INK2))
    for gi, g in enumerate(groups):                     # group bar + name above its time panels
        a = fig.layout[f"xaxis{gi * nT + 1 if gi * nT else ''}"].domain[0]
        b = fig.layout[f"xaxis{(gi + 1) * nT}"].domain[1]
        fig.add_shape(type="rect", xref="paper", yref="paper", x0=a, x1=b, y0=1 + 28 * px, y1=1 + 33 * px,
                      fillcolor=gc.get(g, "#888"), line=dict(width=0))
        fig.add_annotation(x=(a + b) / 2, xref="paper", y=1 + 36 * px, yref="paper", yanchor="bottom",
                           showarrow=False, text=f"<b>{sg.get(g, g)}</b>", font=dict(size=13, color=gc.get(g, INK)))
    fig.add_annotation(x=0.5, xref="paper", y=-34 * px, yref="paper", yanchor="top", showarrow=False,
                       text=f"Mean isotopologue fraction (%) ± {unc}", font=dict(size=11, color=INK))
    for i in [i for i in range(1, len(rpw)) if rpw[i] != rpw[i - 1]]:
        fig.add_shape(type="line", xref="paper", yref="y", x0=0, x1=1, y0=i - 0.5, y1=i - 0.5,
                      line=dict(color="#cfcdd8", width=0.8, dash="dot"))
    for p in pws:
        fig.add_scatter(x=[None], y=[None], mode="markers", marker=dict(symbol="square", size=11, color=pcol[p]),
                        name=p, showlegend=True)
    fig.update_layout(
        template="plotly_white", font=FONT, height=H, bargap=0.22,
        title=dict(text="Label arrival — every group and time point", x=0, y=1 - 8 * px, yanchor="top",
                   font=dict(size=15)),
        margin=dict(l=10, r=20, t=120, b=95),
        legend=dict(orientation="h", y=-60 * px, yanchor="top", x=0, font=dict(size=11), title_text="Pathway: "))
    return fig


def race_grid_mpl(R: pd.DataFrame, groups: list, times: list, unc: str = "SEM", xlab: str = "Fraction (%)",
                  width: float = 12.0, height: float | None = None, short_group: dict | None = None,
                  font: float = 7.0):
    """
    Publication figure laid out in inches so no text can collide: label column measured from the longest label,
    fixed bands for title / group headers / time labels / axis label / legend, and ≥ 0.15 in per row.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    labels, rpw = _order(R, times)
    pws = list(dict.fromkeys(rpw))
    pcol = {p: PATH_COLORS[i % len(PATH_COLORS)] for i, p in enumerate(pws)}
    cols = [(g, h) for g in groups for h in times]
    ec = "sem" if unc == "SEM" else "ci95"
    vmax = float(np.nanmax(R["mean"] + R[ec].fillna(0))) if len(R) else 1.0
    upper, ticks = _nice_ticks(vmax)
    n = len(labels)
    row_in = max(0.15, font * 1.75 / 72)
    lab_in = max(1.0, max(len(x) for x in labels) * font * 0.56 / 72 + 0.15) if labels else 1.0
    T_TITLE, T_GROUP, T_TIME = 0.42, 0.34, 0.22
    B_AX, B_XLAB = 0.24, 0.26
    n_leg_rows = int(np.ceil(len(pws) / max(1, min(len(pws), 5))))
    B_LEG = 0.22 * n_leg_rows + 0.12
    need_h = T_TITLE + T_GROUP + T_TIME + n * row_in + B_AX + B_XLAB + B_LEG
    H = max(float(height or 0), need_h)
    row_in = (H - (need_h - n * row_in)) / max(n, 1)
    W = max(float(width), lab_in + 1.1 * len(cols) + 0.4)
    gap_t, gap_g = 0.08, 0.28
    panel_w = (W - lab_in - 0.25 - gap_t * (len(cols) - len(groups)) - gap_g * (len(groups) - 1)) / len(cols)
    fig = plt.figure(figsize=(W, H))
    fx, fy = (lambda x: x / W), (lambda y: 1 - y / H)            # y measured from the top
    top_axes = T_TITLE + T_GROUP + T_TIME
    hgt = n * row_in
    gc = group_colors(groups)
    x = lab_in
    pos = {}
    for gi, g in enumerate(groups):
        x0g = x
        for ti, h in enumerate(times):
            pos[(g, h)] = x
            x += panel_w + (gap_t if ti < len(times) - 1 else 0)
        x1g = x
        fig.patches.append(Rectangle((fx(x0g), fy(T_TITLE + T_GROUP * 0.62)), fx(x1g - x0g), 0.05 / H,
                                     transform=fig.transFigure, color=gc.get(g, "#888"), figure=fig))
        fig.text(fx((x0g + x1g) / 2), fy(T_TITLE + T_GROUP * 0.62) + 0.07 / H, (short_group or {}).get(g, g),
                 ha="center", va="bottom", fontsize=font + 2.5, fontweight="bold", color=gc.get(g, INK))
        x += gap_g
    y = np.arange(n)
    first = None
    for k, (g, h) in enumerate(cols):
        ax = fig.add_axes([fx(pos[(g, h)]), fy(top_axes + hgt), fx(panel_w), hgt / H], sharey=first)
        first = first or ax
        d = R[(R.Group == g) & np.isclose(R.Hours, h)].set_index("Label").reindex(labels)
        ax.barh(y, d["mean"].fillna(0), xerr=d[ec].fillna(0), color=[pcol[p] for p in rpw], height=0.74,
                error_kw=dict(lw=0.45, capsize=1.2, capthick=0.45, ecolor="#4a4952"), zorder=2)
        ax.set_xlim(0, upper); ax.set_ylim(n - 0.5, -0.5)
        ax.set_xticks(ticks); ax.set_xticklabels([_fmt(t) for t in ticks], fontsize=font - 0.5)
        ax.tick_params(axis="x", length=2, pad=1.5)
        ax.grid(axis="x", color=GRID, lw=0.4, zorder=0)
        for i in range(1, n):
            if rpw[i] != rpw[i - 1]:
                ax.axhline(i - 0.5, color="#bdbbc9", lw=0.4, ls=(0, (2, 2)), zorder=1)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_linewidth(0.5); ax.spines[sp].set_color("#9b9aa3")
        ax.text(0.5, 1 + (T_TIME * 0.35) / hgt, f"{h:g} h", transform=ax.transAxes, ha="center", va="bottom",
                fontsize=font + 0.5, color=INK2)
        if k == 0:
            ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=font)
            ax.tick_params(axis="y", length=0, pad=3)
        else:
            ax.tick_params(axis="y", left=False, labelleft=False)
    fig.text(fx(lab_in + (W - lab_in - 0.25) / 2), fy(top_axes + hgt + B_AX + B_XLAB * 0.5),
             f"Mean isotopologue fraction (%) ± {unc}", ha="center", va="center", fontsize=font + 1)
    handles = [Rectangle((0, 0), 1, 1, color=pcol[p]) for p in pws]
    fig.legend(handles, pws, loc="upper center", ncol=min(len(pws), 5), fontsize=font, frameon=False,
               bbox_to_anchor=(fx(lab_in + (W - lab_in - 0.25) / 2), fy(top_axes + hgt + B_AX + B_XLAB)),
               handlelength=1.0, columnspacing=1.4)
    fig.text(fx(0.08), fy(0.08), "Label arrival — every group and time point", ha="left", va="top",
             fontsize=font + 4, fontweight="bold")
    return fig
