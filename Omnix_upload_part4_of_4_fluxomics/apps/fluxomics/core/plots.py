"""
core.plots — figures. Every figure exists as an interactive Plotly version (for the app) and, where it is meant
for publication, as a matplotlib version (export to PNG / TIFF / PDF / SVG / JPEG at any size and DPI).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .style import (INK, INK2, GRID, MUTED, M0_GRAY, SEQUENTIAL, DIVERGING, group_colors, iso_colors)
from . import pathways as pw

LINESTYLES = ["-", "--", ":", "-."]
DASHES = ["solid", "dash", "dot", "dashdot"]
MARKERS = ["o", "s", "^", "D"]
PLOTLY_FONT = dict(family="Inter, Segoe UI, Helvetica, Arial, sans-serif", size=12, color=INK)


def _plotly_layout(fig, title="", height=None):
    fig.update_layout(template="plotly_white", font=PLOTLY_FONT, title=dict(text=title, x=0, font=dict(size=15)),
                      margin=dict(l=60, r=20, t=60 if title else 30, b=50), legend=dict(orientation="h", y=-0.15),
                      height=height, hoverlabel=dict(font_size=12))
    return fig


def seq_hex(v, vmax=100.0):
    import matplotlib.colors as mc
    if v is None or not np.isfinite(v):
        return "#ffffff"
    return mc.to_hex(SEQUENTIAL(float(min(max(v / vmax, 0.0), 1.0))))


def plotly_seq_scale():
    return [[t, seq_hex(t * 100)] for t in np.linspace(0, 1, 9)]


def plotly_div_scale():
    import matplotlib.colors as mc
    return [[t, mc.to_hex(DIVERGING(t))] for t in np.linspace(0, 1, 11)]


# ---------------------------------------------------------------------------
# Isotopologue distribution
# ---------------------------------------------------------------------------
def mid_table(F: pd.DataFrame, D: pd.DataFrame, compound: str) -> pd.DataFrame:
    a = D[D.Include.astype(bool) & (D.SampleType == "13C labeled")][["Sample", "Group", "Hours", "Condition"]]
    X = F[F.Compound == compound].merge(a, on="Sample")
    return X.groupby(["Condition", "Group", "Hours", "i"], sort=False).Fraction.agg(["mean", "std", "count"]).reset_index()


def mid_bars_plotly(F, D, compound, groups, as_percent=True, show_m0=True):
    T = mid_table(F, D, compound)
    if T.empty:
        return go.Figure()
    n = int(T.i.max())
    cols = iso_colors(n)
    conds = list(dict.fromkeys(T.sort_values(["Group", "Hours"], key=lambda s: s.map(
        {g: k for k, g in enumerate(groups)}) if s.name == "Group" else s).Condition))
    fig = go.Figure()
    k = 100 if as_percent else 1
    # x labels without overlap: two levels (time under each bar, group name centred under its bars) when every
    # group name fits under its bars; otherwise one vertical label per bar (never overlaps, however many samples)
    ci = T.drop_duplicates("Condition").set_index("Condition").reindex(conds)
    tl = [f"{h:g} h" if pd.notna(h) else str(c) for c, h in zip(conds, ci.Hours)]
    gl = [str(g) for g in ci.Group]
    n_per = pd.Series(gl).value_counts()
    two_level = all(len(g) * 7.5 <= n_per[g] * 42 for g in n_per.index) and len(set(zip(gl, tl))) == len(conds)
    xv = [gl, tl] if two_level else conds
    present = sorted(int(x) for x in T.i.unique())          # only isotopologues that exist (no empty layers)
    for i in [x for x in present if show_m0 or x > 0]:
        sub = T[T.i == i].set_index("Condition").reindex(conds)
        sub["mean"] = sub["mean"].fillna(0.0)
        fig.add_bar(x=xv, y=sub["mean"] * k, name=f"M+{i}", marker_color=cols[i], text=conds,
                    textposition="none",
                    error_y=None, customdata=np.c_[sub["std"] * k, sub["count"]],
                    hovertemplate="%{text}<br>M+" + str(i) + ": %{y:.2f}" + ("%" if as_percent else "") +
                                  " (SD %{customdata[0]:.2f}, n=%{customdata[1]})<extra></extra>")
    fig.update_layout(barmode="stack", yaxis_title="Fraction (%)" if as_percent else "Fraction",
                      bargap=0.25)
    if not show_m0:
        fig.update_layout(yaxis_title="Labeled fraction (%) — M+0 hidden" if as_percent else "Labeled fraction")
    max_t = max(len(t) for t in (tl if two_level else conds))
    fig = _plotly_layout(fig, f"{compound} — isotopologue distribution (group mean)", 430 + (0 if two_level else
                                                                                            int(max_t * 6.5)))
    fig.update_xaxes(automargin=True, tickangle=0 if two_level else -90, tickfont=dict(size=11),
                     showdividers=True, dividercolor="#bbbbbb")
    # legend in its own column on the right (top of the stack first), clear of the axis labels
    fig.update_layout(legend=dict(orientation="v", x=1.01, xanchor="left", y=1, yanchor="top",
                                  traceorder="reversed", title=dict(text="Isotopologue")),
                      margin=dict(r=110, b=40))
    return fig


def mid_bars_mpl(F, D, compound, groups, show_m0=True, width=6.5, height=3.6):
    """
    Publication MID bar chart, laid out in inches so no text overlaps at any size:
    title band (top) · plot · legend column (right) · time labels under each bar and the group name centred under a
    bracket (bottom).
    """
    T = mid_table(F, D, compound)
    n = int(T.i.max()) if len(T) else 0
    cols = iso_colors(n)
    order = {g: k for k, g in enumerate(groups)}
    T = T.assign(_g=T.Group.map(order).fillna(len(groups))).sort_values(["_g", "Hours"])
    conds = list(dict.fromkeys(T.Condition))
    cinfo = T.drop_duplicates("Condition").set_index("Condition")[["Group", "Hours"]].reindex(conds)
    present = sorted(int(x) for x in T.i.unique()) if len(T) else []
    shown = [x for x in present if show_m0 or x > 0]
    # ---- bar positions: groups separated by a gap
    xs, x, prev = [], 0.0, None
    for c in conds:
        g = cinfo.loc[c, "Group"]
        if prev is not None and g != prev:
            x += 0.6
        xs.append(x); x += 1.0; prev = g
    xs = np.array(xs)
    # ---- font sizes (pt) and band sizes (in)
    fs_t, fs_ax, fs_tick, fs_leg, fs_grp = 10.0, 8.5, 7.5, 7.5, 8.0
    ch = lambda fs: 0.56 * fs / 72                     # mean character width (in)
    tlabels = [f"{h:g} h" if pd.notna(h) else str(c) for c, h in zip(conds, cinfo.Hours)]
    n_leg = len(shown)
    leg_cols = 1 if n_leg <= 14 else 2
    leg_rows = int(np.ceil(n_leg / leg_cols)) if n_leg else 0
    W_leg = leg_cols * (0.22 + ch(fs_leg) * 5 + 0.12) + 0.15
    W_left = 0.18 + fs_ax / 72 * 1.4 + ch(fs_tick) * 3 + 0.12
    W_right = 0.12 + W_leg
    H_title = 0.16 + fs_t / 72 * 1.5
    W_ax = max(width - W_left - W_right, 1.2)
    slot = W_ax / max(xs[-1] + 1 if len(xs) else 1, 1)
    # group names centred under brackets only when each fits under its own bars; otherwise one vertical
    # "group · time" label per bar (never overlaps, however many samples)
    gcount = cinfo.Group.value_counts()
    brackets = all(len(str(g)) * ch(fs_grp) * 1.08 <= gcount[g] * slot - 0.04 for g in gcount.index)
    if not brackets:
        tlabels = [str(c) for c in conds]
    rot = (not brackets) or max((len(t) for t in tlabels), default=0) * ch(fs_tick) > 0.92 * slot
    H_tick = (max((len(t) for t in tlabels), default=1) * ch(fs_tick) + 0.1) if rot else fs_tick / 72 * 1.6
    H_bottom = 0.06 + H_tick + ((0.10 + fs_grp / 72 * 1.7 + 0.08) if brackets else 0.08)
    H_ax = max(height - H_title - H_bottom, 1.0)
    Wt, Ht = W_left + W_ax + W_right, H_title + H_ax + H_bottom
    fig = plt.figure(figsize=(Wt, Ht))
    ax = fig.add_axes([W_left / Wt, H_bottom / Ht, W_ax / Wt, H_ax / Ht])
    bottom = np.zeros(len(conds))
    for i in shown:
        v = T[T.i == i].set_index("Condition").reindex(conds)["mean"].fillna(0).to_numpy() * 100
        ax.bar(xs, v, bottom=bottom, color=cols[i], width=0.74, edgecolor="white", lw=0.5, label=f"M+{i}")
        bottom += v
    ax.set_xlim(xs.min() - 0.6 if len(xs) else -0.6, xs.max() + 0.6 if len(xs) else 0.6)
    ax.set_ylim(0, max(100.0 if show_m0 else float(np.nanmax(bottom)) * 1.08 if len(bottom) else 1.0, 1e-6))
    ax.set_xticks(xs)
    ax.set_xticklabels(tlabels, fontsize=fs_tick, rotation=90 if rot else 0, ha="center", color="black")
    ax.tick_params(axis="x", length=0, pad=3)
    ax.tick_params(axis="y", labelsize=fs_tick, length=2.5, colors="black")
    ax.set_ylabel("Fraction (%)" if show_m0 else "Labeled fraction (%)", fontsize=fs_ax, color="black")
    ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    # ---- group brackets and names (axes-x / figure-inch y)
    trans = ax.get_xaxis_transform()
    yb = -(0.06 + H_tick + 0.05) / H_ax
    for g in (dict.fromkeys(cinfo.Group) if brackets else []):
        gx = xs[(cinfo.Group == g).to_numpy()]
        ax.plot([gx.min() - 0.37, gx.max() + 0.37], [yb, yb], transform=trans, color="black", lw=0.7,
                clip_on=False)
        ax.text((gx.min() + gx.max()) / 2, yb - 0.05 / H_ax, str(g), transform=trans, ha="center", va="top",
                fontsize=fs_grp, fontweight="bold", color="black", clip_on=False)
    # ---- title band and legend column
    fig.text(W_left / Wt, 1 - 0.1 / Ht, compound, fontsize=fs_t, fontweight="bold", ha="left", va="top",
             color="black")
    if n_leg:
        h, l = ax.get_legend_handles_labels()
        h, l = h[::-1], l[::-1]                       # top of the stack first
        fig.legend(h, l, loc="upper left", bbox_to_anchor=((W_left + W_ax + 0.15) / Wt, (H_bottom + H_ax) / Ht),
                   ncol=leg_cols, frameon=False, fontsize=fs_leg, handlelength=1.0, handleheight=0.9,
                   borderaxespad=0, labelspacing=0.35, columnspacing=0.8, title="Isotopologue",
                   title_fontsize=fs_leg, alignment="left", labelcolor="black")
    return fig


# ---------------------------------------------------------------------------
# Statistics figures
# ---------------------------------------------------------------------------
def volcano_plotly(T, x="log2FC", p="p", q="FDR (BH)", label="Feature", fdr=0.05, fc=1.0, title="", xlab=None):
    T = T.dropna(subset=[x, p]).copy()
    T["-log10 p"] = -np.log10(T[p].clip(lower=1e-300))
    sig_up = (T[q] < fdr) & (T[x] >= fc if fc else T[x] > 0)
    sig_dn = (T[q] < fdr) & (T[x] <= -fc if fc else T[x] < 0)
    T["cls"] = np.where(sig_up, "Higher in B", np.where(sig_dn, "Lower in B", "Not significant"))
    colors = {"Higher in B": "#eb6834", "Lower in B": "#2a78d6", "Not significant": "#c4c3cb"}
    fig = go.Figure()
    for c in ["Not significant", "Lower in B", "Higher in B"]:
        s = T[T.cls == c]
        fig.add_scatter(x=s[x], y=s["-log10 p"], mode="markers", name=c, marker=dict(color=colors[c], size=8,
                        line=dict(color="white", width=0.6)), text=s[label],
                        customdata=np.c_[s[q]], hovertemplate="%{text}<br>effect %{x:.3f}<br>p %{y:.2f} (−log10)"
                                                              "<br>FDR %{customdata[0]:.3g}<extra></extra>")
    if fc:
        for v in (-fc, fc):
            fig.add_vline(x=v, line=dict(color="#bbb", dash="dot", width=1))
    fig.update_layout(xaxis_title=xlab or x, yaxis_title="−log10 p")
    return _plotly_layout(fig, title, 480)


def volcano_mpl(T, x="log2FC", p="p", q="FDR (BH)", label="Feature", fdr=0.05, fc=1.0, title="", xlab=None,
                n_labels=12, width=5.5, height=4.5):
    T = T.dropna(subset=[x, p]).copy()
    y = -np.log10(T[p].clip(lower=1e-300))
    up = (T[q] < fdr) & (T[x] >= fc if fc else T[x] > 0)
    dn = (T[q] < fdr) & (T[x] <= -fc if fc else T[x] < 0)
    fig, ax = plt.subplots(figsize=(width, height))
    ax.scatter(T[x][~(up | dn)], y[~(up | dn)], s=12, color="#c4c3cb", lw=0)
    ax.scatter(T[x][dn], y[dn], s=18, color="#2a78d6", lw=0)
    ax.scatter(T[x][up], y[up], s=18, color="#eb6834", lw=0)
    top = T.assign(_y=y)[(up | dn)].sort_values("_y", ascending=False).head(n_labels)
    for _, r in top.iterrows():
        ax.annotate(str(r[label])[:28], (r[x], r["_y"]), xytext=(3, 2), textcoords="offset points", fontsize=6.5)
    if fc:
        for v in (-fc, fc):
            ax.axvline(v, color="#bbb", ls=":", lw=0.8)
    ax.set_xlabel(xlab or x); ax.set_ylabel("−log10 p"); ax.set_title(title, loc="left")
    fig.tight_layout()
    return fig


def strip_plotly(values: pd.Series, D: pd.DataFrame, groups, title="", ylab=""):
    a = D[D.Include.astype(bool) & (D.SampleType == "13C labeled")]
    gc = group_colors(groups)
    fig = go.Figure()
    for g in groups:
        s = a[a.Group == g]
        fig.add_box(x=s.TimeLabel, y=values.reindex(s.Sample), name=g, marker_color=gc[g], boxpoints="all",
                    jitter=0.35, pointpos=0, line=dict(width=1), fillcolor="rgba(0,0,0,0)",
                    text=s.Sample, hovertemplate="%{text}<br>%{y:.3g}<extra>" + g + "</extra>")
    fig.update_layout(boxmode="group", yaxis_title=ylab, xaxis_title="Time")
    return _plotly_layout(fig, title, 380)
