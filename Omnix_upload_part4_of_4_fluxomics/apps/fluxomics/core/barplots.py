"""
core.barplots — dynamic bar plots (v2.2, replaces the box plots).

One plot per metabolite: x = isotopologue (M+0 … M+n); within each isotopologue one bar per group × time point
(group hue, lighter = earlier time). Bar height = mean of the biological replicates, error bar = SD or SEM,
dots = individual replicates. Isotopologues are never summed.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .style import INK, INK2, GRID, group_colors

FONT = dict(family="Inter, Segoe UI, Helvetica, Arial, sans-serif", size=12, color=INK)
ERRORS = {"SD": "SD", "SEM": "SEM", "none": "No error bar"}


def condition_colors(arms: list, times: list) -> dict:
    """Group hue, lighter → darker with time."""
    import matplotlib.colors as mc
    gc = group_colors(arms)
    out = {}
    for g in arms:
        base = np.array(mc.to_rgb(gc[g]))
        for k, h in enumerate(times):
            f = 0.55 - 0.55 * (k / max(len(times) - 1, 1))          # 0.55 (light) … 0 (full colour)
            out[(g, h)] = mc.to_hex(base + (1 - base) * f)
    return out


def bar_summary(X: pd.DataFrame, arms: list, times: list, isos: list) -> pd.DataFrame:
    """Mean / SD / SEM / n per isotopologue × group × time (X: Isotopologue, Arm, Hours, TimeLabel, Value)."""
    rows = []
    for g in arms:
        for h in times:
            s = X[(X.Arm == g) & np.isclose(X.Hours, h)]
            if s.empty:
                continue
            for iso in isos:
                v = s[s.Isotopologue == iso].Value.dropna().to_numpy(float)
                n = len(v)
                rows.append(dict(Isotopologue=iso, Group=g, Hours=h, Time=s.TimeLabel.iloc[0], n=n,
                                 Mean=float(v.mean()) if n else np.nan,
                                 SD=float(v.std(ddof=1)) if n > 1 else np.nan,
                                 SEM=float(v.std(ddof=1) / np.sqrt(n)) if n > 1 else np.nan))
    return pd.DataFrame(rows)


def _conds(X, arms, times):
    return [(g, h) for g in arms for h in times if len(X[(X.Arm == g) & np.isclose(X.Hours, h)])]


def iso_bar_plotly(X: pd.DataFrame, arms: list, times: list, isos: list, ylab: str, title: str = "",
                   points: bool = True, error: str = "SD", height: int = 500, xlab: str = "Isotopologue"):
    import plotly.graph_objects as go
    cc = condition_colors(arms, times)
    conds = _conds(X, arms, times)
    T = bar_summary(X, arms, times, isos)
    k_n = max(len(conds), 1)
    w = 0.84 / k_n
    xi = {iso: j for j, iso in enumerate(isos)}
    rng = np.random.default_rng(0)
    fig = go.Figure()
    for k, (g, h) in enumerate(conds):
        off = (k - (k_n - 1) / 2) * w
        d = T[(T.Group == g) & np.isclose(T.Hours, h)]
        lab = f"{g} · {d.Time.iloc[0]}"
        err = d[error].to_numpy(float) if error in ("SD", "SEM") else None
        fig.add_trace(go.Bar(x=[xi[i] + off for i in d.Isotopologue], y=d.Mean, width=w * 0.92, name=lab,
                             legendgroup=lab, marker=dict(color=cc[(g, h)], line=dict(width=0)),
                             error_y=dict(type="data", array=np.nan_to_num(err), visible=True, thickness=1.1,
                                          width=2.5, color="#3a3940") if err is not None else None,
                             customdata=np.column_stack([
                                 d.Isotopologue.astype(str).to_numpy(), d.n.to_numpy(),
                                 np.nan_to_num(d[error].to_numpy(float)) if err is not None else np.zeros(len(d))]),
                             hovertemplate="%{customdata[0]}<br>mean %{y:.3g}"
                                           + (f" ± %{{customdata[2]:.2g}} ({error})" if err is not None else "")
                                           + "<br>n = %{customdata[1]}<extra>" + lab + "</extra>"))
        if points:
            s = X[(X.Arm == g) & np.isclose(X.Hours, h) & X.Isotopologue.isin(isos)].dropna(subset=["Value"])
            xs = np.array([xi[i] + off for i in s.Isotopologue]) + rng.uniform(-w * 0.22, w * 0.22, len(s))
            fig.add_trace(go.Scatter(x=xs, y=s.Value, mode="markers", legendgroup=lab, showlegend=False,
                                     marker=dict(size=5, color="white", line=dict(color="#2b2a30", width=0.9)),
                                     customdata=np.c_[s.Sample, s.Isotopologue],
                                     hovertemplate="%{customdata[0]}<br>%{customdata[1]}: %{y:.3g}<extra></extra>"))
    fig.update_layout(template="plotly_white", font=FONT, height=height, barmode="overlay", bargap=0,
                      title=dict(text=title, x=0, font=dict(size=15)), margin=dict(l=60, r=20, t=60, b=50),
                      legend=dict(orientation="h", y=-0.14, x=0, font=dict(size=11)),
                      yaxis_title=ylab, xaxis_title=xlab)
    fig.update_xaxes(tickmode="array", tickvals=list(range(len(isos))), ticktext=isos,
                     range=[-0.6, len(isos) - 0.4], showgrid=False)
    fig.update_yaxes(gridcolor=GRID, rangemode="tozero")
    return fig


def default_size(n_iso: int, n_cond: int) -> tuple:
    """Export size (in) that gives every bar ≈ 0.15 in, plus room for the axis and the legend column."""
    axes_w = max(2.4, n_iso * (n_cond * 0.15 + 0.3))
    return round(min(axes_w + 2.9, 16.0), 1), 3.8


def iso_bar_mpl(X: pd.DataFrame, arms: list, times: list, isos: list, ylab: str, title: str = "", points=True,
                error: str = "SD", width=None, height=None, xlab: str = ""):
    """
    Publication bar plot with an explicit layout (v2.5.6): title (metabolite, bold) and subtitle centred over the
    plot, legend in its own column on the right, y-axis labels measured — nothing overlaps or squeezes the axes.
    """
    import textwrap
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    cc = condition_colors(arms, times)
    conds = _conds(X, arms, times)
    T = bar_summary(X, arms, times, isos)
    k_n = max(len(conds), 1)
    dw, dh = default_size(len(isos), k_n)
    W, H = float(width or dw), float(height or dh)
    w = 0.84 / k_n
    xi = {iso: j for j, iso in enumerate(isos)}
    rng = np.random.default_rng(0)
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0.1, 0.1, 0.6, 0.6])
    labels = []
    for k, (g, h) in enumerate(conds):
        off = (k - (k_n - 1) / 2) * w
        d = T[(T.Group == g) & np.isclose(T.Hours, h)]
        pos = np.array([xi[i] + off for i in d.Isotopologue])
        err = d[error].to_numpy(float) if error in ("SD", "SEM") else None
        ax.bar(pos, d.Mean, width=w * 0.9, color=cc[(g, h)], zorder=2, lw=0,
               yerr=None if err is None else np.nan_to_num(err),
               error_kw=dict(lw=0.7, capsize=min(2.0, 60 * w / max(len(isos), 1)), capthick=0.7, ecolor="black"))
        labels.append((f"{g} · {d.Time.iloc[0]}", cc[(g, h)]))
        if points:
            s = X[(X.Arm == g) & np.isclose(X.Hours, h) & X.Isotopologue.isin(isos)].dropna(subset=["Value"])
            xs = np.array([xi[i] + off for i in s.Isotopologue]) + rng.uniform(-w * 0.2, w * 0.2, len(s))
            ax.scatter(xs, s.Value, s=7, facecolor="white", edgecolor="black", lw=0.45, zorder=3)
    ax.set_xticks(range(len(isos))); ax.set_xticklabels(isos, fontsize=8, color="black")
    ax.set_xlim(-0.55, len(isos) - 0.45)
    ax.set_ylim(bottom=0)
    ax.set_ylabel(ylab, fontsize=9, color="black"); ax.tick_params(labelsize=8, colors="black", length=2.5)
    if xlab:
        ax.set_xlabel(xlab, fontsize=9, color="black")
    ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color("black"); ax.spines[sp].set_linewidth(0.7)
    leg = fig.legend(handles=[Patch(fc=c, ec="none", label=l) for l, c in labels], fontsize=7.5, frameon=False,
                     loc="upper left", handlelength=1.1, handleheight=0.9, labelspacing=0.35, borderaxespad=0,
                     labelcolor="black")
    cap = {"SD": "mean ± SD", "SEM": "mean ± SEM"}.get(error, "mean") + ("; dots = replicates" if points else "")
    head, _, rest = title.partition(" — ")
    sub = f"{rest}  ({cap})" if rest else f"({cap})"
    max_ch = max(int(W * 72 / (9.5 * 0.56)), 20)
    t1 = fig.text(0.5, 1, "\n".join(textwrap.wrap(head, max_ch)), ha="center", va="top", fontsize=10.5,
                  fontweight="bold", color="black")
    t2 = fig.text(0.5, 1, "\n".join(textwrap.wrap(sub, int(max_ch * 1.2))), ha="center", va="top", fontsize=8.5,
                  color="black")
    # ---- measure the text, then place everything in inches
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    dpi = fig.dpi
    inch = lambda bb: (bb.width / dpi, bb.height / dpi)
    lw_, lh_ = inch(leg.get_window_extent(r))
    t1h, t2h = inch(t1.get_window_extent(r))[1], inch(t2.get_window_extent(r))[1]
    for _ in range(2):
        ab = ax.get_window_extent(r); tb = ax.get_tightbbox(r)
        left_dec = (ab.x0 - tb.x0) / dpi
        bot_dec = (ab.y0 - tb.y0) / dpi
        top_band = 0.1 + t1h + 0.04 + t2h + 0.18
        x0 = 0.08 + left_dec
        x1 = W - (0.15 + lw_ + 0.2)
        y0 = 0.08 + bot_dec
        y1 = H - top_band
        ax.set_position([x0 / W, y0 / H, max(x1 - x0, 0.8) / W, max(y1 - y0, 0.8) / H])
        fig.canvas.draw()
    leg.set_bbox_to_anchor((x1 / W + 0.15 / W, y1 / H), transform=fig.transFigure)
    xc = (x0 + x1) / 2 / W
    t1.set_position((xc, 1 - 0.1 / H)); t2.set_position((xc, 1 - (0.1 + t1h + 0.04) / H))
    return fig
