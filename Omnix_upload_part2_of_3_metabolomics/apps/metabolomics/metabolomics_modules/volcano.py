"""
volcano.py - Publication-quality volcano plot module.

Supports customizable thresholds, colors, point styling, label modes (auto/manual/
significant-only) with optional repelling, legend placement/format, axis and title
customization, gridlines, threshold-line styling, metabolite highlighting, marker
shapes, background, figure-size presets, journal-style themes, multi-format/multi-DPI
export, and data/settings export for reproducibility.

Not implemented here (would require a different architecture, noted for transparency):
  - True interactive hover/click (would need Plotly/Bokeh instead of static matplotlib)
  - HMDB/KEGG/pathway-based labeling or filtering (no compound-database integration
    in this app's data model — labels/filters work by compound name only)
  - VIP-based filtering (would require a PLS-DA model, not implemented)
  - Freehand in-app annotation (arrows/rectangles/regions) — a fixed set of
    programmatic highlight/annotate options is provided instead
"""

import io
import json
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D

try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False

matplotlib.use("Agg")

# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------
COLORBLIND_PALETTES = {
    "Default (Red/Teal/Gray)": {"up": "#B2182B", "down": "#00A9A5", "ns": "#B0B0B0"},
    "Classic (Red/Blue/Gray)": {"up": "#D62728", "down": "#1F77B4", "ns": "#C7C7C7"},
    "Colorblind-safe (Okabe-Ito)": {"up": "#D55E00", "down": "#0072B2", "ns": "#999999"},
    "Colorblind-safe (Viridis ends)": {"up": "#FDE725", "down": "#440154", "ns": "#AFAFAF"},
    "Monochrome": {"up": "#000000", "down": "#666666", "ns": "#D9D9D9"},
}

MARKER_SHAPES = {"Circle": "o", "Triangle": "^", "Square": "s", "Diamond": "D", "Cross": "x"}
LINE_STYLES = {"Dashed": "--", "Solid": "-", "Dotted": ":", "None": "None"}

FIGURE_SIZE_PRESETS = {
    "Single column (4x4 in)": (4.0, 4.0),
    "Double column (8x6 in)": (8.0, 6.0),
    "Standard (6x5 in)": (6.0, 5.0),
    "Presentation (10x7.5 in)": (10.0, 7.5),
    "Poster (14x10 in)": (14.0, 10.0),
    "Custom": None,
}

# Approximate journal-style presets: font family, grid, spine, and palette choices.
# These are stylistic approximations, not exact reproductions of any journal's
# official figure specification.
THEMES = {
    "Default": {"font": "sans-serif", "grid": "none", "spines_top_right": False, "palette": "Default (Red/Teal/Gray)"},
    "Nature": {"font": "sans-serif", "grid": "none", "spines_top_right": False, "palette": "Classic (Red/Blue/Gray)"},
    "Cell": {"font": "sans-serif", "grid": "none", "spines_top_right": False, "palette": "Default (Red/Teal/Gray)"},
    "Cancer Research": {"font": "serif", "grid": "none", "spines_top_right": False, "palette": "Classic (Red/Blue/Gray)"},
    "Clinical Cancer Research": {"font": "serif", "grid": "none", "spines_top_right": False, "palette": "Classic (Red/Blue/Gray)"},
    "PNAS": {"font": "sans-serif", "grid": "major", "spines_top_right": True, "palette": "Colorblind-safe (Okabe-Ito)"},
}

LEGEND_POSITIONS = {
    "Right": "right", "Left": "left", "Top": "top", "Bottom": "bottom", "Hidden": "none",
    "Right - top": "right", "Right - middle": "right", "Right - bottom": "right",
    "Left - top": "left", "Left - middle": "left", "Left - bottom": "left",
    "Top (above plot)": "top", "Bottom (below plot)": "bottom",
}
# Vertical alignment of a side (right/left) legend beside the plot. Bare "Right"/"Left"
# (older saved settings) map to middle, matching the historical centred placement.
LEGEND_VALIGN = {
    "Right - top": "top", "Right - middle": "middle", "Right - bottom": "bottom",
    "Left - top": "top", "Left - middle": "middle", "Left - bottom": "bottom",
}
#: Choices shown in the app's "Legend position" dropdown (order = display order).
LEGEND_POSITION_CHOICES = ["Right - top", "Right - middle", "Right - bottom",
                           "Left - top", "Left - middle", "Left - bottom",
                           "Top (above plot)", "Bottom (below plot)", "Hidden"]




def wrap_label(text, max_chars: int, max_lines: int = 3) -> str:
    """Wrap a feature name (protein, gene, metabolite, pathway) to lines of at most `max_chars`
    characters so it never runs into a neighbouring panel. Breaks between words first; a word
    longer than a line is split after a natural break character ( / , - _ : ; ) ] ) and only as
    a last resort mid-word. More than `max_lines` lines end in an ellipsis."""
    text = str(text)
    max_chars = max(6, int(max_chars))
    if len(text) <= max_chars:
        return text
    tokens = []
    for wi, word in enumerate(text.split()):
        if len(word) <= max_chars:
            parts = [word]
        else:
            chunks, cur = [], ""
            for ch in word:
                cur += ch
                if ch in "/,-_:;)]":
                    chunks.append(cur)
                    cur = ""
            if cur:
                chunks.append(cur)
            parts, line = [], ""
            for ch_ in chunks:
                if len(line) + len(ch_) <= max_chars:
                    line += ch_
                    continue
                if line:
                    parts.append(line)
                while len(ch_) > max_chars:
                    parts.append(ch_[:max_chars])
                    ch_ = ch_[max_chars:]
                line = ch_
            if line:
                parts.append(line)
        for pi, part in enumerate(parts):
            tokens.append((part, " " if (wi > 0 and pi == 0) else ""))
    lines, cur = [], tokens[0][0]
    for part, joiner in tokens[1:]:
        if len(cur) + len(joiner) + len(part) <= max_chars:
            cur += joiner + part
        else:
            lines.append(cur)
            cur = part
    lines.append(cur)
    if len(lines) > max_lines:
        lines = lines[:max_lines - 1] + [lines[max_lines - 1][: max_chars - 1].rstrip() + "…"]
    return "\n".join(lines)

def _agg_renderer(fig):
    """
    Return a renderer for `fig` via an explicitly-constructed Agg canvas, regardless of whatever
    backend is globally active. `fig.canvas.get_renderer()` only exists on Agg-family canvases --
    on a deployment where something else (Streamlit itself, or another imported library) selected
    a different backend before this module's `matplotlib.use("Agg")` had a chance to take effect,
    the figure's canvas can end up being a backend that has no `get_renderer` at all, raising
    `AttributeError: 'FigureCanvasBase' object has no attribute 'get_renderer'` the moment any
    control on the page triggers a rerun. Wrapping the figure in its own FigureCanvasAgg sidesteps
    whatever the ambient backend is and always works for this measure-only, throwaway figure.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    return canvas.get_renderer()


def _text_width_in(text: str, fontsize: float) -> float:
    tmp_fig = plt.figure()
    renderer = _agg_renderer(tmp_fig)
    t = tmp_fig.text(0, 0, text, fontsize=fontsize)
    width_in = t.get_window_extent(renderer=renderer).width / tmp_fig.dpi
    plt.close(tmp_fig)
    return width_in


def _text_size_in(text: str, fontsize: float, **kw):
    """(width, height) in inches of `text` as actually rendered (multi-line aware)."""
    tmp_fig = plt.figure()
    renderer = _agg_renderer(tmp_fig)
    t = tmp_fig.text(0, 0, text, fontsize=fontsize, **kw)
    bb = t.get_window_extent(renderer=renderer)
    plt.close(tmp_fig)
    return bb.width / tmp_fig.dpi, bb.height / tmp_fig.dpi


def _measure_legend_size_in(legend_handles, title, fontsize=9.5, title_fontsize=11.5, ncol=1):
    tmp_fig = plt.figure()
    leg = tmp_fig.legend(handles=legend_handles, loc="center left", fontsize=fontsize,
                         frameon=False, title=title, title_fontsize=title_fontsize,
                         alignment="left", ncol=ncol)
    renderer = _agg_renderer(tmp_fig)
    bbox = leg.get_window_extent(renderer=renderer)
    size_in = (bbox.width / tmp_fig.dpi, bbox.height / tmp_fig.dpi)
    plt.close(tmp_fig)
    return size_in


def classify_direction(df: pd.DataFrame, fc_col: str, sig_col: str, sig_cutoff: float, fc_threshold: float):
    def _classify(row):
        if pd.isna(row[fc_col]) or pd.isna(row[sig_col]):
            return "NS"
        if row[sig_col] < sig_cutoff and row[fc_col] > fc_threshold:
            return "Up"
        if row[sig_col] < sig_cutoff and row[fc_col] < -fc_threshold:
            return "Down"
        return "NS"
    return df.apply(_classify, axis=1)


def volcano_plot(
    stats_table: pd.DataFrame,
    # --- 1. Thresholds ---
    fc_col: str = "Log2FC",
    use_fc_not_log2: bool = False,          # if True, x-axis shows linear FC (requires Linear_FC column)
    fc_threshold: float = 0.0,              # in log2FC units regardless of display axis
    y_metric: str = "pvalue",               # 'pvalue' or 'fdr'
    sig_cutoff: float = 0.05,
    fdr_cutoff: float = None,               # optional secondary FDR cutoff (combined significance)
    threshold_line_style: str = "Dashed",
    # --- 2. Colors ---
    palette: str = "Default (Red/Teal/Gray)",
    up_color: str = None, down_color: str = None, ns_color: str = None,
    alpha: float = 0.75, edge_color: str = "none", edge_width: float = 0.0,
    # --- 3. Point size ---
    point_size: float = 14, point_shape: str = "Circle",
    # --- 4. Labels ---
    label_mode: str = "top_n",              # 'top_n' | 'significant_only' | 'manual' | 'none'
    top_label_n: int = 10,
    manual_labels: list = None,
    label_font_size: float = 7.5, label_color: str = None, label_bold: bool = False,
    label_italic: bool = False, repel_labels: bool = True,
    # --- 5. Legend ---
    legend_position: str = "Right", legend_format: str = "full",  # 'full' -> "Significantly Upregulated (n=5)", 'short' -> "Up (5)"
    # --- 6. Axes ---
    x_title: str = None, y_title: str = None, axis_font_size: float = 12, axis_bold: bool = False,
    x_limits: tuple = None, y_limits: tuple = None, x_tick_spacing: float = None, y_decimals: int = None,
    # --- 7. Title ---
    title: str = None, subtitle: str = None, title_font_size: float = 14,
    title_bold: bool = True, title_italic: bool = False, title_align: str = "center", hide_title: bool = False,
    # --- 8. Grid ---
    grid_mode: str = "none",                # 'none' | 'major' | 'both'
    grid_style: str = "Dashed", grid_color: str = "#D9D9D9",
    # --- 9. Threshold line style ---
    threshold_line_color: str = "gray", threshold_line_width: float = 0.7,
    # --- 10. Highlight specific metabolites ---
    highlight_names: list = None, highlight_color: str = "#FFD700",
    highlight_size_mult: float = 2.0, highlight_shape: str = "Star",
    # --- 12. Background ---
    background: str = "white",              # 'white' | 'transparent' | 'gray' | hex color
    # --- 13. Figure size ---
    fig_width_in: float = 7.2, fig_height_in: float = 6.4,
    # --- 15. Stats box ---
    show_stats_box: bool = False,
    # --- theme override ---
    theme: str = None,
    # --- FDR column label ---
    fdr_label: str = "FDR",   # column name for the BH FDR-adjusted p-value in stats_table
                              # (e.g. "BH P Value" for untargeted/unbiased assays)
):
    """
    Build a fully customizable, publication-oriented volcano plot.
    Returns (fig, annotated_df) where annotated_df has 'Direction' and the y-metric
    -log10 column added.
    """
    df = stats_table.copy()
    manual_labels = manual_labels or []
    highlight_names = highlight_names or []

    # Theme overrides (only fill in values the caller left at their defaults)
    if theme and theme in THEMES:
        t = THEMES[theme]
        if palette == "Default (Red/Teal/Gray)":
            palette = t["palette"]
        if grid_mode == "none":
            grid_mode = t["grid"]
        plt.rcParams["font.family"] = t["font"]

    colors = COLORBLIND_PALETTES.get(palette, COLORBLIND_PALETTES["Default (Red/Teal/Gray)"])
    up_c = up_color or colors["up"]
    down_c = down_color or colors["down"]
    ns_c = ns_color or colors["ns"]
    label_c_default = {"Up": up_c, "Down": down_c, "NS": "#666666"}

    sig_col = "p-value" if y_metric == "pvalue" else fdr_label
    y_label_default = f"-log$_{{10}}$({'p-value' if y_metric == 'pvalue' else fdr_label})"
    x_col = "Linear_FC" if use_fc_not_log2 else fc_col
    x_label_default = "FC" if use_fc_not_log2 else "log$_2$FC"

    df["neglog10"] = -np.log10(df[sig_col].replace(0, np.nextafter(0, 1)))
    df["Direction"] = classify_direction(df, fc_col, sig_col, sig_cutoff, fc_threshold)
    if fdr_cutoff is not None and fdr_label in df.columns:
        df.loc[df[fdr_label] >= fdr_cutoff, "Direction"] = "NS"

    counts = {d: int((df["Direction"] == d).sum()) for d in ("Down", "NS", "Up")}
    total = len(df)

    def _legend_label(direction_key, full_name):
        if legend_format == "short":
            short = {"Up": "Up", "Down": "Down", "NS": "NS"}[direction_key]
            return f"{short} ({counts[direction_key]})"
        return f"{full_name} (n={counts[direction_key]})"

    # Legend swatches are sized to visually match the actual plotted point size (scatter's `s` is an
    # area in points^2; Line2D's markersize is a diameter in points), clamped to a legible range --
    # so the legend honestly represents how big the dots on the plot really are instead of always
    # showing the same small fixed circle regardless of the point-size setting.
    legend_marker_size = float(np.clip(2.0 * math.sqrt(point_size / math.pi), 6.0, 13.0))
    legend_handles = [
        Line2D([0], [0], marker=MARKER_SHAPES.get(point_shape, "o"), linestyle="",
               markerfacecolor=up_c, markeredgecolor=edge_color if edge_width else "none",
               markersize=legend_marker_size, label=_legend_label("Up", "Significantly Upregulated")),
        Line2D([0], [0], marker=MARKER_SHAPES.get(point_shape, "o"), linestyle="",
               markerfacecolor=ns_c, markeredgecolor=edge_color if edge_width else "none",
               markersize=legend_marker_size, label=_legend_label("NS", "Not significant")),
        Line2D([0], [0], marker=MARKER_SHAPES.get(point_shape, "o"), linestyle="",
               markerfacecolor=down_c, markeredgecolor=edge_color if edge_width else "none",
               markersize=legend_marker_size, label=_legend_label("Down", "Significantly Downregulated")),
    ]
    show_legend = legend_position != "Hidden"
    legend_pos_key = LEGEND_POSITIONS.get(legend_position, "right")

    # -------------------------------------------------------------------
    # Layout: every band around the plot (title, subtitle, legend, statistics,
    # tick numbers, axis titles) gets space reserved from its MEASURED size at the
    # chosen font sizes, so nothing overlaps however large the fonts are set.
    # -------------------------------------------------------------------
    tick_fs = max(8.5, axis_font_size * 0.78)
    axis_w = "bold" if axis_bold else "normal"
    line_h = lambda f: f / 72.0 * 1.3
    ytick_sample = "0." + "0" * y_decimals if y_decimals else "00.0"
    ytick_w = max(_text_size_in(ytick_sample, tick_fs)[0], _text_size_in("-00", tick_fs)[0])
    left_margin_in = 0.15 + line_h(axis_font_size) + 0.12 + ytick_w + 0.1
    x_decor_in = line_h(tick_fs) + 0.08 + line_h(axis_font_size) + 0.06   # tick numbers + x-axis title
    right_margin_in = 0.3
    bottom_pad_in = 0.2
    gap_in = 0.3

    # Title band (title + optional subtitle), measured.
    title_h_in = subtitle_h_in = 0.0
    if not hide_title:
        _t = title or f"Volcano Plot ({'p-value' if y_metric == 'pvalue' else fdr_label} based)"
        title_h_in = _text_size_in(_t, title_font_size, fontweight="bold" if title_bold else "normal")[1]
        if subtitle:
            subtitle_h_in = _text_size_in(subtitle, title_font_size * 0.7)[1]
    subtitle_gap_in = 0.08 if subtitle_h_in else 0.0
    title_band_in = (0.12 + subtitle_h_in + subtitle_gap_in + title_h_in + 0.15) if title_h_in else 0.25

    legend_w_in = legend_h_in = 0.0
    if show_legend:
        ncol = 3 if legend_pos_key in ("top", "bottom") else 1
        legend_w_in, legend_h_in = _measure_legend_size_in(legend_handles, "Direction", ncol=ncol)
        legend_w_in += 0.25
        legend_h_in += 0.2

    # Statistics summary: plain text (no box), placed OUTSIDE the plot so it can
    # never cover points or their labels. It sits in the side column (under the
    # legend when the legend is left/right), otherwise in its own right-hand column.
    stats_fs = max(8.0, tick_fs * 0.95)
    stats_text = None
    stats_w_in = stats_h_in = 0.0
    if show_stats_box:
        stats_text = (f"Total: {total}\nUp: {counts['Up']}\nDown: {counts['Down']}\n"
                      f"FC cutoff: ±{fc_threshold}\n{fdr_label if fdr_cutoff else sig_col} cutoff: "
                      f"{fdr_cutoff if fdr_cutoff else sig_cutoff}")
        stats_w_in, stats_h_in = _text_size_in(stats_text, stats_fs, linespacing=1.35)
        stats_w_in += 0.1
    side_col_w_in = 0.0
    if legend_pos_key in ("right", "left"):
        side_col_w_in = max(legend_w_in, stats_w_in)
    elif stats_text:
        side_col_w_in = stats_w_in

    legend_band_in = (legend_h_in + gap_in) if legend_pos_key in ("top", "bottom") else 0.0
    stats_on_right = legend_pos_key != "left"

    # A top/bottom legend is centred on the plot; if it is wider than the plot,
    # widen the side margins so it never runs off the figure edge.
    if legend_pos_key in ("top", "bottom"):
        overhang = max(0.0, (legend_w_in - fig_width_in) / 2 + 0.1)
        left_margin_in = max(left_margin_in, overhang)
        right_margin_in = max(right_margin_in, overhang - (gap_in + side_col_w_in if side_col_w_in else 0.0))
    W = left_margin_in + fig_width_in + right_margin_in
    if side_col_w_in:
        W += gap_in + side_col_w_in
    H = title_band_in + legend_band_in + fig_height_in + x_decor_in + bottom_pad_in
    ax_bottom_in = bottom_pad_in + x_decor_in + (legend_band_in if legend_pos_key == "bottom" else 0.0)
    # Left legend: [legend column][gap][y-axis title + numbers][plot]
    ax_left_in = left_margin_in + ((0.15 + side_col_w_in + gap_in) if legend_pos_key == "left" else 0.0)
    if legend_pos_key == "left":
        W += 0.15

    fig = plt.figure(figsize=(W, H))
    if background == "transparent":
        fig.patch.set_alpha(0.0)
    elif background == "gray":
        fig.patch.set_facecolor("#EBEBEB")
    elif background not in ("white", None):
        fig.patch.set_facecolor(background)

    ax_left = ax_left_in / W
    ax_bottom = ax_bottom_in / H
    ax_top_in = ax_bottom_in + fig_height_in
    ax = fig.add_axes([ax_left, ax_bottom, fig_width_in / W, fig_height_in / H])
    if background == "gray":
        ax.set_facecolor("#F5F5F5")
    elif background not in ("white", "transparent", None):
        ax.set_facecolor(background)

    marker = MARKER_SHAPES.get(point_shape, "o")
    for direction, color in [("NS", ns_c), ("Down", down_c), ("Up", up_c)]:
        sub = df[df["Direction"] == direction]
        ax.scatter(sub[x_col], sub["neglog10"], c=color, s=point_size, alpha=alpha,
                   marker=marker, edgecolors=edge_color if edge_width else "none",
                   linewidths=edge_width)

    # Highlight specific metabolites on top of everything else
    if highlight_names:
        hl = df[df.index.isin(highlight_names)]
        hl_marker = MARKER_SHAPES.get(highlight_shape, "*") if highlight_shape in MARKER_SHAPES else "*"
        ax.scatter(hl[x_col], hl["neglog10"], c=highlight_color, s=point_size * highlight_size_mult,
                   marker=hl_marker, edgecolors="black", linewidths=0.8, zorder=5)

    # Threshold lines
    if threshold_line_style != "None":
        ls = LINE_STYLES.get(threshold_line_style, "--")
        # The FC threshold vertical line(s) are only meaningful once a non-zero
        # cutoff is set -- at fc_threshold == 0 both lines would sit on top of
        # each other at x=0, which is visual noise rather than a real cutoff.
        if not use_fc_not_log2 and fc_threshold != 0:
            ax.axvline(fc_threshold, color=threshold_line_color, linestyle=ls, lw=threshold_line_width)
            ax.axvline(-fc_threshold, color=threshold_line_color, linestyle=ls, lw=threshold_line_width)
        ax.axhline(-np.log10(sig_cutoff), color=threshold_line_color, linestyle=ls, lw=threshold_line_width)

    # Grid
    if grid_mode != "none":
        gs = LINE_STYLES.get(grid_style, "--")
        ax.grid(True, which="major", linestyle=gs, color=grid_color, linewidth=0.6)
        if grid_mode == "both":
            ax.minorticks_on()
            ax.grid(True, which="minor", linestyle=gs, color=grid_color, linewidth=0.3, alpha=0.5)
    else:
        ax.grid(False)

    # Labels
    texts = []
    if label_mode != "none":
        if label_mode == "manual":
            to_label = df[df.index.isin(manual_labels)]
        elif label_mode == "significant_only":
            to_label = df[df["Direction"] != "NS"]
        else:  # top_n
            to_label = df[df["Direction"] != "NS"].sort_values(sig_col).head(top_label_n)

        fontweight = "bold" if label_bold else "normal"
        fontstyle = "italic" if label_italic else "normal"
        for name, row in to_label.iterrows():
            color = label_color or label_c_default.get(row["Direction"], "#333333")
            texts.append(ax.text(row[x_col], row["neglog10"], wrap_label(name, 26, max_lines=2), fontsize=label_font_size, linespacing=0.95,
                                  color=color, fontweight=fontweight, fontstyle=fontstyle))

        if repel_labels and texts and HAS_ADJUST_TEXT:
            try:
                if not y_limits:                      # a little headroom so top labels can spread out
                    lo_, hi_ = ax.get_ylim()
                    ax.set_ylim(lo_, hi_ + 0.10 * (hi_ - lo_))
                adjust_text(texts, ax=ax, expand=(1.2, 1.6), force_text=(0.5, 0.8), force_static=(0.2, 0.4),
                            arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))
            except Exception:
                pass  # fall back to unadjusted labels rather than failing the whole plot

    # Axes
    ax.set_xlabel(x_title or x_label_default, fontsize=axis_font_size,
                  fontweight="bold" if axis_bold else "normal")
    ax.set_ylabel(y_title or y_label_default, fontsize=axis_font_size,
                  fontweight="bold" if axis_bold else "normal")
    # Tick numerals scale with the axis-title font size (kept a notch smaller, per typical figure
    # convention) instead of sitting at matplotlib's fixed default regardless of that setting --
    # otherwise a bumped-up axis-title font looks mismatched against small tick numbers.
    ax.tick_params(axis="both", labelsize=max(8.5, axis_font_size * 0.78))
    if x_limits:
        ax.set_xlim(x_limits)
    if y_limits:
        ax.set_ylim(y_limits)
    if x_tick_spacing:
        xmin, xmax = ax.get_xlim()
        ax.set_xticks(np.arange(np.floor(xmin / x_tick_spacing) * x_tick_spacing,
                                 np.ceil(xmax / x_tick_spacing) * x_tick_spacing + 1e-9, x_tick_spacing))
    if y_decimals is not None:
        ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter(f"%.{y_decimals}f"))

    theme_spines = THEMES.get(theme, {}).get("spines_top_right", False) if theme else False
    if not theme_spines:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Title / subtitle: the subtitle sits directly above the plot (or above a top
    # legend) and the title above that, each with its own measured band.
    band_base_pts = (legend_band_in if legend_pos_key == "top" else 0.0) * 72
    if not hide_title:
        final_title = title or f"Volcano Plot ({'p-value' if y_metric == 'pvalue' else fdr_label} based)"
        sub_x = 0.5 if title_align == "center" else (0.0 if title_align == "left" else 1.0)
        if subtitle:
            ax.annotate(subtitle, xy=(sub_x, 1.0), xycoords="axes fraction",
                        xytext=(0, band_base_pts + 0.12 * 72), textcoords="offset points",
                        fontsize=title_font_size * 0.7, ha=title_align, va="bottom", color="#555555",
                        annotation_clip=False)
        ax.set_title(final_title, fontsize=title_font_size, loc=title_align,
                     pad=band_base_pts + (0.12 + subtitle_h_in + subtitle_gap_in) * 72 + 2,
                     fontweight="bold" if title_bold else "normal",
                     fontstyle="italic" if title_italic else "normal")

    # Legend placement (inside the space reserved above)
    side_x_in = (ax_left_in + fig_width_in + gap_in) if legend_pos_key != "left" else 0.15
    side_top_in = ax_top_in
    if show_legend and legend_pos_key in ("right", "left"):
        # The side column (legend, plus the statistics summary under it when shown) is
        # aligned against the plot's top, vertical middle, or bottom, as the user chose.
        valign = LEGEND_VALIGN.get(legend_position, "middle")
        block_h_in = legend_h_in + ((0.15 + stats_h_in) if stats_text else 0.0)
        if valign == "bottom":
            side_top_in = ax_bottom_in + block_h_in
        elif valign == "middle":
            side_top_in = ax_bottom_in + (fig_height_in + block_h_in) / 2
        side_top_in = min(side_top_in, ax_top_in) if block_h_in <= fig_height_in else ax_top_in
    if show_legend:
        leg = None
        common = dict(bbox_transform=fig.transFigure, fontsize=9.5, frameon=False,
                      title="Direction", title_fontsize=11.5)
        if legend_pos_key in ("right", "left"):
            leg = fig.legend(handles=legend_handles, loc="upper left",
                             bbox_to_anchor=(side_x_in / W, side_top_in / H), alignment="left", **common)
            side_top_in -= legend_h_in + 0.15
        elif legend_pos_key == "top":
            leg = fig.legend(handles=legend_handles, loc="lower center", ncol=3,
                             bbox_to_anchor=(ax_left + (fig_width_in / W) / 2, (ax_top_in + 0.12) / H), **common)
        elif legend_pos_key == "bottom":
            leg = fig.legend(handles=legend_handles, loc="lower center", ncol=3,
                             bbox_to_anchor=(ax_left + (fig_width_in / W) / 2, bottom_pad_in / H), **common)
        # A bold legend title reads as more intentional/designed than the default (plain-weight,
        # same-size-as-labels) title matplotlib gives it.
        if leg is not None and leg.get_title() is not None:
            leg.get_title().set_fontweight("bold")

    # Statistics summary (plain text, no box) in the side column
    if stats_text:
        fig.text((side_x_in + 0.12) / W, side_top_in / H, stats_text, fontsize=stats_fs,
                 va="top", ha="left", color="#333333", linespacing=1.35)

    return fig, df


def get_direction_table(annotated_df: pd.DataFrame, direction: str = None) -> pd.DataFrame:
    """direction: 'Up', 'Down', 'NS', or None for all significant (Up+Down)."""
    if direction is None:
        return annotated_df[annotated_df["Direction"] != "NS"]
    return annotated_df[annotated_df["Direction"] == direction]


def export_settings_json(settings: dict) -> str:
    """Serialize the plot settings dict to a pretty JSON string for reproducibility."""
    return json.dumps(settings, indent=2, default=str)


def top_biomarker_labels(stats_table: pd.DataFrame, sig_col: str = "FDR", n: int = 20):
    return stats_table.sort_values(sig_col).head(n)


def export_figure(fig, fmt: str = "png", dpi: int = 300) -> bytes:
    """
    Export a matplotlib figure to bytes. fmt: 'png', 'pdf', 'svg', 'eps', 'jpeg'/'jpg', or 'tiff'.
    dpi applies to raster formats (300/600/1200 supported; ignored for vector pdf/svg/eps).
    """
    fmt = fmt.lower()
    buf = io.BytesIO()
    save_kwargs = {"format": "png" if fmt == "png" else fmt, "bbox_inches": "tight"}
    if fmt in ("jpeg", "jpg", "tiff", "png"):
        save_kwargs["dpi"] = dpi
    if fmt in ("jpeg", "jpg"):
        save_kwargs["facecolor"] = fig.get_facecolor() if fig.get_alpha() != 0 else "white"
    fig.savefig(buf, **save_kwargs)
    buf.seek(0)
    return buf.read()
