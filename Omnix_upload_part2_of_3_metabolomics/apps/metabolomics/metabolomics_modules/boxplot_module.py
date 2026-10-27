"""
boxplot_module.py - "Boxplot of Metabolites" module.

Lets users compare one or more metabolites across any chosen subset of sample
groups, with statistics (p-value/FDR) computed on log2-transformed data (never
raw peak areas, per this app's statistical policy) and displayed on each panel.
"""

import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.figure

from metabolomics_modules import stats_analysis

matplotlib.use("Agg")

GROUP_PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2", "#937860", "#DA8BC3", "#8C8C8C"]


def compute_stats_for_metabolites(log2_df: pd.DataFrame, meta: pd.DataFrame,
                                   metabolites: list, groups: list, group_col: str = "Group",
                                   fdr_label: str = "FDR") -> pd.DataFrame:
    """
    Compute p-value/FDR for the given metabolites across the given groups (from
    `group_col` -- any categorical metadata column), using log2-transformed data.
    If exactly 2 groups are given, runs a Welch's t-test; if 3+, runs one-way
    ANOVA. FDR is corrected across only the selected metabolite subset shown here
    (not the full feature panel) — noted in the UI to avoid confusion with a
    panel-wide FDR from the Statistics tab.

    fdr_label : str, default "FDR"
        Column name for the returned Benjamini-Hochberg FDR-adjusted p-value
        (e.g. "BH P Value" for untargeted/unbiased assays).
    """
    sub = log2_df.loc[log2_df.index.intersection(metabolites)]
    sample_cols = [c for c in sub.columns if meta.loc[c, group_col] in groups]
    sub = sub[sample_cols]

    if len(groups) == 2:
        a_samples = [c for c in sample_cols if meta.loc[c, group_col] == groups[0]]
        b_samples = [c for c in sample_cols if meta.loc[c, group_col] == groups[1]]
        result = stats_analysis.two_group_test(sub, a_samples, b_samples, method="ttest",
                                                 fdr_label=fdr_label)
        return result
    else:
        group_map = meta.loc[sample_cols, group_col]
        anova_table, _ = stats_analysis.anova_test(sub, group_map, posthoc="tukey",
                                                     fdr_label=fdr_label)
        return anova_table.rename(columns={"ANOVA p-value": "p-value"})


class _TextMeasurer:
    """Measures rendered text size (inches) with one reusable off-screen Agg canvas.
    Used to size panels from the actual label text instead of guessing, so nothing
    collides at any font size."""

    def __init__(self, font_family: str = "sans-serif"):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        self._fig = matplotlib.figure.Figure()
        self._canvas = FigureCanvasAgg(self._fig)
        self._r = self._canvas.get_renderer()
        self._family = font_family

    def size(self, text: str, fontsize: float, bold: bool = False):
        t = self._fig.text(0, 0, str(text), fontsize=fontsize, fontfamily=self._family,
                           fontweight="bold" if bold else "normal")
        bb = t.get_window_extent(renderer=self._r)
        t.remove()
        return bb.width / self._fig.dpi, bb.height / self._fig.dpi

    def width(self, text, fontsize, bold=False):
        return self.size(text, fontsize, bold)[0]


def _wrap_to_width(text: str, max_w_in: float, fontsize: float, meas: "_TextMeasurer",
                   bold: bool = False, max_lines: int = 3, split_words: bool = True) -> str:
    """Wrap `text` so every line fits within max_w_in (measured, not guessed).

    Breaks between words first. A single word longer than a line -- common for lipid and
    metabolite names such as PS(20:5(5Z,8Z,11Z,14Z,17Z)/20:3(8Z,11Z,14Z)) -- is split after
    a natural break character ( / , - _ : ; ) ] ), and only as a last resort between two
    characters, so a name never runs into the neighbouring panel. If more than `max_lines`
    lines would be needed, the last line ends with an ellipsis."""
    text = str(text)
    if not text.strip():
        return text
    fits = lambda s: meas.width(s, fontsize, bold) <= max_w_in

    def split_word(word):
        """Pieces of one over-long word, each fitting on a line."""
        chunks, cur = [], ""
        for ch in word:
            cur += ch
            if ch in "/,-_:;)]":
                chunks.append(cur)
                cur = ""
        if cur:
            chunks.append(cur)
        pieces, line_ = [], ""
        for ch_ in chunks:
            if fits(line_ + ch_):
                line_ += ch_
                continue
            if line_:
                pieces.append(line_)
                line_ = ""
            while ch_ and not fits(ch_):          # a chunk longer than a line: split by characters
                k = len(ch_)
                while k > 1 and not fits(ch_[:k]):
                    k -= 1
                pieces.append(ch_[:k])
                ch_ = ch_[k:]
            line_ = ch_
        if line_:
            pieces.append(line_)
        return pieces

    # tokens: (piece, joiner) -- joiner is " " between words, "" between pieces of one word
    tokens = []
    for wi, word in enumerate(text.split()):
        parts = [word] if (fits(word) or not split_words) else split_word(word)
        for pi, part in enumerate(parts):
            tokens.append((part, " " if (wi > 0 and pi == 0) else ""))
    lines, cur = [], tokens[0][0]
    for part, joiner in tokens[1:]:
        trial = cur + joiner + part
        if fits(trial):
            cur = trial
        else:
            lines.append(cur)
            cur = part
    lines.append(cur)
    if len(lines) > max_lines and not split_words:     # caller checks the width and rotates instead
        lines = lines[:max_lines - 1] + [" ".join(lines[max_lines - 1:])]
    elif len(lines) > max_lines:
        last = lines[max_lines - 1]
        while last and not fits(last + "…"):
            last = last[:-1]
        lines = lines[:max_lines - 1] + [last.rstrip() + "…"]
    return "\n".join(lines)


def boxplot_metabolites(
    log2_df: pd.DataFrame, meta: pd.DataFrame, metabolites: list, groups: list,
    stats_table: pd.DataFrame = None,
    fig_width_in: float = None, fig_height_in: float = None,
    font_size: float = 10, font_family: str = "sans-serif",
    group_colors: dict = None, show_points: bool = True,
    show_mean: bool = False, show_median: bool = True,
    ncols: int = None, group_col: str = "Group",
    fdr_label: str = "FDR",
):
    """
    One panel per metabolite, boxes grouped by the selected sample groups (from
    `group_col` -- any categorical metadata column: Group/Diagnosis, Gender,
    Treatment, Ethnicity, etc.). Returns the matplotlib figure.
    """
    n = len(metabolites)
    ncols = int(ncols or min(3, n))
    ncols = max(1, min(ncols, n))
    nrows = int(np.ceil(n / ncols))
    ng = max(len(groups), 1)
    group_colors = group_colors or {g: GROUP_PALETTE[i % len(GROUP_PALETTE)] for i, g in enumerate(groups)}

    # ------------------------------------------------------------------
    # Layout is computed from MEASURED text sizes, so every element (title,
    # p-value lines, y-axis label, group names) gets its own reserved space at any
    # font size, instead of sharing space and colliding as the font grows.
    # ------------------------------------------------------------------
    meas = _TextMeasurer(font_family)
    fs = float(font_size)
    fs_title, fs_tick, fs_stat = fs + 1, max(6.0, fs - 1), max(6.0, fs - 2)
    line = lambda f: f / 72.0 * 1.25          # one text line incl. leading (inches)
    scale = max(1.0, fs / 10.0)

    # --- statistics text per metabolite (no box; one line if it fits, else two) ---
    stat_parts = {}
    if stats_table is not None:
        for met in metabolites:
            if met in stats_table.index:
                row = stats_table.loc[met]
                p = row.get("p-value", np.nan)
                fdr = row.get(fdr_label, np.nan)
                parts = [f"p = {p:.3g}"] + ([f"{fdr_label} = {fdr:.3g}"] if not pd.isna(fdr) else [])
                stat_parts[met] = parts
    stat_w_2line = max((meas.width(pt, fs_stat) for parts in stat_parts.values() for pt in parts), default=0.0)

    # --- horizontal budget ---
    ytick_w = max(meas.width(s_, fs_tick) for s_ in ("00.0", "-00.0"))
    y_area = line(fs) + ytick_w + 0.16                     # y-label + tick labels + tick marks
    right_pad = 0.18 * scale                               # breathing room between panels
    slot_rot_min = line(fs_tick) * 1.3 / np.sin(np.deg2rad(45))   # min tick spacing for 45° labels
    default_cell_w = 2.9 * scale
    cell_w = (fig_width_in / ncols) if fig_width_in else default_cell_w
    ax_w = cell_w - y_area - right_pad
    ax_w = max(ax_w, 1.3 * scale, ng * slot_rot_min, stat_w_2line + 0.1)

    # --- group (x tick) labels: horizontal (wrapped to 2 lines) if they fit, else 45° ---
    slot = ax_w / ng
    wrapped = [_wrap_to_width(g, slot * 0.88, fs_tick, meas, max_lines=2, split_words=False) for g in groups]
    fits_h = all(max(meas.width(l_, fs_tick) for l_ in w.split("\n")) <= slot * 0.88 for w in wrapped)
    if fits_h:
        tick_labels = wrapped
        rotate = False
        n_lines = max(w.count("\n") + 1 for w in wrapped)
        x_depth = n_lines * line(fs_tick) + 0.08
    else:
        tick_labels = [str(g) for g in groups]
        rotate = True
        max_w = max(meas.width(g, fs_tick) for g in tick_labels)
        # Extra tick pad for angled labels: the first label runs down-left toward the
        # y-axis numbers, so it must start far enough below the axis to clear the
        # lowest y tick label (half a text line tall) in the bottom-left corner.
        x_pad_in = line(fs_tick) * 0.6
        x_depth = (max_w + line(fs_tick)) * np.sin(np.deg2rad(45)) + x_pad_in + 0.10

    # --- titles wrapped to the panel width; heights of title + stats band ---
    titles = {m: _wrap_to_width(m, ax_w + y_area * 0.5, fs_title, meas, bold=True, max_lines=3)
              for m in metabolites}
    title_lines = max(t.count("\n") + 1 for t in titles.values())
    title_h = title_lines * line(fs_title)
    one_line_ok = all(meas.width("   ".join(pts), fs_stat) <= ax_w for pts in stat_parts.values())
    stat_lines = 0 if not stat_parts else (1 if one_line_ok else max(len(p_) for p_ in stat_parts.values()))
    stat_h = stat_lines * line(fs_stat)
    gap_axes_stat = 0.07 * scale                           # axes top -> p-value text
    gap_stat_title = 0.12 * scale if stat_lines else 0.0   # p-value text -> metabolite name
    top_band = gap_axes_stat + stat_h + gap_stat_title + title_h + 0.08 * scale

    # --- vertical budget ---
    default_cell_h = 2.7 * scale
    cell_h = (fig_height_in / nrows) if fig_height_in else default_cell_h
    bottom_band = x_depth + 0.10 * scale
    ax_h = max(cell_h - top_band - bottom_band, 1.2 * scale)

    cell_w = y_area + ax_w + right_pad
    cell_h = top_band + ax_h + bottom_band
    fig_w = max(fig_width_in or 0, cell_w * ncols)
    fig_h = max(fig_height_in or 0, cell_h * nrows)
    # If the user asked for a LARGER figure than needed, give the extra room to the plots.
    ax_w += (fig_w - cell_w * ncols) / ncols
    ax_h += (fig_h - cell_h * nrows) / nrows
    cell_w, cell_h = fig_w / ncols, fig_h / nrows

    fig = plt.figure(figsize=(fig_w, fig_h))
    rng = np.random.default_rng(0)

    for idx, met in enumerate(metabolites):
        r, c = divmod(idx, ncols)
        left = c * cell_w + y_area
        bottom = fig_h - (r + 1) * cell_h + bottom_band
        ax = fig.add_axes([left / fig_w, bottom / fig_h, ax_w / fig_w, ax_h / fig_h])

        data_per_group = []
        for g in groups:
            cols = [s_ for s_ in log2_df.columns if meta.loc[s_, group_col] == g]
            vals = log2_df.loc[met, cols].dropna().values if met in log2_df.index else np.array([])
            data_per_group.append(vals)

        bp = ax.boxplot(data_per_group, tick_labels=tick_labels, patch_artist=True,
                        showmeans=show_mean, meanline=show_mean, showfliers=not show_points)
        for patch, g in zip(bp["boxes"], groups):
            patch.set_facecolor(group_colors.get(g, "#4C72B0"))
            patch.set_alpha(0.7)
        for line_ in bp["medians"]:
            if show_median:
                line_.set_color("black")
                line_.set_linewidth(1.2)
            else:
                line_.set_visible(False)

        if show_points:
            for i, gdata in enumerate(data_per_group):
                if len(gdata):
                    x = rng.normal(i + 1, 0.045, size=len(gdata))
                    ax.scatter(x, gdata, color="black", alpha=0.55, s=14 * scale, zorder=3, edgecolors="none")

        # Small symmetric margin so points never touch the frame; the statistics no
        # longer live inside the plotting area, so no data-space headroom is needed.
        nonempty = [g for g in data_per_group if len(g)]
        if nonempty:
            all_vals = np.concatenate(nonempty)
            lo, hi = float(np.nanmin(all_vals)), float(np.nanmax(all_vals))
            rng_ = (hi - lo) or 1.0
            ax.set_ylim(lo - 0.07 * rng_, hi + 0.07 * rng_)

        ax.set_ylabel("Relative Abundance", fontsize=fs, fontfamily=font_family, labelpad=4)
        ax.tick_params(labelsize=fs_tick)
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_fontfamily(font_family)
        if rotate:
            ax.tick_params(axis="x", pad=3.5 + x_pad_in * 72)
            for lbl in ax.get_xticklabels():
                lbl.set_rotation(45)
                lbl.set_ha("right")
                lbl.set_rotation_mode("anchor")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Statistics: plain text (no box) in its own band just above the plot, with
        # clear space between it and the metabolite name above it.
        if met in stat_parts:
            txt = ("   ".join(stat_parts[met]) if stat_lines == 1 else "\n".join(stat_parts[met]))
            ax.annotate(txt, xy=(0.5, 1.0), xycoords="axes fraction",
                        xytext=(0, gap_axes_stat * 72), textcoords="offset points",
                        ha="center", va="bottom", fontsize=fs_stat, fontfamily=font_family,
                        color="#333333", linespacing=1.15, annotation_clip=False)
        title_pad_pts = (gap_axes_stat + stat_h + gap_stat_title) * 72
        ax.set_title(titles[met], fontsize=fs_title, fontweight="bold", loc="center",
                     fontfamily=font_family, pad=title_pad_pts, linespacing=1.1)

    return fig


def export_figure(fig, fmt: str = "png", dpi: int = 300) -> bytes:
    """Export a matplotlib figure to bytes. fmt: 'png', 'pdf', 'svg', 'jpeg'/'jpg', or 'tiff'."""
    fmt = fmt.lower()
    buf = io.BytesIO()
    save_kwargs = {"format": "png" if fmt == "png" else fmt, "bbox_inches": "tight"}
    if fmt in ("jpeg", "jpg", "tiff", "png"):
        save_kwargs["dpi"] = dpi
    if fmt in ("jpeg", "jpg"):
        save_kwargs["facecolor"] = "white"
    fig.savefig(buf, **save_kwargs)
    buf.seek(0)
    return buf.read()
