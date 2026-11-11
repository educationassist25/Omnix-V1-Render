"""
Sample names on figure axes: one rule for every Omnix platform.

* Every sample is labelled (labels are never thinned) and the full name is shown by default.
* The label font follows the room each sample gets on the axis, between 6 and 9 pt. When even 6 pt
  would not fit, the figure grows along that axis, so neighbouring names never overlap.
* The figure also grows across the axis by the length of the longest name, so long names do not
  squeeze the plot itself.
* The Omnix setting "Sample names in figures" (session key ``omnix_sample_label_max``; 0 = full
  names, the default) can shorten long names in the middle, e.g. ``HealthyControl_1…RunA``. Tables and
  downloaded data always keep the full names.

The same file is kept in each platform's module folder (here: core/), so the platforms behave identically.
"""
from __future__ import annotations

import math

MIN_FS = 6.0                 # smallest sample-label font (pt)
MAX_FS = 9.0                 # largest sample-label font (pt)
SPACING = 1.35               # label pitch / font size: room for the glyph height plus a gap
BASE_ROOM_IN = 0.45          # label length the default figure margins already absorb (in)
SETTING_KEY = "omnix_sample_label_max"


def max_chars() -> int:
    """The user's setting: 0 (full names) or the longest label length in characters."""
    try:
        import streamlit as st
        return int(st.session_state.get(SETTING_KEY, 0) or 0)
    except Exception:
        return 0


def shorten(name, limit: int | None = None) -> str:
    """Name for display: unchanged when it fits, otherwise the start and end joined by an ellipsis."""
    s = str(name)
    n = max_chars() if limit is None else int(limit or 0)
    if n <= 0 or len(s) <= n:
        return s
    keep = max(n - 1, 2)
    head = math.ceil(keep * 0.6)
    tail = keep - head
    return s[:head] + "…" + (s[-tail:] if tail > 0 else "")


def display(names) -> list[str]:
    """Display labels for a sequence of sample names; duplicates after shortening get a number."""
    out, seen = [], {}
    for x in names:
        d = shorten(x)
        if d in seen:
            seen[d] += 1
            d = f"{d} ({seen[d]})"
        else:
            seen[d] = 1
        out.append(d)
    return out


def text_width_in(text: str, fontsize: float) -> float:
    """Rendered width of `text` at `fontsize` (inches), from the font outline."""
    try:
        from matplotlib.textpath import TextPath
        return TextPath((0, 0), str(text), size=fontsize).get_extents().width / 72.0
    except Exception:
        return len(str(text)) * fontsize * 0.6 / 72.0


def set_sample_ticks(ax, names, axis: str = "x", positions=None, rotation: float | None = None,
                     fontsize: float | None = None, grow: bool = True, max_fs: float = MAX_FS):
    """Label every sample on `axis` ('x' or 'y') of `ax` without overlaps.

    positions: tick positions (default 0..n-1, as for imshow; boxplots use 1..n).
    Returns the font size used. Call before fig.tight_layout()."""
    fig = ax.figure
    shown = display(names)
    n = max(len(shown), 1)
    pos = list(positions) if positions is not None else list(range(len(shown)))
    W, H = fig.get_size_inches()
    box = ax.get_position()
    along_in = box.width * W if axis == "x" else box.height * H
    room_pt = along_in * 72.0 / n
    fs = fontsize if fontsize else min(max_fs, room_pt / SPACING)
    if grow and fs < MIN_FS:
        extra = MIN_FS * SPACING * n / 72.0 - along_in
        if axis == "x":
            W += extra / max(box.width, 0.1)
        else:
            H += extra / max(box.height, 0.1)
        fs = MIN_FS
    fs = max(fs, MIN_FS) if not fontsize else fs
    longest = max((text_width_in(s, fs) for s in shown), default=0.0)
    grown = getattr(fig, "_omx_label_room", set())
    if grow and axis not in grown and longest > BASE_ROOM_IN:
        if axis == "x":
            H += longest - BASE_ROOM_IN
        else:
            W += longest - BASE_ROOM_IN
        grown.add(axis)
        fig._omx_label_room = grown
    if grow:
        fig.set_size_inches(min(W, 60.0), min(H, 60.0), forward=False)
    if axis == "x":
        ax.set_xticks(pos)
        ax.set_xticklabels(shown, rotation=90 if rotation is None else rotation, fontsize=fs)
    else:
        ax.set_yticks(pos)
        ax.set_yticklabels(shown, rotation=0 if rotation is None else rotation, fontsize=fs)
    return fs


def column_label_fontsize(font_size: float, cell_in: float) -> float:
    """Heatmap column (sample) label size: 75% of the base font, never taller than a column."""
    return max(4.0, min(0.75 * font_size, cell_in * 72.0 / SPACING))


def correlation_matrix_figure(corr, title: str, cbar_label: str = "Pearson r", cmap: str = "viridis",
                              min_side_in: float = 3.0):
    """Sample-by-sample correlation matrix with every sample labelled on both axes.

    The square matrix is laid out in inches: its side grows with the number of samples (never less than
    6 pt per label) and the margins take the longest name, so labels neither overlap nor shrink the
    matrix. Same look in every platform: viridis, colour bar on the right, bold 11 pt title."""
    import matplotlib.pyplot as plt
    shown = display(corr.columns)
    n = max(len(shown), 1)
    side = max(min_side_in, MIN_FS * SPACING * n / 72.0)
    fs = max(MIN_FS, min(MAX_FS, side * 72.0 / n / SPACING))
    lab = max((text_width_in(s, fs) for s in shown), default=0.3)
    left, bottom = lab + 0.2, lab + 0.2
    right, top = 1.15, 0.5
    title_w = text_width_in(title, 11) + 0.3
    W = max(left + side + right, title_w)
    H = bottom + side + top
    fig = plt.figure(figsize=(W, H))
    x0 = left + max(0.0, (W - (left + side + right)) / 2)
    ax = fig.add_axes([x0 / W, bottom / H, side / W, side / H])
    cax = fig.add_axes([(x0 + side + 0.18) / W, (bottom + 0.1 * side) / H, 0.16 / W, 0.8 * side / H])
    vals = corr.values
    im = ax.imshow(vals, cmap=cmap, vmin=float(vals[~(vals != vals)].min()) if vals.size else 0, vmax=1,
                   interpolation="nearest")
    ax.set_xticks(range(len(shown)))
    ax.set_xticklabels(shown, rotation=90, fontsize=fs)
    ax.set_yticks(range(len(shown)))
    ax.set_yticklabels(shown, fontsize=fs)
    ax.set_title(title, fontsize=11, fontweight="bold")
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(cbar_label, fontsize=9)
    cb.ax.tick_params(labelsize=8)
    return fig
