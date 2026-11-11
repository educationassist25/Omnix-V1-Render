"""
One text colour for every figure in Omnix: pure black.

The platforms were written separately and use many greys for figure text (#333, #444, #555, slate blue-greys,
coloured labels). In Omnix every word in a figure — titles, axis labels, tick labels, legends, colour-bar labels,
annotations — is drawn in pure black (#000000), on screen and in every downloaded file.

* Matplotlib: applied just before a figure is drawn or saved (Figure.savefig, which st.pyplot and every download
  use), so it covers all platforms without touching their plotting code.
* Plotly (interactive charts): applied when the chart is sent to the page.

Text that must stay light — white labels on dark heatmap cells or on coloured boxes — is left as it is, and
figures with a dark background (Integrix's dark theme) are not changed.
"""
from __future__ import annotations

BLACK = "#000000"
_installed = False


def _lum(color) -> float | None:
    """Relative brightness 0 (black) .. 1 (white), or None if the colour cannot be read."""
    try:
        from matplotlib import colors as mcolors
        r, g, b, a = mcolors.to_rgba(color)
        if a == 0:
            return None
        return 0.2126 * r + 0.7152 * g + 0.0722 * b
    except Exception:
        return None


def _light_background(fig) -> bool:
    lum = _lum(fig.get_facecolor())
    return lum is None or lum > 0.5


def blacken_mpl(fig):
    """Every visible text in a light-background matplotlib figure becomes pure black (light text stays)."""
    try:
        if not _light_background(fig):
            return
        from matplotlib.text import Text
        for t in fig.findobj(Text):
            lum = _lum(t.get_color())
            if lum is not None and lum > 0.82:       # white / near-white text on a dark cell or box
                continue
            t.set_color(BLACK)
    except Exception:
        pass                                          # styling must never stop a figure from being shown


def _rgba(c):
    """(r, g, b, a) in 0..1 for hex, named or rgb()/rgba() colours; None if unreadable."""
    try:
        c = str(c).strip()
        if c.startswith(("rgba(", "rgb(")):
            v = [float(x) for x in c[c.index("(") + 1:c.index(")")].split(",")]
            return (v[0] / 255, v[1] / 255, v[2] / 255, v[3] if len(v) > 3 else 1.0)
        from matplotlib import colors as mcolors
        return mcolors.to_rgba(c)
    except Exception:
        return None


def _plotly_light(fig) -> bool:
    """True when the chart is drawn for a light page. A transparent background takes the page's colour, so the
    theme is read from the text colour instead (light text = dark theme)."""
    try:
        for bg in (fig.layout.paper_bgcolor, fig.layout.plot_bgcolor):
            v = _rgba(bg) if bg else None
            if v is not None and v[3] > 0.05:
                return 0.2126 * v[0] + 0.7152 * v[1] + 0.0722 * v[2] > 0.5
        fc = _rgba(fig.layout.font.color) if fig.layout.font and fig.layout.font.color else None
        return fc is None or (0.2126 * fc[0] + 0.7152 * fc[1] + 0.0722 * fc[2]) < 0.6
    except Exception:
        return True


_OWN_SCALE = {"heatmap", "heatmapgl", "contour", "histogram2d", "histogram2dcontour", "surface", "choropleth",
              "choroplethmapbox", "choroplethmap", "densitymapbox", "densitymap", "volume", "isosurface", "cone",
              "streamtube"}


def blacken_plotly(fig):
    """Plotly figure: titles, axes, legend, colour bars and annotations in pure black (light background only)."""
    try:
        if not _plotly_light(fig):
            return
        fig.update_layout(font_color=BLACK, title_font_color=BLACK, legend_font_color=BLACK,
                          legend_title_font_color=BLACK)
        fig.update_xaxes(tickfont_color=BLACK, title_font_color=BLACK)
        fig.update_yaxes(tickfont_color=BLACK, title_font_color=BLACK)
        for ann in fig.layout.annotations or ():
            c = _rgba(ann.font.color) if ann.font and ann.font.color else None
            if c is not None and 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2] > 0.82:
                continue
            ann.font.color = BLACK
        # Colour bars: restyle only a colour bar the chart really shows. Writing colour-bar settings onto a trace
        # or layout that has none makes Plotly draw an extra, empty colour bar for it, and those bars stack on
        # top of each other at the right of the chart.
        for name, ca in fig.layout.to_plotly_json().items():
            if name.startswith("coloraxis") and isinstance(ca, dict) and ca.get("showscale", True) is not False:
                fig.layout[name].colorbar.tickfont.color = BLACK
                fig.layout[name].colorbar.title.font.color = BLACK

        def _tr(t):
            try:
                j = t.to_plotly_json()
                # traces that draw their own colour bar unless showscale is False
                if (t.type in _OWN_SCALE and not j.get("coloraxis") and j.get("showscale", True) is not False):
                    t.colorbar.tickfont.color = BLACK
                    t.colorbar.title.font.color = BLACK
                mk = j.get("marker")
                # marker colour bars are off unless showscale is True
                if isinstance(mk, dict) and mk.get("showscale") is True and not mk.get("coloraxis"):
                    t.marker.colorbar.tickfont.color = BLACK
                    t.marker.colorbar.title.font.color = BLACK
            except Exception:
                pass
        fig.for_each_trace(_tr)
    except Exception:
        pass


def install():
    """Patch matplotlib and st.plotly_chart once per server process."""
    global _installed
    if _installed:
        return
    _installed = True
    try:
        from matplotlib.figure import Figure
        _savefig = Figure.savefig

        def savefig(self, *args, **kwargs):
            blacken_mpl(self)
            return _savefig(self, *args, **kwargs)
        savefig.__wrapped__ = _savefig
        Figure.savefig = savefig
    except Exception:
        pass
    try:
        import functools
        import streamlit as st
        from streamlit.delta_generator import DeltaGenerator
        _pc = DeltaGenerator.plotly_chart

        @functools.wraps(_pc)
        def plotly_chart(self, figure_or_data, *args, **kwargs):
            if hasattr(figure_or_data, "update_layout"):
                blacken_plotly(figure_or_data)
            return _pc(self, figure_or_data, *args, **kwargs)
        DeltaGenerator.plotly_chart = plotly_chart
        main_dg = getattr(st, "_main", None)
        if main_dg is not None:
            st.plotly_chart = functools.wraps(st.plotly_chart)(
                lambda figure_or_data, *a, **k: plotly_chart(main_dg, figure_or_data, *a, **k))
    except Exception:
        pass
