# Omnix: Integrix chart error, overlapping size legend, black figure text (2026-11-11)

## StreamlitDuplicateElementId in Integrix (Results)
Choosing 15 pathways on the Pathways tab made its canonical-pathway chart identical to the one in Overview. On
Streamlit 1.52 both tabs are drawn in the same run, and two identical charts got the same automatic ID.
Every Integrix chart and table now has its own key. New test: `test_same_chart_twice_on_one_page`.
Files: `apps/integrix/app.py` (part 5), `tests/test_integrix_portal.py` (part 1).

## Size legend circles overlapping (Integrix pathway activity dot plot)
Each legend row is now as tall as the largest circle, so the q 0.05 / 0.001 / 1e-6 circles never touch. The class
labels next to the dots are black (were teal).
File: `apps/integrix/trace_core/pubfig.py` (part 5).

## All figure text in pure black, every platform
New `omnix_portal/figure_style.py`, switched on in `portal.py`:
- Matplotlib figures (all platforms, on screen and in every download): titles, axis labels, tick labels, legends,
  colour bars and annotations are drawn in #000000.
- Interactive (Plotly) charts in Integrix and Fluxomics: the same.
- Text that has to stay light (white labels on dark heatmap cells or coloured boxes) and Integrix's dark theme are
  left unchanged.
Files: `omnix_portal/figure_style.py` (new), `omnix_portal/portal.py` (part 1).

## Checks
- Portal tests: 77 passed, 1 skipped on Streamlit 1.52.2 and on 1.67.
- Integrix tests: 49 passed.
- Browser: the Proteomics QC page shows all figure text in black.
