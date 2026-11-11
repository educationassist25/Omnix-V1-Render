# Omnix: final review (2026-11-10)

This update also contains the 2026-11-09 speed and sign-in fixes (`requirements.txt`, `render_start.sh`,
`omnix_portal/account.py`): upload everything below together.

## 1. Proteomics: Network Analysis removed
- The **Network Analysis → P-P Interaction** page and its code (`proteomics_modules/ppi.py`) are removed.
- Removed from the home card, Help, the Proteomics tutorial and the Proteomics README.
- STRING *enrichment* stays on the pathway page.
- The page called "GSEA" in Proteomics and Transcriptomics is now **Pathway Enrichment**, because it runs GSEA,
  ORA and STRING. The tutorials and Help use the same name.
- An old link to a removed page opens the first page instead of failing (all three platforms).

## 2. Long sample names in figures
New `sample_labels.py`, the same file in Metabolomics/Lipidomics, Proteomics, Transcriptomics and Fluxomics:
- **Every sample is labelled**: the Transcriptomics distribution and correlation plots used to label only every
  2nd or 3rd sample.
- **Full names by default.** The font follows the space per sample (6–9 pt, never the 4–5 pt it was).
- **No squeezed or overlapping plots**: the figure grows to fit the labels.
- **Correlation matrices**: one layout in all platforms. The matrix grows with the number of samples and the
  margins take the longest name.
- **Heatmaps**: sample labels are 75% of the base font (was 60%). The colour-bar label reads
  "Z-score (row-scaled)" in all three platforms.
- **Fluxomics**: per-sample plots use the same rule. The legend no longer covers the title. The on-screen plot
  no longer cuts long names.
- **New setting in the side panel of every platform, "Sample names in figures"**: Full names (default) or
  shortened to 30, 20 or 12 characters (`HealthyControl_1…RunA`). Tables and downloads always keep full names.
- **Other overlap fixes:**
  - Proteomics and Transcriptomics PCA loading plots move protein and gene names apart, as Metabolomics already
    did.
  - Log-axis minor labels no longer hit the "Before" panel title in the normalization figure.

## 3. Downloads: one standard everywhere
New `figure_export.py`, the same file in Metabolomics/Lipidomics, Proteomics and Transcriptomics.

| | Before | Now (all six platforms) |
|---|---|---|
| Formats | PNG, JPEG, TIFF, SVG, PDF (EPS on some figures) | PNG, TIFF, JPEG, SVG, PDF, EPS, or all six as a ZIP |
| Resolution | 150 / 300 / 600, default 300 (some 300–1200) | 300 / 600 / 1200 DPI, default 600 |
| Tables | mostly CSV only | CSV **and** Excel for every table |

- **TIFF** is a flat RGB file with LZW compression and the DPI tag set.
- **SVG, PDF and EPS** keep editable text.
- **Very large raster exports** are capped at 60 million pixels (server memory); the panel says when that
  happens.
- The panel shows the file name, size in inches and pixel size.
- Fluxomics and Integrix use the same formats and DPI; they keep their width, height and font controls.

## 4. Texts
- **FAQ:**
  - Upgrading no longer promises "you pay only the difference"; plans are changed from the Account page.
  - Extra Lab users: "Contact us to add users ($390 per user per year)".
  - "Without a subscription" now includes Integrix.
- **Tutorials:**
  - The download and sample-name standard is described in all six tutorials.
  - Old DPI and format lists are removed.
  - The Integrix error message is quoted exactly.
- **Home cards:** Proteomics and Transcriptomics list "Pathway enrichment (GSEA, ORA, STRING)".
- **Data preparation:** the lipid and protein annotation examples say how many rows they cover.
- **Checked:**
  - every quoted error message and file name in the tutorials exists in the apps;
  - every example matrix matches its metadata;
  - prices are the same on the site, README, billing doc and email.

## 5. Other fixes
- The demo-mode banner text no longer touches its border inside some platforms.
- **Tests now pass on Streamlit 1.52.2**, the version Render uses: `tests/conftest.py` works around two limits of
  that version's test tool.
- One Lab-team test is skipped on 1.52.2 because its test tool cannot submit that form; the same steps were
  checked in a browser.

## Checks
- **Portal tests:** 74 passed and 1 skipped on Streamlit 1.52.2; 75 passed on 1.67.
- **Platform tests:** Metabolomics 88, Fluxomics 67, Integrix 49 and the Proteomics pipeline: all passed.
- **Long-name test:** every demo was run with long sample names, checking 54 figures. No text overlaps. The
  two remaining flags are tick labels outside the plot area that are not drawn.
- **Browser check** (Streamlit 1.52.2): every page, the Fluxomics and Proteomics demos and the Lab team flow, with
  no errors.

## Files
**Part 1, core:**
- `omnix_portal/`: `config.py`, `portal.py`, `help_content.py`, `demo_files.py`, `account.py`, and all six
  `tutorials/*.md`
- `requirements.txt`, `render_start.sh`
- `tests/`: `test_omnix.py`, `test_integrix_portal.py`, new `conftest.py`

**Part 2, Metabolomics:**
- `app.py`
- `metabolomics_modules/`: new `sample_labels.py`, new `figure_export.py`; changed `qc.py`, `heatmap_module.py`,
  `normalization.py`, `dataset_manager.py`, `utils.py`, `pathway_analysis.py`, `rewiring_ui.py`

**Part 3, Proteomics:**
- `app.py`, `README.md`, `test_pipeline.py`
- `proteomics_modules/`: new `sample_labels.py`, new `figure_export.py`; changed `qc.py`, `heatmap_module.py`,
  `normalization.py`, `pca_module.py`, `utils.py`; **deleted** `ppi.py`

**Part 3, Transcriptomics:**
- `app.py`
- `transcriptomics_modules/`: new `sample_labels.py`, new `figure_export.py`; changed `qc.py`,
  `filtering_rnaseq.py`, `normalization_deseq2.py`, `heatmap_module.py`, `pca_module.py`, `utils.py`

**Part 4, Fluxomics:**
- `core/`: new `sample_labels.py`; changed `style.py`, `heatmap.py`, `plots.py`
- `ui/`: `common.py`, `step4_normalize.py`, `step8_heatmap.py`
- `tests/test_core.py`

**Part 5, Integrix:** `app.py` (60-megapixel export cap).
