# Fluxomics added to Omnix (Render)

## Upload to GitHub
Upload the new **part 4** like the other parts, so that the platforms sit in `apps/`:

```
omnix_app.py, render.yaml, render_start.sh, requirements.txt, omnix_portal/, tools/, tests/, docs/   (part 1)
apps/metabolomics/                                                                                    (part 2)
apps/proteomics/, apps/transcriptomics/                                                               (part 3)
apps/fluxomics/                                                                                       (part 4, NEW)
```
Then replace the changed part-1 files listed below. Render redeploys by itself after the push. `requirements.txt` now
includes `plotly`, which Fluxomics needs.

## Changed files (part 1, core)
| File | Change |
|---|---|
| `omnix_portal/config.py` | New **Fluxomics** platform (menu tab, home card, chips, highlights). Plan texts say "four" instead of "three". Prices are unchanged. |
| `omnix_portal/portal.py` | Fluxomics tab after Transcriptomics in the header, side panel, home (4 cards), footer, Help menu, `?page=fluxomics` and `?page=help-fluxomics`. Plan and compare texts ("all four", "1 of 4"). Data preparation tab. The platform finder accepts a platform's own code folder (`marker`). Tighter header spacing so the extra tab fits on one row. Locked uploaders look disabled. |
| `omnix_portal/access.py`, `keys.py`, `billing.py` | Platform list taken from `config.PLATFORMS`, so subscriptions, access keys, Stripe plans and the Owner console include Fluxomics. Single-Omics can be Fluxomics; Multi-Omics, Lab and Enterprise cover all four. |
| `omnix_portal/account_page.py` | "All four platforms"; the header gap is adjusted. |
| `omnix_portal/help_content.py` | Fluxomics quick start, input-file specifications and its own upload checklist. |
| `omnix_portal/demo_files.py` | Fluxomics example files on Data preparation: generic tables and metadata as Excel workbooks, plus the original Compound Discoverer templates. |
| `omnix_portal/tutorials/fluxomics.md` | **New** tutorial (Help → About Fluxomics). |
| `tools/make_access_key.py` | `--platforms fluxomics` |
| `requirements.txt` | `plotly` added |
| `tests/test_omnix.py` | Fluxomics tests: layout, full workflow, Compound Discoverer demo, data isolation, Reset keeps access, Single-Omics key, help page, Excel examples |
| `README.md`, `docs/DEPLOY_ON_RENDER.md`, `omnix_app.py`, `isolation.py` | Text only |

## Fluxomics (part 4)
`apps/fluxomics` is Omnix Fluxomics 2.5.11. **Only `app.py` (the shell) changed.** `core/`, `ui/`, `data/` and
`example_data/` are untouched, so every calculation, figure and step is unchanged.
- Same layout as Metabolomics. The side panel lists the workflow groups: Setup, Normalization, Isotope Labeling,
  Statistical Analysis, Exploratory Analysis, Pathway Analysis and Export & Report. The page shows the group's steps
  as a row of buttons (the selected step is orange, ✓ = done). Action buttons are teal and download buttons orange.
- The Dataset card, the "Intensity dataset for downstream steps" switch and Session → Reset stay in the side panel.
- Inside Omnix, Fluxomics leaves the page title to the portal. **Reset all steps** clears only Fluxomics, so the
  visitor's Omnix access and page are kept.
- It still runs on its own: `streamlit run apps/fluxomics/app.py`.

## Tests
`python -m pytest tests -q`: 49 passed. `(cd apps/fluxomics && python -m pytest tests -q)`: 67 passed, 1 skipped.
