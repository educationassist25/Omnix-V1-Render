# Omnix: Integrix added (2026-11-06)

Integrix v1.8.3 (tri-omics integration: transcriptomics · proteomics · metabolomics/lipidomics) is now the sixth page
of the Omnix portal. The five data platforms' analysis code is unchanged.

## Upload

| Repository folder | What changed |
|---|---|
| `Omnix_upload_part1_of_3_core/` | portal files below (replace the folder's files with these) |
| `Omnix_upload_part5_integrix/apps/integrix/` | **new**: the Integrix app (engine, knowledge base, demo data, tests) |
| parts 2, 3, 4 | unchanged |

On Render nothing else changes: `render.yaml` and `render_start.sh` stay as they are. The next deploy installs
`mofapy2` and `h5py` from `requirements.txt`.

## Files

**New**
- `omnix_portal/integrix_bridge.py`: the **Send to Integrix** / **Download for Integrix (.zip)** box in the side
  panel of Metabolomics, Lipidomics, Proteomics and Transcriptomics (after Normalization), and the session inbox.
- `omnix_portal/tutorials/integrix.md`: Help → Integrix tutorial.
- `tests/test_integrix_portal.py`: 13 tests (registration, plans, bridge round trip, deep links, demo lock, plan
  unlock, session hand-over with a full run, Proteomics → Send to Integrix).
- `apps/integrix/` (part 5).

**Changed**
- `omnix_portal/config.py`: Integrix platform entry, `DATA_PLATFORMS` / `INTEGRATION_PLATFORMS`, plan feature lines,
  subscription intro, a new FAQ ("What is Integrix, and which plans include it?").
- `omnix_portal/portal.py`: Integrix in the header menu and side menu, the home page's "Integrate across
  platforms" card, the Integrix box under each module's side panel, a Data preparation tab, `OMNIX_STOP` so an
  app can end its page early without hiding the access panel, and a tighter header menu so all items and Login
  stay on one line.
- `omnix_portal/access.py`: a key listing all five data platforms also covers Integrix.
- `omnix_portal/billing.py`, `omnix_portal/account_page.py`: Single-Omics and Duo-Omics choose only among the five
  data platforms; "All five platforms + Integrix" wording.
- `omnix_portal/help_content.py`, `omnix_portal/demo_files.py`: Integrix quick start, file formats, checklist and
  example files.
- `requirements.txt`: `mofapy2>=0.7`, `h5py>=3.9`.
- `tools/make_access_key.py`: help text mentions Integrix.
- `tests/test_omnix.py`: Data preparation download count; the Excel-reader test skips Integrix (it has no Excel
  reader).
- `README.md`: platform table, deep links, plans, layout, tests.

**Integrix app changes for the portal** (`apps/integrix/`)
- `app.py`: portal mode (no `set_page_config`); new data source **From this Omnix session** (the default once two
  or more datasets have been sent); demo users cannot use own data, sent data or exports; ends early with
  `OMNIX_STOP` inside the portal; keeps up to 3 finished analyses on the server in the portal (8 standalone).
- `trace_core/ui.py`: Integrix's card style no longer applies to the portal's header containers.
- `docs/OMNIX_INTEGRATION.md`, `README.md`: describe the integration as built.

## Plans

| Plan | Integrix |
|---|---|
| Single-Omics, Duo-Omics | demo only |
| Multi-Omics, Lab, Enterprise | own data, sent data and exports |

## Checks

- Portal tests: 74 passed (61 existing + 13 new Integrix tests). Integrix's own tests: 49 passed.
- Browser check (1280, 1500 and 1920 px): home, Integrix, Help → Integrix and Data preparation render without
  errors; the header menu stays on one line.
- Memory: one full Integrix demo run (all engines) peaks at about 850 MB in the server process. On Render's
  2 GB Standard plan that allows about two simultaneous Integrix runs; move to a larger plan if several people will
  integrate at the same time.
