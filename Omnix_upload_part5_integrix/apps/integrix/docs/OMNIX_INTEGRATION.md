# Integrix inside Omnix

Integrix runs as the sixth page of the Omnix portal (`?page=integrix`) and still runs on its own
(`streamlit run apps/integrix/app.py`). No module's analysis code was changed to add it.

## Layout

```text
Omnix_upload_part5_integrix/
  apps/integrix/
    app.py          the app; inside the portal it runs with OMNIX=True (no set_page_config, portal header)
    trace_core/     the engine
    kb/             knowledge base (Human-GEM, pathways, mappings)
    demo_data/      normalized and raw demo studies
    assets/ docs/ tests/
```

The portal finds it through `config.PLATFORMS["integrix"]` (folder `integrix`, marker `trace_core`,
`"integration": True`). `config.DATA_PLATFORMS` lists the five data platforms and `config.INTEGRATION_PLATFORMS`
lists Integrix, so the home page shows it in its own "Integrate across platforms" row.

## Plans

| Plan | Integrix |
|---|---|
| Single-Omics, Duo-Omics | demo only (Integrix cannot be chosen as one of the plan's platforms) |
| Multi-Omics, Lab, Enterprise | own data |

An access key without a platform list, with `all`, or listing all five data platforms covers Integrix.
Without a covering plan, every upload and the "From this Omnix session" source are locked, as in the other platforms.

## From the modules to Integrix (`omnix_portal/integrix_bridge.py`)

After **Normalization** in Metabolomics, Lipidomics, Proteomics or Transcriptomics, the portal adds an
**Integrix** box to the module's side panel:

- **Send to Integrix** puts the module's export in the session inbox (`omnix_integrix_inbox`; the `omnix_` prefix
  keeps it while the visitor moves between platforms). In Integrix: **Data → From this Omnix session**.
- **Download for Integrix (.zip)** gives the same export as a file, for **Data → Upload Omnix exports** later.

The export follows [`OMNIX_EXPORT_CONTRACT.md`](OMNIX_EXPORT_CONTRACT.md) 1.0: the module's `log2_data` with its
identifier columns (Metabolite/Lipid + HMDB/KEGG/ChEBI from the row annotations; Protein + Gene; Gene), the
metadata with QC injections removed, and a sidecar whose normalization, imputation and batch-correction fields are
taken from the module's processing notes. Integrix uses the matrix as it is and cites those steps in its Methods.

Metabolomics and Lipidomics both fill the metabolite layer; when both are sent, the second is unticked by default.

## Portal-specific behaviour

- The page ends early through `OMNIX_STOP` (raised instead of `st.stop()`), so the portal's access panel still draws.
- Finished analyses are kept on the server: up to 3 inside the portal (8 standalone).
- `mofapy2` and `h5py` are in the portal `requirements.txt`.

## Tests

```bash
python -m pytest tests/test_integrix_portal.py -q      # registration, plans, bridge, session hand-over
(cd apps/integrix && python -m pytest tests -q)         # the engine
```

## Optional knowledge-base add-ons (any internet-connected computer)

```bash
python -m trace_core.kb_build --skip-core --kegg --collectri --string --out apps/integrix/kb
```

Commit the new parquet files and the updated `kb_manifest.json`.
