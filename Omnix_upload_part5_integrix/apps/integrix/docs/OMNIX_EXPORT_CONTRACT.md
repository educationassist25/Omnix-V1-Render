# Omnix export contract 1.0 — normalized data for Integrix

Each Omnix module (Transcriptomics / BulkRNAAI Pro, Proteomics / ProteoAI Pro, Metabolomics and Lipidomics /
MetaboAI Pro) exports the matrix it already uses for its own statistics. Integrix uses that matrix as provided and
skips every step the module has already done. Its work starts at the per-layer statistics and the integration.

## What a module writes

Either **one zip per module**:

```text
<module>_integrix_export.zip
  matrix.csv      normalized matrix: identifier column(s) first, then one column per sample
  metadata.csv    sample metadata (Sample + the module's grouping columns)
  info.json       the sidecar below
```

or **one combined zip** for the whole study:

```text
transcriptomics.csv   transcriptomics_omnix.json
proteomics.csv        proteomics_omnix.json
metabolomics.csv      metabolomics_omnix.json
metadata.csv
```

## The sidecar

```json
{
  "omnix_export": "1.0",
  "module": "MetaboAI Pro",
  "module_version": "2.3.0",
  "layer": "metabolomics",
  "matrix": "metabolomics.csv",
  "id_columns": ["Metabolite", "HMDB"],
  "scale": "log2",
  "normalization": "PQN, log2",
  "imputation": "none",
  "batch_corrected": false,
  "batch_method": null,
  "filtering": "features missing in > 50% of samples removed",
  "created": "2026-11-05T09:00:00"
}
```

| Key | Required | Meaning | What Integrix does with it |
|---|---|---|---|
| `layer` | yes | transcriptomics, proteomics, metabolomics or lipidomics (lipidomics joins the metabolite layer) | puts the matrix in that layer |
| `scale` | yes | `log2` for normalized exports | no log transform; warns if the values look raw |
| `id_columns` | recommended | identifier columns, first = main ID | used for ID mapping (gene symbol, UniProt, HMDB…) |
| `matrix` | optional | file name of the matrix in the same zip | otherwise `matrix.csv` or the file named after the layer |
| `normalization` | recommended | free text | cited in the Methods |
| `imputation` | recommended | `none`, `kNN`, `half-minimum`… | Integrix imputes only gaps that are left |
| `batch_corrected`, `batch_method` | recommended | whether the module removed batch effects | Integrix's batch correction is switched off for that layer |
| `module`, `module_version`, `created`, `filtering` | recommended | provenance | Methods text and run manifest |

## Rules for the matrix

- Values already normalized, on log2 scale. Leave missing values empty if the module did not impute.
- The same sample names in every module's matrix and metadata (capitals and surrounding spaces are ignored).
- One row per feature; duplicated identifiers are averaged by Integrix.

## In Integrix

In the standalone app, upload normalized tables with **Data → Upload my data → Normalized (log2)**: no log
transform is applied and the Methods state that the tables were supplied normalized. Inside the Omnix platform,
`trace_core/omnix_io.py` reads the zips and sidecars directly (one zip per module or one combined zip), lists each
layer's module, normalization, imputation and batch correction, and cites them in the Methods. The **Omnix
normalized** demo's download is a complete example of this format.
