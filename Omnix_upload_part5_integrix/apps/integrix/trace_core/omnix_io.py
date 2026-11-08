"""omnix_io.py - read normalized exports from the Omnix modules (the "Omnix export contract").

Each Omnix module (transcriptomics, proteomics, metabolomics / lipidomics) exports the matrix it already uses for
its own statistics, plus a small JSON sidecar that says how that matrix was made. Integrix then skips every step
the module has already done (log2 transform, normalization, imputation, batch correction) and starts at the
per-layer statistics and the integration.

A sidecar (``<layer>_omnix.json``, or ``info.json`` inside a single-module zip)::

    {
      "omnix_export": "1.0",
      "module": "MetaboAI Pro",            "module_version": "2.3.0",
      "layer": "metabolomics",             # transcriptomics | proteomics | metabolomics | lipidomics
      "matrix": "metabolomics.csv",        # optional: the matrix file in the same zip
      "id_columns": ["Metabolite", "HMDB"],
      "scale": "log2",                     # log2 | raw
      "normalization": "PQN",              # free text, cited in the Methods
      "imputation": "none",                # none | kNN | half-minimum | ...
      "batch_corrected": false,  "batch_method": null,
      "filtering": "features missing in > 50% of samples removed",
      "created": "2026-11-05T09:00:00"
    }

A zip may hold one module's export (``matrix.csv`` + ``metadata.csv`` + ``info.json``) or a combined bundle
(``transcriptomics.csv``, ``proteomics.csv``, ``metabolomics.csv``, ``metadata.csv`` and one sidecar per layer).
Several zips can be uploaded together.
"""
from __future__ import annotations

import io
import json
import zipfile

import numpy as np
import pandas as pd

CONTRACT_VERSION = "1.0"
LAYER_KEY = {"transcriptomics": "rna", "rna": "rna", "rna-seq": "rna", "mrna": "rna",
             "proteomics": "prot", "protein": "prot", "prot": "prot",
             "metabolomics": "met", "lipidomics": "met", "metabolite": "met", "met": "met"}
LAYER_FILE = {"rna": "transcriptomics", "prot": "proteomics", "met": "metabolomics"}
REQUIRED = ("layer", "scale")


def _read_csv(name: str, data: bytes) -> pd.DataFrame:
    sep = "\t" if name.lower().endswith((".tsv", ".txt")) else ","
    if name.lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(io.BytesIO(data))
    return pd.read_csv(io.BytesIO(data), sep=sep)


def _guess_layer(name: str) -> str | None:
    n = name.lower().rsplit("/", 1)[-1]
    if "metadata" in n or n.startswith("meta") or "sample" in n:
        return None
    for word, k in (("transcript", "rna"), ("rna", "rna"), ("gene", "rna"), ("prote", "prot"),
                    ("metabol", "met"), ("lipid", "met")):
        if word in n:
            return k
    return None


def looks_log(df: pd.DataFrame) -> bool:
    """True when the numeric values look like log-scale data (max below 100, as for log2 intensities)."""
    num = df.select_dtypes("number").values
    v = num[np.isfinite(num)]
    return bool(v.size and np.nanmax(v) < 100)


def check_info(info: dict) -> list[str]:
    """Problems with a sidecar: missing required keys or an unknown layer."""
    out = [f"missing '{k}'" for k in REQUIRED if k not in info]
    if "layer" in info and str(info["layer"]).lower() not in LAYER_KEY:
        out.append(f"unknown layer '{info['layer']}'")
    return out


def describe(info: dict) -> str:
    """One line for the screen and the Methods: module, normalization, imputation, batch correction."""
    if info.get("declared"):
        return "declared normalized (log2) by the user"
    ver = f" {info['module_version']}" if info.get("module_version") else ""
    how = str(info.get("normalization") or "normalized")
    scale = str(info.get("scale", "log2"))
    if scale.lower() not in how.lower():
        how += f", {scale}"
    imp = str(info.get("imputation") or "none")
    bits = [how, "no imputation" if imp.lower() == "none" else f"imputed ({imp})"]
    if info.get("batch_corrected"):
        bits.append(f"batch-corrected ({info.get('batch_method') or 'yes'})")
    return f"{info.get('module', 'Omnix module')}{ver}: " + "; ".join(bits)


def read_exports(files: list[tuple[str, bytes]]) -> tuple[dict, pd.DataFrame | None, dict, list[str]]:
    """Read one or more Omnix export zips (or loose files). Returns (raw tables by layer, metadata,
    sidecars by layer, messages). Messages explain anything that was guessed or could not be used."""
    entries: list[tuple[str, bytes, str]] = []           # (name, bytes, group id = source zip)
    for fname, data in files:
        if fname.lower().endswith(".zip"):
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                for n in z.namelist():
                    if not n.endswith("/") and not n.split("/")[-1].startswith((".", "__MACOSX")) \
                            and "__MACOSX" not in n:
                        entries.append((n, z.read(n), fname))
        else:
            entries.append((fname, data, fname))
    raw, infos, msgs, metas = {}, {}, [], []
    sidecars = [(n, d, g) for n, d, g in entries if n.lower().endswith(".json") and "ground_truth" not in n.lower()]
    tables = [(n, d, g) for n, d, g in entries if n.lower().endswith((".csv", ".tsv", ".txt", ".xlsx"))
              and not n.lower().endswith("readme.txt")]
    for n, d, g in tables:
        if "metadata" in n.lower().rsplit("/", 1)[-1] or n.lower().rsplit("/", 1)[-1].startswith("samples"):
            metas.append((n, _read_csv(n, d)))
    used = set()
    for n, d, g in sidecars:
        try:
            info = json.loads(d.decode("utf-8"))
        except Exception:
            msgs.append(f"{n}: not a valid JSON sidecar, ignored")
            continue
        if not isinstance(info, dict) or "layer" not in info:
            continue
        bad = check_info(info)
        if bad:
            msgs.append(f"{n}: " + "; ".join(bad))
            continue
        k = LAYER_KEY[str(info["layer"]).lower()]
        want = info.get("matrix")
        cand = [t for t in tables if t[2] == g and t[0] not in [m[0] for m in metas]]
        pick = None
        if want:
            pick = next((t for t in cand if t[0].rsplit("/", 1)[-1] == want), None)
        if pick is None:
            pick = next((t for t in cand if _guess_layer(t[0]) == k), None)
        if pick is None:
            pick = next((t for t in cand if t[0].rsplit("/", 1)[-1].lower().startswith("matrix")), None)
        if pick is None:
            msgs.append(f"{n}: no matrix file found for {LAYER_FILE[k]}")
            continue
        if k in raw:
            msgs.append(f"{LAYER_FILE[k]}: more than one export; the first was used")
            continue
        raw[k] = _read_csv(pick[0], pick[1])
        info = dict(info, matrix=pick[0].rsplit("/", 1)[-1], source=g)
        infos[k] = info
        used.add(pick[0])
        if str(info.get("scale", "")).lower() == "log2" and not looks_log(raw[k]):
            msgs.append(f"{LAYER_FILE[k]}: the sidecar says log2 but the values look like raw intensities "
                        f"(maximum ≥ 100). Check the export.")
    # tables without a sidecar: recognised by name, treated as not normalized
    for n, d, g in tables:
        k = _guess_layer(n)
        if k and n not in used and k not in raw:
            raw[k] = _read_csv(n, d)
            msgs.append(f"{n}: no Omnix sidecar, so it is treated as raw data (Integrix normalizes it)")
    meta = metas[0][1] if metas else None
    if len(metas) > 1:
        ids = [set(m[1].iloc[:, 0].astype(str).str.strip().str.lower()) for m in metas]
        if any(i != ids[0] for i in ids[1:]):
            msgs.append("The metadata files list different samples; the first one was used: " + metas[0][0])
    if meta is None and raw:
        msgs.append("No metadata.csv found in the exports: add the sample metadata.")
    bc = {k: bool(v.get("batch_corrected")) for k, v in infos.items()}
    if len(set(bc.values())) > 1:
        msgs.append("Batch correction differs between layers: "
                    + ", ".join(f"{LAYER_FILE[k]} {'corrected' if v else 'not corrected'}" for k, v in bc.items())
                    + ". Integrix can correct the others (Study design).")
    return raw, meta, infos, msgs


def sidecar(layer: str, module: str, normalization: str, id_columns: list, matrix: str, imputation="none",
            batch_corrected=False, batch_method=None, filtering="", module_version="", created="") -> dict:
    return {"omnix_export": CONTRACT_VERSION, "module": module, "module_version": module_version, "layer": layer,
            "matrix": matrix, "id_columns": id_columns, "scale": "log2", "normalization": normalization,
            "imputation": imputation, "batch_corrected": batch_corrected, "batch_method": batch_method,
            "filtering": filtering, "created": created}
