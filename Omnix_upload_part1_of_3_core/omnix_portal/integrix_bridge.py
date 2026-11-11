"""
Integrix bridge: hands a module's normalized data to Omnix Integrix (tri-omics integration).

The four quantitative modules (Metabolomics, Lipidomics, Proteomics, Transcriptomics) keep their final normalized
matrix in the session as ``log2_data`` (features x samples), the sample metadata as ``meta`` (indexed by sample)
and their steps in ``processing_notes``. After a module has run, the portal draws a small Integrix panel in its side
panel (see ``portal.run_app``) with:

* **Send to Integrix**: puts the module's export in this session's Integrix inbox
  (``st.session_state["omnix_integrix_inbox"]``; the ``omnix_`` prefix keeps it when the visitor switches apps),
  so Integrix can use it at once (Data > From this Omnix session).
* **Download for Integrix**: the same export as a zip, to keep or to upload later (Data > Upload Omnix exports).

The export follows the Omnix export contract 1.0 (apps/integrix/docs/OMNIX_EXPORT_CONTRACT.md): the matrix with
its identifier columns, the metadata (QC injections removed) and a JSON sidecar saying how the matrix was made.
No module's analysis code is changed: the panel only reads what the module already keeps in the session.
"""
from __future__ import annotations

import datetime as dt
import io
import json
import zipfile

import pandas as pd
import streamlit as st

INBOX = "omnix_integrix_inbox"
MODULES = {   # platform key -> (contract layer, module name for the sidecar and Methods)
    "metabolomics": ("metabolomics", "Omnix Metabolomics"),
    "lipidomics": ("lipidomics", "Omnix Lipidomics"),
    "proteomics": ("proteomics", "Omnix Proteomics"),
    "transcriptomics": ("transcriptomics", "Omnix Transcriptomics"),
}
_ID_COLS = {"hmdb": "HMDB", "hmdb_id": "HMDB", "kegg": "KEGG", "kegg_id": "KEGG", "chebi": "ChEBI", "chebi_id": "ChEBI",
            "pubchem": "PubChem", "pubchemcid": "PubChem", "uniprot": "UniProt", "gene": "Gene", "gene_symbol": "Gene"}


def _truthy(v) -> bool:
    return str(v).strip().lower() in ("true", "1", "yes", "y", "qc")


def _notes(state) -> list:
    n = state.get("processing_notes") or []
    return [str(x) for x in n] if isinstance(n, (list, tuple)) else [str(n)]


def build(platform: str, state=None):
    """(matrix with identifier columns first, metadata without QC injections, sidecar dict) from the module's
    session, or None when the module has no normalized matrix with at least 3 samples yet."""
    state = st.session_state if state is None else state
    X = state.get("log2_data")
    meta = state.get("meta")
    if not isinstance(X, pd.DataFrame) or X.empty or not isinstance(meta, pd.DataFrame):
        return None
    layer, module = MODULES[platform]
    meta = meta.copy()
    if "Sample" in meta.columns:
        meta = meta.set_index("Sample")
    meta.index = meta.index.astype(str)
    samples = [str(c) for c in X.columns if str(c) in meta.index]
    if "IsQC" in meta.columns:
        samples = [s for s in samples if not _truthy(meta.at[s, "IsQC"])]
    if len(samples) < 3:
        return None
    X = X.copy()
    X.columns = [str(c) for c in X.columns]
    X = X[samples]
    ids = pd.DataFrame(index=X.index)
    if layer in ("metabolomics", "lipidomics"):
        main = "Lipid" if layer == "lipidomics" else "Metabolite"
        ids[main] = X.index.astype(str)
        ann = state.get("row_annotations")
        if isinstance(ann, pd.DataFrame):
            for c in ann.columns:
                std = _ID_COLS.get(str(c).strip().lower())
                if std and std not in ids.columns:
                    ids[std] = ann[c].reindex(X.index).values
    elif layer == "proteomics":
        ids["Protein"] = X.index.astype(str)
        gm = state.get("gene_map")
        if isinstance(gm, pd.Series):
            gm = gm.to_dict()
        if isinstance(gm, dict) and gm:
            ids["Gene"] = [gm.get(p, "") for p in X.index]
    else:
        ids["Gene"] = X.index.astype(str)
    matrix = pd.concat([ids.reset_index(drop=True), X.reset_index(drop=True)], axis=1)
    md = meta.loc[samples].drop(columns=[c for c in ("IsQC",) if c in meta.columns])
    md = md.reset_index().rename(columns={md.index.name or "index": "Sample"})
    notes = _notes(state)
    low = [n.lower() for n in notes]
    norm = "; ".join(n for n, l in zip(notes, low) if any(w in l for w in ("normaliz", "log2", "vst", "cpm", "istd",
                                                                           "iqr", "pqn", "median"))) or "see module notes"
    imp = "; ".join(n for n, l in zip(notes, low) if "imput" in l) or "none"
    batch = any(("combat" in l or "batch correct" in l or "batch-correct" in l or "removebatch" in l) for l in low)
    info = {"omnix_export": "1.0", "module": module, "module_version": "", "layer": layer,
            "matrix": f"{layer}.csv", "id_columns": list(ids.columns), "scale": "log2",
            "normalization": norm[:600], "imputation": imp[:300], "batch_corrected": batch, "batch_method": None,
            "filtering": "", "created": dt.datetime.now().isoformat(timespec="seconds")}
    info["processing_notes"] = notes[-12:]
    return matrix, md, info


def to_zip(matrix: pd.DataFrame, meta: pd.DataFrame, info: dict) -> bytes:
    buf = io.BytesIO()
    layer = info["layer"]
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{layer}.csv", matrix.to_csv(index=False))
        z.writestr("metadata.csv", meta.to_csv(index=False))
        z.writestr(f"{layer}_omnix.json", json.dumps(info, indent=2, default=str))
        z.writestr("README.txt",
                   f"Omnix export for Integrix (contract {info.get('omnix_export', '1.0')})\n"
                   f"{info['module']}: {matrix.shape[0]} features x {meta.shape[0]} samples, normalized, log2.\n"
                   "Use it in Omnix Integrix: Data > Upload Omnix exports (drop this zip and the other modules' zips).\n")
    return buf.getvalue()


def send(platform: str):
    out = build(platform)
    if out is None:
        st.session_state["omnix_integrix_msg"] = (platform, "error")
        return
    matrix, meta, info = out
    inbox = st.session_state.get(INBOX) or {}
    other = {"metabolomics": "lipidomics", "lipidomics": "metabolomics"}.get(platform)
    inbox[platform] = {
        "platform_name": MODULES[platform][1].replace("Omnix ", ""), "layer": info["layer"],
        "zip": to_zip(matrix, meta, info), "zip_name": f"{platform}_integrix_export.zip",
        "n_features": int(matrix.shape[0]), "n_samples": int(meta.shape[0]),
        "sent": dt.datetime.now().strftime("%H:%M"),
        "default_use": not (other in inbox and inbox[other].get("default_use", True)),
    }
    st.session_state[INBOX] = inbox
    st.session_state["omnix_integrix_msg"] = (platform, "sent")


def sidebar_panel(platform: str, go=None):
    """The Integrix box in a module's side panel."""
    if platform not in MODULES:
        return
    X = st.session_state.get("log2_data")
    demo_only = st.session_state.get("omnix_demo_only", False)
    with st.sidebar.container(key="omnix_integrix_panel", border=True):
        st.markdown("<div class='omx-ix-t'>Integrix · tri-omics integration</div>", unsafe_allow_html=True)
        if not isinstance(X, pd.DataFrame) or X.empty:
            st.caption("After **Normalization**, send this dataset to Integrix to integrate it with the other "
                       "omics layers of the same samples.")
            return
        if demo_only:
            st.caption("Sending your own normalized data to Integrix is part of a subscription.")
            return
        inbox = st.session_state.get(INBOX) or {}
        st.button("Send to Integrix", key="omnix_ix_send", on_click=send, args=(platform,), icon=":material/send:",
                  use_container_width=True, type="primary",
                  help="Hands this module's normalized matrix, sample metadata and processing notes to Integrix in "
                       "this session. Integrix uses the matrix as it is (no second normalization).")
        state = st.session_state
        snap = {k: state.get(k) for k in ("log2_data", "meta", "processing_notes", "row_annotations", "gene_map")}

        def _zip():
            out = build(platform, snap)
            return to_zip(*out) if out else b""
        st.download_button("Download for Integrix (.zip)", _zip, file_name=f"{platform}_integrix_export.zip",
                           mime="application/zip", key="omnix_ix_dl", icon=":material/download:",
                           use_container_width=True,
                           help="The same export as a zip (matrix, metadata, how it was made), to upload later in "
                                "Integrix: Data > Upload Omnix exports.")
        msg = st.session_state.get("omnix_integrix_msg")
        if msg and msg[0] == platform:
            if msg[1] == "error":
                st.caption(":material/error: No normalized matrix with at least 3 non-QC samples to send.")
        if inbox:
            names = ", ".join(v["platform_name"] for v in inbox.values())
            st.caption(f":material/check_circle: In Integrix now: {names}.")
            if go is not None:
                st.button("Open Integrix  →", key="omnix_ix_open", on_click=go, args=("integrix",),
                          use_container_width=True)
