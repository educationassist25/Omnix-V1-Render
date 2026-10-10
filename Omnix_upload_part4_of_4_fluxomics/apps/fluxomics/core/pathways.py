"""
core.pathways — Step 6: metabolite annotation, pathway assignment and pathway diagrams.

Mapping confidence (a metabolite is assigned to pathways only at High or Medium):
  High    exact name / synonym match AND identical molecular formula (or a user-confirmed ID)
  Medium  match after removing stereo descriptors (D-/L-/DL-/(+)/(−)) with identical formula —
          LC-MS does not resolve enantiomers, so the cellular (usually L-) form is assumed
  Low     possible identity (isomer with the same formula, or a known mis-annotation) — NOT assigned,
          shown as a suggestion for the user to accept or reject
  Class   compound-class rule (e.g. acyl-carnitines -> fatty acid metabolism) — shown as Low, not assigned
          unless the user accepts it
  Unmapped no library match

Pathway drawings (positions, reactions) live in core.pathway_templates; metabolites are matched to them by KEGG ID.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"
STEREO = re.compile(r"^(l|d|dl|r|s|rs|\(\+\)|\(-\)|\(r\)|\(s\)|\(rs\))$", re.I)
CONF_ORDER = {"High": 0, "Medium": 1, "Low": 2, "Class rule": 3, "Unmapped": 4}


def load_library() -> pd.DataFrame:
    L = pd.read_csv(DATA / "metabolite_library.csv", dtype=str).fillna("")
    return L


def load_pathway_names() -> dict:
    P = pd.read_csv(DATA / "pathways.csv", dtype=str)
    return dict(zip(P.Code, P.Pathway))


def load_suggestions() -> dict:
    S = pd.read_csv(DATA / "annotation_suggestions.csv", dtype=str).fillna("")
    return {r.DataName: (r.LibraryID, r.Reason) for r in S.itertuples()}


def _norm(name: str, strip_stereo: bool = False) -> str:
    s = str(name).lower().replace("′", "'").replace("’", "'")
    toks = [t for t in re.split(r"[\s,\-_]+", s) if t]
    if strip_stereo:
        toks = [t for t in toks if not STEREO.match(t)]
        toks = [re.sub(r"^\(?[+-]\)?", "", t) for t in toks]
    return " ".join(t for t in toks if t)


def _formula_key(f: str) -> str:
    d = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", str(f).replace(" ", "")):
        d[el] = d.get(el, 0) + int(n or 1)
    return "".join(f"{k}{d[k]}" for k in sorted(d))


CLASS_RULES = [  # (regex on data name, pathway code, description)
    (r"carnitine", "FA", "acyl-carnitine (carnitine shuttle / β-oxidation)"),
    (r"^(gamma|γ)-glutamyl", "RED", "γ-glutamyl dipeptide (γ-glutamyl cycle)"),
    (r"glutamyl.*glutam|glutamylcysteinyl", "RED", "γ-glutamyl peptide"),
    (r"(alanyl|glycyl|leucyl|isoleucyl|valyl|aspartyl|lysyl|glutaminyl|prolyl|methionyl|asparaginyl|seryl|"
     r"threonyl|phenylalanyl|tyrosyl|tryptophyl|histidyl)\w*|^(ala|gly|leu|ile|val|asp|lys|gln|glu|pro|met|asn|ser|"
     r"thr|phe|tyr|trp|his|arg|cys)-(ala|gly|leu|ile|val|asp|lys|gln|glu|pro|met|asn|ser|thr|phe|tyr|trp|his|arg|cys)$",
     "", "dipeptide (protein turnover; no pathway assigned)"),
    (r"glycero.*phosph|phosph.*glycer", "PL", "glycerophospho- compound (phospholipid turnover)"),
]


def map_metabolites(data: pd.DataFrame, library: pd.DataFrame | None = None,
                    overrides: pd.DataFrame | None = None) -> pd.DataFrame:
    """
    data: Compound, Formula, Method.  Returns the annotation table (one row per data metabolite).
    overrides: rows with Compound + any of LibraryID / KEGG / HMDB / PubChemCID / ChEBI / Pathways / Confidence,
               applied last (Confidence 'User-confirmed').
    """
    L = library if library is not None else load_library()
    pw = load_pathway_names()
    sugg = load_suggestions()
    L = L.assign(_f=L.Formula.map(_formula_key))
    exact, loose = {}, {}
    for _, r in L.iterrows():
        for nm in [r.StandardName] + [s for s in r.Synonyms.split("|") if s]:
            exact.setdefault(_norm(nm), []).append(r.LibraryID)
            loose.setdefault(_norm(nm, True), []).append(r.LibraryID)
    by_id = L.set_index("LibraryID")
    rows = []
    for _, d in data.iterrows():
        fk = _formula_key(d.Formula)
        conf, lid, status, note = "Unmapped", "", "Not in library", ""
        cand = [i for i in exact.get(_norm(d.Compound), []) if by_id.loc[i, "_f"] == fk]
        if cand:
            conf, lid, status = "High", cand[0], "Exact name + formula"
        else:
            cand = [i for i in loose.get(_norm(d.Compound, True), []) if by_id.loc[i, "_f"] == fk]
            if cand:
                conf, lid, status = "Medium", cand[0], "Name match ignoring stereo prefix (D/L not resolved by LC-MS)"
            elif d.Compound in sugg:
                conf, lid = "Low", sugg[d.Compound][0]
                status, note = "Suggested (isomer / probable mis-annotation) — not assigned", sugg[d.Compound][1]
            else:
                iso = L.loc[L._f == fk, "LibraryID"].tolist()
                rule = next(((code, desc) for rx, code, desc in CLASS_RULES if re.search(rx, d.Compound, re.I)), None)
                if exact.get(_norm(d.Compound)) or loose.get(_norm(d.Compound, True)):
                    conf, status = "Low", "Name matches but formula differs — not assigned"
                elif rule:
                    conf, status, note = "Class rule", f"Compound class: {rule[1]}", rule[0]
                elif iso:
                    conf, status = "Low", "Same formula as " + ", ".join(by_id.loc[i, "StandardName"] for i in iso[:3]) + \
                                   " — different name, not assigned"
        rec = {"Compound": d.Compound, "Method": d.get("Method", ""), "Formula": d.Formula, "LibraryID": lid,
               "StandardName": "", "Label": "", "KEGG": "", "HMDB": "", "PubChemCID": "", "ChEBI": "",
               "Pathways": "", "Confidence": conf, "MappingStatus": status, "Note": note, "IDSource": ""}
        rec["Candidate"] = ""
        if lid:
            e = by_id.loc[lid]
            rec["CandidatePathways"] = "; ".join(pw.get(c, c) for c in e.Pathways.split("|") if c)
            if conf in ("High", "Medium"):
                rec.update(StandardName=e.StandardName, Label=e.Label, KEGG=e.KEGG, HMDB=e.HMDB,
                           PubChemCID=e.PubChemCID, ChEBI=e.ChEBI, IDSource=e.IDSource, Note=(note or e.Note))
                rec["Pathways"] = rec["CandidatePathways"]
            else:   # suggestion only: identifiers are NOT attached until the user accepts
                rec["Candidate"] = e.StandardName
                rec["Note"] = note or e.Note
        elif conf == "Class rule" and note:
            rec["CandidatePathways"] = pw.get(note, note)
            rec["Note"] = ""
        else:
            rec["CandidatePathways"] = ""
        rows.append(rec)
    A = pd.DataFrame(rows)
    if overrides is not None and len(overrides):
        A = apply_overrides(A, overrides)
    A["_o"] = A.Confidence.map(CONF_ORDER).fillna(9)
    return A.sort_values(["_o", "Compound"]).drop(columns="_o").reset_index(drop=True)


def apply_overrides(A: pd.DataFrame, O: pd.DataFrame) -> pd.DataFrame:
    A = A.set_index("Compound")
    for _, o in O.iterrows():
        if o.Compound not in A.index:
            continue
        for c in ["StandardName", "Label", "KEGG", "HMDB", "PubChemCID", "ChEBI", "Pathways", "LibraryID"]:
            if c in o and isinstance(o[c], str) and o[c].strip() != "":
                A.loc[o.Compound, c] = o[c].strip()
        A.loc[o.Compound, "Confidence"] = "User-confirmed"
        A.loc[o.Compound, "MappingStatus"] = "Edited by user"
    return A.reset_index()


def accept_suggestions(A: pd.DataFrame, compounds: list, library: pd.DataFrame | None = None) -> pd.DataFrame:
    """User accepts Low / Class-rule candidates: identifiers and candidate pathways become assigned."""
    A = A.copy()
    L = (library if library is not None else load_library()).set_index("LibraryID")
    m = A.Compound.isin(compounds)
    for i in A.index[m]:
        lid = A.at[i, "LibraryID"]
        if lid and lid in L.index:
            e = L.loc[lid]
            for c in ["StandardName", "Label", "KEGG", "HMDB", "PubChemCID", "ChEBI", "IDSource"]:
                A.at[i, c] = e[c]
    A.loc[m, "Pathways"] = A.loc[m, "CandidatePathways"]
    A.loc[m, "Confidence"] = "User-confirmed"
    A.loc[m, "MappingStatus"] = A.loc[m, "MappingStatus"].str.replace(" — not assigned", "") + " (accepted by user)"
    return A


def pathway_members(A: pd.DataFrame) -> dict:
    """{pathway name: [compounds]} for assigned metabolites only."""
    out = {}
    for _, r in A.iterrows():
        for p in [x.strip() for x in str(r.Pathways).split(";") if x.strip()]:
            out.setdefault(p, []).append(r.Compound)
    return out


# ---------------------------------------------------------------------------
# Online verification (runs on the user's computer; needs internet)
# ---------------------------------------------------------------------------
def online_available(timeout: float = 3.0) -> bool:
    """Quick reachability check for KEGG REST (used before automatic verification)."""
    try:
        import requests
        return requests.get("https://rest.kegg.jp/info/kegg", timeout=timeout).ok
    except Exception:
        return False


def verify_online(A: pd.DataFrame, cache_file: Path | None = None, timeout: float = 8.0, progress=None,
                  compounds=None) -> pd.DataFrame:
    """
    For each mapped metabolite: KEGG REST (name, formula, ChEBI/PubChem links) and PubChem PUG REST
    (CID by name, HMDB/ChEBI/KEGG from synonyms). Fills blanks, flags mismatches. Successful lookups are cached
    (failed ones are retried next time). compounds: limit the check to these data compound names.
    """
    import requests
    cache_file = cache_file or (DATA / "id_cache.json")
    cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}
    A = A.copy()
    A["Verification"] = A.get("Verification", "")
    sel = A.LibraryID.astype(str) != ""
    if compounds is not None:
        sel &= A.Compound.isin(set(compounds))
    todo = A[sel].index.tolist()
    for n, i in enumerate(todo):
        r = A.loc[i]
        key = r.StandardName or r.Compound
        if progress:
            progress((n + 1) / max(len(todo), 1), key)
        info = cache.get(key)
        if info is None:
            info = {}
            try:
                if r.KEGG:
                    t = requests.get(f"https://rest.kegg.jp/get/{r.KEGG}", timeout=timeout).text
                    m = re.search(r"^FORMULA\s+(\S+)", t, re.M)
                    info["kegg_formula"] = m.group(1) if m else ""
                    m = re.search(r"ChEBI:\s+(\d+)", t)
                    info["kegg_chebi"] = m.group(1) if m else ""
                    m = re.search(r"^NAME\s+(.+?);?$", t, re.M)
                    info["kegg_name"] = m.group(1) if m else ""
                q = requests.utils.quote(key)
                cid = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/cids/TXT",
                                   timeout=timeout).text.split()
                cid = cid[0] if cid and cid[0].isdigit() else (r.PubChemCID or "")
                info["cid"] = cid
                if cid:
                    syn = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON",
                                       timeout=timeout).json()
                    s = syn["InformationList"]["Information"][0].get("Synonym", [])
                    info["hmdb"] = next((x for x in s if re.fullmatch(r"HMDB\d{7}", x)), "")
                    info["chebi"] = next((x.split(":")[1] for x in s if re.fullmatch(r"CHEBI:\d+", x)), "")
                    info["kegg"] = next((x for x in s if re.fullmatch(r"C\d{5}", x)), "")
                time.sleep(0.25)
            except Exception as e:
                info["error"] = f"{type(e).__name__}"
            if not info.get("error"):
                cache[key] = info
        notes = []
        if info.get("error"):
            notes.append("lookup failed (" + info["error"] + ")")
        else:
            if info.get("kegg_formula") and _formula_key(info["kegg_formula"]) != _formula_key(r.Formula):
                notes.append(f"KEGG formula {info['kegg_formula']} ≠ data formula")
            for col, k in (("PubChemCID", "cid"), ("HMDB", "hmdb"), ("ChEBI", "chebi")):
                v = info.get(k, "")
                if v and not r[col]:
                    A.loc[i, col] = v; notes.append(f"{col} filled")
                elif v and r[col] and str(v) != str(r[col]):
                    notes.append(f"{col} differs online ({v})")
            if info.get("kegg") and r.KEGG and info["kegg"] != r.KEGG:
                notes.append(f"KEGG differs online ({info['kegg']})")
        A.loc[i, "Verification"] = "; ".join(notes) or "verified"
        A.loc[i, "IDSource"] = "verified online (KEGG/PubChem)" if not notes or all("filled" in x for x in notes) \
            else "checked online — see Verification"
    try:
        cache_file.write_text(json.dumps(cache, indent=1))
    except OSError:
        pass
    return A
