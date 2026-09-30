"""
rewiring_reactions.py - Reaction and pathway knowledge for the Correlation Rewiring Map.

Two jobs:

1. Map measured features to metabolites and pathways, using the SAME identifier
   standardization and metabolite-set libraries as the Pathway Analysis (MSEA) page
   (metabolomics_modules/pathway_libraries.py + pathway_analysis.standardize_for_option_a), so a
   metabolite sits in the same KEGG / SMPDB pathway on both pages.

2. Link measured metabolites through enzymatic reactions, using the Human-GEM reaction-
   pair library built by build_reaction_library.py (metabolomics_modules/data/humangem_*.parquet):
     * direct (1 step): one reaction interconverts the two metabolites;
     * 2 steps: A -> X -> B through ONE unmeasured intermediate X, following reaction
       directions (a route, not two co-products of one hydrolysis), with X not a hub
       (route degree <= hub_cap; acetyl-CoA, carnitine, UDP-glucuronate, ... excluded), not a
       peptide, not reached through peptide-hydrolysis / protein-degradation reactions, and with
       the two steps coming from different reactions (co-products of one hydrolysis are not a route).
   The links are rebuilt for whichever metabolites a dataset measures.

Everything here is data lookup; statistics stay in rewiring.py.
"""

from __future__ import annotations

import collections
import json
import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
PAIRS_PATH = os.path.join(_DATA, "humangem_reaction_pairs.parquet")
METS_PATH = os.path.join(_DATA, "humangem_metabolites.parquet")
INFO_PATH = os.path.join(_DATA, "humangem_library_info.json")

XWALK_PATH = os.path.join(_DATA, "hmdb_kegg_crosswalk.parquet")
HMDB_NAMES_PATH = os.path.join(_DATA, "hmdb_metabolite_gene_mapping.parquet")

ROUTE_EXCLUDED_SUBSYSTEMS = ("Peptide metabolism", "Protein degradation")

# Abbreviations and report spellings that neither the MSEA reference nor HMDB's official
# names cover. Values are HMDB accessions (checked against HMDB).
NAME_SYNONYMS = {
    "gaba": "HMDB0000112", "gammaaminobutyrate": "HMDB0000112", "4aminobutyrate": "HMDB0000112",
    "4aminobutanoate": "HMDB0000112",
    "gsh": "HMDB0000125", "reducedglutathione": "HMDB0000125",
    "gssg": "HMDB0003337", "glutathionedisulfide": "HMDB0003337", "oxidizedglutathione": "HMDB0003337",
    "sam": "HMDB0001185", "adomet": "HMDB0001185", "sadenosyllmethionine": "HMDB0001185",
    "sah": "HMDB0000939", "adohcy": "HMDB0000939", "sadenosyllhomocysteine": "HMDB0000939",
    "akg": "HMDB0000208", "2og": "HMDB0000208", "2oxoglutarate": "HMDB0000208",
    "oaa": "HMDB0000223", "pep": "HMDB0000263", "dhap": "HMDB0001473", "g6p": "HMDB0001401",
    "f6p": "HMDB0000124", "fbp": "HMDB0001058", "f16bp": "HMDB0001058", "3pg": "HMDB0000807",
    "2pg": "HMDB0000362", "g3p": "HMDB0001112",
    "5htp": "HMDB0000472", "5hiaa": "HMDB0000763", "5hydroxyindoleacetate": "HMDB0000763",
    "nfk": "HMDB0001200", "nformylkynurenine": "HMDB0001200", "formylkynurenine": "HMDB0001200",
    "lformylkynurenine": "HMDB0001200",
    "3hk": "HMDB0000732", "3hydroxylkynurenine": "HMDB0000732",
    "3haa": "HMDB0001476", "3hydroxyanthranilate": "HMDB0001476",
    "kyna": "HMDB0000715", "kynurenate": "HMDB0000715", "quin": "HMDB0000232", "quinolinicacid": "HMDB0000232",
    "anthranilate": "HMDB0001123", "xanthurenate": "HMDB0000881",
    "urocanate": "HMDB0000301", "transurocanate": "HMDB0000301",
    "nag": "HMDB0001138", "nacetylglutamate": "HMDB0001138", "nacetyllglutamate": "HMDB0001138",
    "6phosphogluconolactone": "HMDB0001127", "6phosphodglucono15lactone": "HMDB0001127",
    "3methyl2oxobutanoate": "HMDB0000019", "alphaketoisovalerate": "HMDB0000019", "ketoisovalerate": "HMDB0000019",
    "4methyl2oxopentanoate": "HMDB0000695", "alphaketoisocaproate": "HMDB0000695", "ketoisocaproate": "HMDB0000695",
    "cystathionine": "HMDB0000099", "lcystathionine": "HMDB0000099",
    "uricacid": "HMDB0000289", "urate": "HMDB0000289",
}
DEFAULT_HUB_CAP = 45

# KEGG draws its maps with anomer-specific compound IDs (glycolysis uses beta-D-fructose-1,6-bisphosphate
# C05378, alpha/beta-D-glucose-6-phosphate C00668/C01172, ...). A metabolomics measurement reports the
# generic compound (C00354, C00092), which would then miss glycolysis entirely. Membership in a KEGG
# set is therefore tested for the whole anomer group. (Checked against the bundled KEGG library.)
KEGG_ANOMER_GROUPS = [
    {"C00031", "C00267", "C00221"},   # D-glucose, alpha-, beta-
    {"C00092", "C00668", "C01172"},   # D-glucose 6-phosphate, alpha-, beta-
    {"C00085", "C05345"},             # D-fructose 6-phosphate, beta-
    {"C00354", "C05378"},             # D-fructose 1,6-bisphosphate, beta-
    {"C00095", "C02336"},             # D-fructose, beta-
    {"C00124", "C00984", "C00962"},   # D-galactose, alpha-, beta-
    {"C00159", "C00936"},             # D-mannose, alpha-
    {"C00103", "C00663"},             # D-glucose 1-phosphate, beta-
]
_ANOMER = {k: g for g in KEGG_ANOMER_GROUPS for k in g}

LIB_ANNOTATION = "annotation"          # a column of the row-annotation table
LIB_CLUSTER = "cluster"                # correlation clusters
LIB_GEM = "Human-GEM subsystems"       # subsystems of the reaction library itself


# ---------------------------------------------------------------------------
# Name normalization (shared with build_reaction_library.py)
# ---------------------------------------------------------------------------
def name_keys(name) -> list:
    """Normalized keys for a metabolite name: as given, without stereo prefix, '-ic acid' <->
    '-ate', and without / only a trailing parenthetical ('Glutathione (GSH)')."""
    import re
    text = str(name or "").strip()
    if not text:
        return []
    forms = {text}
    m = re.match(r"^(.*?)\s*\(([^()]+)\)\s*$", text)
    if m:
        forms.update({m.group(1), m.group(2)})
    for f in list(forms):
        forms.add(re.sub(r"^(l|d|dl|sn)-", "", f, flags=re.I))
    for f in list(forms):   # internal stereo marker: 5-Hydroxy-L-tryptophan -> 5-Hydroxytryptophan
        forms.add(re.sub(r"(?<=[a-z0-9])-(l|d)-(?=[a-z])", "", f, flags=re.I))
    for f in list(forms):
        if re.search(r"ic acid$", f, re.I):
            forms.add(re.sub(r"ic acid$", "ate", f, flags=re.I))
        elif re.search(r"ate$", f, re.I):
            forms.add(re.sub(r"ate$", "ic acid", f, flags=re.I))
    return [k for k in dict.fromkeys(re.sub(r"[^a-z0-9]", "", f.lower()) for f in forms) if k]


# ---------------------------------------------------------------------------
# Reaction library
# ---------------------------------------------------------------------------
class ReactionLibrary:
    def __init__(self, pairs: pd.DataFrame, mets: pd.DataFrame, info: dict, hub_cap: int = DEFAULT_HUB_CAP):
        self.pairs = pairs.reset_index(drop=True)
        self.mets = mets.set_index("met_id")
        self.info = info
        self.hub_cap = hub_cap
        from metabolomics_modules import pathway_reference_db as ref
        self._norm = ref.normalize_name
        self.by_hmdb, self.by_kegg, self.by_name = (collections.defaultdict(list) for _ in range(3))
        for mid, r in self.mets.iterrows():
            if r.hmdb_id:
                self.by_hmdb[r.hmdb_id].append(mid)
            if r.kegg_id:
                self.by_kegg[r.kegg_id].append(mid)
            for v in self.name_variants(r["name"]):
                self.by_name[v].append(mid)
        # HMDB official names (bundled with MetaboAI Pro) -> HMDB accession, and the app's
        # HMDB <-> KEGG crosswalk as extra KEGG aliases (e.g. C00095 D-fructose vs Human-GEM's
        # beta-D-fructose C02336 share HMDB0000660)
        self.hmdb_by_name = {}
        self.kegg_of_hmdb = {}
        try:
            xw = pd.read_parquet(XWALK_PATH)
            for h, k in zip(xw["HMDB_ID"], xw["KEGG_ID"]):
                if isinstance(k, str) and k.startswith("C"):
                    self.kegg_of_hmdb.setdefault(h, k)
            for h, k, n in zip(xw["HMDB_ID"], xw["KEGG_ID"], xw["HMDB_Name"]):
                # alias only when HMDB's own name confirms the Human-GEM metabolite (the
                # crosswalk is many-to-one in places, e.g. C00117 is listed under 3 HMDB IDs)
                if h in self.by_hmdb and isinstance(k, str) and k.startswith("C") and k not in self.by_kegg:
                    nk = set(name_keys(n))
                    ok = [mid for mid in self.by_hmdb[h] if nk & set(name_keys(self.mets.at[mid, "name"]))]
                    if ok:
                        self.by_kegg[k].extend(ok)
                for v in self.name_variants(n):
                    self.hmdb_by_name.setdefault(v, h)
            gm = pd.read_parquet(HMDB_NAMES_PATH, columns=["HMDB_ID", "Metabolite"]).drop_duplicates("HMDB_ID")
            for h, n in zip(gm["HMDB_ID"], gm["Metabolite"]):
                for v in self.name_variants(n):
                    self.hmdb_by_name.setdefault(v, h)
        except Exception:
            pass
        for k, h in NAME_SYNONYMS.items():
            self.hmdb_by_name[k] = h
        self.no_route = set(self.mets.index[self.mets["is_peptide"]]) if "is_peptide" in self.mets else set()
        self.edge = {}
        self.degree = collections.Counter()
        self.route_degree = collections.Counter()
        self.succ = collections.defaultdict(set)        # routing graph (no peptide/protein edges)
        for k, r in enumerate(self.pairs.itertuples()):
            a, b = r.met_a, r.met_b
            self.edge[(a, b)] = self.edge[(b, a)] = k
            self.degree[a] += 1
            self.degree[b] += 1
            if any(s in (r.subsystems or "") for s in ROUTE_EXCLUDED_SUBSYSTEMS) and not self._has_other_subsystem(r.subsystems):
                continue
            self.route_degree[a] += 1
            self.route_degree[b] += 1
            if r.direction in ("a>b", "both"):
                self.succ[a].add(b)
            if r.direction in ("b>a", "both"):
                self.succ[b].add(a)

    @staticmethod
    def _has_other_subsystem(subs: str) -> bool:
        parts = [s for s in (subs or "").split(";") if s]
        return any(p not in ROUTE_EXCLUDED_SUBSYSTEMS for p in parts)

    # -- identifiers -------------------------------------------------------
    def _best(self, cands):
        cands = list(dict.fromkeys(cands))
        if not cands:
            return None
        return sorted(cands, key=lambda m: (-self.degree[m], m))[0]

    def name_variants(self, name) -> list:
        return name_keys(name)

    def resolve(self, hmdb: str = "", kegg: str = "", name: str = ""):
        """(met_id, how) for one metabolite: HMDB, KEGG, HMDB official name / synonym, then
        Human-GEM name."""
        if hmdb and self.by_hmdb.get(hmdb):
            return self._best(self.by_hmdb[hmdb]), "HMDB"
        if kegg and self.by_kegg.get(kegg):
            return self._best(self.by_kegg[kegg]), "KEGG"
        keys = self.name_variants(name)
        for key in keys:
            h = self.hmdb_by_name.get(key)
            if h and self.by_hmdb.get(h):
                return self._best(self.by_hmdb[h]), "HMDB name"
        for h in [hmdb] + [self.hmdb_by_name.get(k) for k in keys]:
            k = self.kegg_of_hmdb.get(h or "")
            if k and self.by_kegg.get(k):
                return self._best(self.by_kegg[k]), "HMDB name → KEGG"
        for key in keys:
            if self.by_name.get(key):
                return self._best(self.by_name[key]), "Human-GEM name"
        return None, ""

    def converts(self, m, n) -> bool:
        """True if some reaction in the library converts metabolite m into n."""
        k = self.edge.get((m, n))
        if k is None:
            return False
        r = self.pairs.iloc[k]
        return r.direction == "both" or (r.direction == "a>b") == (r.met_a == m)

    def met_name(self, mid) -> str:
        return str(self.mets["name"].get(mid, mid)) if mid else ""

    # -- links between measured metabolites --------------------------------
    def links(self, met_ids: list, max_steps: int = 2) -> dict:
        """met_ids: one Human-GEM id (or None) per feature. Returns {(i, j): link} for i < j,
        link = {steps, direction ('i>j' | 'j>i' | 'both'), via (id or None), reactions,
        reaction_names, subsystems, ec}."""
        idx_of = collections.defaultdict(list)
        for i, m in enumerate(met_ids):
            if m:
                idx_of[m].append(i)
        measured = set(idx_of)
        out = {}

        def put(i, j, rec):
            key = (min(i, j), max(i, j))
            if key in out and out[key]["steps"] < rec["steps"]:
                return
            if key in out and out[key]["steps"] == rec["steps"]:
                d0, d1 = out[key]["direction"], rec["direction"]
                if d0 != d1:
                    out[key]["direction"] = "both"
                return
            out[key] = rec

        # direct
        for m in measured:
            for n in measured:
                if m >= n or (m, n) not in self.edge:
                    continue
                r = self.pairs.iloc[self.edge[(m, n)]]
                m_to_n = self.converts(m, n)
                n_to_m = self.converts(n, m)
                for i in idx_of[m]:
                    for j in idx_of[n]:
                        if i == j:
                            continue
                        a, b = (i, j) if i < j else (j, i)       # a < b; feature a carries m iff a == i
                        a_to_b = m_to_n if a == i else n_to_m
                        b_to_a = n_to_m if a == i else m_to_n
                        d = "both" if (a_to_b and b_to_a) else ("i>j" if a_to_b else "j>i")
                        put(a, b, {"steps": 1, "direction": d, "via": None, "via_name": "",
                                   "reactions": r.reactions, "reaction_names": r.reaction_names,
                                   "subsystems": r.subsystems, "ec": r.ec})
        if max_steps < 2:
            return out
        # two steps: m -> x -> n, x unmeasured, not a hub
        for m in measured:
            for x in self.succ.get(m, ()):
                if x in measured or self.route_degree[x] > self.hub_cap or x in self.no_route:
                    continue
                for n in self.succ.get(x, ()):
                    if n == m or n not in measured or (m, n) in self.edge:
                        continue
                    r1 = self.pairs.iloc[self.edge[(m, x)]]
                    r2 = self.pairs.iloc[self.edge[(x, n)]]
                    # both steps from the same reaction = two co-products (or co-substrates) of one
                    # conversion, e.g. a peptide or glutathione conjugate hydrolysed into glutamate AND
                    # glycine -- not a route from one to the other
                    if set(r1.reactions.split(";")) & set(r2.reactions.split(";")):
                        continue
                    for i in idx_of[m]:
                        for j in idx_of[n]:
                            if i == j:
                                continue
                            a, b = (i, j) if i < j else (j, i)
                            d = "i>j" if i < j else "j>i"   # route runs from feature i's metabolite
                            put(a, b, {"steps": 2, "direction": d, "via": x, "via_name": self.met_name(x),
                                       "reactions": f"{r1.reactions} → {r2.reactions}",
                                       "reaction_names": f"{r1.reaction_names.split(' | ')[0]} → "
                                                         f"{r2.reaction_names.split(' | ')[0]}",
                                       "subsystems": ";".join(sorted(set(filter(None, (r1.subsystems + ';' + r2.subsystems).split(';'))))),
                                       "ec": ";".join(sorted(set(filter(None, (r1.ec + ';' + r2.ec).split(';')))))})
        return out

    def subsystem_sets(self) -> dict:
        """Human-GEM subsystem -> set of met_ids (from the reaction pairs)."""
        sets = collections.defaultdict(set)
        for r in self.pairs.itertuples():
            for s in (r.subsystems or "").split(";"):
                if s and s not in ROUTE_EXCLUDED_SUBSYSTEMS and s not in ("Isolated", "Miscellaneous"):
                    sets[s].update((r.met_a, r.met_b))
        return dict(sets)


_LIB = {}


def reaction_library(hub_cap: int = DEFAULT_HUB_CAP) -> ReactionLibrary:
    if hub_cap not in _LIB:
        with open(INFO_PATH, encoding="utf-8") as fh:
            info = json.load(fh)
        _LIB[hub_cap] = ReactionLibrary(pd.read_parquet(PAIRS_PATH), pd.read_parquet(METS_PATH), info, hub_cap)
    return _LIB[hub_cap]


# ---------------------------------------------------------------------------
# Feature -> metabolite -> pathway mapping
# ---------------------------------------------------------------------------
@dataclass
class FeatureMap:
    table: pd.DataFrame                    # one row per feature (ID mapping + GEM + pathway)
    pathway: list                          # primary pathway per feature (ring grouping)
    sets: dict                             # set name -> sorted feature indices (>= 1 member)
    library_name: str
    notes: list = field(default_factory=list)


def _primary_sets(n: int, sets_by_name: dict, set_size: dict) -> list:
    """Ring pathway per feature: among the sets containing it (preferring sets with >= 3, then >= 2
    measured members -- 3 is also the pathway-test minimum), the one with the highest measured^2 / set size -- well covered and specific, so TCA
    intermediates sit in 'Citrate cycle' rather than the broad 'Glyoxylate and dicarboxylate
    metabolism' map, and fructose-1,6-bisphosphate in glycolysis."""
    member_of = collections.defaultdict(list)
    for sname, idx in sets_by_name.items():
        for i in idx:
            member_of[i].append(sname)
    out = []
    for i in range(n):
        cands = member_of.get(i, [])
        if not cands:
            out.append(None)
            continue
        multi = ([c for c in cands if len(sets_by_name[c]) >= 3] or [c for c in cands if len(sets_by_name[c]) >= 2]
                 or cands)
        multi.sort(key=lambda c: (-(len(sets_by_name[c]) ** 2) / max(1, set_size.get(c, len(sets_by_name[c]))),
                                  -len(sets_by_name[c]), c))
        out.append(multi[0])
    return out


# ---------------------------------------------------------------------------
# Cleaning of metabolite names before database matching (Network Analysis only; the MSEA page's own
# matcher is unchanged). Display names are never altered; only the query used for matching.
# ---------------------------------------------------------------------------
_DASHES = dict.fromkeys(map(ord, "\u2010\u2011\u2012\u2013\u2014\u2015\u2212\ufe63\uff0d"), "-")
_PRIMES = dict.fromkeys(map(ord, "\u2032\u2033\u2035\u2018\u2019\u201b`\u00b4"), "'")
_SPACES = dict.fromkeys(map(ord, "\u00a0\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\t"), " ")
_ZERO_WIDTH = dict.fromkeys(map(ord, "\u200b\u200c\u200d\ufeff\u00ad"), None)
# instrument / export annotations that are not part of the compound name
_ANNOTATION_PATTERNS = [
    r"[\s_]*\[\s*\d*M\s*[+-][^\]]*\]\s*\d*[+-]*",        # [M+H]+, [M-H]-, [2M+Na]+
    r"[\s_]*\(\s*\d*M\s*[+-][^)]*\)\s*\d*[+-]*",          # (M+H)+
    r"[\s_]*\[\s*(?:HILIC|RP|C18|pos|neg|ESI)[^\]]*\]",     # [HILIC-pos], [RP neg]
    r"[\s_]*\(\s*(?:pos|neg|positive|negative|ESI[+-]?)\s*\)",
    r"[\s_]+(?:pos|neg|ESI[+-]?)$",
    r"[\s_@]+RT[\s=_:]*\d+(?:\.\d+)?\s*(?:min)?$",         # _RT5.21, @RT 5.2 min
    r"\s*@\s*\d+(?:\.\d+)?\s*(?:min)?$",                    # @5.21
    r"[\s*\u2020\u2021#?]+$",                                  # trailing footnote marks: * † ‡ # ?
]


def _greek_to_latin(text: str) -> str:
    """β-Alanine -> beta-Alanine (the Greek letter is part of the name; deleting it would turn
    β-alanine into alanine)."""
    import unicodedata
    out = []
    for ch in text:
        if "\u0370" <= ch <= "\u03ff":
            nm = unicodedata.name(ch, "")
            if nm.startswith("GREEK SMALL LETTER ") or nm.startswith("GREEK CAPITAL LETTER "):
                word = nm.split("LETTER ", 1)[1].split()[0].lower()
                word = {"lamda": "lambda", "final": "sigma"}.get(word, word)
                out.append(word)
                continue
        out.append(ch)
    return "".join(out)


def clean_query_name(name) -> tuple:
    """(cleaned name, list of alternatives, list of what was changed) for database matching.
    Alternatives are the parts of an isomer group such as 'Leucine/Isoleucine' or
    'Citrate; Isocitrate' (tried in order); otherwise just the cleaned name."""
    import re
    import unicodedata
    raw = str(name or "")
    changes = []
    t = unicodedata.normalize("NFKC", raw)            # superscripts, full-width, µ -> μ, ligatures
    t = t.translate(_ZERO_WIDTH).translate(_SPACES).translate(_DASHES).translate(_PRIMES)
    if t != raw:
        changes.append("normalized Unicode characters (dashes, primes, spaces, superscripts)")
    g = _greek_to_latin(t)
    if g != t:
        changes.append("Greek letters spelled out")
        t = g
    for pat in _ANNOTATION_PATTERNS:
        u = re.sub(pat, "", t, flags=re.I)
        if u != t and u.strip():
            changes.append("removed annotation '" + t[len(u):].strip() + "'" if t.startswith(u) else "removed annotation")
            t = u
    t = re.sub(r"\s+", " ", t).strip(" _-;,")
    # isomer groups: split on / ; | or ' or ' outside parentheses, when every part is a name
    alts = [t]
    depth, parts, cur = 0, [], ""
    i = 0
    while i < len(t):
        ch = t[i]
        depth += ch in "([{"
        depth -= ch in ")]}"
        if depth == 0 and (ch in "/;|" or t[i:i + 4].lower() == " or "):
            parts.append(cur)
            cur = ""
            i += 4 if t[i:i + 4].lower() == " or " else 1
            continue
        cur += ch
        i += 1
    parts.append(cur)
    parts = [x.strip(" _-,") for x in parts]
    if len(parts) > 1 and all(len(re.sub(r"[^A-Za-z]", "", x)) >= 3 for x in parts):
        alts = parts
        changes.append("isomer group: " + " | ".join(parts))
    # last-resort forms for two unambiguous locant conventions that references often omit:
    # 2'-deoxy nucleosides (2'-Deoxyguanosine = Deoxyguanosine) and N,N(,N)-methyl names
    # (N,N-Dimethylglycine = Dimethylglycine). Tried only after the full name.
    extra = []
    for a in alts:
        for pat in (r"^2'-(?=deoxy)", r"^N,N(?:,N)?-(?=(?:di|tri)methyl)"):
            b = re.sub(pat, "", a, flags=re.I)
            if b != a and b not in alts + extra:
                extra.append(b)
    if extra:
        alts = alts + extra
        changes.append("fallback without locant: " + " | ".join(extra))
    return t, alts, changes


def _display(name) -> str:
    text = str(name)
    return text.split("::", 1)[1] if "::" in text else text


def map_features(features: list, row_annotations: pd.DataFrame | None = None, library: str = "KEGG",
                 annotation_column: str | None = None) -> FeatureMap:
    """
    Standardize feature IDs exactly as MSEA does, attach the Human-GEM metabolite, and
    assign pathways from `library`:
      * "KEGG" / "SMPDB" / "LIPID MAPS" : MSEA libraries (pathway_libraries.get_library)
      * LIB_GEM                         : Human-GEM subsystems
      * LIB_ANNOTATION                  : `annotation_column` of row_annotations
      * LIB_CLUSTER                     : no sets (caller clusters)
    Each feature's primary (ring) pathway is the set containing it with the most measured
    members (ties: smaller library set, then name); features in no set get "Not in library".
    """
    from metabolomics_modules import pathway_analysis as pa
    from metabolomics_modules import pathway_libraries as libs

    notes = []
    feats = [str(f) for f in features]

    def _std(names, ann):
        try:
            return pa.standardize_for_option_a(names, ann).reset_index(drop=True)
        except Exception:  # never let ID standardization stop the analysis
            return pd.DataFrame({"Metabolite": names, "Matched_Name": "", "HMDB_ID": "", "KEGG_ID": "", "ChEBI_ID": "",
                                 "PubChem_CID": "", "Match_Method": libs.MATCH_UNMAPPED, "Matched_On": ""})

    # 1) IDs from the annotation columns and matches on the name as given (keyed by the original name)
    idm = _std(feats, row_annotations)
    # 2) names with special characters, instrument tags or isomer groups: match the cleaned name instead
    #    (unless an annotation-column ID already identified the metabolite)
    ann_methods = {libs.MATCH_ANN_HMDB, libs.MATCH_ANN_KEGG, libs.MATCH_ANN_CHEBI, libs.MATCH_ANN_PUBCHEM,
                   libs.MATCH_ID_IN_NAME}
    cleaned = [clean_query_name(_display(f)) for f in feats]
    queries = [c[1] for c in cleaned]
    idm["Name_Cleaning"] = ["; ".join(c[2]) for c in cleaned]
    redo = [k for k, (f, c) in enumerate(zip(feats, cleaned))
            if c[2] and idm.at[k, "Match_Method"] not in ann_methods]
    if redo:
        flat = [(k, a) for k in redo for a in cleaned[k][1]]
        alt = _std([a for _, a in flat], None)
        cols = [c for c in idm.columns if c not in ("Metabolite", "Name_Cleaning")]
        for k in redo:
            rows = [n for n, (kk, _) in enumerate(flat) if kk == k]
            hit = next((n for n in rows if alt.at[n, "Match_Method"] != libs.MATCH_UNMAPPED), None)
            had_greek = "Greek letters spelled out" in cleaned[k][2]
            if hit is not None:
                for c in cols:
                    if c in alt.columns:
                        idm.at[k, c] = alt.at[hit, c]
                if len(rows) > 1:
                    idm.at[k, "Name_Cleaning"] += f" -> matched '{flat[hit][1]}'"
            elif had_greek:
                # the name as given matched only because the Greek letter was dropped; do not keep it
                for c in ("Matched_Name", "HMDB_ID", "KEGG_ID", "ChEBI_ID", "PubChem_CID", "Matched_On"):
                    if c in idm.columns:
                        idm.at[k, c] = ""
                idm.at[k, "Match_Method"] = libs.MATCH_UNMAPPED
    rl = reaction_library()
    gem, how = [], []
    for k, (r, f) in enumerate(zip(idm.itertuples(), feats)):
        mid, h = rl.resolve(r.HMDB_ID, r.KEGG_ID, r.Matched_Name or queries[k][0])
        if mid is None:
            for q in queries[k]:
                mid, h = rl.resolve(name=q)
                if mid is not None:
                    break
        gem.append(mid)
        how.append(h)
    idm["GEM_ID"] = [g or "" for g in gem]
    idm["GEM_Name"] = [rl.met_name(g) for g in gem]
    idm["GEM_Match"] = how
    # Names the MSEA reference does not know (3-Hydroxykynurenine, Quinolinate, "Glutathione (GSH)")
    # but the reaction library resolved through HMDB's own names: carry the IDs over so the
    # metabolite can also be placed in the MSEA libraries.
    for k in range(len(idm)):
        if idm.at[k, "Match_Method"] == libs.MATCH_UNMAPPED and gem[k] and how[k] in ("HMDB name", "HMDB name → KEGG",
                                                                                    "Human-GEM name"):
            rec = rl.mets.loc[gem[k]]
            if rec.hmdb_id or rec.kegg_id:
                idm.at[k, "HMDB_ID"] = rec.hmdb_id or ""
                idm.at[k, "KEGG_ID"] = rec.kegg_id or ""
                idm.at[k, "Matched_Name"] = rec["name"]
                idm.at[k, "Match_Method"] = "Reaction-library name match"
                idm.at[k, "Matched_On"] = feats[k]

    sets_by_name = {}
    set_size = {}
    if library in (LIB_CLUSTER,):
        pass
    elif library == LIB_ANNOTATION:
        if row_annotations is None or annotation_column not in getattr(row_annotations, "columns", []):
            raise ValueError(f"Annotation column '{annotation_column}' not found.")
        col = row_annotations[annotation_column]
        lut = {str(k): v for k, v in col.items() if pd.notna(v) and str(v).strip()}
        for i, f in enumerate(feats):
            v = lut.get(f, lut.get(_display(f)))
            if v is not None:
                sets_by_name.setdefault(str(v), []).append(i)
        set_size = {k: len(v) for k, v in sets_by_name.items()}
    elif library == LIB_GEM:
        sub = rl.subsystem_sets()
        for s, members in sub.items():
            idx = [i for i, g in enumerate(gem) if g and g in members]
            if idx:
                sets_by_name[s] = idx
            set_size[s] = len(members)
    else:
        lib = libs.get_library(library)
        keys = libs.mapping_keys(idm, lib.key_type, lib)
        expanded = [(_ANOMER.get(k, {k}) if k else set()) for k in keys]
        for s, members in lib.sets.items():
            idx = [i for i, ks in enumerate(expanded) if ks & members]
            if idx:
                sets_by_name[s] = idx
            set_size[s] = len(members)

    primary = _primary_sets(len(feats), sets_by_name, set_size)
    fallback_n = 0
    if library not in (LIB_CLUSTER, LIB_ANNOTATION, LIB_GEM) and any(p is None and g for p, g in zip(primary, gem)):
        # metabolites in no set of the chosen library (most acylcarnitines are in no KEGG map):
        # place them by their Human-GEM subsystem on the ring rather than in one catch-all group
        sub = rl.subsystem_sets()
        gsets, gsize = {}, {}
        for sname, members in sub.items():
            idx = [i for i, g in enumerate(gem) if g and g in members]
            if idx:
                gsets[sname] = idx
            gsize[sname] = len(members)
        gprim = _primary_sets(len(feats), gsets, gsize)
        for i in range(len(feats)):
            if primary[i] is None and gprim[i] is not None:
                primary[i] = f"{gprim[i]} (Human-GEM)"
                fallback_n += 1
    if fallback_n:
        notes.append(f"{fallback_n} feature(s) in no '{library}' set were placed by their Human-GEM subsystem on the "
                     "ring (marked '(Human-GEM)'); the pathway-level test uses only the chosen library's sets.")
    idm["Pathway"] = [p or "" for p in primary]
    idm["All_Pathways"] = ["; ".join(sorted(s for s, idx in sets_by_name.items() if i in idx)) for i in range(len(feats))]
    n_mapped = int((idm["Match_Method"] != libs.MATCH_UNMAPPED).sum())
    n_gem = int((idm["GEM_ID"] != "").sum())
    notes.append(f"{n_mapped} of {len(feats)} features matched to a reference metabolite (same ID standardization as "
                 f"MSEA); {n_gem} matched to the Human-GEM reaction network.")
    if library not in (LIB_CLUSTER,):
        n_in = sum(p is not None for p in primary)
        notes.append(f"{n_in} of {len(feats)} features belong to at least one '{library}' set.")
    return FeatureMap(table=idm, pathway=primary, sets={k: sorted(v) for k, v in sets_by_name.items()},
                      library_name=library, notes=notes)
