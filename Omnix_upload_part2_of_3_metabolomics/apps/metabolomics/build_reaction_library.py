"""
build_reaction_library.py - Build the metabolite reaction-pair library used by the
Correlation Rewiring Map (metabolomics_modules/rewiring_reactions.py).

Source
------
Human-GEM (SysBioChalmers), the community genome-scale metabolic model of a generic human
cell, licensed CC BY 4.0. https://github.com/SysBioChalmers/Human-GEM
  Human1: Robinson JL et al. An atlas of human metabolism. Sci. Signal. 13, eaaz1482 (2020).
  Human2: Luo J et al. Reconstruction of human metabolic models with large language models.
          PNAS 123:e2516511123 (2026).

What "reaction pair" means here
-------------------------------
A pair (substrate -> product) of metabolites that one enzymatic reaction interconverts,
i.e. the pairs whose levels a metabolomics experiment can meaningfully relate
(glutamine -> glutamate, tryptophan -> kynurenine). Derivation, per reaction:

 1. Collapse compartments (MAM01975c and MAM01975m are both "glutamine") and sum the
    stoichiometry. Metabolites with net zero (moved between compartments) drop out, so
    pure transport reactions produce no pairs.
 2. Remove currency / cofactor species (H+, H2O, ATP/ADP/AMP, NAD(P)(H), FAD(H2), CoA,
    Pi/PPi, O2, CO2, NH3, H2O2, nucleotide carriers UDP/GDP/CMP..., PAP/PAPS, sulfate,
    inorganic ions, electron-carrier proteins, ...) and generic/macromolecular species
    (pools, tRNAs, proteins or anything with > MAX_CARBONS carbons, "R"-group formulas). The list is CURRENCY_NAMES below;
    it is applied only to decide pairs, never to hide a measured metabolite.
 3. Pair the remaining substrates and products:
      * one substrate or one product: every substrate is paired with every product
        (A + B -> C condensations and C -> A + B cleavages: both are carbon-carrying);
      * otherwise: optimal one-to-one assignment of substrates to products maximizing the
        conserved-atom fraction (maximum common substructure / larger molecule, RDKit FMCS on
        Human-GEM SMILES; carbon-count ratio when a SMILES is missing); leftover species join
        their best partner if they share at least half their atoms (MIN_EXTRA_SIMILARITY); pairs
        scoring < MIN_SIMILARITY are dropped. Acyl-CoA thioesters are compared by their acyl
        group only, so CoA-transferases pair acetoacetate -> acetoacetyl-CoA, not the two CoA
        esters with each other. This is how
        transaminations resolve into aspartate -> oxaloacetate and 2-oxoglutarate ->
        glutamate (not aspartate -> glutamate), and methyl transfers into betaine ->
        dimethylglycine (not betaine -> methionine).
 3b. Contract routes through cofactor-bound acyl intermediates (hydroxyethyl-ThPP,
     S-acyl-dihydrolipoamides) so the PDH / OGDH / BCKDH complexes give pyruvate ->
     acetyl-CoA, 2-oxoglutarate -> succinyl-CoA, ... as single conversions.
 3c. Transketolase (EC 2.2.1.1) and transaldolase (EC 2.2.1.2) are paired by the carbon unit
     they move (C2 / C3 from a ketose donor to an aldose acceptor), not by look-alike matching.
 4. Deduplicate to unordered metabolite pairs, keeping every supporting reaction, the
    direction(s) seen, reversibility, subsystem(s), EC numbers and gene count.

Identifiers: each Human-GEM metabolite carries HMDB / KEGG / ChEBI / PubChem IDs from
Human-GEM's metabolites.tsv; missing HMDB/KEGG links are filled from MetaboAI Pro's
own HMDB<->KEGG crosswalk (metabolomics_modules/data/hmdb_kegg_crosswalk.parquet).

Outputs (metabolomics_modules/data/)
-----------------------
  humangem_reaction_pairs.parquet   one row per unordered metabolite pair
  humangem_metabolites.parquet      one row per metabolite (compartment-free)
  humangem_library_info.json        source version, counts, build settings

Usage
-----
    python build_reaction_library.py [path/to/Human-GEM checkout or folder with
                                      Human-GEM.yml + metabolites.tsv]
Without a path the two files are downloaded from GitHub (main branch). Needs PyYAML,
RDKit (pip install rdkit) and pandas; the app itself only reads the parquet outputs.
"""

from __future__ import annotations

import collections
import json
import os
import re
import sys
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "metabolomics_modules", "data")
sys.path.insert(0, HERE)
RAW_BASE = "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/"
MIN_SIMILARITY = 0.10
MIN_EXTRA_SIMILARITY = 0.5   # a species left over after the one-to-one assignment joins a partner only if this similar
CARBON_WEIGHT = 0.2        # weight of carbon-count agreement added to the MCS fraction
MCS_TIMEOUT = 2           # seconds per maximum-common-substructure search
MAX_CARBONS = 150          # above this a species is a protein/macromolecule, not a metabolite

# Currency / cofactor species, matched on Human-GEM metabolite names (case-insensitive,
# exact). Kept deliberately explicit so the list can be reviewed.
CURRENCY_NAMES = {
    # protons, water, gases, reactive oxygen
    "h+", "h2o", "o2", "co2", "h2o2", "nh3", "hco3-", "superoxide", "no", "co", "h2s",
    # energy / phosphate
    "atp", "adp", "amp", "pi", "ppi", "gtp", "gdp", "gmp", "utp", "udp", "ump", "ctp", "cdp", "cmp",
    "itp", "idp", "dattp", "datp", "dadp", "dgtp", "dgdp", "dctp", "dcdp",
    "dttp", "dtdp", "dutp", "dudp",
    # redox cofactors
    "nad+", "nadh", "nadp+", "nadph", "fad", "fadh2", "fmn", "fmnh2",
    "ubiquinol", "ubiquinone", "ferricytochrome c", "ferrocytochrome c", "fe2+", "fe3+",
    "thioredoxin", "thioredoxin disulfide", "glutaredoxin", "glutaredoxin disulfide",
    "reduced acceptor", "acceptor", "electron-transfer flavoprotein", "reduced electron-transfer flavoprotein",
    "ferredoxin", "reduced ferredoxin", "oxidized ferredoxin",
    # carriers (thiamin-PP and lipoamide carry acyl groups inside the PDH / OGDH / BCKDH complexes)
    "thiamin-pp", "lipoamide", "dihydrolipoamide",
    "coa", "pap", "paps", "sulfate", "sulfite", "acyl-carrier protein", "acp",
    # inorganic ions
    "na+", "k+", "ca2+", "mg2+", "zn2+", "cu2+", "cu+", "cl-", "i-", "mn2+", "co2+", "cd2+",
    "hg2+", "pb2+", "ni2+", "fluoride", "iodide", "chloride", "selenite", "phosphate",
    # spellings actually used in Human-GEM 2.x (checked against the model's names)
    "nh4+", "li+", "o2-", "cytochrome-c", "apocytochrome-c", "ferricytochrome b5", "ferrocytochrome b5",
    "oxidized thioredoxin", "mitothioredoxin", "mitooxidized thioredoxin",
    "oxidized adrenal ferredoxin", "reduced adrenal ferredoxin", "2fe2s iron-sulfur cluster",
}
# Cofactors that are currency only when they act as a cofactor, i.e. when a member of the same
# family is on the other side of the reaction (NAD+ -> NADH, ATP -> ADP). Otherwise they are
# real substrates/products and are kept: NAD+ made by NAD synthetase, NAD+ consumed by PARPs,
# AMP -> IMP, ATP -> cAMP.
CONTEXT_CURRENCY_FAMILIES = [
    {"nad+", "nadh"}, {"nadp+", "nadph"}, {"fad", "fadh2"}, {"fmn", "fmnh2"},
    {"atp", "adp", "amp", "datp", "dadp", "dattp"}, {"gtp", "gdp", "gmp", "dgtp", "dgdp"},
    {"utp", "udp", "ump", "dutp", "dudp"}, {"ctp", "cdp", "cmp", "dctp", "dcdp"}, {"itp", "idp"},
    {"dttp", "dtdp"},
]
_FAMILY_OF = {n: i for i, fam in enumerate(CONTEXT_CURRENCY_FAMILIES) for n in fam}
# Nucleotide-activated donors (UDP-glucose, GDP-mannose, CDP-choline, CMP-sialate, ADP-ribose ...)
# count as partners of their nucleotide family: in UDP-glucuronate + X -> X-glucuronide + UDP the
# released UDP is a carrier, not a product of X.
_ACTIVATED_PREFIX = {"udp": "udp", "gdp": "gdp", "cdp": "cdp", "cmp": "cmp", "dtdp": "dtdp", "adp": "adp",
                     "tdp": "dtdp"}


def _family_of_name(name: str):
    n = name.strip().lower()
    if n in _FAMILY_OF:
        return _FAMILY_OF[n], True             # (family, is the cofactor itself)
    m = re.match(r"^(udp|gdp|cdp|cmp|dtdp|tdp|adp)-", n)
    if m:
        return _FAMILY_OF[_ACTIVATED_PREFIX[m.group(1)]], False
    return None, False

# generic / macromolecular species: excluded from pairing (they are not measured and
# would create spurious two-step links)
GENERIC_NAME_RE = re.compile(
    r"\bpool\b|activated methyl group|trna|\bprotein\b|biomass|lipid droplet|\bLD\b|apoprotein|\bdna\b|\brna\b|"
    r"glycogenin|\bR-total\b|\bfatty acid\b.*\bpool|cofactor|vitamin.*pool|nucleotide pool", re.I)
GENERIC_FORMULA_RE = re.compile(r"R|X|Z|\(|^$")
# acyl groups bound to a cofactor inside a multi-enzyme complex (hydroxyethyl-ThPP,
# S-acetyldihydrolipoamide, ...): pairs are contracted through them, so pyruvate -> acetyl-CoA
# and 2-oxoglutarate -> succinyl-CoA come out as one conversion, as in KEGG.
BOUND_INTERMEDIATE_RE = re.compile(r"thiamine?[- ]?diphosphate|thpp|dihydrolipoamide", re.I)


def _download(name: str, dest: str):
    url = RAW_BASE + name
    print(f"downloading {url}")
    urllib.request.urlretrieve(url, dest)


def load_humangem(src: str | None):
    import yaml
    folder = src or os.path.join(HERE, ".humangem_cache")
    if src and os.path.isdir(os.path.join(src, "model")):
        folder = os.path.join(src, "model")
    os.makedirs(folder, exist_ok=True)
    yml = os.path.join(folder, "Human-GEM.yml")
    tsv = os.path.join(folder, "metabolites.tsv")
    for f, p in (("Human-GEM.yml", yml), ("metabolites.tsv", tsv)):
        if not os.path.exists(p):
            _download(f, p)
    base_loader = getattr(yaml, "CSafeLoader", yaml.SafeLoader)

    class loader(base_loader):  # noqa: N801
        pass

    # YAML 1.1 resolves bare NO / YES / ON / OFF to booleans, which would turn the metabolite
    # "NO" (nitric oxide) into False. Resolve booleans only from true/false.
    loader.yaml_implicit_resolvers = {k: [(tag, rx) for tag, rx in v if tag != "tag:yaml.org,2002:bool"]
                                      for k, v in base_loader.yaml_implicit_resolvers.items()}
    loader.add_implicit_resolver("tag:yaml.org,2002:bool", re.compile(r"^(?:true|True|TRUE|false|False|FALSE)$"),
                                 list("tTfF"))
    with open(yml, encoding="utf-8") as fh:
        doc = dict(yaml.load(fh, Loader=loader))
    meta = dict(doc["metaData"])
    mets = [dict(m) for m in doc["metabolites"]]
    rxns = []
    for r in doc["reactions"]:
        r = dict(r)
        r["metabolites"] = dict(r.get("metabolites") or {})
        r["annotation"] = dict(r.get("annotation") or {})
        rxns.append(r)
    ann = pd.read_csv(tsv, sep="\t", dtype=str).fillna("")
    return meta, mets, rxns, ann


def _norm_hmdb(x: str) -> str:
    m = re.match(r"HMDB0*(\d+)$", str(x).strip())
    return f"HMDB{int(m.group(1)):07d}" if m else ""


def _first(x: str) -> str:
    return str(x).split(";")[0].strip()


def metabolite_table(mets, ann, crosswalk: dict | None, hmdb_names: dict | None = None):
    base = {}
    for m in mets:
        mid = m["id"][:-1]
        if mid not in base:
            base[mid] = {"met_id": mid, "name": m.get("name", ""), "formula": m.get("formula", "") or ""}
    a = ann.drop_duplicates("metsNoComp").set_index("metsNoComp")
    rows = []
    kegg_to_hmdbs = collections.defaultdict(set)
    if crosswalk:
        for h, k in crosswalk.items():
            kegg_to_hmdbs[k].add(h)
    hmdb_names = hmdb_names or {}
    for mid, rec in base.items():
        r = a.loc[mid] if mid in a.index else None
        hmdb = _norm_hmdb(_first(r["metHMDBID"])) if r is not None else ""
        kegg_all = [k for k in (r["metKEGGID"].split(";") if r is not None else []) if k.startswith("C")]
        kegg = kegg_all[0] if kegg_all else ""
        filled = ""
        name = rec["name"]
        # fill a missing ID from the crosswalk only when it is unambiguous AND HMDB's own name
        # for that accession matches the Human-GEM name
        if crosswalk:
            if not kegg and hmdb and hmdb in crosswalk and _names_agree(name, hmdb_names.get(hmdb, ())):
                kegg, filled = crosswalk[hmdb], "kegg_from_crosswalk"
            elif (not hmdb and kegg and len(kegg_to_hmdbs.get(kegg, ())) == 1):
                h = next(iter(kegg_to_hmdbs[kegg]))
                if _names_agree(name, hmdb_names.get(h, ())):
                    hmdb, filled = h, "hmdb_from_crosswalk"
        is_currency = name.strip().lower() in CURRENCY_NAMES
        is_generic = (bool(GENERIC_NAME_RE.search(name)) or bool(GENERIC_FORMULA_RE.search(rec["formula"]))
                      or _carbons(rec["formula"]) > MAX_CARBONS)   # proteins (albumin, apoB100, ...)
        rows.append({
            **rec, "hmdb_id": hmdb, "kegg_id": kegg,
            "chebi_id": _first(r["metChEBIID"]) if r is not None else "",
            "pubchem_cid": _first(r["metPubChemID"]) if r is not None else "",
            "smiles": _first(r["metSmiles"]) if r is not None else "",
            "id_filled": filled, "is_currency": is_currency, "is_generic": is_generic,
        })
    df = pd.DataFrame(rows)
    return _resolve_hmdb_conflicts(df, hmdb_names)


def _names_agree(gem_name, hmdb_name_list) -> bool:
    from metabolomics_modules.rewiring_reactions import name_keys
    g = set(name_keys(gem_name))
    return any(g & set(name_keys(n)) for n in hmdb_name_list)


def _resolve_hmdb_conflicts(df, hmdb_names):
    """Human-GEM occasionally gives one HMDB accession to two different metabolites (e.g.
    HMDB0000573, elaidic acid, on both elaidate and (9E)-octadecenoyl-CoA). The accession stays
    with the metabolite whose name matches HMDB's own name for it and is removed from the
    others; if no name matches, nothing is changed."""
    df = df.copy()
    df["id_conflict"] = ""
    for h, g in df[df.hmdb_id != ""].groupby("hmdb_id"):
        if len(g) < 2:
            continue
        owners = [i for i, r in g.iterrows() if _names_agree(r["name"], hmdb_names.get(h, ()))]
        if len(owners) != 1:
            continue
        for idx in g.index:
            if idx != owners[0]:
                df.at[idx, "hmdb_id"] = ""
                df.at[idx, "id_conflict"] = f"{h} removed: it belongs to {df.at[owners[0], 'name']}"
    return df


def _carbons(formula: str) -> int:
    m = re.search(r"C(\d*)(?![a-z])", formula or "")
    if not m:
        return 0
    return int(m.group(1) or 1)


class Similarity:
    """Pairing score for a substrate/product pair, in [0, 1]:

        ( exact MCS fraction + skeleton MCS fraction ) / 2  +  CARBON_WEIGHT x carbon-count ratio
        ------------------------------------------------------------------------------------
                                     1 + CARBON_WEIGHT

    * exact MCS fraction: maximum common substructure (RDKit FMCS; elements, bond orders and
      chirality must match) / heavy atoms of the larger molecule -- conserved atoms (chirality
      keeps galactose-1-P with UDP-galactose, not with its epimer glucose-1-P);
    * skeleton MCS fraction: the same, but N, O and S count as one kind of atom and any bond
      order matches -- a conserved skeleton whose reaction centre changed. This is what pairs
      serine -> hydroxypyruvate and pyruvate -> alanine in serine:pyruvate transamination,
      where exact MCS ties;
    * carbon-count ratio: tie-breaker in favour of the conserved carbon skeleton.
    Carriers are compared without the carrier: acyl-CoA thioesters by their acyl group
    (CoA-transferases then pair acetoacetate -> acetoacetyl-CoA), nucleoside-diphosphate sugars
    by their sugar phosphate (UDP-glucose -> glucose-1-phosphate). Without a SMILES the carbon-
    count ratio (x 0.8) is used, with CoA's 21 carbons discounted for CoA esters."""

    def __init__(self, met_df):
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        self._Chem = Chem
        self._pant = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
        self._ndp = Chem.MolFromSmarts("P(~O)(~O)(~O)OP(~O)(~O)OC[C;R]1O[C;R](n)[C;R][C;R]1")
        self.mol, self.c = {}, {}
        for r in met_df.itertuples():
            c = _carbons(r.formula)
            if re.search(r"coa\b", str(r.name), re.I) and c > 21:
                c -= 21                                   # acyl carbons of an acyl-CoA
            self.c[r.met_id] = c
            if r.smiles:
                mol = Chem.MolFromSmiles(r.smiles)
                if mol is not None and mol.GetNumHeavyAtoms() > 0:
                    self.mol[r.met_id] = self._strip_carrier(mol)
        self._cache = {}

    def _fragment_with(self, mol, cut_a, cut_b, keep_atom):
        Chem = self._Chem
        em = Chem.RWMol(mol)
        em.RemoveBond(cut_a, cut_b)
        m2 = em.GetMol()
        for f, ids in zip(Chem.GetMolFrags(m2, asMols=True, sanitizeFrags=False), Chem.GetMolFrags(m2)):
            if keep_atom in ids:
                return f
        return mol

    def _strip_carrier(self, mol):
        hit = mol.GetSubstructMatch(self._pant)
        if hit:                                           # acyl-CoA: keep acyl + S
            s_idx, ch2 = hit[0], hit[1]
            if any(n.GetIdx() not in hit for n in mol.GetAtomWithIdx(s_idx).GetNeighbors()):
                return self._fragment_with(mol, s_idx, ch2, s_idx)
            return mol
        hit = mol.GetSubstructMatch(self._ndp)
        if hit:                                           # NDP-sugar: keep sugar-1-phosphate
            p1, bridge, p2 = hit[0], hit[4], hit[5]
            return self._fragment_with(mol, bridge, p2, p1)
        return mol

    def _mcs(self, ma, mb, skeleton: bool) -> float:
        from rdkit.Chem import rdFMCS
        if skeleton:
            ma, mb = self._skeleton(ma), self._skeleton(mb)
            res = rdFMCS.FindMCS([ma, mb], atomCompare=rdFMCS.AtomCompare.CompareIsotopes,
                                 bondCompare=rdFMCS.BondCompare.CompareAny, ringMatchesRingOnly=True,
                                 timeout=MCS_TIMEOUT)
        else:
            res = rdFMCS.FindMCS([ma, mb], atomCompare=rdFMCS.AtomCompare.CompareElements,
                                 bondCompare=rdFMCS.BondCompare.CompareOrder, ringMatchesRingOnly=True,
                                 matchChiralTag=True, timeout=MCS_TIMEOUT)
        return res.numAtoms / max(ma.GetNumHeavyAtoms(), mb.GetNumHeavyAtoms())

    def _skeleton(self, mol):
        m = self._Chem.Mol(mol)
        for a in m.GetAtoms():                            # label: carbon 1, N/O/S 2, others 3 + Z
            z = a.GetAtomicNum()
            a.SetIsotope(1 if z == 6 else 2 if z in (7, 8, 16) else 3 + z)
        return m

    def __call__(self, a, b) -> tuple[float, str]:
        key = (a, b) if a < b else (b, a)
        if key in self._cache:
            return self._cache[key]
        ca, cb = self.c.get(a, 0), self.c.get(b, 0)
        cr = min(ca, cb) / max(ca, cb) if ca and cb else 0.0
        if a in self.mol and b in self.mol:
            ma, mb = self.mol[a], self.mol[b]
            frac = 0.5 * (self._mcs(ma, mb, False) + self._mcs(ma, mb, True))
            out = ((frac + CARBON_WEIGHT * cr) / (1 + CARBON_WEIGHT), "mcs")
        elif cr:
            out = (0.8 * cr, "carbon_ratio")
        else:
            out = (0.0, "none")
        self._cache[key] = out
        return out


# EC -> carbons moved per turnover, for group transfers whose atom mapping look-alike matching
# gets wrong (the products resemble the "other" substrate). Transketolase moves a C2 ketol
# unit (X5P + E4P -> G3P + F6P), transaldolase a C3 dihydroxyacetone unit (S7P + G3P -> E4P + F6P).
CARBON_TRANSFER_EC = {"2.2.1.1": 2, "2.2.1.2": 3}


def _carbon_shift_rule(r, S, P, sim):
    ecs = set(r.get("annotation", {}).get("ec-code", []) or [])
    k = next((v for e, v in CARBON_TRANSFER_EC.items() if e in ecs), None)
    if k is None or len(S) != 2 or len(P) != 2:
        return None
    c = sim.c
    out = []
    for s_ in S:
        match = [p for p in P if abs(c.get(s_, 0) - c.get(p, 0)) == k]
        if len(match) != 1:
            return None
        out.append((s_, match[0]))
    return out if len({p for _, p in out}) == 2 else None


PEPTIDE_BOND_SMARTS = "[NX3][CX4][CX3](=O)[NX3][CX4][CX3](=O)"   # two consecutive residues
PEPTIDE_NAME_RE = re.compile(r"peptide|angiotensin|bradykinin|kinin\b|tensin\b|enkephalin|endorphin|"
                             r"(yl-){2,}|glycyl|alanyl|leucyl|valyl|seryl|prolyl|histidyl|lysyl|arginyl|"
                             r"tyrosyl|tryptophanyl|phenylalanyl|glutamyl-(?!.*(cysteine|glycine)$)", re.I)


def flag_peptides(met_df):
    """Peptides (dipeptides and longer, e.g. kinetensin, angiotensins, Leu-Trp) are real metabolites
    but must not serve as two-step route intermediates: hydrolysis of one peptide releases several
    amino acids, which would otherwise look 'two reactions apart'. Detected by a two-residue
    backbone in the SMILES or by name."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    patt = Chem.MolFromSmarts(PEPTIDE_BOND_SMARTS)
    flags = []
    for r in met_df.itertuples():
        hit = bool(PEPTIDE_NAME_RE.search(str(r.name)))
        if not hit and r.smiles:
            mol = Chem.MolFromSmiles(r.smiles)
            hit = mol is not None and mol.HasSubstructMatch(patt)
        flags.append(hit)
    out = met_df.copy()
    out["is_peptide"] = flags
    return out


def derive_pairs(rxns, met_df):
    fam_info = {mid: _family_of_name(n) for mid, n in zip(met_df.met_id, met_df["name"])}
    family = {mid: f for mid, (f, _) in fam_info.items() if f is not None}          # cofactors + activated donors
    cofactor = {mid for mid, (f, own) in fam_info.items() if f is not None and own}  # candidates for removal
    excluded_always = set(met_df.loc[(met_df.is_currency & ~met_df.met_id.isin(cofactor)) | met_df.is_generic, "met_id"])
    sim = Similarity(met_df)
    stats = collections.Counter()
    records = []
    for r in rxns:
        net = collections.defaultdict(float)
        for m, c in r["metabolites"].items():
            net[m[:-1]] += float(c)
        net = {m: c for m, c in net.items() if abs(c) > 1e-9}
        if not net:
            stats["transport_only"] += 1
            continue
        fam_s = {family[m] for m, c in net.items() if c < 0 and m in family}
        fam_p = {family[m] for m, c in net.items() if c > 0 and m in family}
        excluded = set(excluded_always)
        for m, c in net.items():                       # cofactor role: same family on both sides
            if m in cofactor and family[m] in (fam_p if c < 0 else fam_s):
                excluded.add(m)
        S = [m for m, c in net.items() if c < 0 and m not in excluded]
        P = [m for m, c in net.items() if c > 0 and m not in excluded]
        if not S or not P:
            stats["no_pair_after_filter"] += 1
            continue
        stats["paired_reactions"] += 1
        pairs = {}
        if len(S) == 1 or len(P) == 1:
            for s in S:
                for p in P:
                    if s != p:
                        pairs[(s, p)] = sim(s, p) + ("single_side",)
        elif _carbon_shift_rule(r, S, P, sim) is not None:
            # transketolase / transaldolase: a C2 / C3 unit moves from a ketose to an aldose
            for (a, b) in _carbon_shift_rule(r, S, P, sim):
                pairs[(a, b)] = sim(a, b) + ("carbon_transfer_rule",)
        else:
            # optimal one-to-one assignment (max total score), then any substrate/product left
            # over (unequal counts) joins its best partner
            from scipy.optimize import linear_sum_assignment
            M = {(s, p): sim(s, p) if s != p else (-1.0, "same") for s in S for p in P}
            W = np.array([[M[(s, p)][0] for p in P] for s in S])
            rows, cols = linear_sum_assignment(-W)
            done_s, done_p = set(), set()
            for i, j in zip(rows, cols):
                if W[i, j] >= MIN_SIMILARITY:
                    pairs[(S[i], P[j])] = M[(S[i], P[j])] + ("assignment",)
                    done_s.add(i)
                    done_p.add(j)
            for i in range(len(S)):
                if i not in done_s:
                    j = int(np.argmax(W[i]))
                    if W[i, j] >= MIN_EXTRA_SIMILARITY:
                        pairs[(S[i], P[j])] = M[(S[i], P[j])] + ("assignment_extra",)
            for j in range(len(P)):
                if j not in done_p:
                    i = int(np.argmax(W[:, j]))
                    if W[i, j] >= MIN_EXTRA_SIMILARITY:
                        pairs[(S[i], P[j])] = M[(S[i], P[j])] + ("assignment_extra",)
        rev = float(r.get("lower_bound", 0)) < 0
        ec = ";".join(r["annotation"].get("ec-code", []) or [])
        genes = r.get("gene_reaction_rule", "") or ""
        n_genes = len(set(re.findall(r"ENSG\d+", genes)))
        for (s, p), (score, how, rule) in pairs.items():
            records.append({
                "substrate": s, "product": p, "rxn_id": r["id"], "rxn_name": r.get("name", ""),
                "subsystem": ";".join(r.get("subsystem") or []), "ec": ec, "n_genes": n_genes,
                "reversible": rev, "similarity": score, "similarity_method": how, "rule": rule,
            })
    return pd.DataFrame(records), stats


def contract_bound_intermediates(rec: pd.DataFrame, met_df: pd.DataFrame, max_depth: int = 4):
    """Add u -> v for every directed route u -> b1 -> ... -> bk -> v whose inner nodes are all
    cofactor-bound intermediates (BOUND_INTERMEDIATE_RE) and u, v are not. Pairs that touch a
    bound intermediate are then dropped (such species are never measured)."""
    bound = set(met_df.loc[met_df["name"].str.contains(BOUND_INTERMEDIATE_RE) & ~met_df.is_currency, "met_id"])
    if not bound:
        return rec, 0
    succ = collections.defaultdict(list)
    for r in rec.itertuples():
        succ[r.substrate].append(r)
        if r.reversible:
            succ[r.product].append(r._replace(substrate=r.product, product=r.substrate))
    new = []
    starts = {r.substrate for r in rec.itertuples() if r.substrate not in bound and r.product in bound}
    for u in starts:
        stack = [(e.product, [e]) for e in succ[u] if e.product in bound]
        while stack:
            node, path = stack.pop()
            for e in succ.get(node, []):
                if e.product in (p.substrate for p in path) or e.product == u:
                    continue
                if e.product in bound:
                    if len(path) < max_depth:
                        stack.append((e.product, path + [e]))
                    continue
                chain = path + [e]
                new.append({
                    "substrate": u, "product": e.product, "rxn_id": "+".join(c.rxn_id for c in chain),
                    "rxn_name": " → ".join(c.rxn_name for c in chain) + " (via enzyme-bound intermediates)",
                    "subsystem": chain[-1].subsystem, "ec": ";".join(sorted({x for c in chain for x in c.ec.split(";") if x})),
                    "n_genes": max(c.n_genes for c in chain), "reversible": all(c.reversible for c in chain),
                    "similarity": min(c.similarity for c in chain), "similarity_method": "chain",
                    "rule": "bound_intermediate",
                })
    keep = rec[~rec.substrate.isin(bound) & ~rec["product"].isin(bound)]
    extra = pd.DataFrame(new).drop_duplicates(["substrate", "product"]) if new else pd.DataFrame(columns=rec.columns)
    return pd.concat([keep, extra], ignore_index=True), len(extra)


def collapse_pairs(rec: pd.DataFrame) -> pd.DataFrame:
    rec = rec.copy()
    rec["a"] = np.where(rec.substrate < rec["product"], rec.substrate, rec["product"])
    rec["b"] = np.where(rec.substrate < rec["product"], rec["product"], rec.substrate)
    rec["fwd"] = rec.substrate == rec.a            # a -> b as written
    out = []
    for (a, b), g in rec.groupby(["a", "b"], sort=True):
        dirs = set()
        for f, rv in zip(g.fwd, g.reversible):
            dirs.add("a>b" if f else "b>a")
            if rv:
                dirs.add("b>a" if f else "a>b")
        direction = "both" if len(dirs) == 2 else dirs.pop()
        out.append({
            "met_a": a, "met_b": b, "direction": direction,
            "n_reactions": int(g.rxn_id.nunique()),
            "reactions": ";".join(sorted(g.rxn_id.unique())),
            "reaction_names": " | ".join(list(dict.fromkeys(g.rxn_name))[:5]),
            "subsystems": ";".join(sorted({s for x in g.subsystem for s in x.split(";") if s})),
            "ec": ";".join(sorted({e for x in g.ec for e in x.split(";") if e})),
            "max_genes": int(g.n_genes.max()),
            "similarity": float(g.similarity.max()),
            "similarity_method": g.similarity_method.iloc[int(np.argmax(g.similarity.values))],
            "rules": ";".join(sorted(set(g.rule))),
        })
    return pd.DataFrame(out)


def main(src=None):
    meta, mets, rxns, ann = load_humangem(src)
    cw_path = os.path.join(OUT_DIR, "hmdb_kegg_crosswalk.parquet")
    crosswalk = None
    if os.path.exists(cw_path):
        cw = pd.read_parquet(cw_path)
        hc = next(c for c in cw.columns if "hmdb" in c.lower())
        kc = next(c for c in cw.columns if "kegg" in c.lower())
        crosswalk = {_norm_hmdb(h): k for h, k in zip(cw[hc], cw[kc]) if _norm_hmdb(h) and str(k).startswith("C")}
    hmdb_names = collections.defaultdict(list)
    if os.path.exists(cw_path):
        for h, n in zip(cw[hc], cw["HMDB_Name"] if "HMDB_Name" in cw else []):
            hmdb_names[_norm_hmdb(h)].append(n)
    gm_path = os.path.join(OUT_DIR, "hmdb_metabolite_gene_mapping.parquet")
    if os.path.exists(gm_path):
        gm = pd.read_parquet(gm_path, columns=["HMDB_ID", "Metabolite"]).drop_duplicates()
        for h, n in zip(gm["HMDB_ID"], gm["Metabolite"]):
            hmdb_names[_norm_hmdb(h)].append(n)
    met_df = metabolite_table(mets, ann, crosswalk, hmdb_names)
    met_df = flag_peptides(met_df)
    n_conflicts = int((met_df.id_conflict != "").sum())
    rec, stats = derive_pairs(rxns, met_df)
    rec, n_contracted = contract_bound_intermediates(rec, met_df)
    stats["pairs_via_bound_intermediates"] = n_contracted
    pairs = collapse_pairs(rec)
    used = set(pairs.met_a) | set(pairs.met_b)
    met_out = met_df[met_df.met_id.isin(used) | (met_df.hmdb_id != "") | (met_df.kegg_id != "")].drop(columns=["smiles"])
    os.makedirs(OUT_DIR, exist_ok=True)
    pairs.to_parquet(os.path.join(OUT_DIR, "humangem_reaction_pairs.parquet"), index=False)
    met_out.to_parquet(os.path.join(OUT_DIR, "humangem_metabolites.parquet"), index=False)
    info = {
        "source": "Human-GEM (SysBioChalmers), CC BY 4.0",
        "source_url": meta.get("sourceUrl", "https://github.com/SysBioChalmers/Human-GEM"),
        "version": meta.get("version"), "model_date": meta.get("date"),
        "citation": "Robinson JL et al. Sci. Signal. 13, eaaz1482 (2020); Luo J et al. PNAS 123:e2516511123 (2026)",
        "n_reactions_in_model": len(rxns), "reaction_stats": dict(stats),
        "n_pairs": int(len(pairs)), "n_metabolites_in_pairs": int(len(used)),
        "n_metabolites_with_hmdb": int((met_out.hmdb_id != "").sum()),
        "n_metabolites_with_kegg": int((met_out.kegg_id != "").sum()),
        "min_similarity": MIN_SIMILARITY, "pairing_score": "mean of exact and skeleton maximum-common-substructure fractions + carbon-count tie-breaker (RDKit FMCS)",
        "hmdb_conflicts_resolved": n_conflicts, "n_peptides_flagged": int(met_df.is_peptide.sum()), "currency_names": sorted(CURRENCY_NAMES),
    }
    with open(os.path.join(OUT_DIR, "humangem_library_info.json"), "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=2)
    print(json.dumps({k: v for k, v in info.items() if k != "currency_names"}, indent=2))
    return pairs, met_out, info


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
