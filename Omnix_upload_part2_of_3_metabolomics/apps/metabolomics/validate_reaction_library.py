"""
validate_reaction_library.py - Checks the Human-GEM reaction-pair library used by the
Rewiring Map (metabolomics_modules/data/humangem_*.parquet, built by build_reaction_library.py).

    python validate_reaction_library.py            # prints a markdown report, exit 1 on failure

Checks
------
1. Gold standard: textbook substrate -> product conversions (glycolysis, TCA, glutamine/
   glutamate, urea cycle, one-carbon, tryptophan/kynurenine, purines, ...). Metabolite
   names are resolved exactly the way user data is (MSEA ID standardization, then Human-GEM),
   so this also tests the name -> ID -> network path. Each must be linked within the
   expected number of steps, in the right direction.
2. Negative controls: pairs with no direct enzymatic link must not be linked directly,
   including transamination partners (aspartate -/-> glutamate) that naive pairing gets wrong.
3. Currency / generic species never appear in a pair.
4. ID agreement: HMDB <-> KEGG pairs in the library vs MetaboAI Pro's own reference
   database and HMDB-KEGG crosswalk.
5. Coverage of the bundled demo datasets.
"""

from __future__ import annotations

import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from metabolomics_modules import rewiring_reactions as rr  # noqa: E402

# (substrate, product, max steps). Names as a user would write them.
GOLD = [
    # glycolysis
    ("Glucose", "Glucose-6-phosphate", 1), ("Glucose-6-phosphate", "Fructose-6-phosphate", 1),
    ("Fructose-6-phosphate", "Fructose-1,6-bisphosphate", 1), ("Fructose-1,6-bisphosphate", "Dihydroxyacetone phosphate", 1),
    ("Fructose-1,6-bisphosphate", "Glyceraldehyde-3-phosphate", 1), ("Dihydroxyacetone phosphate", "Glyceraldehyde-3-phosphate", 1),
    ("3-Phosphoglycerate", "2-Phosphoglycerate", 1), ("2-Phosphoglycerate", "Phosphoenolpyruvate", 1),
    ("Phosphoenolpyruvate", "Pyruvate", 1), ("Pyruvate", "Lactate", 1), ("Pyruvate", "Alanine", 1),
    ("Pyruvate", "Acetyl-CoA", 1), ("Pyruvate", "Oxaloacetate", 1),
    # TCA
    ("Oxaloacetate", "Citrate", 1), ("Acetyl-CoA", "Citrate", 1), ("Citrate", "cis-Aconitate", 1),
    ("cis-Aconitate", "Isocitrate", 1), ("Citrate", "Isocitrate", 2), ("Isocitrate", "alpha-Ketoglutarate", 1),
    ("alpha-Ketoglutarate", "Succinyl-CoA", 1), ("Succinyl-CoA", "Succinate", 1), ("Succinate", "Fumarate", 1),
    ("Fumarate", "Malate", 1), ("Malate", "Oxaloacetate", 1),
    # glutamine / amino acids
    ("Glutamine", "Glutamate", 1), ("Glutamate", "alpha-Ketoglutarate", 1), ("Glutamine", "alpha-Ketoglutarate", 2),
    ("Aspartate", "Oxaloacetate", 1), ("Aspartate", "Asparagine", 1), ("Glutamate", "GABA", 1),
    ("Glutamate", "N-Acetylglutamate", 1), ("Proline", "Glutamate", 2),   # glutamate -> proline is 4 steps (P5CS x2, spontaneous, PYCR) ("Serine", "Glycine", 1),
    ("Phenylalanine", "Tyrosine", 1), ("Tyrosine", "L-DOPA", 1), ("L-DOPA", "Dopamine", 1),
    ("Histidine", "Urocanic acid", 1), ("Histidine", "Histamine", 1),
    ("Leucine", "4-Methyl-2-oxopentanoate", 1), ("Valine", "3-Methyl-2-oxobutanoate", 1),
    # urea cycle / polyamines / creatine
    ("Arginine", "Ornithine", 1), ("Ornithine", "Citrulline", 1), ("Citrulline", "Argininosuccinate", 1),
    ("Argininosuccinate", "Arginine", 1), ("Arginine", "Citrulline", 1), ("Ornithine", "Putrescine", 1),
    ("Putrescine", "Spermidine", 1), ("Spermidine", "Spermine", 1), ("Guanidinoacetate", "Creatine", 1),
    ("Glycine", "Guanidinoacetate", 1),   # AGAT: the carbon skeleton comes from glycine (KEGG RCLASS main pair)
    ("Aspartate", "Asparagine", 1), ("Glutamine", "Glutamate", 1),  # asparagine synthetase main pairs
    # one-carbon / sulfur
    ("Methionine", "S-Adenosylmethionine", 1), ("S-Adenosylmethionine", "S-Adenosylhomocysteine", 1),
    ("S-Adenosylhomocysteine", "Homocysteine", 1), ("Homocysteine", "Methionine", 1),
    ("Homocysteine", "Cystathionine", 1), ("Cystathionine", "Cysteine", 1), ("Glutathione", "Glutathione disulfide", 1),
    ("Choline", "Phosphocholine", 1), ("Choline", "Betaine", 2), ("Betaine", "Dimethylglycine", 1),
    ("Hypotaurine", "Taurine", 1),
    # tryptophan / kynurenine / NAD
    ("Tryptophan", "N-Formylkynurenine", 1), ("N-Formylkynurenine", "Kynurenine", 1), ("Tryptophan", "Kynurenine", 2),
    ("Kynurenine", "Kynurenic acid", 1), ("Kynurenine", "3-Hydroxykynurenine", 1), ("Kynurenine", "Anthranilic acid", 1),
    ("3-Hydroxykynurenine", "3-Hydroxyanthranilic acid", 1), ("3-Hydroxykynurenine", "Xanthurenic acid", 1),
    ("Tryptophan", "5-Hydroxytryptophan", 1), ("5-Hydroxytryptophan", "Serotonin", 1), ("Tryptophan", "Serotonin", 2),
    ("Serotonin", "5-Hydroxyindoleacetic acid", 2), ("Tryptophan", "Tryptamine", 1),
    # purines / pyrimidines
    ("Adenosine", "Inosine", 1), ("Inosine", "Hypoxanthine", 1), ("Hypoxanthine", "Xanthine", 1), ("Xanthine", "Uric acid", 1),
    ("Guanosine", "Guanine", 1), ("Guanine", "Xanthine", 1), ("Uridine", "Uracil", 1), ("Cytidine", "Uridine", 1),
    # lipids / carnitine / others
    ("Carnitine", "Acetylcarnitine", 1), ("Glycerol", "Glycerol-3-phosphate", 1),
    ("Glucose-6-phosphate", "6-Phosphogluconolactone", 1), ("Ribulose-5-phosphate", "Ribose-5-phosphate", 1),
    ("Cholesterol", "7alpha-Hydroxycholesterol", 1), ("Acetoacetate", "3-Hydroxybutyrate", 1),
]

# pairs that must NOT be directly linked (no single enzyme interconverts them)
NEGATIVE_DIRECT = [
    ("Aspartate", "Glutamate"),          # transamination partners, not substrate/product
    ("Alanine", "Glutamate"), ("Leucine", "Glutamate"), ("Pyruvate", "Glutamate"),
    ("Glucose", "Taurine"), ("Citrate", "Tryptophan"), ("Lactate", "Kynurenine"), ("Glutamine", "Serotonin"),
    ("Hypoxanthine", "Taurine"), ("Palmitic acid", "Hypoxanthine"), ("Tryptophan", "Leucine"),
    ("Alanine", "Glycine"), ("Serine", "Leucine"), ("Glutamine", "Citrate"),
    ("Glutamine", "Asparagine"),         # asparagine synthetase: amide transfer, not a skeleton conversion
    ("Tyrosine", "Glutamate"), ("Phenylalanine", "Glutamate"),
]
# pairs that must not be linked even within 2 steps (would only connect through hubs or peptides)
NEGATIVE_ANY = [
    ("Tryptophan", "Leucine"), ("Tryptophan", "Phenylalanine"), ("Valine", "Histidine"), ("Glucose", "Taurine"),
    ("Hypoxanthine", "Taurine"), ("Glutamine", "Serotonin"), ("Lactate", "Kynurenine"),
]


def _resolve(names):
    fm = rr.map_features(names, None, library=rr.LIB_CLUSTER)
    return dict(zip(names, fm.table["GEM_ID"])), fm.table


def check_gold(lib):
    names = sorted({n for s, p, _ in GOLD for n in (s, p)})
    gem, table = _resolve(names)
    rows, unresolved = [], [n for n in names if not gem[n]]
    for s, p, k in GOLD:
        ms, mp = gem[s], gem[p]
        if not ms or not mp:
            rows.append((s, p, k, None, "unresolved", ""))
            continue
        L = lib.links([ms, mp], max_steps=2).get((0, 1))
        if L is None:
            rows.append((s, p, k, None, "MISSING", ""))
            continue
        forward = L["direction"] in ("i>j", "both")
        ok = L["steps"] <= k and forward
        why = "ok" if ok else ("wrong direction" if not forward else f"found in {L['steps']} steps")
        rows.append((s, p, k, L["steps"], why, (L["via_name"] or L["reaction_names"].split(" | ")[0])[:60]))
    return pd.DataFrame(rows, columns=["substrate", "product", "max_steps", "steps", "result", "via / reaction"]), unresolved, table


def check_negatives(lib):
    out = []
    for s, p in NEGATIVE_DIRECT:
        gem, _ = _resolve([s, p])
        L = lib.links([gem[s], gem[p]], 1).get((0, 1)) if gem[s] and gem[p] else None
        out.append((s, p, "direct", "ok" if L is None else f"LINKED by {L['reaction_names'][:60]}"))
    for s, p in NEGATIVE_ANY:
        gem, _ = _resolve([s, p])
        L = lib.links([gem[s], gem[p]], 2).get((0, 1)) if gem[s] and gem[p] else None
        out.append((s, p, "≤2 steps", "ok" if L is None else f"LINKED ({L['steps']} steps via {L['via_name'] or '-'})"))
    return pd.DataFrame(out, columns=["a", "b", "scope", "result"])


def check_currency(lib):
    """Always-currency and generic species must never be in a pair. Context cofactors (NAD+,
    ATP, ... - see CONTEXT_CURRENCY_FAMILIES in the build script) may appear only as genuine
    substrates/products; they are reported with their degree and must not be usable as
    two-step intermediates (route degree above the hub cap) unless their degree is small."""
    import build_reaction_library as B
    m = lib.mets
    ctx = {n for fam in B.CONTEXT_CURRENCY_FAMILIES for n in fam}
    lower = m["name"].str.strip().str.lower()
    always = set(m.index[(m.is_currency & ~lower.isin(ctx)) | m.is_generic])
    used = set(lib.pairs.met_a) | set(lib.pairs.met_b)
    bad = sorted(m.loc[list(always & used), "name"])
    ctx_used = [(m.at[i, "name"], lib.degree[i], lib.route_degree[i]) for i in used if lower.get(i) in ctx]
    return bad, sorted(ctx_used, key=lambda x: -x[1])


def check_ids(lib):
    from metabolomics_modules import pathway_reference_db as ref
    m = lib.mets
    agree = disagree = absent = 0
    bad = []
    for rec in ref.METABOLITES:
        h, k = rec["hmdb_id"], rec["kegg_id"]
        if not h or not k:
            continue
        g = m[m.hmdb_id == h]
        if g.empty:
            absent += 1
            continue
        if k in set(g.kegg_id) or (g.kegg_id == "").all():
            agree += 1
        else:
            disagree += 1
            bad.append((rec["name"], h, k, ";".join(sorted(set(g.kegg_id))), ";".join(g["name"])))
    return agree, disagree, absent, bad


def check_coverage():
    rows = []
    for label, path in [("Targeted demo", "sample_peak_area_matrix_targeted.csv"),
                        ("Untargeted demo", "sample_peak_area_matrix_untargeted.csv"),
                        ("Rich Clinical demo", "sample_peak_area_matrix_richdemo.csv"),
                        ("Rewiring demo", "sample_peak_area_matrix_rewiring_demo.csv")]:
        p = os.path.join(HERE, path)
        if not os.path.exists(p):
            continue
        names = [n for n in pd.read_csv(p).iloc[:, 0].astype(str) if "ISTD" not in n]   # not biology
        ann_path = p.replace("sample_peak_area_matrix_", "metabolite_row_annotations_")
        ann = pd.read_csv(ann_path).set_index("Metabolite") if os.path.exists(ann_path) else None
        fm = rr.map_features(names, ann, library="KEGG")
        gem_ids = [g or None for g in fm.table["GEM_ID"]]
        links = rr.reaction_library().links(gem_ids, 2)
        n1 = sum(1 for v in links.values() if v["steps"] == 1)
        rows.append((label, len(names), int((fm.table["Match_Method"] != "Unmapped").sum()),
                     int((fm.table["GEM_ID"] != "").sum()), sum(p is not None for p in fm.pathway), n1, len(links) - n1))
    return pd.DataFrame(rows, columns=["dataset", "features", "ID-mapped", "in Human-GEM", "in a KEGG set",
                                       "direct reaction pairs", "2-step pairs"])


def main() -> int:
    lib = rr.reaction_library()
    info = lib.info
    print(f"# Reaction library validation\n\nSource: {info['source']} v{info['version']} ({info['model_date']}); "
          f"{info['n_pairs']} metabolite pairs from {info['reaction_stats']['paired_reactions']} reactions, "
          f"{info['n_metabolites_in_pairs']} metabolites.\n")
    gold, unresolved, table = check_gold(lib)
    n_ok = int((gold.result == "ok").sum())
    print(f"## 1. Gold-standard conversions: {n_ok}/{len(gold)} correct\n")
    print(gold.to_markdown(index=False))
    if unresolved:
        print(f"\nNames not resolved: {', '.join(unresolved)}")
    neg = check_negatives(lib)
    n_neg = int((neg.result == "ok").sum())
    print(f"\n## 2. Negative controls: {n_neg}/{len(neg)} correctly not linked\n")
    print(neg.to_markdown(index=False))
    cur, ctx = check_currency(lib)
    print(f"\n## 3. Always-currency / generic species in pairs: {len(cur)}" + (f" ({', '.join(cur[:20])})" if cur else ""))
    print("\nContext cofactors kept where they are real substrates/products (name: pairs, route degree; "
          f"two-step intermediates need route degree ≤ {lib.hub_cap}):\n")
    print(", ".join(f"{n}: {d}, {r}" for n, d, r in ctx))
    agree, disagree, absent, bad = check_ids(lib)
    print(f"\n## 4. HMDB↔KEGG agreement with MetaboAI Pro's reference: {agree} agree, {disagree} disagree, "
          f"{absent} reference HMDB IDs not in Human-GEM")
    for b in bad[:25]:
        print(f"- {b[0]} ({b[1]}): reference KEGG {b[2]}, Human-GEM {b[3]} [{b[4]}]")
    print("\n## 5. Coverage of bundled datasets\n")
    print(check_coverage().to_markdown(index=False))
    failed = (n_ok < len(gold)) or (n_neg < len(neg)) or bool(cur)
    print(f"\n**{'FAIL' if failed else 'PASS'}**")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
