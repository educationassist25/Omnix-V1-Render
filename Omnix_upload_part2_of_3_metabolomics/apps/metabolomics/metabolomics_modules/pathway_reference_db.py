"""
pathway_reference_db.py - Embedded, offline metabolite/pathway/gene reference
database for MetaboAI Pro's Pathway Analysis tab.

Contents
--------
1. PATHWAYS
   Curated metabolic pathway definitions. The ten "primary" pathways use the
   exact category names of the Rich Clinical Demo's row-annotation `Pathway`
   column (so the reference is internally consistent with the demo data); the
   remaining "sub-pathways" are finer, KEGG-style pathways that a metabolite
   can additionally belong to (e.g. succinate -> TCA Cycle + Oxidative
   Phosphorylation + Alanine, Aspartate and Glutamate Metabolism). Each has a
   stable internal ID (MPA-###) and, for transparency, the closest KEGG human
   map. These are CURATED approximations of KEGG/SMPDB pathway membership --
   not a verbatim copy of either database.

2. METABOLITES
   One record per metabolite: display name, HMDB accession (HMDB + 7 digits),
   primary pathway, additional pathways, and synonyms used for name matching of
   user-uploaded data. Metabolite -> gene associations are NOT stored here any
   more: the earlier hand-curated gene lists (a few "representative" genes per
   metabolite) were replaced by the real HMDB export (858,077 associations,
   22,849 metabolites) in metabolomics_modules/hmdb_gene_mapping.py. `demo=True`
   marks the 200 metabolites used by the Rich Clinical Demo dataset;
   `demo=False` entries are extra, reference-only metabolites so that
   (a) common metabolites in users' own data can still be mapped and
   (b) the "full reference" background option is genuinely larger than the
   demo's detected set.

   HMDB accessions were checked against hmdb.ca metabocards during curation and
   have since been cross-checked against the real HMDB export: all 214 of the
   234 accessions present in the export carry the same compound (identical name
   or a genuine synonym, e.g. Glycodeoxycholic acid = "Deoxycholic acid glycine
   conjugate"). 3-Hydroxybutyric acid's HMDB0000357 is an HMDB SECONDARY
   accession of HMDB0000011 (hmdb.ca); Option B resolves it. The other 19 absent
   accessions are primary HMDB records with no protein associations (mostly
   acylcarnitines, pregnenolone sulfate, tauro-deoxy/urso-deoxycholic acid).

2b. _XREF / BY_KEGG / BY_CHEBI / BY_PUBCHEM
   KEGG Compound, ChEBI and PubChem CID cross-references per HMDB accession
   (provenance and verification level documented at the table itself). Used only by
   Option A's ID standardization and metabolite-set libraries (pathway_libraries.py).

3. GENE_SETS
   Curated human metabolic gene sets (current HGNC symbols), defined
   INDEPENDENTLY of the HMDB metabolite->gene associations (whole-pathway enzyme
   lists, KEGG-style), offered as one gene-set collection for gene-based ORA.

Everything here is plain Python data -- no network access is ever required.
"""

import re

# ---------------------------------------------------------------------------
# 1. Pathways
# ---------------------------------------------------------------------------
# key: (pathway name, internal id, closest KEGG human map / note, tier)
_PW = {
    # --- primary categories (match the Rich Clinical Demo annotation file) ---
    "AA":   ("Amino Acid Metabolism", "MPA-001", "KEGG class 09105 (amino acid metabolism, several maps)", "primary"),
    "TCA":  ("TCA Cycle", "MPA-002", "hsa00020 Citrate cycle (TCA cycle)", "primary"),
    "GLY":  ("Glycolysis", "MPA-003", "hsa00010 Glycolysis / Gluconeogenesis", "primary"),
    "LIP":  ("Lipid Metabolism", "MPA-004", "hsa00564 Glycerophospholipid + hsa00561 Glycerolipid metabolism", "primary"),
    "FAO":  ("Fatty Acid Oxidation", "MPA-005", "hsa00071 Fatty acid degradation (+ carnitine shuttle)", "primary"),
    "BA":   ("Bile Acid Metabolism", "MPA-006", "hsa00120 + hsa00121 (primary + secondary bile acids)", "primary"),
    "SPH":  ("Sphingolipid Metabolism", "MPA-007", "hsa00600 Sphingolipid metabolism", "primary"),
    "PUR":  ("Purine Metabolism", "MPA-008", "hsa00230 Purine metabolism", "primary"),
    "NUC":  ("Nucleotide Metabolism", "MPA-009", "hsa00240 Pyrimidine metabolism (purines are MPA-008)", "primary"),
    "STE":  ("Steroid Hormone Biosynthesis", "MPA-010", "hsa00140 Steroid hormone biosynthesis", "primary"),
    # --- sub-pathways ---
    "PPP":  ("Pentose Phosphate Pathway", "MPA-011", "hsa00030 Pentose phosphate pathway", "sub"),
    "GAL":  ("Galactose Metabolism", "MPA-012", "hsa00052 Galactose metabolism", "sub"),
    "FRU":  ("Fructose and Mannose Metabolism", "MPA-013", "hsa00051 Fructose and mannose metabolism", "sub"),
    "ADG":  ("Alanine, Aspartate and Glutamate Metabolism", "MPA-014", "hsa00250", "sub"),
    "GST":  ("Glycine, Serine and Threonine Metabolism", "MPA-015", "hsa00260", "sub"),
    "BCAA": ("Valine, Leucine and Isoleucine Degradation", "MPA-016", "hsa00280", "sub"),
    "UREA": ("Urea Cycle", "MPA-017", "hsa00220 Arginine biosynthesis (urea cycle)", "sub"),
    "KET":  ("Ketone Body Metabolism", "MPA-018", "hsa00072 Synthesis and degradation of ketone bodies", "sub"),
    "FAS":  ("Fatty Acid Biosynthesis", "MPA-019", "hsa00061 Fatty acid biosynthesis", "sub"),
    "UFA":  ("Biosynthesis of Unsaturated Fatty Acids", "MPA-020", "hsa01040", "sub"),
    "PBA":  ("Primary Bile Acid Biosynthesis", "MPA-021", "hsa00120", "sub"),
    "SBA":  ("Secondary Bile Acid Metabolism", "MPA-022", "hsa00121 (gut-microbial; no human gene set)", "sub"),
    "OXP":  ("Oxidative Phosphorylation", "MPA-023", "hsa00190", "sub"),
    "GC":   ("Glucocorticoid & Mineralocorticoid Biosynthesis", "MPA-024", "subset of hsa00140", "sub"),
    "AE":   ("Androgen & Estrogen Metabolism", "MPA-025", "subset of hsa00140", "sub"),
    "CER":  ("Ceramide De Novo Synthesis", "MPA-026", "subset of hsa00600", "sub"),
}

PATHWAYS = {v[0]: {"id": v[1], "kegg_reference": v[2], "tier": v[3]} for v in _PW.values()}
PRIMARY_PATHWAYS = [v[0] for v in _PW.values() if v[3] == "primary"]


def _p(*keys):
    return [_PW[k][0] for k in keys]


# ---------------------------------------------------------------------------
# 2. Metabolites
# ---------------------------------------------------------------------------
# Tuple layout:
#   (name, HMDB ID, primary key, [extra pathway keys], lipid_class, demo, [synonyms])
# (The former hand-curated [genes] element was removed: metabolite->gene associations now come
#  from the real HMDB export via metabolomics_modules/hmdb_gene_mapping.py.)
# lipid_class=True marks lipid species (acyl chains / sterol-like lipids typically
# reported by a lipidomics assay) -- used only to assign demo names sensibly to
# the demo's lipidomics-method rows.
_M = [
    # ======================= TCA Cycle (18 demo) =======================
    ("Citric acid", "HMDB0000094", "TCA", [], False, True, []),
    ("cis-Aconitic acid", "HMDB0000072", "TCA", [], False, True, []),
    ("Isocitric acid", "HMDB0000193", "TCA", [], False, True, []),
    ("Oxalosuccinic acid", "HMDB0003974", "TCA", [], False, True, []),
    ("Oxoglutaric acid", "HMDB0000208", "TCA", ["ADG"], False, True,
     ["alpha-Ketoglutaric acid", "alpha-Ketoglutarate", "2-Oxoglutarate", "2-Oxoglutaric acid", "AKG", "2-Ketoglutarate"]),
    ("Succinyl-CoA", "HMDB0001022", "TCA", ["BCAA"], False, True, []),
    ("Succinic acid", "HMDB0000254", "TCA", ["OXP", "ADG"], False, True, []),
    ("Fumaric acid", "HMDB0000134", "TCA", ["OXP", "UREA", "ADG"], False, True, []),
    ("L-Malic acid", "HMDB0000156", "TCA", [], False, True, ["Malic acid", "Malate"]),
    ("Oxaloacetic acid", "HMDB0000223", "TCA", ["ADG"], False, True,
     ["Oxalacetic acid", "Oxaloacetate", "OAA"]),
    ("Acetyl-CoA", "HMDB0001206", "TCA", ["FAO", "KET", "FAS"], False, True, ["Acetyl coenzyme A"]),
    ("Coenzyme A", "HMDB0001423", "TCA", [], False, True, ["CoA"]),
    ("NAD", "HMDB0000902", "TCA", ["OXP"], False, True, ["NAD+", "Nicotinamide adenine dinucleotide"]),
    ("NADH", "HMDB0001487", "TCA", ["OXP"], False, True, []),
    ("FAD", "HMDB0001248", "TCA", ["OXP"], False, True, ["Flavin adenine dinucleotide"]),
    ("FADH2", "HMDB0001197", "TCA", ["OXP"], False, True, ["FADH"]),
    ("Itaconic acid", "HMDB0002092", "TCA", [], False, True, []),
    ("D-2-Hydroxyglutaric acid", "HMDB0000606", "TCA", [], False, True, ["D-2-Hydroxyglutarate", "2-Hydroxyglutarate", "2-HG"]),

    # ======================= Glycolysis (24 demo) =======================
    ("D-Glucose", "HMDB0000122", "GLY", ["GAL"], False, True, ["Glucose"]),
    ("Glucose 6-phosphate", "HMDB0001401", "GLY", ["PPP", "GAL"], False, True, ["G6P"]),
    ("Glucose 1-phosphate", "HMDB0001586", "GLY", ["GAL"], False, True, ["G1P"]),
    ("Fructose 6-phosphate", "HMDB0000124", "GLY", ["PPP", "FRU"], False, True, ["F6P"]),
    ("Fructose 1,6-bisphosphate", "HMDB0001058", "GLY", ["FRU"], False, True, ["F16BP", "Fructose-1,6-diphosphate"]),
    ("Fructose 2,6-bisphosphate", "HMDB0001047", "GLY", ["FRU"], False, True, ["D-Fructose 2,6-bisphosphate"]),
    ("Dihydroxyacetone phosphate", "HMDB0001473", "GLY", ["FRU"], False, True, ["DHAP"]),
    ("Glyceraldehyde 3-phosphate", "HMDB0001112", "GLY", ["PPP", "FRU"], False, True, ["D-Glyceraldehyde 3-phosphate", "GAP"]),
    ("1,3-Bisphosphoglyceric acid", "HMDB0001270", "GLY", [], False, True, ["Glyceric acid 1,3-biphosphate", "1,3-BPG"]),
    ("3-Phosphoglyceric acid", "HMDB0000807", "GLY", ["GST"], False, True, ["3-Phosphoglycerate", "3PG"]),
    ("2-Phosphoglyceric acid", "HMDB0000362", "GLY", [], False, True, ["2-Phosphoglycerate", "2PG"]),
    ("Phosphoenolpyruvic acid", "HMDB0000263", "GLY", [], False, True, ["Phosphoenolpyruvate", "PEP"]),
    ("Pyruvic acid", "HMDB0000243", "GLY", ["ADG"], False, True, []),
    ("L-Lactic acid", "HMDB0000190", "GLY", [], False, True, ["Lactic acid", "Lactate"]),
    ("2,3-Bisphosphoglyceric acid", "HMDB0001294", "GLY", [], False, True, ["2,3-Diphosphoglyceric acid", "2,3-BPG", "2,3-DPG"]),
    ("D-Fructose", "HMDB0000660", "GLY", ["FRU"], False, True, ["Fructose"]),
    ("Glycerol 3-phosphate", "HMDB0000126", "GLY", ["LIP"], False, True, ["sn-Glycerol 3-phosphate", "Glycerophosphate"]),
    ("Glycerol", "HMDB0000131", "GLY", ["LIP", "GAL"], False, True, []),
    ("D-Mannose", "HMDB0000169", "GLY", ["FRU"], False, True, ["Mannose"]),
    ("Mannose 6-phosphate", "HMDB0001078", "GLY", ["FRU"], False, True, []),
    ("D-Galactose", "HMDB0000143", "GLY", ["GAL"], False, True, ["Galactose"]),
    ("Galactose 1-phosphate", "HMDB0000645", "GLY", ["GAL"], False, True, []),
    ("Glyceric acid", "HMDB0000139", "GLY", ["GST"], False, True, ["Glycerate"]),
    ("Methylglyoxal", "HMDB0001167", "GLY", [], False, True, ["Pyruvaldehyde"]),

    # =================== Amino Acid Metabolism (18 demo) ===================
    ("L-Alanine", "HMDB0000161", "AA", ["ADG"], False, True, ["Alanine"]),
    ("L-Arginine", "HMDB0000517", "AA", ["UREA"], False, True, ["Arginine"]),
    ("L-Aspartic acid", "HMDB0000191", "AA", ["ADG", "UREA"], False, True, ["Aspartic acid", "Aspartate"]),
    ("L-Glutamic acid", "HMDB0000148", "AA", ["ADG"], False, True, ["Glutamic acid", "Glutamate"]),
    ("L-Glutamine", "HMDB0000641", "AA", ["ADG"], False, True, ["Glutamine"]),
    ("Glycine", "HMDB0000123", "AA", ["GST"], False, True, []),
    ("L-Histidine", "HMDB0000177", "AA", [], False, True, ["Histidine"]),
    ("L-Isoleucine", "HMDB0000172", "AA", ["BCAA"], False, True, ["Isoleucine"]),
    ("L-Leucine", "HMDB0000687", "AA", ["BCAA"], False, True, ["Leucine"]),
    ("L-Lysine", "HMDB0000182", "AA", [], False, True, ["Lysine"]),
    ("L-Methionine", "HMDB0000696", "AA", [], False, True, ["Methionine"]),
    ("L-Phenylalanine", "HMDB0000159", "AA", [], False, True, ["Phenylalanine"]),
    ("L-Proline", "HMDB0000162", "AA", [], False, True, ["Proline"]),
    ("L-Serine", "HMDB0000187", "AA", ["GST", "SPH", "CER"], False, True, ["Serine"]),
    ("L-Threonine", "HMDB0000167", "AA", ["GST"], False, True, ["Threonine"]),
    ("L-Tryptophan", "HMDB0000929", "AA", [], False, True, ["Tryptophan"]),
    ("L-Tyrosine", "HMDB0000158", "AA", [], False, True, ["Tyrosine"]),
    ("L-Valine", "HMDB0000883", "AA", ["BCAA"], False, True, ["Valine"]),

    # ================ Nucleotide Metabolism (pyrimidines, 11 demo) ================
    ("Uracil", "HMDB0000300", "NUC", [], False, True, []),
    ("Uridine", "HMDB0000296", "NUC", [], False, True, []),
    ("Uridine 5'-monophosphate", "HMDB0000288", "NUC", [], False, True, ["UMP", "Uridine monophosphate"]),
    ("Uridine 5'-diphosphate", "HMDB0000295", "NUC", [], False, True, ["UDP", "Uridine diphosphate"]),
    ("Uridine triphosphate", "HMDB0000285", "NUC", [], False, True, ["UTP"]),
    ("Cytidine", "HMDB0000089", "NUC", [], False, True, []),
    ("Cytidine monophosphate", "HMDB0000095", "NUC", [], False, True, ["CMP"]),
    ("Cytidine triphosphate", "HMDB0000082", "NUC", [], False, True, ["CTP"]),
    ("Thymine", "HMDB0000262", "NUC", [], False, True, []),
    ("Thymidine", "HMDB0000273", "NUC", [], False, True, []),
    ("Orotic acid", "HMDB0000226", "NUC", [], False, True, ["Orotate"]),

    # =================== Purine Metabolism (23 demo) ===================
    ("Adenine", "HMDB0000034", "PUR", [], False, True, []),
    ("Guanine", "HMDB0000132", "PUR", [], False, True, []),
    ("Hypoxanthine", "HMDB0000157", "PUR", [], False, True, []),
    ("Xanthine", "HMDB0000292", "PUR", [], False, True, []),
    ("Uric acid", "HMDB0000289", "PUR", [], False, True, ["Urate"]),
    ("Inosine", "HMDB0000195", "PUR", [], False, True, []),
    ("Adenosine", "HMDB0000050", "PUR", [], False, True, []),
    ("Guanosine", "HMDB0000133", "PUR", [], False, True, []),
    ("Xanthosine", "HMDB0000299", "PUR", [], False, True, []),
    ("Adenosine monophosphate", "HMDB0000045", "PUR", [], False, True, ["AMP"]),
    ("ADP", "HMDB0001341", "PUR", ["OXP"], False, True, ["Adenosine diphosphate"]),
    ("Adenosine triphosphate", "HMDB0000538", "PUR", ["OXP"], False, True, ["ATP"]),
    ("Inosinic acid", "HMDB0000175", "PUR", [], False, True, ["IMP", "Inosine monophosphate"]),
    ("Guanosine monophosphate", "HMDB0001397", "PUR", [], False, True, ["GMP"]),
    ("Guanosine diphosphate", "HMDB0001201", "PUR", [], False, True, ["GDP"]),
    ("Guanosine triphosphate", "HMDB0001273", "PUR", [], False, True, ["GTP"]),
    ("Xanthylic acid", "HMDB0001554", "PUR", [], False, True, ["XMP", "Xanthosine monophosphate"]),
    ("Deoxyadenosine", "HMDB0000101", "PUR", [], False, True, []),
    ("Deoxyguanosine", "HMDB0000085", "PUR", [], False, True, []),
    ("Deoxyinosine", "HMDB0000071", "PUR", [], False, True, []),
    # Allantoin: humans lack urate oxidase (UOX is a pseudogene); it forms
    # non-enzymatically from urate -- deliberately has no gene association.
    ("Allantoin", "HMDB0000462", "PUR", [], False, True, []),
    ("Cyclic AMP", "HMDB0000058", "PUR", [], False, True, ["cAMP"]),
    ("Cyclic GMP", "HMDB0001314", "PUR", [], False, True, ["cGMP"]),

    # ================ Steroid Hormone Biosynthesis (22 demo) ================
    ("Cholesterol", "HMDB0000067", "STE", ["PBA"], True, True, []),
    ("Pregnenolone", "HMDB0000253", "STE", ["GC"], True, True, []),
    ("Pregnenolone sulfate", "HMDB0000774", "STE", [], True, True, []),
    ("17alpha-Hydroxypregnenolone", "HMDB0000363", "STE", ["AE"], True, True, ["17a-Hydroxypregnenolone", "17-Hydroxypregnenolone"]),
    ("Progesterone", "HMDB0001830", "STE", ["GC"], True, True, []),
    ("17alpha-Hydroxyprogesterone", "HMDB0000374", "STE", ["GC"], True, True, ["17-Hydroxyprogesterone", "17-OHP"]),
    ("11-Deoxycorticosterone", "HMDB0000016", "STE", ["GC"], True, True, ["Deoxycorticosterone", "DOC"]),
    ("Corticosterone", "HMDB0001547", "STE", ["GC"], True, True, []),
    ("18-Hydroxycorticosterone", "HMDB0000319", "STE", ["GC"], True, True, []),
    ("Aldosterone", "HMDB0000037", "STE", ["GC"], True, True, []),
    ("11-Deoxycortisol", "HMDB0000015", "STE", ["GC"], True, True, ["Cortexolone"]),
    ("Cortisol", "HMDB0000063", "STE", ["GC"], True, True, ["Hydrocortisone"]),
    ("Cortisone", "HMDB0002802", "STE", ["GC"], True, True, []),
    ("Dehydroepiandrosterone", "HMDB0000077", "STE", ["AE"], True, True, ["DHEA"]),
    ("Dehydroepiandrosterone sulfate", "HMDB0001032", "STE", ["AE"], True, True, ["DHEA-S", "DHEAS"]),
    ("Androstenedione", "HMDB0000053", "STE", ["AE"], True, True, []),
    ("Testosterone", "HMDB0000234", "STE", ["AE"], True, True, []),
    ("Dihydrotestosterone", "HMDB0002961", "STE", ["AE"], True, True, ["DHT"]),
    ("Estrone", "HMDB0000145", "STE", ["AE"], True, True, []),
    ("Estradiol", "HMDB0000151", "STE", ["AE"], True, True, ["17beta-Estradiol"]),
    ("Estriol", "HMDB0000153", "STE", ["AE"], True, True, []),
    ("Androsterone", "HMDB0000031", "STE", ["AE"], True, True, []),

    # =================== Bile Acid Metabolism (19 demo) ===================
    ("Cholic acid", "HMDB0000619", "BA", ["PBA"], True, True, ["Cholate"]),
    ("Chenodeoxycholic acid", "HMDB0000518", "BA", ["PBA"], True, True, ["CDCA"]),
    ("Glycocholic acid", "HMDB0000138", "BA", ["PBA"], True, True, ["GCA"]),
    ("Taurocholic acid", "HMDB0000036", "BA", ["PBA"], True, True, ["Taurocholate"]),
    ("Glycochenodeoxycholic acid", "HMDB0000637", "BA", ["PBA"], True, True, ["GCDCA", "Chenodeoxycholic acid glycine conjugate"]),
    ("Taurochenodeoxycholic acid", "HMDB0000951", "BA", ["PBA"], True, True, ["TCDCA"]),
    ("7alpha-Hydroxycholesterol", "HMDB0001496", "BA", ["PBA"], True, True, []),
    ("7alpha-Hydroxy-4-cholesten-3-one", "HMDB0001993", "BA", ["PBA"], True, True, ["7a-Hydroxy-cholestene-3-one"]),
    ("27-Hydroxycholesterol", "HMDB0002103", "BA", ["PBA"], True, True, ["26-Hydroxycholesterol"]),
    ("Deoxycholic acid", "HMDB0000626", "BA", ["SBA"], True, True, ["DCA"]),
    ("Lithocholic acid", "HMDB0000761", "BA", ["SBA"], True, True, ["LCA"]),
    ("Ursodeoxycholic acid", "HMDB0000946", "BA", ["SBA"], True, True, ["UDCA"]),
    ("Hyodeoxycholic acid", "HMDB0000733", "BA", ["SBA"], True, True, ["HDCA"]),
    ("Glycodeoxycholic acid", "HMDB0000631", "BA", ["SBA"], True, True, ["GDCA"]),
    ("Taurodeoxycholic acid", "HMDB0000896", "BA", ["SBA"], True, True, ["TDCA"]),
    ("Glycolithocholic acid", "HMDB0000698", "BA", ["SBA"], True, True, ["GLCA"]),
    ("Taurolithocholic acid", "HMDB0000722", "BA", ["SBA"], True, True, ["TLCA", "Lithocholyltaurine"]),
    ("Glycoursodeoxycholic acid", "HMDB0000708", "BA", ["SBA"], True, True, ["GUDCA"]),
    ("Tauroursodeoxycholic acid", "HMDB0000874", "BA", ["SBA"], True, True, ["TUDCA"]),

    # ================= Sphingolipid Metabolism (20 demo) =================
    ("3-Ketosphinganine", "HMDB0001480", "SPH", ["CER"], False, True, ["3-Dehydrosphinganine", "KDS"]),
    ("Sphinganine", "HMDB0000269", "SPH", ["CER"], False, True, ["Dihydrosphingosine"]),
    ("Sphingosine", "HMDB0000252", "SPH", [], False, True, []),
    ("Sphingosine 1-phosphate", "HMDB0000277", "SPH", [], False, True, ["S1P"]),
    ("Sphinganine 1-phosphate", "HMDB0001383", "SPH", [], False, True, ["Dihydrosphingosine 1-phosphate"]),
    ("Cer(d18:0/16:0)", "HMDB0011760", "SPH", ["CER"], True, True, ["Dihydroceramide (d18:0/16:0)"]),
    ("Cer(d18:1/16:0)", "HMDB0004949", "SPH", ["CER"], True, True, ["C16 Ceramide"]),
    ("Cer(d18:1/18:0)", "HMDB0004950", "SPH", ["CER"], True, True, ["C18 Ceramide"]),
    ("Cer(d18:1/20:0)", "HMDB0004951", "SPH", ["CER"], True, True, []),
    ("Cer(d18:1/22:0)", "HMDB0004952", "SPH", ["CER"], True, True, []),
    ("Cer(d18:1/24:0)", "HMDB0004956", "SPH", ["CER"], True, True, ["C24 Ceramide"]),
    ("Cer(d18:1/24:1)", "HMDB0004953", "SPH", ["CER"], True, True, ["Cer(d18:1/24:1(15Z))", "Nervonic ceramide"]),
    ("SM(d18:1/16:0)", "HMDB0010169", "SPH", [], True, True, []),
    ("SM(d18:1/18:0)", "HMDB0001348", "SPH", [], True, True, []),
    ("SM(d18:1/18:1)", "HMDB0012101", "SPH", [], True, True, ["SM(d18:1/18:1(9Z))"]),
    ("SM(d18:1/24:0)", "HMDB0011697", "SPH", [], True, True, []),
    ("GlcCer(d18:1/16:0)", "HMDB0004971", "SPH", [], True, True, ["Glucosylceramide (d18:1/16:0)"]),
    ("GlcCer(d18:1/18:0)", "HMDB0004972", "SPH", [], True, True, ["Glucosylceramide (d18:1/18:0)"]),
    ("LacCer(d18:1/16:0)", "HMDB0006750", "SPH", [], True, True, ["Lactosylceramide (d18:1/16:0)"]),
    ("Ganglioside GM3 (d18:1/16:0)", "HMDB0004844", "SPH", [], True, True, ["GM3(d18:1/16:0)"]),

    # ===================== Lipid Metabolism (20 demo) =====================
    ("Choline", "HMDB0000097", "LIP", [], False, True, []),
    ("Phosphocholine", "HMDB0001565", "LIP", [], False, True, ["Phosphorylcholine"]),
    ("CDP-choline", "HMDB0001413", "LIP", [], False, True, ["Citicoline"]),
    ("Glycerophosphocholine", "HMDB0000086", "LIP", [], False, True, ["GPC", "alpha-Glycerophosphocholine"]),
    ("Ethanolamine", "HMDB0000149", "LIP", [], False, True, []),
    ("O-Phosphoethanolamine", "HMDB0000224", "LIP", [], False, True, ["Phosphoethanolamine"]),
    ("CDP-ethanolamine", "HMDB0001564", "LIP", [], False, True, []),
    ("LysoPC(16:0)", "HMDB0010382", "LIP", [], True, True, ["LysoPC(16:0/0:0)", "LPC 16:0"]),
    ("LysoPC(18:0)", "HMDB0010384", "LIP", [], True, True, ["LysoPC(18:0/0:0)", "LPC 18:0"]),
    ("LysoPC(18:1)", "HMDB0002815", "LIP", [], True, True, ["LysoPC(18:1/0:0)", "LysoPC(18:1(9Z)/0:0)", "LPC 18:1"]),
    ("PC(16:0/18:1)", "HMDB0007972", "LIP", [], True, True, ["PC(16:0/18:1(9Z))", "POPC"]),
    ("PC(16:0/18:2)", "HMDB0007973", "LIP", [], True, True, ["PC(16:0/18:2(9Z,12Z))"]),
    ("PC(18:0/18:2)", "HMDB0008039", "LIP", [], True, True, ["PC(18:0/18:2(9Z,12Z))"]),
    ("PE(18:0/20:4)", "HMDB0009003", "LIP", [], True, True, ["PE(18:0/20:4(5Z,8Z,11Z,14Z))"]),
    ("PA(16:0/18:1)", "HMDB0007859", "LIP", [], True, True, ["PA(16:0/18:1(9Z))"]),
    ("DG(16:0/18:1/0:0)", "HMDB0007102", "LIP", [], True, True, ["DG(16:0/18:1(9Z)/0:0)"]),
    ("TG(16:0/18:1/18:1)", "HMDB0005382", "LIP", [], True, True, ["TG(16:0/18:1(9Z)/18:1(9Z))"]),
    ("Arachidonic acid", "HMDB0001043", "LIP", ["UFA"], True, True, ["Arachidonate", "FA(20:4)"]),
    ("Linoleic acid", "HMDB0000673", "LIP", ["UFA"], True, True, ["Linoleate", "FA(18:2)"]),
    ("Docosahexaenoic acid", "HMDB0002183", "LIP", ["UFA"], True, True, ["DHA", "FA(22:6)"]),

    # ================== Fatty Acid Oxidation (25 demo) ==================
    ("L-Carnitine", "HMDB0000062", "FAO", [], False, True, ["Carnitine", "C0"]),
    ("gamma-Butyrobetaine", "HMDB0001161", "FAO", [], False, True, ["Deoxycarnitine", "4-Trimethylammoniobutanoic acid"]),
    ("L-Acetylcarnitine", "HMDB0000201", "FAO", [], False, True, ["Acetylcarnitine", "C2"]),
    ("Propionylcarnitine", "HMDB0000824", "FAO", ["BCAA"], False, True, ["C3"]),
    ("Malonylcarnitine", "HMDB0002095", "FAO", [], False, True, ["C3-DC"]),
    ("Butyrylcarnitine", "HMDB0002013", "FAO", [], False, True, ["C4"]),
    ("3-Hydroxybutyrylcarnitine", "HMDB0013127", "FAO", ["KET"], False, True, ["C4-OH"]),
    ("Isovalerylcarnitine", "HMDB0000688", "FAO", ["BCAA"], False, True, ["C5"]),
    ("Glutarylcarnitine", "HMDB0013130", "FAO", [], False, True, ["C5-DC"]),
    ("Hexanoylcarnitine", "HMDB0000756", "FAO", [], False, True, ["C6"]),
    ("Octanoylcarnitine", "HMDB0000791", "FAO", [], False, True, ["C8"]),
    ("Decanoylcarnitine", "HMDB0000651", "FAO", [], False, True, ["C10"]),
    ("Dodecanoylcarnitine", "HMDB0002250", "FAO", [], True, True, ["C12", "Lauroylcarnitine"]),
    ("Tetradecanoylcarnitine", "HMDB0005066", "FAO", [], True, True, ["C14", "Myristoylcarnitine"]),
    ("Palmitoylcarnitine", "HMDB0000222", "FAO", [], True, True, ["C16", "Hexadecanoylcarnitine"]),
    ("Stearoylcarnitine", "HMDB0000848", "FAO", [], True, True, ["C18", "Octadecanoylcarnitine"]),
    ("Oleoylcarnitine", "HMDB0005065", "FAO", [], True, True, ["C18:1"]),
    ("Linoleyl carnitine", "HMDB0006469", "FAO", [], True, True, ["Linoleylcarnitine", "C18:2"]),
    ("Palmitoyl-CoA", "HMDB0001338", "FAO", ["FAS", "SPH", "CER"], True, True, ["Palmityl-CoA"]),
    ("Malonyl-CoA", "HMDB0001175", "FAO", ["FAS"], False, True, []),
    ("3-Hydroxybutyric acid", "HMDB0000357", "FAO", ["KET"], False, True, ["beta-Hydroxybutyric acid", "beta-Hydroxybutyrate", "3-Hydroxybutyrate", "BHB"]),
    ("Acetoacetic acid", "HMDB0000060", "FAO", ["KET"], False, True, ["Acetoacetate"]),
    ("Palmitic acid", "HMDB0000220", "FAO", ["FAS"], True, True, ["Palmitate", "Hexadecanoic acid", "FA(16:0)"]),
    ("Oleic acid", "HMDB0000207", "FAO", ["UFA"], True, True, ["Oleate", "FA(18:1)"]),
    ("Myristic acid", "HMDB0000806", "FAO", ["FAS"], True, True, ["Myristate", "Tetradecanoic acid", "FA(14:0)"]),

    # ============ Reference-only metabolites (not in the demo dataset) ============
    ("L-Ornithine", "HMDB0000214", "AA", ["UREA"], False, False, ["Ornithine"]),
    ("L-Citrulline", "HMDB0000904", "AA", ["UREA"], False, False, ["Citrulline"]),
    ("L-Asparagine", "HMDB0000168", "AA", ["ADG"], False, False, ["Asparagine"]),
    ("L-Cysteine", "HMDB0000574", "AA", ["GST"], False, False, ["Cysteine"]),
    ("Taurine", "HMDB0000251", "AA", ["PBA"], False, False, []),
    ("Creatine", "HMDB0000064", "AA", [], False, False, []),
    ("Betaine", "HMDB0000043", "AA", ["GST"], False, False, ["Glycine betaine"]),
    ("Kynurenine", "HMDB0000684", "AA", [], False, False, ["L-Kynurenine"]),
    ("Serotonin", "HMDB0000259", "AA", [], False, False, ["5-Hydroxytryptamine"]),
    ("Hydroxyproline", "HMDB0000725", "AA", [], False, False, ["trans-4-Hydroxy-L-proline"]),
    ("Methylmalonic acid", "HMDB0000202", "AA", ["BCAA"], False, False, ["Methylmalonate", "MMA"]),
    ("Ribose 5-phosphate", "HMDB0001548", "GLY", ["PPP"], False, False, ["D-Ribose 5-phosphate"]),
    ("Sedoheptulose 7-phosphate", "HMDB0001068", "GLY", ["PPP"], False, False, []),
    ("Fructose 1-phosphate", "HMDB0001076", "GLY", ["FRU"], False, False, []),
    ("Glucose 1,6-bisphosphate", "HMDB0003514", "GLY", ["GAL"], False, False, ["alpha-D-Glucose 1,6-bisphosphate"]),
    ("Valerylcarnitine", "HMDB0013128", "FAO", [], False, False, []),
    ("Tiglylcarnitine", "HMDB0002366", "FAO", ["BCAA"], False, False, ["C5:1"]),
    ("Hexadecenoylcarnitine", "HMDB0013207", "FAO", [], True, False, ["C16:1"]),
    ("Tetradecenoylcarnitine", "HMDB0002014", "FAO", [], True, False, ["C14:1"]),
    ("Acetoacetyl-CoA", "HMDB0001484", "FAO", ["KET"], False, False, []),
    ("LysoPC(18:2)", "HMDB0010386", "LIP", [], True, False, ["LysoPC(18:2/0:0)", "LPC 18:2"]),
    ("LysoPC(20:4)", "HMDB0010395", "LIP", [], True, False, ["LysoPC(20:4/0:0)", "LPC 20:4"]),
    ("GlcCer(d18:1/24:0)", "HMDB0004978", "SPH", [], True, False, []),
    ("GlcCer(d18:1/22:0)", "HMDB0004974", "SPH", [], True, False, []),
    ("LacCer(d18:1/24:0)", "HMDB0011595", "SPH", [], True, False, []),
    ("Epiandrosterone", "HMDB0000365", "STE", ["AE"], True, False, []),
    ("Androsterone sulfate", "HMDB0002759", "STE", ["AE"], True, False, []),
    ("Estrone sulfate", "HMDB0001425", "STE", ["AE"], True, False, []),
    ("Tetrahydrocortisone", "HMDB0000903", "STE", ["GC"], True, False, []),
    ("Hyocholic acid", "HMDB0000760", "BA", [], True, False, []),
    ("Isoursodeoxycholic acid", "HMDB0000686", "BA", ["SBA"], True, False, []),
    ("Orotidylic acid", "HMDB0000218", "NUC", [], False, False, ["OMP", "Orotidine 5'-monophosphate"]),
    ("Deoxyadenosine monophosphate", "HMDB0000905", "PUR", [], False, False, ["dAMP"]),
    ("Deoxyadenosine triphosphate", "HMDB0001532", "PUR", [], False, False, ["dATP"]),
]

# Formerly: 27 canonical IDs taken from curator knowledge without an online re-check. All 27 have
# since been confirmed against the real HMDB export (same compound under that accession), so the
# set is now empty. Kept as a (backward-compatible) name.
HMDB_IDS_NOT_WEB_VERIFIED = set()

# ---------------------------------------------------------------------------
# 2b. Cross-references (KEGG Compound, ChEBI, PubChem CID) per HMDB accession
# ---------------------------------------------------------------------------
# Used by the Option A "ID standardization" step (multi-ID matching) and by the
# KEGG metabolite-set library. Provenance / confidence:
#   * KEGG Compound IDs: every one of the 189 IDs below was checked against the
#     KEGG REST API (rest.kegg.jp/list/...) and its KEGG name matched the
#     metabolite (four draft IDs were corrected this way, e.g. taurolithocholic
#     acid is C02592, not C03642 = its 3-sulfate). Blank = no species-level KEGG
#     entry (typical for individual lipid species such as Cer(d18:1/16:0)).
#   * ChEBI IDs: for every metabolite that has a KEGG ID, the ChEBI ID is one of
#     the ChEBI entries KEGG itself cross-links (rest.kegg.jp/conv/chebi/...).
#     The four lipid-species ChEBI IDs without a KEGG entry are listed in
#     CHEBI_IDS_NOT_WEB_VERIFIED (curator knowledge only).
#   * PubChem CIDs: 85 were confirmed against the PubChem record title (PubChem
#     compound pages / search-indexed PubChem titles); the remainder (mostly very
#     common, long-standing CIDs for amino acids, nucleotides, steroids, etc.) are
#     from curator knowledge and listed in PUBCHEM_CIDS_NOT_WEB_VERIFIED. Blank =
#     no CID asserted (not confident enough to include).
#   * METLIN IDs are deliberately NOT included: METLIN identifiers are not
#     consistently public/citable, and none could be verified here, so none are
#     asserted rather than risk fabricating them.
_XREF = {
    "HMDB0000094": ("C00158", "CHEBI:30769", "311"),
    "HMDB0000072": ("C00417", "CHEBI:32805", "643757"),
    "HMDB0000193": ("C00311", "CHEBI:30887", "1198"),
    "HMDB0003974": ("C05379", "CHEBI:7815", "972"),
    "HMDB0000208": ("C00026", "CHEBI:30915", "51"),
    "HMDB0001022": ("C00091", "CHEBI:15380", "92133"),
    "HMDB0000254": ("C00042", "CHEBI:15741", "1110"),
    "HMDB0000134": ("C00122", "CHEBI:18012", "444972"),
    "HMDB0000156": ("C00149", "CHEBI:30797", "222656"),
    "HMDB0000223": ("C00036", "CHEBI:30744", "970"),
    "HMDB0001206": ("C00024", "CHEBI:15351", "444493"),
    "HMDB0001423": ("C00010", "CHEBI:15346", "87642"),
    "HMDB0000902": ("C00003", "CHEBI:15846", "5892"),
    "HMDB0001487": ("C00004", "CHEBI:16908", "439153"),
    "HMDB0001248": ("C00016", "CHEBI:16238", "643975"),
    "HMDB0001197": ("C01352", "CHEBI:17877", "446013"),
    "HMDB0002092": ("C00490", "CHEBI:30838", "811"),
    "HMDB0000606": ("C01087", "CHEBI:32796", "439391"),
    "HMDB0000122": ("C00031", "CHEBI:4167", "5793"),
    "HMDB0001401": ("C00092", "CHEBI:4170", "5958"),
    "HMDB0001586": ("C00103", "CHEBI:29042", "439165"),
    "HMDB0000124": ("C00085", "CHEBI:15946", "69507"),
    "HMDB0001058": ("C00354", "CHEBI:37736", "172313"),
    "HMDB0001047": ("C00665", "CHEBI:28602", "105021"),
    "HMDB0001473": ("C00111", "CHEBI:16108", "668"),
    "HMDB0001112": ("C00118", "CHEBI:29052", "729"),
    "HMDB0001270": ("C00236", "CHEBI:16001", "683"),
    "HMDB0000807": ("C00197", "CHEBI:17794", "724"),
    "HMDB0000362": ("C00631", "CHEBI:17835", "59"),
    "HMDB0000263": ("C00074", "CHEBI:18021", "1005"),
    "HMDB0000243": ("C00022", "CHEBI:32816", "1060"),
    "HMDB0000190": ("C00186", "CHEBI:422", "107689"),
    "HMDB0001294": ("C01159", "CHEBI:17720", "186004"),
    "HMDB0000660": ("C00095", "CHEBI:15824", "5984"),
    "HMDB0000126": ("C00093", "CHEBI:15978", "439162"),
    "HMDB0000131": ("C00116", "CHEBI:17754", "753"),
    "HMDB0000169": ("C00159", "CHEBI:4208", "18950"),
    "HMDB0001078": ("C00275", "CHEBI:17369", "65127"),
    "HMDB0000143": ("C00124", "CHEBI:4139", "6036"),
    "HMDB0000645": ("C00446", "CHEBI:17973", "123912"),
    "HMDB0000139": ("C00258", "CHEBI:32398", "752"),
    "HMDB0001167": ("C00546", "CHEBI:17158", "880"),
    "HMDB0000161": ("C00041", "CHEBI:16977", "5950"),
    "HMDB0000517": ("C00062", "CHEBI:16467", "6322"),
    "HMDB0000191": ("C00049", "CHEBI:17053", "5960"),
    "HMDB0000148": ("C00025", "CHEBI:16015", "33032"),
    "HMDB0000641": ("C00064", "CHEBI:18050", "5961"),
    "HMDB0000123": ("C00037", "CHEBI:15428", "750"),
    "HMDB0000177": ("C00135", "CHEBI:15971", "6274"),
    "HMDB0000172": ("C00407", "CHEBI:17191", "6306"),
    "HMDB0000687": ("C00123", "CHEBI:15603", "6106"),
    "HMDB0000182": ("C00047", "CHEBI:18019", "5962"),
    "HMDB0000696": ("C00073", "CHEBI:16643", "6137"),
    "HMDB0000159": ("C00079", "CHEBI:17295", "6140"),
    "HMDB0000162": ("C00148", "CHEBI:17203", "145742"),
    "HMDB0000187": ("C00065", "CHEBI:17115", "5951"),
    "HMDB0000167": ("C00188", "CHEBI:16857", "6288"),
    "HMDB0000929": ("C00078", "CHEBI:16828", "6305"),
    "HMDB0000158": ("C00082", "CHEBI:17895", "6057"),
    "HMDB0000883": ("C00183", "CHEBI:16414", "6287"),
    "HMDB0000300": ("C00106", "CHEBI:17568", "1174"),
    "HMDB0000296": ("C00299", "CHEBI:16704", "6029"),
    "HMDB0000288": ("C00105", "CHEBI:16695", "6030"),
    "HMDB0000295": ("C00015", "CHEBI:17659", "6031"),
    "HMDB0000285": ("C00075", "CHEBI:15713", "6133"),
    "HMDB0000089": ("C00475", "CHEBI:17562", "6175"),
    "HMDB0000095": ("C00055", "CHEBI:17361", "6131"),
    "HMDB0000082": ("C00063", "CHEBI:17677", "6176"),
    "HMDB0000262": ("C00178", "CHEBI:17821", "1135"),
    "HMDB0000273": ("C00214", "CHEBI:17748", "5789"),
    "HMDB0000226": ("C00295", "CHEBI:16742", "967"),
    "HMDB0000034": ("C00147", "CHEBI:16708", "190"),
    "HMDB0000132": ("C00242", "CHEBI:16235", "135398634"),
    "HMDB0000157": ("C00262", "CHEBI:17368", "135398638"),
    "HMDB0000292": ("C00385", "CHEBI:17712", "1188"),
    "HMDB0000289": ("C00366", "CHEBI:17775", "1175"),
    "HMDB0000195": ("C00294", "CHEBI:17596", "135398641"),
    "HMDB0000050": ("C00212", "CHEBI:16335", "60961"),
    "HMDB0000133": ("C00387", "CHEBI:16750", "135398635"),
    "HMDB0000299": ("C01762", "CHEBI:18107", "64959"),
    "HMDB0000045": ("C00020", "CHEBI:16027", "6083"),
    "HMDB0001341": ("C00008", "CHEBI:16761", "6022"),
    "HMDB0000538": ("C00002", "CHEBI:15422", "5957"),
    "HMDB0000175": ("C00130", "CHEBI:17202", "135398640"),
    "HMDB0001397": ("C00144", "CHEBI:17345", "135398631"),
    "HMDB0001201": ("C00035", "CHEBI:17552", "135398619"),
    "HMDB0001273": ("C00044", "CHEBI:15996", "135398633"),
    "HMDB0001554": ("C00655", "CHEBI:15652", "73323"),
    "HMDB0000101": ("C00559", "CHEBI:17256", "13730"),
    "HMDB0000085": ("C00330", "CHEBI:17172", "135398592"),
    "HMDB0000071": ("C05512", "CHEBI:28997", "135398593"),
    "HMDB0000462": ("C01551", "CHEBI:15676", "204"),
    "HMDB0000058": ("C00575", "CHEBI:17489", "6076"),
    "HMDB0001314": ("C00942", "CHEBI:16356", "135398570"),
    "HMDB0000067": ("C00187", "CHEBI:16113", "5997"),
    "HMDB0000253": ("C01953", "CHEBI:16581", "8955"),
    "HMDB0000774": ("C18044", "CHEBI:35420", "105074"),
    "HMDB0000363": ("C05138", "CHEBI:28750", "91451"),
    "HMDB0001830": ("C00410", "CHEBI:17026", "5994"),
    "HMDB0000374": ("C01176", "CHEBI:17252", "6238"),
    "HMDB0000016": ("C03205", "CHEBI:16973", "6166"),
    "HMDB0001547": ("C02140", "CHEBI:16827", "5753"),
    "HMDB0000319": ("C01124", "CHEBI:16485", ""),
    "HMDB0000037": ("C01780", "CHEBI:27584", "5839"),
    "HMDB0000015": ("C05488", "CHEBI:28324", "440707"),
    "HMDB0000063": ("C00735", "CHEBI:17650", "5754"),
    "HMDB0002802": ("C00762", "CHEBI:16962", "222786"),
    "HMDB0000077": ("C01227", "CHEBI:28689", "5881"),
    "HMDB0001032": ("C04555", "CHEBI:16814", "12594"),
    "HMDB0000053": ("C00280", "CHEBI:16422", "6128"),
    "HMDB0000234": ("C00535", "CHEBI:17347", "6013"),
    "HMDB0002961": ("C03917", "CHEBI:16330", "10635"),
    "HMDB0000145": ("C00468", "CHEBI:17263", "5870"),
    "HMDB0000151": ("C00951", "CHEBI:16469", "5757"),
    "HMDB0000153": ("C05141", "CHEBI:27974", "5756"),
    "HMDB0000031": ("C00523", "CHEBI:16032", "5879"),
    "HMDB0000619": ("C00695", "CHEBI:16359", "221493"),
    "HMDB0000518": ("C02528", "CHEBI:16755", "10133"),
    "HMDB0000138": ("C01921", "CHEBI:17687", "10140"),
    "HMDB0000036": ("C05122", "CHEBI:28865", "6675"),
    "HMDB0000637": ("C05466", "CHEBI:36274", "12544"),
    "HMDB0000951": ("C05465", "CHEBI:16525", "387316"),
    "HMDB0001496": ("C03594", "CHEBI:17500", "107722"),
    "HMDB0001993": ("C05455", "CHEBI:17899", ""),
    "HMDB0002103": ("C15610", "CHEBI:17703", ""),
    "HMDB0000626": ("C04483", "CHEBI:28834", "222528"),
    "HMDB0000761": ("C03990", "CHEBI:16325", "9903"),
    "HMDB0000946": ("C07880", "CHEBI:9907", "31401"),
    "HMDB0000733": ("C15517", "CHEBI:52023", "5283820"),
    "HMDB0000631": ("C05464", "CHEBI:27471", "3035026"),
    "HMDB0000896": ("C05463", "CHEBI:36261", "2733768"),
    "HMDB0000698": ("C15557", "CHEBI:37998", "115245"),
    "HMDB0000722": ("C02592", "CHEBI:36259", "439763"),
    "HMDB0000708": ("", "", "12310288"),
    "HMDB0000874": ("C16868", "CHEBI:80774", "9848818"),
    "HMDB0001480": ("C02934", "CHEBI:17862", ""),
    "HMDB0000269": ("C00836", "CHEBI:16566", "91486"),
    "HMDB0000252": ("C00319", "CHEBI:16393", "5280335"),
    "HMDB0000277": ("C06124", "CHEBI:37550", "5283560"),
    "HMDB0001383": ("C01120", "CHEBI:16893", ""),
    "HMDB0011760": ("", "", ""),
    "HMDB0004949": ("", "CHEBI:72959", "5283564"),
    "HMDB0004950": ("", "CHEBI:72961", "5283565"),
    "HMDB0004951": ("", "", "5283566"),
    "HMDB0004952": ("", "", ""),
    "HMDB0004956": ("", "", ""),
    "HMDB0004953": ("", "", ""),
    "HMDB0010169": ("", "", ""),
    "HMDB0001348": ("", "", ""),
    "HMDB0012101": ("", "", ""),
    "HMDB0011697": ("", "", ""),
    "HMDB0004971": ("", "", ""),
    "HMDB0004972": ("", "", ""),
    "HMDB0006750": ("", "", ""),
    "HMDB0004844": ("", "", ""),
    "HMDB0000097": ("C00114", "CHEBI:15354", "305"),
    "HMDB0001565": ("C00588", "CHEBI:18132", "1014"),
    "HMDB0001413": ("C00307", "CHEBI:16436", "13804"),
    "HMDB0000086": ("C00670", "CHEBI:16870", "71920"),
    "HMDB0000149": ("C00189", "CHEBI:16000", "700"),
    "HMDB0000224": ("C00346", "CHEBI:17553", "1015"),
    "HMDB0001564": ("C00570", "CHEBI:16732", "123727"),
    "HMDB0010382": ("", "CHEBI:72998", "460602"),
    "HMDB0010384": ("", "", "497299"),
    "HMDB0002815": ("", "", "16081932"),
    "HMDB0007972": ("", "CHEBI:73001", "5497103"),
    "HMDB0007973": ("", "", ""),
    "HMDB0008039": ("", "", ""),
    "HMDB0009003": ("", "", ""),
    "HMDB0007859": ("", "", ""),
    "HMDB0007102": ("", "", ""),
    "HMDB0005382": ("", "", ""),
    "HMDB0001043": ("C00219", "CHEBI:15843", "444899"),
    "HMDB0000673": ("C01595", "CHEBI:17351", "5280450"),
    "HMDB0002183": ("C06429", "CHEBI:28125", "445580"),
    "HMDB0000062": ("C00318", "CHEBI:16347", "10917"),
    "HMDB0001161": ("C01181", "CHEBI:1941", "725"),
    "HMDB0000201": ("C02571", "CHEBI:57589", "7045767"),
    "HMDB0000824": ("C03017", "CHEBI:28867", "107738"),
    "HMDB0002095": ("", "", ""),
    "HMDB0002013": ("C02862", "CHEBI:7676", "213144"),
    "HMDB0013127": ("", "", ""),
    "HMDB0000688": ("C20826", "", "169235"),
    "HMDB0013130": ("", "", ""),
    "HMDB0000756": ("", "", "6426853"),
    "HMDB0000791": ("C02838", "CHEBI:18102", "11953814"),
    "HMDB0000651": ("C03299", "CHEBI:28717", "11953821"),
    "HMDB0002250": ("", "", ""),
    "HMDB0005066": ("", "", "53477791"),
    "HMDB0000222": ("C02990", "CHEBI:17490", "11953816"),
    "HMDB0000848": ("", "", "52922056"),
    "HMDB0005065": ("", "", "6441392"),
    "HMDB0006469": ("", "", "6450015"),
    "HMDB0001338": ("C00154", "CHEBI:15525", "644109"),
    "HMDB0001175": ("C00083", "CHEBI:15531", "644066"),
    "HMDB0000357": ("C01089", "CHEBI:17066", "92135"),
    "HMDB0000060": ("C00164", "CHEBI:15344", "96"),
    "HMDB0000220": ("C00249", "CHEBI:15756", "985"),
    "HMDB0000207": ("C00712", "CHEBI:16196", "445639"),
    "HMDB0000806": ("C06424", "CHEBI:28875", "11005"),
    "HMDB0000214": ("C00077", "CHEBI:15729", "6262"),
    "HMDB0000904": ("C00327", "CHEBI:16349", "9750"),
    "HMDB0000168": ("C00152", "CHEBI:17196", "6267"),
    "HMDB0000574": ("C00097", "CHEBI:17561", "5862"),
    "HMDB0000251": ("C00245", "CHEBI:15891", "1123"),
    "HMDB0000064": ("C00300", "CHEBI:16919", "586"),
    "HMDB0000043": ("C00719", "CHEBI:17750", "247"),
    "HMDB0000684": ("C00328", "CHEBI:16946", "161166"),
    "HMDB0000259": ("C00780", "CHEBI:28790", "5202"),
    "HMDB0000725": ("C01157", "CHEBI:18095", "5810"),
    "HMDB0000202": ("C02170", "CHEBI:30860", "487"),
    "HMDB0001548": ("C00117", "CHEBI:17797", "77982"),
    "HMDB0001068": ("C05382", "CHEBI:15721", "165007"),
    "HMDB0001076": ("C01094", "CHEBI:18105", ""),
    "HMDB0003514": ("C01231", "CHEBI:18148", "82400"),
    "HMDB0013128": ("", "", ""),
    "HMDB0002366": ("", "", ""),
    "HMDB0013207": ("", "", ""),
    "HMDB0002014": ("", "", ""),
    "HMDB0001484": ("C00332", "CHEBI:15345", "92153"),
    "HMDB0010386": ("", "", ""),
    "HMDB0010395": ("", "", ""),
    "HMDB0004978": ("", "", ""),
    "HMDB0004974": ("", "", ""),
    "HMDB0011595": ("", "", ""),
    "HMDB0000365": ("C07635", "CHEBI:541975", "441302"),
    "HMDB0002759": ("", "", ""),
    "HMDB0001425": ("C02538", "CHEBI:17474", ""),
    "HMDB0000903": ("C05470", "CHEBI:9901", "5866"),
    "HMDB0000760": ("C17649", "CHEBI:81244", "92805"),
    "HMDB0000686": ("C17662", "CHEBI:43419", ""),
    "HMDB0000218": ("C01103", "CHEBI:15842", "160"),
    "HMDB0000905": ("C00360", "CHEBI:17713", "12599"),
    "HMDB0001532": ("C00131", "CHEBI:16284", "15993"),
}

CHEBI_IDS_NOT_WEB_VERIFIED = {"CHEBI:72959", "CHEBI:72961", "CHEBI:72998", "CHEBI:73001"}
PUBCHEM_CIDS_NOT_WEB_VERIFIED = set(['1005', '10133', '1014', '1015', '1060', '10635', '107722', '10917', '11005', '1123', '1135', '1174', '1175', '1188', '1198', '12594', '12599', '13730', '13804', '145742', '15993', '160', '161166', '18950', '190', '204', '221493', '222528', '222786', '247', '305', '31401', '33032', '444899', '445580', '445639', '487', '5202', '5280335', '5280450', '5283560', '5753', '5754', '5756', '5757', '5789', '5810', '5839', '586', '5862', '5870', '5879', '5881', '59', '5950', '5951', '5957', '5960', '5961', '5962', '5994', '5997', '6013', '6022', '6029', '6030', '6031', '6036', '6057', '6076', '6083', '60961', '6106', '6128', '6131', '6133', '6137', '6140', '6166', '6175', '6176', '6238', '6262', '6267', '6274', '6287', '6288', '6305', '6306', '6322', '644066', '644109', '64959', '6675', '700', '71920', '724', '725', '750', '752', '753', '77982', '880', '92135', '96', '967', '9750', '985', '9903'])

METABOLITES = [
    {
        "name": m[0], "hmdb_id": m[1], "primary_pathway": _PW[m[2]][0],
        "pathways": [_PW[m[2]][0]] + [_PW[k][0] for k in m[3]],
        "lipid_class": m[4], "demo": m[5],
        "synonyms": list(m[6]),
        "kegg_id": _XREF.get(m[1], ("", "", ""))[0],
        "chebi_id": _XREF.get(m[1], ("", "", ""))[1],
        "pubchem_cid": _XREF.get(m[1], ("", "", ""))[2],
    }
    for m in _M
]

BY_HMDB = {m["hmdb_id"]: m for m in METABOLITES}
BY_KEGG = {m["kegg_id"]: m for m in METABOLITES if m["kegg_id"]}
BY_CHEBI = {m["chebi_id"]: m for m in METABOLITES if m["chebi_id"]}
BY_PUBCHEM = {m["pubchem_cid"]: m for m in METABOLITES if m["pubchem_cid"]}


# ---------------------------------------------------------------------------
# 3. Gene sets (independent, whole-pathway enzyme lists; HGNC symbols)
# ---------------------------------------------------------------------------
_GS = {
    "TCA": ["CS", "ACO1", "ACO2", "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G", "OGDH", "DLST", "DLD",
            "SUCLG1", "SUCLG2", "SUCLA2", "SDHA", "SDHB", "SDHC", "SDHD", "FH", "MDH1", "MDH2", "PC",
            "PCK1", "PCK2", "PDHA1", "PDHB", "DLAT", "ACLY"],
    "GLY": ["HK1", "HK2", "HK3", "GCK", "GPI", "PFKL", "PFKM", "PFKP", "FBP1", "FBP2", "ALDOA", "ALDOB",
            "ALDOC", "TPI1", "GAPDH", "PGK1", "PGK2", "PGAM1", "PGAM2", "ENO1", "ENO2", "ENO3", "PKM",
            "PKLR", "LDHA", "LDHB", "LDHC", "PCK1", "PCK2", "G6PC1", "PGM1", "PGM2", "BPGM", "MINPP1",
            "PDHA1", "PDHB", "DLAT"],
    "AA": ["GPT", "GPT2", "GOT1", "GOT2", "GLS", "GLS2", "GLUD1", "GLUL", "ASNS", "ASS1", "ASL", "ARG1",
           "ARG2", "OTC", "CPS1", "PAH", "TAT", "HPD", "HGD", "BCAT1", "BCAT2", "BCKDHA", "BCKDHB", "DBT",
           "IVD", "SHMT1", "SHMT2", "PHGDH", "PSAT1", "PSPH", "GLDC", "MAT1A", "MAT2A", "CBS", "CTH",
           "PYCR1", "PRODH", "ALDH18A1", "HAL", "AASS", "TDO2", "IDO1", "SDS", "OAT"],
    "LIP": ["CHKA", "CHKB", "PCYT1A", "PCYT1B", "CHPT1", "CEPT1", "ETNK1", "ETNK2", "PCYT2", "SELENOI",
            "PEMT", "PISD", "PTDSS1", "PTDSS2", "LPCAT1", "LPCAT3", "AGPAT1", "AGPAT2", "GPAM", "GPAT4",
            "LPIN1", "LPIN2", "DGAT1", "DGAT2", "PNPLA2", "LIPE", "MGLL", "LPL", "PLA2G4A", "PLD1", "PLD2",
            "DGKA", "GPD1", "GK", "LCAT"],
    "FAO": ["CPT1A", "CPT1B", "CPT2", "SLC25A20", "CRAT", "CROT", "SLC22A5", "ACADS", "ACADM", "ACADL",
            "ACADVL", "HADHA", "HADHB", "HADH", "ECHS1", "ACAA2", "ETFA", "ETFB", "ETFDH", "ACSL1", "ACOX1",
            "EHHADH", "ACAA1", "DECR1", "ECI1", "MLYCD", "BBOX1", "TMLHE", "HMGCS2", "HMGCL", "BDH1",
            "OXCT1", "ACAT1"],
    "BA": ["CYP7A1", "CYP7B1", "CYP8B1", "CYP27A1", "HSD3B7", "AKR1D1", "AKR1C4", "SLC27A5", "BAAT",
           "AMACR", "ACOX2", "HSD17B4", "SCP2", "ABCB11", "SLC10A1", "SLC10A2", "SLCO1B1", "NR1H4",
           "FGF19", "SULT2A1", "CYP3A4"],
    "SPH": ["SPTLC1", "SPTLC2", "SPTSSA", "KDSR", "CERS1", "CERS2", "CERS3", "CERS4", "CERS5", "CERS6",
            "DEGS1", "DEGS2", "ASAH1", "ASAH2", "ACER1", "ACER2", "ACER3", "SPHK1", "SPHK2", "SGPP1",
            "SGPL1", "SGMS1", "SGMS2", "SMPD1", "SMPD2", "SMPD3", "UGCG", "GBA1", "GBA2", "B4GALT5",
            "B4GALT6", "ST3GAL5", "GLB1", "GALC", "CERK", "CERT1", "UGT8"],
    "PUR": ["PPAT", "GART", "PFAS", "PAICS", "ADSL", "ATIC", "ADSS1", "ADSS2", "IMPDH1", "IMPDH2", "GMPS",
            "HPRT1", "APRT", "PNP", "ADA", "ADK", "AMPD1", "AMPD2", "AMPD3", "XDH", "GDA", "NT5C2", "NT5E",
            "GUK1", "AK1", "AK2", "NME1", "NME2", "PRPS1", "ENTPD1", "ADCY3", "PDE4B", "GMPR"],
    "NUC": ["CAD", "DHODH", "UMPS", "CMPK1", "CTPS1", "CTPS2", "NME1", "NME2", "UCK1", "UCK2", "UPP1",
            "UPP2", "CDA", "DPYD", "DPYS", "UPB1", "TYMS", "TK1", "TK2", "TYMP", "DUT", "DCK", "RRM1",
            "RRM2", "NT5C", "NT5C3A", "DTYMK"],
    "STE": ["STAR", "CYP11A1", "CYP17A1", "HSD3B1", "HSD3B2", "CYP21A2", "CYP11B1", "CYP11B2", "HSD11B1",
            "HSD11B2", "HSD17B1", "HSD17B2", "HSD17B3", "AKR1C1", "AKR1C2", "AKR1C3", "AKR1C4", "SRD5A1",
            "SRD5A2", "CYP19A1", "SULT2A1", "SULT1E1", "STS", "UGT2B7", "UGT2B17", "CYP1B1", "CYP3A4",
            "POR", "CYB5A", "FDX1", "FDXR"],
    "PPP": ["G6PD", "PGLS", "PGD", "RPIA", "RPE", "TKT", "TALDO1", "PRPS1"],
    "GAL": ["GALK1", "GALT", "GALE", "GALM", "PGM1", "UGP2", "LCT", "GLB1", "HK1"],
    "FRU": ["KHK", "ALDOB", "MPI", "PMM1", "PMM2", "PFKFB1", "PFKFB2", "PFKFB3", "PFKFB4", "SORD", "TPI1", "HK1"],
    "ADG": ["GPT", "GPT2", "GOT1", "GOT2", "GLS", "GLS2", "GLUD1", "GLUL", "ASNS", "ASS1", "ASL", "ADSL",
            "ADSS2", "CAD", "GAD1", "GFPT1", "PPAT", "CPS1"],
    "GST": ["PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2", "GLDC", "AMT", "GATM", "SDS", "SRR", "CBS", "CTH",
            "GLYCTK", "GRHPR", "BHMT", "CHDH"],
    "BCAA": ["BCAT1", "BCAT2", "BCKDHA", "BCKDHB", "DBT", "DLD", "IVD", "MCCC1", "MCCC2", "ACADSB", "HIBCH",
             "HIBADH", "PCCA", "PCCB", "MMUT"],
    "UREA": ["CPS1", "OTC", "ASS1", "ASL", "ARG1", "NAGS", "SLC25A15"],
    "KET": ["HMGCS2", "HMGCL", "BDH1", "OXCT1", "ACAT1", "SLC16A1"],
    "FAS": ["ACACA", "ACACB", "FASN", "OXSM", "MCAT", "ELOVL6", "SCD", "ACSL1", "ACSL3"],
    "UFA": ["FADS1", "FADS2", "ELOVL2", "ELOVL5", "SCD", "ACOX1"],
    "PBA": ["CYP7A1", "CYP7B1", "CYP8B1", "CYP27A1", "HSD3B7", "AKR1D1", "AKR1C4", "SLC27A5", "BAAT",
            "AMACR", "ACOX2", "HSD17B4", "SCP2", "CH25H", "CYP39A1"],
    "OXP": ["NDUFS1", "NDUFV1", "NDUFA9", "SDHA", "SDHB", "SDHC", "SDHD", "UQCRC1", "CYC1", "COX4I1",
            "ATP5F1A", "ATP5F1B", "ATP5PO"],
    "GC": ["STAR", "CYP11A1", "HSD3B2", "CYP21A2", "CYP11B1", "CYP11B2", "CYP17A1", "HSD11B1", "HSD11B2",
           "POR", "FDX1", "FDXR"],
    "AE": ["CYP17A1", "HSD3B2", "HSD17B1", "HSD17B2", "HSD17B3", "AKR1C3", "SRD5A1", "SRD5A2", "CYP19A1",
           "SULT2A1", "SULT1E1", "STS", "UGT2B17", "UGT2B7", "CYB5A", "CYP1B1"],
    "CER": ["SPTLC1", "SPTLC2", "SPTSSA", "KDSR", "CERS1", "CERS2", "CERS3", "CERS4", "CERS5", "CERS6", "DEGS1"],
    # SBA intentionally omitted: secondary bile acids are produced by gut microbes,
    # so there is no meaningful human biosynthetic gene set for it.
}

GENE_SETS = {_PW[k][0]: sorted(set(v)) for k, v in _GS.items()}
GENE_SET_COLLECTION_NAME = "Curated Metabolic Pathway Gene Sets (embedded, KEGG-style, human)"


# ---------------------------------------------------------------------------
# Name matching
# ---------------------------------------------------------------------------
_STEREO_PREFIX = re.compile(r"^(l|d|dl|sn)-", re.IGNORECASE)


def normalize_name(name) -> str:
    """Case/punctuation-insensitive key: 'L-Lactic acid' -> 'llacticacid'."""
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def _name_variants(name: str):
    """Automatic variants: with/without a stereo prefix, and '-ic acid' <-> '-ate'."""
    out = {name}
    stripped = _STEREO_PREFIX.sub("", name)
    out.add(stripped)
    for n in list(out):
        m = re.match(r"^(.*)ic acid$", n, flags=re.IGNORECASE)
        if m:
            out.add(m.group(1) + "ate")
    return out


def _build_name_index():
    index = {}
    for m in METABOLITES:
        for label in [m["name"]] + m["synonyms"]:
            for v in _name_variants(label):
                index.setdefault(normalize_name(v), m["hmdb_id"])
        index.setdefault(normalize_name(m["hmdb_id"]), m["hmdb_id"])
    return index


NAME_INDEX = _build_name_index()


def lookup_by_name(name):
    """Return the reference record for a metabolite name/synonym, or None."""
    if name is None:
        return None
    for v in _name_variants(str(name).strip()):
        hid = NAME_INDEX.get(normalize_name(v))
        if hid:
            return BY_HMDB[hid]
    return None


# ---------------------------------------------------------------------------
# Demo-name assignment (used by generate_rich_demo_data.py)
# ---------------------------------------------------------------------------
def demo_names_by_pathway() -> dict:
    """{primary pathway: [(name, lipid_class), ...]} for demo=True metabolites, in curated order."""
    out = {}
    for m in METABOLITES:
        if m["demo"]:
            out.setdefault(m["primary_pathway"], []).append((m["name"], m["lipid_class"]))
    return out


def assign_demo_names(pathways, methods):
    """
    Deterministically assign one real, unique metabolite name to each demo row,
    given that row's (already fixed) Pathway and Method. Rows measured by a
    lipidomics method take lipid-class names first (e.g. ceramides, acyl-
    carnitines, PCs), other rows take polar metabolites first; each queue falls
    back to the other once exhausted. Consumes no random numbers, so it never
    perturbs the demo generator's RNG stream.
    """
    pools = demo_names_by_pathway()
    queues = {pw: {"lipid": [n for n, lip in names if lip], "polar": [n for n, lip in names if not lip]}
              for pw, names in pools.items()}
    assigned = []
    for pw, method in zip(pathways, methods):
        q = queues[pw]
        first, second = ("lipid", "polar") if "Lipid" in str(method) else ("polar", "lipid")
        pool = q[first] if q[first] else q[second]
        if not pool:
            raise ValueError(f"Not enough curated demo metabolite names for pathway '{pw}'.")
        assigned.append(pool.pop(0))
    return assigned


def demo_hmdb_lookup() -> dict:
    return {m["name"]: m["hmdb_id"] for m in METABOLITES if m["demo"]}
