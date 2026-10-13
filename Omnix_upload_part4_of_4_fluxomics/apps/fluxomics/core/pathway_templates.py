"""
core.pathway_templates — hand-drawn pathway layouts (v2.5).

Pathways 1–10 are the layouts of build_pathways.py (Central-carbon pathway SVG templates): positions, reaction
arrows, cofactor labels, links to other pathways, sub-panels and the TCA ring are copied unchanged. Pathways 11 + are
drawn in the same style. Coordinates are template pixels (x right, y DOWN), exactly as in the SVG templates;
core.pathmap expands them so that every metabolite can carry its isotope-incorporation / Z-score card.

Metabolites are matched to the data by KEGG compound ID (name as fallback). Nodes marked verify=True are generic
class intermediates (e.g. 3-hydroxyacyl-CoA) that are matched only to annotations you confirmed yourself.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass
class Template:
    key: str
    title: str
    fill: str
    border: str
    arrow: str
    annot: tuple = ()                     # annotation pathway names whose other detected members are listed below
    nodes: dict = field(default_factory=dict)
    edges: list = field(default_factory=list)
    circles: list = field(default_factory=list)   # (cx, cy, r, [text lines])
    boxes: list = field(default_factory=list)     # (x, y, w, h, title) sub-panels
    notes: list = field(default_factory=list)     # (x, y, text, rotation)

    # ---- builder API (same as build_pathways.py)
    def met(self, nid, label, kegg, x, y, verify=False):
        self.nodes[nid] = dict(kind="met", label=label, kegg=kegg, x=float(x), y=float(y), verify=verify)

    def stub(self, nid, label, x, y):
        self.nodes[nid] = dict(kind="stub", label=label, x=float(x), y=float(y))

    def edge(self, a, b, style="solid", via=(), label=None, color=None, bi=False, side=1):
        self.edges.append(dict(a=a, b=b, style=style, via=[(float(x), float(y)) for x, y in via], label=label,
                               color=color, bi=bi, side=side))

    def arc(self, a, b, cx, cy, r, a1, a2, bi=False):
        self.edges.append(dict(a=a, b=b, arc=(cx, cy, r, a1, a2), bi=bi, style="solid", color=None, label=None,
                               via=[], side=1))

    def circle(self, cx, cy, r, lines=()):
        self.circles.append((cx, cy, r, list(lines)))

    def box(self, x, y, w, h, title):
        self.boxes.append((x, y, w, h, title))

    def note(self, x, y, text, rot=0):
        self.notes.append((x, y, text, rot))


T: dict[str, Template] = {}


def _add(t: Template) -> Template:
    T[t.key] = t
    return t


# =============================================================================================================
# 1–10: layouts from build_pathways.py
# =============================================================================================================
# 1. Glycolysis
g = _add(Template("glycolysis", "Glycolysis", "#FDEEF1", "#F4BCC9", "#D6336C", ("Glycolysis",)))
X = 290
g.met("glc", "Glucose", "C00031", X, 95)
g.met("g6p", "Glucose-6P", "C00092", X, 175)
g.met("f6p", "Fructose-6P", "C00085", X, 255)
g.met("fbp", "Fructose-1,6BP", "C00354", X, 335)
g.met("dhap", "DHAP", "C00111", 120, 415)
g.met("gap", "Glyceraldehyde-3P", "C00118", X, 415)
g.met("bpg", "1,3-BPG", "C00236", X, 495)
g.met("pg3", "3-PG", "C00197", X, 575)
g.met("pg2", "2-PG", "C00631", X, 655)
g.met("pep", "PEP", "C00074", X, 735)
g.met("pyr", "Pyruvate", "C00022", X, 820)
g.met("lac", "Lactate", "C00186", 110, 820)
g.met("ser", "Serine", "C00065", 110, 575)
g.stub("ppp_out", "→ PPP", 470, 175)
g.stub("ppp_in", "← PPP", 470, 335)
g.stub("tca", "→ Acetyl-CoA / TCA", 455, 820)
for a, b, lab in [("glc", "g6p", "ATP"), ("g6p", "f6p", None), ("f6p", "fbp", "ATP"), ("fbp", "gap", None),
                  ("gap", "bpg", "NADH"), ("bpg", "pg3", "ATP"), ("pg3", "pg2", None), ("pg2", "pep", None),
                  ("pep", "pyr", "ATP")]:
    g.edge(a, b, label=lab, bi=b in ("f6p", "gap", "bpg", "pg3", "pg2", "pep"))
g.edge("fbp", "dhap", bi=True)
g.edge("dhap", "gap", bi=True)
g.edge("pyr", "lac", bi=True, label="NADH")
g.edge("pg3", "ser", label="serine synthesis", style="dashed")
g.edge("g6p", "ppp_out", style="dashed", color="#1C7ED6")
g.edge("ppp_in", "f6p", style="dashed", color="#1C7ED6", via=[(470, 255)])
g.edge("ppp_in", "gap", style="dashed", color="#1C7ED6", via=[(470, 415)])
g.edge("pyr", "tca", style="dashed", color="#1864AB")

# 2. Pentose phosphate pathway
p = _add(Template("ppp", "Pentose Phosphate Pathway (PPP)", "#EAF4FD", "#A9CFF0", "#1C7ED6",
                  ("Pentose phosphate pathway",)))
p.met("g6p", "Glucose-6P", "C00092", 110, 120)
p.met("pgl", "6-P-Gluconolactone", "C01236", 330, 120)
p.met("pg6", "6-Phosphogluconate", "C00345", 560, 120)
p.met("ru5p", "Ribulose-5P", "C00199", 760, 120)
p.met("r5p", "Ribose-5P", "C00117", 610, 250)
p.met("xu5p", "Xylulose-5P", "C00231", 760, 340)
p.met("s7p", "Sedoheptulose-7P", "C05382", 400, 340)
p.met("gap", "Glyceraldehyde-3P", "C00118", 610, 440)
p.met("e4p", "Erythrose-4P", "C00279", 150, 340)
p.met("f6p", "Fructose-6P", "C00085", 380, 440)
p.stub("nuc", "→ Nucleotide synthesis (PRPP)", 380, 250)
p.stub("gly", "→ Glycolysis", 150, 440)
p.edge("g6p", "pgl", label="NADPH")
p.edge("pgl", "pg6")
p.edge("pg6", "ru5p", label="NADPH, CO₂")
p.edge("ru5p", "r5p", bi=True)
p.edge("ru5p", "xu5p", bi=True)
p.edge("r5p", "s7p", bi=True, label="TKT", side=-1)
p.edge("xu5p", "gap", bi=True)
p.edge("s7p", "e4p", bi=True, label="TALDO")
p.edge("s7p", "f6p", bi=True, style="dashed")
p.edge("r5p", "nuc", style="dashed")
p.edge("f6p", "gly", style="dashed", color="#D6336C")

# 3. TCA cycle
t = _add(Template("tca", "TCA Cycle", "#EEF7FE", "#A5D0F3", "#1864AB", ("TCA cycle",)))
CX, CY, R = 400, 370, 200
t.circle(CX, CY, R + 38, ["TCA", "Cycle"])
ring = [("cit", "Citrate", "C00158", 90), ("cisaco", "cis-Aconitate", "C00417", 50),
        ("icit", "Isocitrate", "C00311", 10), ("akg", "α-Ketoglutarate", "C00026", -30),
        ("sca", "Succinyl-CoA", "C00091", -75), ("suc", "Succinate", "C00042", -120),
        ("fum", "Fumarate", "C00122", -165), ("mal", "Malate", "C00149", 160),
        ("oaa", "Oxaloacetate", "C00036", 125)]
for k, lab, kegg, ang in ring:
    t.met(k, lab, kegg, CX + R * math.cos(math.radians(ang)), CY - R * math.sin(math.radians(ang)))
fix = [90, 50, 10, -30, -75, -120, -165, -200, -235, -270]
bis = {"cisaco", "icit", "fum", "mal", "oaa"}
for i in range(len(ring)):
    a, b = ring[i][0], ring[(i + 1) % len(ring)][0]
    t.arc(a, b, CX, CY, R, fix[i], fix[i + 1], bi=(b in bis and b != "cit"))
t.met("pyr", "Pyruvate", "C00022", 110, 85)
t.met("acoa", "Acetyl-CoA", "C00024", 300, 85)
t.met("glu", "Glutamate", "C00025", 800, 499)
t.met("asp", "Aspartate", "C00049", 20, 206)
t.stub("fas", "→ Fatty acid synthesis", 610, 85)
t.stub("bcaa", "← Val / Ile", 720, 610)
t.edge("pyr", "acoa", label="PDH")
t.edge("acoa", "cit", via=[(380, 85)], color="#E8590C")
t.edge("cit", "fas", style="dashed", color="#E8590C", via=[(470, 85)])
t.edge("glu", "akg", bi=True, color="#C92A2A")
t.edge("asp", "oaa", bi=True, color="#C92A2A")
t.edge("bcaa", "sca", style="dashed", color="#3B5BDB")

# 4. One-carbon / folate & methionine cycle
o = _add(Template("one_carbon", "One-Carbon / Folate & Methionine Metabolism", "#F4EFFC", "#C9B6EE", "#7048E8",
                  ("One-carbon metabolism",)))
o.met("ser", "Serine", "C00065", 90, 130)
o.met("gly", "Glycine", "C00037", 250, 130)
o.met("ch2", "5,10-CH₂-THF", "C00143", 450, 130)
o.met("ch3", "5-CH₃-THF", "C00440", 680, 130)
o.met("thf", "THF", "C00101", 560, 280)
o.met("cho", "10-Formyl-THF", "C00234", 400, 400)
o.met("for", "Formate", "C00058", 200, 400)
o.met("met", "Methionine", "C00073", 860, 130)
o.met("sam", "SAM", "C00019", 1050, 130)
o.met("sah", "SAH", "C00021", 1050, 300)
o.met("hcy", "Homocysteine", "C00155", 860, 300)
o.met("cth", "Cystathionine", "C02291", 860, 430)
o.met("cys", "Cysteine", "C00097", 1050, 430)
o.stub("dtmp", "→ dTMP", 450, 30)
o.stub("pur", "→ Purines", 620, 400)
o.edge("ser", "gly", bi=True, label="SHMT")
o.edge("thf", "ch2", style="dashed", label="+ C1 from Ser", side=-1)
o.edge("gly", "ch2", style="dashed", label="GCS", via=[(250, 200), (400, 200)])
o.edge("ch2", "ch3", label="MTHFR")
o.edge("ch3", "thf", via=[(680, 280)])
o.edge("ch2", "dtmp", style="dashed")
o.edge("for", "cho", label="+ THF")
o.edge("thf", "cho", bi=True, style="dashed")
o.edge("cho", "pur", style="dashed")
o.edge("ch3", "met", style="dashed", label="CH₃ to Hcy", color="#9C36B5")
o.edge("met", "sam", label="ATP")
o.edge("sam", "sah", label="methylation")
o.edge("sah", "hcy", bi=True)
o.edge("hcy", "met", label="MTR", side=-1)
o.edge("hcy", "cth", label="+ Ser")
o.edge("cth", "cys")

# 5. Purine synthesis (de novo)
u = _add(Template("purine", "Purine Synthesis (de novo)", "#EBFAF3", "#A6E3C8", "#0B7285",
                  ("Nucleotide biosynthesis",)))
row1 = [("r5p", "Ribose-5P", "C00117"), ("prpp", "PRPP", "C00119"), ("pra", "PRA", "C03090"),
        ("gar", "GAR", "C03838"), ("fgar", "FGAR", "C04376"), ("fgam", "FGAM", "C04640")]
row2 = [("air", "AIR", "C03373"), ("cair", "CAIR", "C04751"), ("saicar", "SAICAR", "C04823"),
        ("aicar", "AICAR", "C04677"), ("faicar", "FAICAR", "C04734"), ("imp", "IMP", "C00130")]
xs = [110, 270, 430, 590, 750, 900]
for (k, lab, c), x in zip(row1, xs):
    u.met(k, lab, c, x, 120)
for (k, lab, c), x in zip(row2, xs[::-1]):
    u.met(k, lab, c, x, 260)
labs1 = [None, "Gln", "Gly", "10-CHO-THF", "Gln"]
for i in range(5):
    u.edge(row1[i][0], row1[i + 1][0], label=labs1[i])
u.edge("fgam", "air")
labs2 = ["CO₂", "Asp", None, "10-CHO-THF", None]
for i in range(5):
    u.edge(row2[i][0], row2[i + 1][0], label=labs2[i])
u.met("samp", "Adenylosuccinate", "C03794", 110, 390)
u.met("amp", "AMP", "C00020", 110, 500)
u.met("xmp", "XMP", "C00655", 380, 390)
u.met("gmp", "GMP", "C00144", 380, 500)
u.stub("ppp", "← PPP", 110, 30)
u.edge("ppp", "r5p", style="dashed", color="#1C7ED6")
u.edge("imp", "samp", label="Asp")
u.edge("samp", "amp")
u.edge("imp", "xmp", via=[(200, 330), (380, 330)], label="NAD⁺", side=-1)
u.edge("xmp", "gmp", label="Gln")

# 6. Pyrimidine synthesis (de novo)
y = _add(Template("pyrimidine", "Pyrimidine Synthesis (de novo)", "#F1EEFC", "#C5BCEF", "#6741D9"))
y.met("gln", "Glutamine", "C00064", 100, 130)
y.met("cp", "Carbamoyl-P", "C00169", 280, 130)
y.met("ca", "Carbamoyl-Asp", "C00438", 470, 130)
y.met("dho", "Dihydroorotate", "C00337", 660, 130)
y.met("oro", "Orotate", "C00295", 840, 130)
y.met("asp", "Aspartate", "C00049", 470, 230)
y.met("omp", "OMP", "C01103", 840, 260)
y.met("ump", "UMP", "C00105", 660, 260)
y.met("udp", "UDP", "C00015", 470, 340)
y.met("utp", "UTP", "C00075", 280, 260)
y.met("ctp", "CTP", "C00063", 100, 260)
y.met("dump", "dUMP", "C00365", 660, 380)
y.met("dtmp", "dTMP", "C00364", 840, 380)
y.edge("gln", "cp", label="HCO₃⁻, 2 ATP")
y.edge("cp", "ca")
y.edge("asp", "ca", style="dashed")
y.edge("ca", "dho", bi=True)
y.edge("dho", "oro", label="DHODH")
y.edge("oro", "omp", label="PRPP")
y.edge("omp", "ump", label="CO₂")
y.edge("ump", "udp", via=[(660, 340)])
y.edge("udp", "utp", via=[(280, 340)])
y.edge("utp", "ctp", label="Gln")
y.edge("udp", "dump", style="dashed", via=[(470, 380)])
y.edge("dump", "dtmp", label="5,10-CH₂-THF")

# 7. Glutathione
h = _add(Template("glutathione", "Glutathione Metabolism", "#FDEEEE", "#F3B5B5", "#E03131",
                  ("Redox metabolism",)))
h.met("glu", "Glutamate", "C00025", 120, 110)
h.met("cys", "Cysteine", "C00097", 360, 110)
h.met("gec", "γ-Glu-Cys", "C00669", 240, 210)
h.met("gly", "Glycine", "C00037", 470, 210)
h.met("gsh", "GSH", "C00051", 240, 320)
h.met("gssg", "GSSG", "C00127", 500, 320)
h.met("nadph", "NADPH", "C00005", 500, 420)
h.stub("ts", "← Transsulfuration", 570, 110)
h.edge("glu", "gec", label="GCL", side=-1)
h.edge("cys", "gec")
h.edge("gec", "gsh", label="GSS", side=-1)
h.edge("gly", "gsh", style="dashed", via=[(470, 265), (240, 265)])
h.edge("gsh", "gssg", label="GPX · ROS")
h.edge("gssg", "gsh", label="GSR", via=[(500, 365), (240, 365)], color="#2B8A3E", side=-1)
h.edge("nadph", "gssg", style="dashed", color="#2B8A3E", label="from PPP")
h.edge("ts", "cys", style="dashed", color="#7048E8")

# 8. Glutamine / glutamate
q = _add(Template("glutamine", "Glutamine / Glutamate Metabolism", "#FDEDF1", "#F4BCC9", "#C92A2A",
                  ("Glutamine metabolism", "Glutamate metabolism")))
q.met("gln", "Glutamine", "C00064", 110, 130)
q.met("glu", "Glutamate", "C00025", 330, 130)
q.met("akg", "α-Ketoglutarate", "C00026", 560, 130)
q.met("hg2", "2-Hydroxyglutarate", "C02630", 560, 270)
q.met("pro", "Proline", "C00148", 330, 270)
q.stub("tca", "→ TCA", 740, 130)
q.stub("gsh", "→ GSH", 110, 270)
q.edge("gln", "glu", label="GLS")
q.edge("glu", "akg", bi=True, label="GLUD / transaminases")
q.edge("akg", "tca", style="dashed", color="#1864AB")
q.edge("akg", "hg2", style="dashed", label="IDH-mut / LDH")
q.edge("glu", "pro", bi=True, label="via P5C")
q.edge("glu", "gsh", style="dashed", color="#E03131")

# 9. Fatty acid metabolism
f = _add(Template("fatty_acid", "Fatty Acid Metabolism", "#FFF9DB", "#F2D27A", "#E8590C",
                  ("Fatty acid metabolism",)))
f.box(30, 60, 300, 550, "Fatty Acid Synthesis")
f.box(520, 60, 350, 550, "β-Oxidation")
f.met("cit", "Citrate", "C00158", 180, 135)
f.met("acoa", "Acetyl-CoA", "C00024", 180, 235)
f.met("mcoa", "Malonyl-CoA", "C00083", 180, 335)
f.met("pal", "Palmitate", "C00249", 180, 445)
f.met("lcfa", "Long-chain fatty acid", "C00638", 690, 135, verify=True)
f.met("acyl", "Acyl-CoA", "C00040", 690, 215)
f.met("enoyl", "trans-2-Enoyl-CoA", "C00658", 690, 295, verify=True)
f.met("hacyl", "3-Hydroxyacyl-CoA", "C01144", 690, 375, verify=True)
f.met("kacyl", "3-Ketoacyl-CoA", "C00264", 690, 455, verify=True)
f.met("acoa2", "Acetyl-CoA", "C00024", 690, 545)
f.edge("cit", "acoa", label="ACLY")
f.edge("acoa", "mcoa", label="ACC")
f.edge("mcoa", "pal", label="FASN · NADPH")
f.edge("pal", "lcfa", via=[(425, 445), (425, 135)], style="dashed", color="#A0711C")
f.edge("lcfa", "acyl")
f.edge("acyl", "enoyl", label="FAD")
f.edge("enoyl", "hacyl")
f.edge("hacyl", "kacyl", label="NAD⁺")
f.edge("kacyl", "acoa2", label="CoA")
f.edge("kacyl", "acyl", style="dashed", via=[(830, 455), (830, 215)])
f.note(846, 335, "Acyl-CoA (n−2) · next cycle", 90)

# 10. Amino acid links to TCA
m = _add(Template("aa_to_tca", "Amino Acid Links to TCA", "#EEF0FC", "#B7BDF0", "#3B5BDB",
                  ("Amino acid metabolism",)))
AA, MID, TC = 110, 330, 560
m.met("ala", "Alanine", "C00041", AA, 100); m.met("pyr", "Pyruvate", "C00022", TC, 100)
m.met("asn", "Asparagine", "C00152", AA, 180); m.met("asp", "Aspartate", "C00049", MID, 180)
m.met("oaa", "Oxaloacetate", "C00036", TC, 180)
m.met("gln", "Glutamine", "C00064", AA, 260); m.met("pro", "Proline", "C00148", AA, 330)
m.met("arg", "Arginine", "C00062", AA, 400); m.met("glu", "Glutamate", "C00025", MID, 330)
m.met("akg", "α-Ketoglutarate", "C00026", TC, 330)
m.met("val", "Valine", "C00183", AA, 480); m.met("ile", "Isoleucine", "C00407", AA, 550)
m.met("ppc", "Propionyl-CoA", "C00100", MID, 515); m.met("sca", "Succinyl-CoA", "C00091", TC, 515)
m.met("leu", "Leucine", "C00123", AA, 620); m.met("acoa", "Acetyl-CoA", "C00024", TC, 620)
m.met("phe", "Phenylalanine", "C00079", AA, 710); m.met("tyr", "Tyrosine", "C00082", MID, 710)
m.met("fum", "Fumarate", "C00122", TC, 710)
for a, b in [("ala", "pyr"), ("asn", "asp"), ("asp", "oaa"), ("gln", "glu"), ("pro", "glu"), ("arg", "glu"),
             ("glu", "akg"), ("val", "ppc"), ("ile", "ppc"), ("ppc", "sca"), ("leu", "acoa"),
             ("phe", "tyr"), ("tyr", "fum")]:
    m.edge(a, b, bi=(a, b) in {("ala", "pyr"), ("asp", "oaa"), ("glu", "akg")})
m.edge("ile", "acoa", style="dashed")

# =============================================================================================================
# 11 +: additional pathways in the same style
# =============================================================================================================
# 11. Serine / glycine synthesis
s = _add(Template("serine", "Serine & Glycine Synthesis", "#FFF4E6", "#FFD8A8", "#D9480F"))
s.stub("gly_in", "← Glycolysis", 110, 40)
s.met("pg3", "3-PG", "C00197", 110, 130)
s.met("php", "3-Phosphohydroxypyruvate", "C03232", 380, 130)
s.met("pser", "Phosphoserine", "C01005", 650, 130)
s.met("ser", "Serine", "C00065", 900, 130)
s.met("gly", "Glycine", "C00037", 900, 270)
s.met("cys", "Cysteine", "C00097", 650, 270)
s.stub("oc", "→ One-carbon (5,10-CH₂-THF)", 380, 270)
s.stub("gsh", "→ GSH / purines", 900, 380)
s.edge("gly_in", "pg3", style="dashed", color="#D6336C")
s.edge("pg3", "php", label="PHGDH · NAD⁺")
s.edge("php", "pser", label="PSAT1 · Glu→αKG")
s.edge("pser", "ser", label="PSPH")
s.edge("ser", "gly", bi=True, label="SHMT1/2")
s.edge("ser", "cys", style="dashed", label="transsulfuration", via=[(780, 200), (650, 200)])
s.edge("gly", "oc", style="dashed", via=[(780, 270)], label="GCS")
s.edge("gly", "gsh", style="dashed")

# 12. Urea cycle & polyamines
w = _add(Template("urea", "Urea Cycle & Polyamines", "#F3F9EC", "#C0DFA0", "#2B8A3E",
                  ("Urea cycle & polyamines",)))
w.met("cp", "Carbamoyl-P", "C00169", 120, 130)
w.met("citr", "Citrulline", "C00327", 400, 130)
w.met("argsuc", "Argininosuccinate", "C03406", 680, 130)
w.met("asp", "Aspartate", "C00049", 940, 130)
w.met("fum", "Fumarate", "C00122", 940, 240)
w.met("arg", "Arginine", "C00062", 680, 340)
w.met("urea", "Urea", "C00086", 540, 240)
w.met("orn", "Ornithine", "C00077", 400, 340)
w.met("put", "Putrescine", "C00134", 400, 470)
w.met("acput", "N-Acetylputrescine", "C02714", 120, 470)
w.met("spd", "Spermidine", "C00315", 680, 470)
w.met("spm", "Spermine", "C00750", 940, 470)
w.stub("nh4", "← NH₄⁺ + HCO₃⁻ (CPS1)", 120, 40)
w.stub("p5c", "↔ Proline / Glutamate (P5C)", 120, 340)
w.stub("tca", "→ TCA", 940, 340)
w.edge("nh4", "cp", style="dashed")
w.edge("cp", "citr", label="OTC")
w.edge("orn", "citr")
w.edge("citr", "argsuc", label="ASS1 · ATP")
w.edge("asp", "argsuc")
w.edge("argsuc", "arg", label="ASL")
w.edge("argsuc", "fum", via=[(800, 240)])
w.edge("fum", "tca", style="dashed", color="#1864AB")
w.edge("arg", "orn", label="ARG1", side=-1)
w.edge("arg", "urea", style="dashed")
w.edge("orn", "p5c", bi=True, style="dashed", label="OAT")
w.edge("orn", "put", label="ODC1")
w.edge("put", "spd", label="SRM · dcSAM")
w.edge("spd", "spm", label="SMS")
w.edge("put", "acput", style="dashed", label="SAT1")

# 13. Creatine
c = _add(Template("creatine", "Creatine Metabolism", "#F8F0FC", "#E1C2F0", "#9C36B5"))
c.met("arg", "Arginine", "C00062", 110, 130)
c.met("gly", "Glycine", "C00037", 110, 280)
c.met("orn", "Ornithine", "C00077", 380, 50)
c.met("gaa", "Guanidinoacetate", "C00581", 380, 205)
c.met("cr", "Creatine", "C00300", 660, 205)
c.met("pcr", "Phosphocreatine", "C02305", 660, 340)
c.met("crn", "Creatinine", "C00791", 940, 205)
c.edge("arg", "gaa", label="AGAT")
c.edge("gly", "gaa")
c.edge("arg", "orn", style="dashed")
c.edge("gaa", "cr", label="GAMT · SAM→SAH")
c.edge("cr", "pcr", bi=True, label="CK · ATP")
c.edge("cr", "crn", style="dashed", label="non-enzymatic")

# 14. Hexosamine & nucleotide sugars
x_ = _add(Template("hexosamine", "Hexosamine & Nucleotide Sugars", "#FFF0F6", "#FCC2D7", "#C2255C",
                   ("Hexosamine & nucleotide sugars",)))
x_.met("f6p", "Fructose-6P", "C00085", 120, 130)
x_.met("gln6p", "Glucosamine-6P", "C00352", 360, 130)
x_.met("glcnac6p", "GlcNAc-6P", "C00357", 600, 130)
x_.met("glcnac1p", "GlcNAc-1P", "C04256", 840, 130)
x_.met("udpglcnac", "UDP-GlcNAc", "C00043", 840, 260)
x_.met("udpgalnac", "UDP-GalNAc", "C00203", 600, 260)
x_.met("neu5ac", "Neu5Ac", "C00270|C19910", 840, 390)
x_.met("cmpneu", "CMP-Neu5Ac", "C00128", 600, 390)
x_.met("g6p", "Glucose-6P", "C00092", 120, 260)
x_.met("g1p", "Glucose-1P", "C00103", 120, 390)
x_.met("udpglc", "UDP-Glucose", "C00029", 360, 390)
x_.met("udpgal", "UDP-Galactose", "C00052", 360, 510)
x_.met("udpglca", "UDP-Glucuronate", "C00167", 600, 510)
x_.met("udpxyl", "UDP-Xylose", "C00190", 840, 510)
x_.edge("f6p", "gln6p", label="GFAT · Gln")
x_.edge("gln6p", "glcnac6p", label="Acetyl-CoA")
x_.edge("glcnac6p", "glcnac1p")
x_.edge("glcnac1p", "udpglcnac", label="UTP")
x_.edge("udpglcnac", "udpgalnac", bi=True, label="GALE")
x_.edge("udpglcnac", "neu5ac", style="dashed", label="GNE · NANS")
x_.edge("neu5ac", "cmpneu", label="CTP")
x_.edge("g6p", "f6p", bi=True)
x_.edge("g6p", "g1p", bi=True, label="PGM")
x_.edge("g1p", "udpglc", label="UTP")
x_.edge("udpglc", "udpgal", bi=True, label="GALE")
x_.edge("udpglc", "udpglca", label="UGDH · 2 NAD⁺")
x_.edge("udpglca", "udpxyl", label="UXS1")

# 15. Phospholipid head groups (Kennedy pathway)
k_ = _add(Template("phospholipid", "Phospholipid Head Groups (Kennedy Pathway)", "#E6FCF5", "#96F2D7", "#087F5B",
                   ("Phospholipid head groups",)))
k_.met("cho", "Choline", "C00114", 110, 130)
k_.met("pcho", "Phosphocholine", "C00588", 360, 130)
k_.met("cdpcho", "CDP-Choline", "C00307", 610, 130)
k_.met("pc", "Phosphatidylcholine", "C00157", 860, 130)
k_.met("gpc", "Glycerophosphocholine", "C00670", 1110, 130)
k_.met("bet", "Betaine", "C00719", 110, 260)
k_.met("dmg", "Dimethylglycine", "C01026", 360, 260)
k_.met("sar", "Sarcosine", "C00213", 610, 260)
k_.met("etn", "Ethanolamine", "C00189", 110, 390)
k_.met("petn", "Phosphoethanolamine", "C00346", 360, 390)
k_.met("cdpetn", "CDP-Ethanolamine", "C00570", 610, 390)
k_.met("pe", "Phosphatidylethanolamine", "C00350", 860, 390)
k_.met("gpe", "Glycerophosphoethanolamine", "C01233", 1110, 390)
k_.edge("cho", "pcho", label="CHK")
k_.edge("pcho", "cdpcho", label="PCYT1 · CTP")
k_.edge("cdpcho", "pc", label="CHPT1 · DAG")
k_.edge("pc", "gpc", style="dashed", label="PLA2 / PLB")
k_.edge("cho", "bet", label="CHDH · BADH")
k_.edge("bet", "dmg", label="BHMT · Hcy→Met")
k_.edge("dmg", "sar")
k_.edge("etn", "petn", label="ETNK")
k_.edge("petn", "cdpetn", label="PCYT2 · CTP")
k_.edge("cdpetn", "pe", label="EPT1 · DAG")
k_.edge("pe", "gpe", style="dashed", label="PLA2 / PLB")
k_.edge("pe", "pc", style="dashed", label="PEMT · 3 SAM", side=-1)

# 16. Glycerolipid synthesis
gl = _add(Template("glycerolipid", "Glycerol-3-Phosphate & Glycerolipids", "#FFF9DB", "#F2D27A", "#A0711C"))
gl.stub("gly", "← Glycolysis", 110, 40)
gl.met("dhap", "DHAP", "C00111", 110, 140)
gl.met("g3p", "Glycerol-3P", "C00093", 380, 140)
gl.met("glyc", "Glycerol", "C00116", 380, 270)
gl.met("lpa", "Lysophosphatidate", "C00681", 650, 140)
gl.met("pa", "Phosphatidate", "C00416", 900, 140)
gl.met("dag", "Diacylglycerol", "C00641", 900, 270)
gl.met("tag", "Triacylglycerol", "C00422", 900, 400)
gl.stub("pl", "→ PC / PE (Kennedy)", 650, 270)
gl.edge("gly", "dhap", style="dashed", color="#D6336C")
gl.edge("dhap", "g3p", bi=True, label="GPD1 · NADH")
gl.edge("glyc", "g3p", label="GK · ATP")
gl.edge("g3p", "lpa", label="GPAT · acyl-CoA")
gl.edge("lpa", "pa", label="AGPAT")
gl.edge("pa", "dag", label="LIPIN")
gl.edge("dag", "tag", label="DGAT")
gl.edge("dag", "pl", style="dashed")

# 17. Carnitine shuttle & acylcarnitines
cn = _add(Template("carnitine", "Carnitine Shuttle & Acylcarnitines", "#FFF4E6", "#FFD8A8", "#E8590C",
                   ("Fatty acid metabolism",)))
cn.box(560, 60, 300, 420, "Mitochondrion")
cn.box(930, 60, 260, 520, "Acylcarnitines")
cn.met("pal", "Palmitate", "C00249", 110, 130)
cn.met("palcoa", "Palmitoyl-CoA", "C00154", 110, 280)
cn.met("car", "Carnitine", "C00318", 380, 130)
cn.met("palcar", "Palmitoylcarnitine", "C02990", 380, 280)
cn.met("palcoa_m", "Palmitoyl-CoA", "C00154", 700, 280)
cn.met("acoa", "Acetyl-CoA", "C00024", 700, 410)
cn.met("accar", "Acetylcarnitine", "C02571", 380, 410)
cn.met("c3", "Propionylcarnitine", "C03017", 1060, 130)
cn.met("c4", "Butyrylcarnitine", "C02862", 1060, 240)
cn.met("c8", "Octanoylcarnitine", "C02838", 1060, 350)
cn.met("c10", "Decanoylcarnitine", "C03299", 1060, 460)
cn.stub("tca", "→ TCA", 700, 130)
cn.edge("pal", "palcoa", label="ACSL · ATP")
cn.edge("palcoa", "palcar", label="CPT1")
cn.edge("car", "palcar", style="dashed")
cn.edge("palcar", "palcoa_m", label="CACT · CPT2")
cn.edge("palcoa_m", "acoa", label="β-oxidation", style="dashed")
cn.edge("acoa", "accar", bi=True, label="CrAT")
cn.edge("acoa", "tca", style="dashed", color="#1864AB", via=[(820, 410), (820, 130)])

# 18. Ketone bodies
kb = _add(Template("ketone", "Ketone Body Metabolism", "#EDF2FF", "#BAC8FF", "#3B5BDB"))
kb.met("acoa", "Acetyl-CoA", "C00024", 110, 130)
kb.met("acacoa", "Acetoacetyl-CoA", "C00332", 380, 130)
kb.met("hmg", "HMG-CoA", "C00356", 660, 130)
kb.met("acac", "Acetoacetate", "C00164", 660, 280)
kb.met("bhb", "3-Hydroxybutyrate", "C01089", 380, 280)
kb.met("acet", "Acetone", "C00207", 660, 410)
kb.met("sca", "Succinyl-CoA", "C00091", 110, 280)
kb.stub("leu", "← Leucine", 920, 130)
kb.edge("acoa", "acacoa", bi=True, label="ACAT1")
kb.edge("acacoa", "hmg", label="HMGCS2")
kb.edge("hmg", "acac", label="HMGCL")
kb.edge("acac", "bhb", bi=True, label="BDH1 · NADH")
kb.edge("acac", "acet", style="dashed", label="spontaneous")
kb.edge("sca", "acacoa", style="dashed", label="OXCT1 (ketolysis)", side=-1)
kb.edge("leu", "hmg", style="dashed", color="#2B8A3E")

# 19. Branched-chain amino acid catabolism
b_ = _add(Template("bcaa", "Branched-Chain Amino Acid Catabolism", "#EBFBEE", "#B2F2BB", "#2B8A3E"))
b_.met("leu", "Leucine", "C00123", 110, 120)
b_.met("kic", "4-Methyl-2-oxopentanoate", "C00233", 370, 120)
b_.met("ivcoa", "Isovaleryl-CoA", "C02939", 630, 120)
b_.met("hmg", "HMG-CoA", "C00356", 880, 120)
b_.met("acac", "Acetoacetate", "C00164", 1130, 120)
b_.met("ile", "Isoleucine", "C00407", 110, 250)
b_.met("kmv", "3-Methyl-2-oxopentanoate", "C00671", 370, 250)
b_.met("mbcoa", "2-Methylbutanoyl-CoA", "C01033", 630, 250)
b_.met("acoa", "Acetyl-CoA", "C00024", 1130, 250)
b_.met("val", "Valine", "C00183", 110, 380)
b_.met("kiv", "3-Methyl-2-oxobutanoate", "C00141", 370, 380)
b_.met("ibcoa", "Isobutyryl-CoA", "C00630", 630, 380)
b_.met("ppcoa", "Propionyl-CoA", "C00100", 880, 380)
b_.met("mma", "Methylmalonate", "C02170", 630, 500)
b_.met("mmcoa", "Methylmalonyl-CoA", "C00683", 880, 500)
b_.met("sca", "Succinyl-CoA", "C00091", 1130, 500)
b_.stub("tca", "→ TCA", 1130, 380)
for a, b, lab in [("leu", "kic", "BCAT"), ("ile", "kmv", "BCAT"), ("val", "kiv", "BCAT")]:
    b_.edge(a, b, bi=True, label=lab)
for a, b in [("kic", "ivcoa"), ("kmv", "mbcoa"), ("kiv", "ibcoa")]:
    b_.edge(a, b, label="BCKDH")
b_.edge("ivcoa", "hmg", style="dashed", label="MCC · AUH")
b_.edge("hmg", "acac", label="HMGCL")
b_.edge("hmg", "acoa")
b_.edge("mbcoa", "acoa", style="dashed", via=[(880, 250)])
b_.edge("mbcoa", "ppcoa", style="dashed")
b_.edge("ibcoa", "ppcoa", style="dashed", label="via 3-HIB")
b_.edge("ppcoa", "mmcoa", label="PCC · biotin")
b_.edge("mmcoa", "sca", label="MUT · B12")
b_.edge("mmcoa", "mma", style="dashed")
b_.edge("acoa", "tca", style="dashed", color="#1864AB")
b_.edge("sca", "tca", style="dashed", color="#1864AB")

# 20. GABA shunt & N-acetylaspartate
gb = _add(Template("gaba", "GABA Shunt & N-Acetylaspartate", "#FDEDF1", "#F4BCC9", "#C92A2A"))
gb.met("gln", "Glutamine", "C00064", 110, 130)
gb.met("glu", "Glutamate", "C00025", 380, 130)
gb.met("akg", "α-Ketoglutarate", "C00026", 650, 130)
gb.met("suc", "Succinate", "C00042", 920, 260)
gb.met("gaba", "GABA", "C00334", 380, 260)
gb.met("ssa", "Succinate semialdehyde", "C00232", 650, 260)
gb.met("ghb", "4-Hydroxybutyrate", "C00989", 920, 390)
gb.met("pglu", "Pyroglutamate", "C01879", 110, 260)
gb.met("asp", "Aspartate", "C00049", 110, 390)
gb.met("naa", "N-Acetylaspartate", "C01042", 380, 390)
gb.met("naag", "N-Acetylaspartylglutamate", "C12270", 650, 390)
gb.edge("gln", "glu", label="GLS")
gb.edge("glu", "akg", bi=True, label="GLUD / transaminases")
gb.edge("akg", "suc", style="dashed", color="#1864AB", label="TCA", via=[(920, 130)])
gb.edge("glu", "gaba", label="GAD1/2")
gb.edge("gaba", "ssa", label="ABAT")
gb.edge("ssa", "suc", label="ALDH5A1")
gb.edge("ssa", "ghb", style="dashed", label="AKR7A2", via=[(650, 330), (920, 330)])
gb.edge("glu", "pglu", bi=True, style="dashed")
gb.edge("asp", "naa", label="NAT8L · Acetyl-CoA")
gb.edge("naa", "naag", label="+ Glu")

# 21. NAD+ synthesis & kynurenine
n_ = _add(Template("nad", "NAD⁺ Synthesis & Kynurenine Pathway", "#E7F5FF", "#A5D8FF", "#1971C2",
                   ("Cofactors & vitamins",)))
n_.met("trp", "Tryptophan", "C00078", 110, 120)
n_.met("nfk", "N-Formylkynurenine", "C02700", 360, 120)
n_.met("kyn", "Kynurenine", "C00328", 610, 120)
n_.met("hk", "3-Hydroxykynurenine", "C03227", 860, 120)
n_.met("kyna", "Kynurenic acid", "C01717", 610, 240)
n_.met("haa", "3-Hydroxyanthranilate", "C00632", 860, 240)
n_.met("quin", "Quinolinate", "C03722", 860, 360)
n_.met("na", "Nicotinate", "C00253", 360, 240)
n_.met("namn", "NaMN", "C01185", 610, 360)
n_.met("naad", "NaAD", "C00857", 360, 360)
n_.met("nad", "NAD⁺", "C00003", 110, 360)
n_.met("nadp", "NADP⁺", "C00006", 110, 240)
n_.met("nmn", "NMN", "C00455", 360, 480)
n_.met("nam", "Nicotinamide", "C00153", 610, 480)
n_.met("mnam", "1-Methylnicotinamide", "C02918", 860, 480)
n_.met("nr", "Nicotinamide riboside", "C03150", 360, 590)
n_.edge("trp", "nfk", label="IDO1 / TDO2")
n_.edge("nfk", "kyn")
n_.edge("kyn", "hk", label="KMO")
n_.edge("kyn", "kyna", style="dashed", label="KAT")
n_.edge("hk", "haa", label="KYNU")
n_.edge("haa", "quin", label="HAAO")
n_.edge("quin", "namn", label="QPRT")
n_.edge("na", "namn", label="NAPRT", via=[(360, 300), (610, 300)])
n_.edge("namn", "naad", label="NMNAT")
n_.edge("naad", "nad", label="NADSYN · Gln")
n_.edge("nad", "nadp", label="NADK")
n_.edge("nam", "nmn", label="NAMPT")
n_.edge("nmn", "nad", label="NMNAT", via=[(110, 480)])
n_.edge("nr", "nmn", label="NRK")
n_.edge("nam", "mnam", label="NNMT · SAM")

# 22. Purine salvage & degradation
ps = _add(Template("purine_salvage", "Purine Salvage & Degradation", "#EBFAF3", "#A6E3C8", "#0B7285",
                   ("Purine metabolism",)))
ps.met("atp", "ATP", "C00002", 110, 120)
ps.met("adp", "ADP", "C00008", 330, 120)
ps.met("amp", "AMP", "C00020", 550, 120)
ps.met("ado", "Adenosine", "C00212", 770, 120)
ps.met("ade", "Adenine", "C00147", 990, 120)
ps.met("camp", "cAMP", "C00575", 110, 250)
ps.met("imp", "IMP", "C00130", 550, 250)
ps.met("ino", "Inosine", "C00294", 770, 250)
ps.met("hyp", "Hypoxanthine", "C00262", 990, 250)
ps.met("xan", "Xanthine", "C00385", 1210, 250)
ps.met("urate", "Urate", "C00366", 1210, 120)
ps.met("gmp", "GMP", "C00144", 550, 380)
ps.met("guo", "Guanosine", "C00387", 770, 380)
ps.met("gua", "Guanine", "C00242", 990, 380)
ps.stub("denovo", "← de novo synthesis", 330, 250)
ps.edge("atp", "adp", bi=True)
ps.edge("adp", "amp", bi=True)
ps.edge("atp", "camp", label="AC")
ps.edge("amp", "ado", label="NT5E")
ps.edge("amp", "imp", label="AMPD")
ps.edge("ado", "ino", label="ADA")
ps.edge("imp", "ino", label="NT5C2")
ps.edge("ino", "hyp", bi=True, label="PNP")
ps.edge("hyp", "xan", label="XDH")
ps.edge("xan", "urate", label="XDH")
ps.edge("imp", "gmp", style="dashed", label="IMPDH · GMPS")
ps.edge("gmp", "guo")
ps.edge("guo", "gua", bi=True, label="PNP")
ps.edge("gua", "xan", label="GDA")
ps.edge("denovo", "imp", style="dashed")
ps.edge("ade", "amp", style="dashed", label="APRT (salvage)", via=[(990, 45), (550, 45)], color="#2B8A3E")
ps.edge("hyp", "imp", style="dashed", label="HPRT (salvage)", via=[(990, 315), (550, 315)], color="#2B8A3E")
ps.edge("gua", "gmp", style="dashed", label="HPRT (salvage)", via=[(990, 450), (550, 450)], color="#2B8A3E")

# 23. Pyrimidine salvage & degradation
py = _add(Template("pyrimidine_salvage", "Pyrimidine Salvage & Degradation", "#F1EEFC", "#C5BCEF", "#6741D9",
                   ("Pyrimidine metabolism",)))
py.met("utp", "UTP", "C00075", 110, 120)
py.met("udp", "UDP", "C00015", 330, 120)
py.met("ump", "UMP", "C00105", 550, 120)
py.met("urd", "Uridine", "C00299", 770, 120)
py.met("ura", "Uracil", "C00106", 990, 120)
py.met("dhu", "Dihydrouracil", "C00429", 1210, 120)
py.met("ncba", "N-Carbamoyl-β-alanine", "C02642", 1210, 250)
py.met("bala", "β-Alanine", "C00099", 1210, 380)
py.met("ctp", "CTP", "C00063", 110, 250)
py.met("cmp", "CMP", "C00055", 550, 250)
py.met("cyd", "Cytidine", "C00475", 770, 250)
py.met("dump", "dUMP", "C00365", 330, 380)
py.met("dtmp", "dTMP", "C00364", 550, 380)
py.met("thd", "Thymidine", "C00214", 770, 380)
py.met("thy", "Thymine", "C00178", 990, 380)
py.stub("bai", "→ β-Aminoisobutyrate", 990, 490)
py.edge("utp", "udp", bi=True)
py.edge("udp", "ump", bi=True)
py.edge("utp", "ctp", label="CTPS · Gln")
py.edge("ump", "urd", label="NT5C")
py.edge("urd", "ura", bi=True, label="UPP")
py.edge("ura", "dhu", label="DPYD · NADPH")
py.edge("dhu", "ncba", label="DPYS")
py.edge("ncba", "bala", label="UPB1")
py.edge("cmp", "cyd")
py.edge("cyd", "urd", label="CDA")
py.edge("dump", "dtmp", label="TYMS")
py.edge("dtmp", "thd")
py.edge("thd", "thy", bi=True, label="TYMP")
py.edge("thy", "bai", style="dashed", label="DPYD")

# 24. Phenylalanine / tyrosine
ar = _add(Template("aromatic", "Phenylalanine & Tyrosine Metabolism", "#FFF9DB", "#FFE066", "#E67700"))
ar.met("phe", "Phenylalanine", "C00079", 110, 130)
ar.met("tyr", "Tyrosine", "C00082", 380, 130)
ar.met("hpp", "4-Hydroxyphenylpyruvate", "C01179", 650, 130)
ar.met("hga", "Homogentisate", "C00544", 920, 130)
ar.met("fum", "Fumarate", "C00122", 920, 270)
ar.met("acac", "Acetoacetate", "C00164", 920, 400)
ar.met("ppy", "Phenylpyruvate", "C00166", 110, 270)
ar.met("plac", "Phenyllactate", "C05607", 110, 400)
ar.met("dopa", "L-DOPA", "C00355", 380, 270)
ar.met("da", "Dopamine", "C03758", 380, 400)
ar.edge("phe", "tyr", label="PAH · BH4")
ar.edge("tyr", "hpp", bi=True, label="TAT")
ar.edge("hpp", "hga", label="HPD")
ar.edge("hga", "fum", label="HGD · FAH")
ar.edge("hga", "acac", style="dashed", via=[(1040, 130), (1040, 400)])
ar.edge("phe", "ppy", style="dashed", label="transaminase")
ar.edge("ppy", "plac")
ar.edge("tyr", "dopa", label="TH")
ar.edge("dopa", "da", label="DDC")

# 25. Histidine & β-alanine
hs = _add(Template("histidine", "Histidine Metabolism", "#F8F0FC", "#E1C2F0", "#862E9C"))
hs.met("his", "Histidine", "C00135", 110, 130)
hs.met("uro", "Urocanate", "C00785", 380, 130)
hs.met("ipr", "Imidazolone-propanoate", "C03680", 650, 130)
hs.met("figlu", "N-Formiminoglutamate", "C00439", 920, 130)
hs.met("glu", "Glutamate", "C00025", 920, 270)
hs.met("hist", "Histamine", "C00388", 110, 270)
hs.met("carn", "Carnosine", "C00386", 380, 270)
hs.met("bala", "β-Alanine", "C00099", 380, 400)
hs.stub("thf", "→ 5,10-CH=THF (one-carbon)", 650, 270)
hs.edge("his", "uro", label="HAL")
hs.edge("uro", "ipr", label="UROC1")
hs.edge("ipr", "figlu", label="AMDHD1")
hs.edge("figlu", "glu", label="FTCD · THF")
hs.edge("figlu", "thf", style="dashed")
hs.edge("his", "hist", label="HDC")
hs.edge("his", "carn", style="dashed", label="CARNS1")
hs.edge("bala", "carn")

ORDER = list(T)


def options() -> list:
    """(key, title) in display order."""
    return [(k, T[k].title) for k in ORDER]
