"""
core.interpret — rule-based interpretation summary.

The text is assembled from the statistical results with controlled wording. It keeps five kinds of statements
apart and never claims that a pathway is activated / inhibited or that flux increased because a fraction rose:
  1 observed peak-intensity (pool size) changes
  2 observed isotope-incorporation changes
  3 isotopologue-distribution differences (with standard route interpretations, stated as "consistent with")
  4 model-derived APPARENT labeling kinetics (only from supported fits)
  5 hypotheses (explicitly labeled)
plus experimental limitations.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

FORBIDDEN = [r"\bactivat(ed|ion|es)\b", r"\binhibit(ed|ion|s)\b", r"\bincreased flux\b", r"\bdecreased flux\b",
             r"\bhigher flux\b", r"\blower flux\b", r"\bupregulat\w*\b", r"\bdownregulat\w*\b", r"\bflux (is|was) (increased|elevated|reduced)\b"]

# standard route meaning of isotopologues for [U-13C]-glucose tracing (library id -> {i: meaning})
ROUTES = {
    "cit": {2: "acetyl-CoA entry (PDH, first TCA turn)", 3: "pyruvate-carboxylase (anaplerotic) entry",
            4: "label retained over more than one TCA turn", 5: "PC + PDH combined / multiple turns"},
    "akg": {2: "first-turn PDH-derived carbon", 3: "pyruvate-carboxylase-derived carbon", 4: "multiple TCA turns"},
    "glu": {2: "first-turn PDH-derived carbon (via α-KG)", 3: "pyruvate-carboxylase-derived carbon",
            4: "multiple TCA turns", 5: "fully labeled α-KG backbone"},
    "gln": {2: "glutamate M+2 (PDH-derived) carried into glutamine", 3: "PC-derived carbon"},
    "mal": {2: "oxidative TCA (first turn)", 3: "pyruvate carboxylase / reductive entry", 4: "multiple turns"},
    "asp": {2: "oxidative TCA (first turn)", 3: "pyruvate carboxylase / reductive entry", 4: "multiple turns"},
    "fum": {2: "oxidative TCA (first turn)", 3: "pyruvate carboxylase entry"},
    "suc": {2: "oxidative TCA (first turn)"},
    "lac": {3: "glycolytic pyruvate"}, "pyr": {3: "glycolytic pyruvate"}, "pep": {3: "glycolysis"},
    "ser": {3: "de novo serine synthesis from glycolytic 3-phosphoglycerate"}, "gly": {2: "serine-derived glycine"},
    "ump": {5: "ribose from glucose (PPP)", 3: "base labeled via aspartate", 8: "ribose + labeled base"},
    "cmp": {5: "ribose from glucose (PPP)", 3: "base labeled via aspartate"},
    "udp": {5: "ribose from glucose (PPP)", 3: "base labeled via aspartate"},
    "amp": {5: "ribose from glucose (PPP)"}, "adp": {5: "ribose from glucose (PPP)"}, "imp": {5: "ribose from glucose (PPP)"},
    "atp": {5: "ribose from glucose (PPP)"}, "nad": {5: "one ribose labeled", 10: "both riboses labeled"},
    "udpglcnac": {2: "acetyl group", 5: "ribose (UTP)", 6: "glucosamine (hexosamine pathway)",
                  11: "glucosamine + ribose", 13: "glucosamine + ribose + acetyl"},
    "udpgal": {6: "hexose moiety", 5: "ribose (UTP)", 11: "hexose + ribose"},
    "gsh": {2: "glutamate- or glycine-derived carbon"}, "gssg": {2: "glutamate- or glycine-derived carbon"},
    "c2car": {2: "glucose-derived acetyl group"},
}


def lint(text: str) -> list:
    """Return forbidden causal phrases present in a (possibly user-edited) interpretation."""
    hits = []
    for rx in FORBIDDEN:
        for m in re.finditer(rx, text, re.I):
            hits.append(m.group(0))
    return sorted(set(hits))


def _fmt_p(q):
    return "FDR < 0.001" if q < 0.001 else f"FDR = {q:.3f}"


def _pw_tag(A, comp):
    if A is None or comp not in set(A.Compound):
        return ""
    p = A.set_index("Compound").loc[comp, "Pathways"]
    return f" [{p}]" if isinstance(p, str) and p else ""


def _label(A, comp):
    if A is None or comp not in set(A.Compound):
        return comp
    r = A.set_index("Compound").loc[comp]
    ok = r.Confidence in ("High", "Medium", "User-confirmed")
    l = r.Label
    return l if ok and isinstance(l, str) and l and l.lower() != comp.lower() else comp


def build_summary(*, group_a: str, group_b: str, intensity: pd.DataFrame | None, incorporation: pd.DataFrame | None,
                  channels: pd.DataFrame | None, fits: pd.DataFrame | None, A: pd.DataFrame | None,
                  design_info: dict, fdr: float = 0.05, min_lfc: float = 0.58, min_diff: float = 2.0,
                  top_n: int = 8) -> str:
    L = []
    w = L.append
    w(f"# Interpretation summary — {group_b} vs {group_a}")
    w("")
    w(f"*Direction: positive values mean higher in **{group_b}** than in **{group_a}**. "
      f"Significance threshold: FDR < {fdr}. Effects are reported per time point.*")
    w("")
    w("> This summary is generated automatically from the statistics. It separates what was **observed** "
      "(pool sizes, isotope incorporation, isotopologue patterns) from what is **model-derived** (apparent labeling "
      "kinetics) and from **hypotheses**. Isotope-labeling data alone do not establish pathway flux.")
    w("")
    # 1 intensity
    w("## 1. Observed changes in peak intensity (metabolite abundance)")
    if intensity is None or intensity.empty:
        w("Not computed.")
    else:
        sig = intensity[(intensity["FDR (BH)"] < fdr) & (intensity.log2FC.abs() >= min_lfc)]
        if sig.empty:
            w(f"No metabolite pool differed at FDR < {fdr} with |log2FC| ≥ {min_lfc:g}.")
        for t, d in sig.groupby("Time (h)"):
            up = d[d.log2FC > 0].sort_values("FDR (BH)").head(top_n)
            dn = d[d.log2FC < 0].sort_values("FDR (BH)").head(top_n)
            w(f"**{t:g} h:** {len(d[d.log2FC > 0])} higher and {len(d[d.log2FC < 0])} lower in {group_b}.")
            for lab, part in (("Higher", up), ("Lower", dn)):
                if len(part):
                    items = "; ".join(f"{_label(A, r.Feature)} (log2FC {r.log2FC:+.2f}, {_fmt_p(r['FDR (BH)'])}){_pw_tag(A, r.Feature)}"
                                      for _, r in part.iterrows())
                    w(f"- {lab}: {items}")
        w("")
        w("*Pool-size changes describe how much of a metabolite is present; they do not by themselves show "
          "whether production or consumption changed.*")
    w("")
    # 2 incorporation
    w("## 2. Observed changes in isotope incorporation")
    if incorporation is None or incorporation.empty:
        w("Not computed.")
    else:
        dcol = "Difference (B−A)"
        sig = incorporation[(incorporation["FDR (BH)"] < fdr) & (incorporation[dcol].abs() >= min_diff)]
        if sig.empty:
            w(f"No metabolite showed a difference in total incorporation at FDR < {fdr} and ≥ {min_diff:g} percentage points.")
        for t, d in sig.groupby("Time (h)"):
            w(f"**{t:g} h:** {int((d[dcol] > 0).sum())} metabolites with more and {int((d[dcol] < 0).sum())} with "
              f"less label in {group_b}.")
            for lab, part in (("More label", d[d[dcol] > 0]), ("Less label", d[d[dcol] < 0])):
                part = part.reindex(part[dcol].abs().sort_values(ascending=False).index).head(top_n)
                if len(part):
                    items = "; ".join(f"{_label(A, r.Feature)} ({r[dcol]:+.1f} points, {_fmt_p(r['FDR (BH)'])}){_pw_tag(A, r.Feature)}"
                                      for _, r in part.iterrows())
                    w(f"- {lab}: {items}")
        w("")
        w("*Incorporation (%) is the share of molecules carrying at least one label; it is not the share of carbon "
          "atoms labeled.*")
    w("")
    # 3 isotopologues
    w("## 3. Differences in isotopologue distribution")
    if channels is None or channels.empty:
        w("Not computed.")
    else:
        dcol = "Difference (B−A)"
        sig = channels[(channels["FDR (BH)"] < fdr) & (channels[dcol].abs() * 100 >= min_diff)].copy()
        if sig.empty:
            w(f"No isotopologue fraction differed at FDR < {fdr} by ≥ {min_diff:g} percentage points.")
        else:
            sig[["Compound", "Channel"]] = sig.Feature.str.split(" \\| ", n=1, expand=True)
            lib = dict(zip(A.Compound, A.LibraryID)) if A is not None else {}
            for t, d in sig.groupby("Time (h)"):
                w(f"**{t:g} h:**")
                for comp, dd in list(d.groupby("Compound"))[:top_n * 2]:
                    parts = []
                    for _, r in dd.sort_values("Channel").iterrows():
                        i = int(r.Channel.replace("M+", ""))
                        meaning = ROUTES.get(lib.get(comp, ""), {}).get(i)
                        parts.append(f"{r.Channel} {r[dcol] * 100:+.1f} points ({_fmt_p(r['FDR (BH)'])})"
                                     + (f" — consistent with {meaning}" if meaning else ""))
                    w(f"- {_label(A, comp)}: " + "; ".join(parts))
            w("")
            w("*Route meanings are the standard interpretation for uniformly labeled glucose and assume natural-"
              "abundance-corrected data. A shifted isotopologue pattern shows that the relative contribution of "
              "labeled routes differs; it does not quantify flux through them.*")
    w("")
    # 4 kinetics
    w("## 4. Model-derived apparent labeling kinetics")
    if fits is None or fits.empty:
        w("Not computed (requires ≥ 3 time points).")
    else:
        ok = fits[fits.status == "Fitted"]
        nf = int((fits.status == "Not fitted").sum())
        w(f"First-order fits: {len(ok)} supported, {int((fits.status == 'Fitted with warnings').sum())} fitted with "
          f"warnings, {nf} not fitted (mostly because labeling did not rise steadily over time).")
        both = ok.pivot_table(index="Compound", columns="Group", values=["t_half", "t_half_lo", "t_half_hi"], aggfunc="first")
        lines = []
        if not both.empty and group_a in both["t_half"].columns and group_b in both["t_half"].columns:
            for c, r in both.iterrows():
                a, b = r[("t_half", group_a)], r[("t_half", group_b)]
                if np.isfinite(a) and np.isfinite(b):
                    sep = r[("t_half_hi", group_b)] < r[("t_half_lo", group_a)] or r[("t_half_hi", group_a)] < r[("t_half_lo", group_b)]
                    if sep:
                        lines.append(f"- {_label(A, c)}: apparent t½ {b:.1f} h ({group_b}) vs {a:.1f} h ({group_a}); 95 % CIs do not overlap.")
        if lines:
            w("Metabolites whose apparent labeling half-time differs between groups (non-overlapping CIs):")
            L.extend(lines)
        else:
            w("No metabolite had supported fits in both groups with non-overlapping half-time confidence intervals.")
        w("")
        w("*Apparent kinetics describe how fast label appears in a pool. They depend on pool size, exchange and "
          "precursor labeling and are not absolute pathway fluxes.*")
    w("")
    # 5 pathway concordance + hypotheses
    w("## 5. Pathway-level pattern")
    hyp = []
    if incorporation is not None and not incorporation.empty and A is not None:
        dcol = "Difference (B−A)"
        sig = incorporation[(incorporation["FDR (BH)"] < fdr) & (incorporation[dcol].abs() >= min_diff)]
        sig = sig.merge(A[["Compound", "Pathways"]], left_on="Feature", right_on="Compound", how="left")
        rows = []
        for _, r in sig.iterrows():
            for p in [x.strip() for x in str(r.Pathways).split(";") if x.strip() and x.strip() != "nan"]:
                rows.append((p, r["Time (h)"], np.sign(r[dcol])))
        if rows:
            P = pd.DataFrame(rows, columns=["Pathway", "Time", "Dir"])
            for p, d in P.groupby("Pathway"):
                up, dn = int((d.Dir > 0).sum()), int((d.Dir < 0).sum())
                w(f"- **{p}:** {up} metabolite×time results with more label and {dn} with less label in {group_b}.")
                if up >= 3 and up >= 3 * max(dn, 1):
                    hyp.append(f"Several {p.lower()} metabolites carry more glucose-derived label in {group_b}. This is "
                               f"consistent with a greater relative contribution of glucose carbon to these pools, which "
                               f"could reflect faster synthesis, slower turnover of pre-existing unlabeled pools, or "
                               f"more highly labeled precursors.")
                elif dn >= 3 and dn >= 3 * max(up, 1):
                    hyp.append(f"Several {p.lower()} metabolites carry less glucose-derived label in {group_b}, consistent "
                               f"with a smaller relative contribution of glucose carbon (e.g. more input from other "
                               f"substrates, slower synthesis, or larger pre-existing pools).")
        else:
            w("No pathway had significant incorporation differences.")
    else:
        w("Pathway assignments or incorporation statistics are not available.")
    w("")
    w("## 6. Hypotheses (to be tested)")
    if hyp:
        for h in hyp:
            w(f"- {h}")
        w("- Distinguishing these explanations would need, for example, flux modeling (e.g. INCA / 13C-MFA), "
          "additional tracers ([1,2-13C2]-glucose, [U-13C5]-glutamine), finer early time points, or enzyme-level assays.")
    else:
        w("- No pathway-level pattern was strong enough to propose a hypothesis automatically.")
    w("")
    w("## 7. Experimental limitations")
    for x in design_info.get("limitations", []):
        w(f"- {x}")
    w("- Metabolite identities are based on accurate mass and retention time (no MS/MS confirmation implied); "
      "isomers with the same formula may not be resolved.")
    w("- Statistics were run per metabolite with Benjamini–Hochberg correction; small groups give low power and "
      "imprecise variance estimates.")
    return "\n".join(L)


def design_limitations(D_an: pd.DataFrame, correction_mode: str, n_times: int, unmapped: int) -> list:
    out = []
    n = D_an.groupby(["Group", "Hours"]).size()
    if len(n) and n.min() < 4:
        out.append(f"Small groups (n = {n.min()}{'' if n.min() == n.max() else '–' + str(n.max())} per group and time point).")
    if n_times < 4:
        out.append(f"Only {n_times} labeling time point(s): early kinetics and plateau are poorly defined.")
    out.append({"already": "Isotopologue areas were used as supplied (natural-abundance correction assumed already done "
                           "by the vendor software).",
                "apply": "Natural-abundance correction was applied in this app (matrix method).",
                "none": "No natural-abundance correction: M+1/M+2 contain natural-isotope signal; mean label number "
                        "is overestimated."}[correction_mode])
    if unmapped:
        out.append(f"{unmapped} metabolites could not be mapped confidently to pathways and are excluded from "
                   "pathway-level statements.")
    return out


def markdown_to_html(md: str) -> str:
    """Small Markdown converter for the report (headings, bold, italics, lists, quotes)."""
    import html
    out, in_list = [], False
    for line in md.splitlines():
        s = html.escape(line)
        s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s)
        s = re.sub(r"\*(.+?)\*", r"<i>\1</i>", s)
        if line.startswith("- "):
            if not in_list:
                out.append("<ul>"); in_list = True
            out.append(f"<li>{s[2:]}</li>"); continue
        if in_list:
            out.append("</ul>"); in_list = False
        if line.startswith("# "):
            out.append(f"<h2>{s[2:]}</h2>")
        elif line.startswith("## "):
            out.append(f"<h3>{s[3:]}</h3>")
        elif line.startswith("&gt; "):
            out.append(f"<blockquote>{s[5:]}</blockquote>")
        elif line.strip():
            out.append(f"<p>{s}</p>")
    if in_list:
        out.append("</ul>")
    return "\n".join(out)
