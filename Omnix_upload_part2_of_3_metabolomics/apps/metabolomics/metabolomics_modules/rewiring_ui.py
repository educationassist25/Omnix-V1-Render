"""
rewiring_ui.py - Streamlit UI for the Correlation Rewiring Map.

Used in two places:
  * rewiring_app.py        standalone app (`streamlit run rewiring_app.py`) with its own
                           data loading (simulated demo, MetaboAI Pro demo, or upload)
  * app.py                 "Network Analysis -> Correlation Rewiring Map" page, which
                           analyses the pipeline's normalized log2 data

All statistics live in metabolomics_modules/rewiring.py; this file only collects settings, runs
the engine, and renders results (an interactive HTML view + native tables/downloads).
"""

from __future__ import annotations

import datetime as _dt
import json
import math
import os
import re

import pandas as pd

from metabolomics_modules import rewiring as rw
from metabolomics_modules import rewiring_figures as rf
from metabolomics_modules import rewiring_reactions as rr
from metabolomics_modules import utils

_ASSET = os.path.join(os.path.dirname(__file__), "rewiring_assets", "rewiring_view.html")
LIBRARY_OPTIONS = {
    "KEGG (as in MSEA)": "KEGG",
    "SMPDB (as in MSEA)": "SMPDB",
    "Human-GEM subsystems": rr.LIB_GEM,
    "LIPID MAPS (as in MSEA)": "LIPID MAPS",
}
CUTOFFS = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.25, 0.5, 1.0]
PW_CUTOFFS = [0.001, 0.01, 0.05, 0.1, 0.2, 0.25, 0.5, 1.0]
SCOPES = {
    "All metabolite pairs": "all",
    "Reaction-linked pairs (≤ 2 steps)": "reaction",
    "Direct reactions only": "direct",
}
DEMO_LABEL = "Simulated rewiring demo (Control vs Tumor, 30 + 30 samples)"


# ---------------------------------------------------------------------------
# HTML view
# ---------------------------------------------------------------------------
def _embed_json(obj) -> str:
    """JSON for a <script type="application/json"> block; '</' and '<!--' are escaped so
    a metabolite name can never close the script element."""
    text = json.dumps(obj, ensure_ascii=False, allow_nan=False)
    return text.replace("</", "<\\/").replace("<!--", "<\\!--")


def render_view_html(result: rw.RewiringResult, *, title: str | None = None, theme: str | None = None,
                     show_table: bool = True, generated: str | None = None, max_cards: int = 6,
                     pw_top: int | None = None, pw_metric: str = "p", pw_cutoff: float = 0.05) -> str:
    payload = result.to_payload(max_cards=max_cards)
    payload["options"] = {"title": title, "theme": theme, "show_table": show_table, "generated": generated,
                          "pw_top": pw_top, "pw_metric": pw_metric, "pw_cutoff": pw_cutoff}
    with open(_ASSET, encoding="utf-8") as fh:
        template = fh.read()
    return template.replace("/*__REWIRING_PAYLOAD__*/", _embed_json(payload))


def _safe_tag(s) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_") or "Data"


def _theme(st) -> str | None:
    try:
        t = st.context.theme.type
        return t if t in ("light", "dark") else None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Demo helpers
# ---------------------------------------------------------------------------
def load_simulated_demo():
    """Simulated demo in MetaboAI Pro's own formats: (peak_df, meta indexed by Sample,
    annotations indexed by Metabolite)."""
    d = rw.simulate_rewiring_demo()
    meta = d["meta"].set_index("Sample")
    ann = d["annotations"].set_index("Metabolite")
    return d["peak_df"], meta, ann


def demo_csv_bytes() -> dict:
    d = rw.simulate_rewiring_demo()
    return {
        "sample_peak_area_matrix_rewiring_demo.csv": d["peak_df"].to_csv().encode("utf-8-sig"),
        "sample_metadata_rewiring_demo.csv": d["meta"].to_csv(index=False).encode("utf-8-sig"),
        "metabolite_row_annotations_rewiring_demo.csv": d["annotations"].to_csv(index=False).encode("utf-8-sig"),
    }


# ---------------------------------------------------------------------------
# Workspace: settings -> run -> results
# ---------------------------------------------------------------------------
def render_workspace(st, data: pd.DataFrame, meta: pd.DataFrame, annotations: pd.DataFrame | None, *,
                     key_prefix: str, source_label: str, data_is_log2: bool,
                     default_median_center: bool = False):
    """
    data        : features x samples (raw peak areas or log2 values)
    meta        : sample metadata indexed by sample name
    annotations : optional features x annotation-columns table indexed by metabolite name
    data_is_log2: True when data already went through MetaboAI Pro normalization + log2
    """
    k = lambda name: f"{key_prefix}_{name}"   # noqa: E731
    samples = [c for c in data.columns if c in meta.index]
    is_qc = meta["IsQC"] if "IsQC" in meta.columns else pd.Series(False, index=meta.index)
    bio = [s for s in samples if not bool(is_qc.get(s, False))]
    cat_cols = [c for c in utils.get_categorical_metadata_columns(meta, bio) if c != "Batch"]
    if not cat_cols:
        st.error("The metadata has no categorical column with at least two groups to compare.")
        return

    st.markdown("##### Comparison")
    c1, c2, c3 = st.columns(3)
    group_col = c1.selectbox("Grouping variable", cat_cols,
                             index=cat_cols.index("Group") if "Group" in cat_cols else 0, key=k("group_col"))
    levels = [str(v) for v in pd.unique(meta.loc[bio, group_col].dropna().astype(str))]
    if len(levels) < 2:
        st.warning(f"'{group_col}' has fewer than two groups among biological samples.")
        return
    a_default = next((i for i, v in enumerate(levels) if re.search(r"control|healthy|normal|untreated|baseline", v, re.I)), 0)
    group_a = c2.selectbox("Reference group (A)", levels, index=a_default, key=k("group_a"))
    b_choices = [v for v in levels if v != group_a]
    group_b = c3.selectbox("Comparison group (B)", b_choices, index=0, key=k("group_b"))
    counts = meta.loc[bio, group_col].astype(str).value_counts()
    n_a, n_b = int(counts.get(group_a, 0)), int(counts.get(group_b, 0))
    st.caption(f"{n_a} samples in **{group_a}**, {n_b} in **{group_b}** (QC samples excluded).")
    n_min = min(n_a, n_b)
    if n_min < 3:
        st.error(f"'{group_a if n_a < 3 else group_b}' has {n_min} biological sample(s); a correlation needs at least 3.")
    elif n_min <= 4:
        splits = math.comb(n_a + n_b, n_a)
        floor = rw.exact_p_floor(n_a, n_b)
        st.caption(
            f"ℹ️ Small groups: all {splits} possible label splits are enumerated (exact permutation test), so the "
            f"smallest achievable p-value is {floor:.3g}"
            + (" and nothing can reach FDR q-value < 0.05. Results are a ranked, exploratory readout (useful for choosing "
               "what to validate and for planning a larger study)." if floor >= 0.05 else
               ". Only very strong, pathway-wide rewiring can reach significance; treat results as exploratory.")
        )
    elif n_min < 15:
        st.caption(f"ℹ️ {n_min} samples in the smaller group: single pairs need large changes (roughly |Δr| > 0.8) to "
                   "reach significance; the pathway-level test pools pairs and has more power.")

    c4, c6 = st.columns(2)
    lib_opts = list(LIBRARY_OPTIONS)
    pw_col = c4.selectbox("Pathway library (ring grouping + pathway test)", lib_opts, index=0, key=k("pw_col"),
                          help="KEGG and SMPDB are the same libraries and ID standardization as the Pathway "
                               "Analysis (MSEA) page. Human-GEM subsystems come from the reaction library.")
    batch_col = "None"          # no batch stratification on this page
    method = c6.radio("Correlation", ["Pearson", "Spearman"], horizontal=True, key=k("method"))

    c7, c8, c9 = st.columns(3)
    scope_label = c7.radio("Pairs to test", list(SCOPES), key=k("scope"),
                           help="Reaction-linked scopes test only metabolite pairs connected by an enzymatic reaction "
                                "in the Human-GEM library (directly, or through one unmeasured intermediate). Far "
                                "fewer tests than all pairs, so more power, and every hit has a mechanism.")
    cutoff_on = "fdr"
    fdr = c8.select_slider("Rewired-pair FDR q-value cutoff (≤)", options=CUTOFFS, value=0.05, key=k("fdr"),
                           help="False discovery rate (permutation FDR or Benjamini–Hochberg, see Advanced settings) "
                                "across all tested pairs.")
    min_delta = c9.slider("Minimum |Δr|", 0.0, 1.2, 0.45, 0.05, key=k("dmin"))
    if fdr > 0.25:
        c8.caption("⚠️ A high FDR q-value cutoff lets through many false positives. Use for exploration.")

    with st.expander("Advanced settings"):
        a1, a2, a3 = st.columns(3)
        max_features = a1.number_input("Max features analysed (most variable kept)", 20, 1000, 300, 10, key=k("maxf"))
        max_missing = a2.slider("Drop features missing in more than this share of a group", 0.0, 0.8, 0.3, 0.05, key=k("maxmiss"))
        n_perm = a3.number_input("Label permutations (0 = skip)", 0, 5000, 500, 50, key=k("nperm"),
                                 help="Used for the permutation FDR, the pathway-level test and the global test. When "
                                      "the groups are small enough that there are fewer distinct label splits than "
                                      "this, every split is enumerated (exact test).")
        a4, a5, a6 = st.columns(3)
        n_boot = a4.number_input("Bootstrap resamples per module", 20, 1000, 200, 20, key=k("nboot"))
        if data_is_log2:
            log_choice = "no"
            a5.caption("Input is already log2 (from Normalization), so no transform is applied.")
        else:
            log_choice = {"Auto-detect": "auto", "Yes": "yes", "No (already log scale)": "no"}[
                a5.selectbox("Log2 transform", ["Auto-detect", "Yes", "No (already log scale)"], key=k("log"))]
        median_center = a6.checkbox("Median-centre each sample (removes dilution/loading effects)",
                                    value=default_median_center, key=k("medc"))
        seed = a4.number_input("Random seed", 0, 10_000, 0, 1, key=k("seed"))
        fdr_label = a5.radio("FDR method", ["Permutation FDR (recommended)", "Benjamini–Hochberg (analytic)"],
                             key=k("fdrm"),
                             help="The analytic Fisher-z p-values are too optimistic in the far tail when groups are "
                                  "small, so BH can over-call pairs. The permutation FDR uses the data's own null.")
        fdr_method = "permutation" if fdr_label.startswith("Permutation") else "bh"
        if fdr_method == "permutation" and int(n_perm) == 0:
            st.caption("Permutations are set to 0, so Benjamini–Hochberg will be used and the pathway-level test "
                       "is skipped.")

    run = st.button("Run Rewiring Analysis", type="primary", icon="🕸️", key=k("run"))
    state_key = k("result")
    scope = SCOPES[scope_label]
    settings = dict(group_col=group_col, group_a=group_a, group_b=group_b, pw_col=pw_col, batch_col=batch_col, scope=scope,
                    cutoff_on=cutoff_on,
                    method=method.lower(), fdr=float(fdr), min_delta=float(min_delta), max_features=int(max_features),
                    max_missing=float(max_missing), n_perm=int(n_perm), n_boot=int(n_boot), log=log_choice,
                    median_center=bool(median_center), seed=int(seed), fdr_method=fdr_method, source=source_label,
                    shape=tuple(data.shape))
    if run:
        library, ann_col = LIBRARY_OPTIONS[pw_col], None
        try:
            with st.spinner("Mapping metabolites to pathways and reactions, computing per-group correlations, "
                            "permutations and robustness checks…"):
                prep = rw.prepare_data(
                    data, meta, group_col, group_a, group_b,
                    pathway_library=library, row_annotations=annotations, annotation_column=ann_col,
                    batch_col=None if batch_col == "None" else batch_col,
                    log_transform=log_choice, max_features=int(max_features), max_missing_frac=float(max_missing),
                    sample_median_center=bool(median_center), seed=int(seed),
                )
                cfg = rw.RewiringConfig(method=method.lower(), fdr=float(fdr), cutoff_on=cutoff_on,
                                        min_delta=float(min_delta),
                                        fdr_method=fdr_method, pair_scope=scope,
                                        n_permutations=int(n_perm), n_bootstrap=int(n_boot), seed=int(seed))
                result = rw.run_rewiring(prep, cfg, run_permutation=int(n_perm) > 0)
            st.session_state[state_key] = {"result": result, "settings": settings}
        except rw.RewiringInputError as err:
            st.session_state[state_key] = None
            st.error(str(err))
            return

    stored = st.session_state.get(state_key)
    if not stored:
        st.info("Choose the two groups to compare and click **Run Rewiring Analysis**.")
        return
    if stored["settings"] != settings:
        st.caption("⚠️ Settings changed since the last run. Click **Run Rewiring Analysis** to update the results below.")
    render_results(st, stored["result"], key_prefix=key_prefix, source_label=stored["settings"]["source"])


def show_html(st, html: str, fallback_height: int = 1720):
    """Embed the self-contained view. st.iframe (Streamlit >= 1.50) sizes itself to the
    content; older versions fall back to components.html with a fixed height."""
    if hasattr(st, "iframe"):
        st.iframe(html, height="content")
    else:
        import streamlit.components.v1 as components
        components.html(html, height=fallback_height, scrolling=True)


COLUMN_LABELS = {
    "p_value": "p-value (Fisher z)", "q_bh": "FDR q-value (BH)", "p_perm": "p-value (permutation)",
    "q_perm": "FDR q-value (permutation)", "q_value": "FDR q-value (used)", "p_used": "p-value (used)",
    "ratio_p": "ratio p-value (Welch t)", "ratio_q": "ratio FDR q-value", "pair_q": "pair FDR q-value",
    "delta_r": "Δr", "z": "z", "tested": "tested",
}


PATHWAY_LABELS = {
    "pathway": "Pathway", "n_measured": "Measured metabolites", "n_pairs": "Pairs", "members": "Members",
    "statistic": "Statistic (mean z²)", "mean_abs_delta_r": "Mean |Δr|", "mean_abs_r_a": "Mean |r| A",
    "mean_abs_r_b": "Mean |r| B", "direction": "Direction", "n_lost": "Lost pairs",
    "n_gained": "Gained pairs", "n_flipped": "Flipped pairs", "n_reaction_pairs": "Reaction-linked pairs",
    "top_pairs": "Top pairs", "p_value": "p-value (permutation)", "q_value": "FDR q-value (BH)",
}


def labelled(df: pd.DataFrame, kind: str = "pairs", groups: tuple | None = None) -> pd.DataFrame:
    """Column names written out for display and CSV downloads on this page (p-value, FDR q-value, ...)."""
    if df is None or len(df) == 0:
        return df
    if kind == "pathways":
        out = df.rename(columns=PATHWAY_LABELS)
        if groups:
            out = out.rename(columns={"Mean |r| A": f"Mean |r| {groups[0]}", "Mean |r| B": f"Mean |r| {groups[1]}"})
        return out
    if kind == "features":
        return df.rename(columns={"p_value": "p-value (Welch t)", "q_value": "FDR q-value (BH)"})
    return df.rename(columns=COLUMN_LABELS)


def render_results(st, result: rw.RewiringResult, *, key_prefix: str, source_label: str):
    ga, gb = result.prep.group_labels
    comparison = f"{_safe_tag(ga)}_vs_{_safe_tag(gb)}"
    c = result.counts()
    S = result.summary()
    n_rew = c["lost"] + c["gained"] + c["flipped"]
    perm = result.permutation
    scope_txt = {"all": "all pairs", "reaction": "reaction-linked", "direct": "direct reactions"}[S["pair_scope"]]
    msg = (f"{n_rew} rewired pairs ({c['lost']} lost, {c['gained']} gained, {c['flipped']} flipped) among "
           f"{S['n_tested_pairs']} tested pairs ({scope_txt}) of {len(result.prep.features)} metabolites; "
           f"{S['reaction_links']['rewired_with_reaction']} of the rewired pairs are linked by a reaction.")
    if S["n_pathways_tested"]:
        msg += (f" Pathway level: {S['n_pathways_significant']} of {S['n_pathways_tested']} "
                f"{S['pathway_library']} pathways rewired (FDR q-value < 0.05).")
    if perm:
        kind = "exact, all label splits" if perm.get("exact") else "label permutations"
        msg += f" Global test p-value = {perm['p_value']:.3g} ({perm['n_permutations']} {kind})."
    if result.cfg.cutoff_on == "p":
        msg += (f" Pairs called on {'permutation ' if result.fdr_used == 'permutation' else ''}p-value ≤ {result.cfg.fdr:g} "
                "(no multiple-testing correction).")
    else:
        msg += (f" Pairs called on the FDR q-value ({'permutation FDR' if result.fdr_used == 'permutation' else 'Benjamini–Hochberg'})"
                f" ≤ {result.cfg.fdr:g}.")
    (st.success if (n_rew or S["n_pathways_significant"]) else st.info)(msg)

    # network display: which pathways the ring shows (also used for the figure downloads)
    has_pw = result.pathway_tests is not None and not result.pathway_tests.empty
    if has_pw:
        n1, n2, n3 = st.columns(3)
        top_label = n1.selectbox("Network: pathways shown", list(rf.TOP_CHOICES), index=1, key=f"{key_prefix}_net_top")
        metric_label = n2.radio("Rank pathways by", ["p-value", "FDR q-value"], horizontal=True,
                                key=f"{key_prefix}_net_metric")
        pw_cut = n3.select_slider("Pathway cutoff (≤)", options=PW_CUTOFFS, value=0.05, key=f"{key_prefix}_net_cut")
        pw_top, pw_metric = rf.TOP_CHOICES[top_label], ("q" if metric_label == "FDR q-value" else "p")
    else:
        pw_top, pw_metric, pw_cut = 0, "q", 1.0
        st.caption("The pathway-level test was not run (no permutations), so the network shows all metabolites.")

    html = render_view_html(result, theme=_theme(st), pw_top=pw_top, pw_metric=pw_metric, pw_cutoff=pw_cut)
    show_html(st, html)

    st.markdown("##### Downloads")
    d1, d2, d3, d4 = st.columns(4)
    d1.download_button("Rewired pairs (CSV)", utils.to_download_bytes_csv(labelled(result.rewired_pairs())),
                       file_name=f"{comparison}_Rewired_Pairs.csv", mime="text/csv", key=f"{key_prefix}_dl_pairs")
    all_pairs = result.pairs.drop(columns=["i", "j"]).rename(columns={"r_a": f"r_{ga}", "r_b": f"r_{gb}"})
    d2.download_button("All pairs (CSV)", utils.to_download_bytes_csv(labelled(all_pairs)),
                       file_name=f"{comparison}_All_Pair_Correlations.csv", mime="text/csv", key=f"{key_prefix}_dl_all")
    d3.download_button("Hypothesis evidence (JSON)", result.evidence_json().encode("utf-8"),
                       file_name=f"{comparison}_Rewiring_Hypotheses.json", mime="application/json",
                       key=f"{key_prefix}_dl_json")
    report = render_view_html(result, title=f"{ga} vs {gb}", theme="light",
                              generated=f"{source_label} · {_dt.date.today().isoformat()}",
                              pw_top=pw_top, pw_metric=pw_metric, pw_cutoff=pw_cut)
    d4.download_button("Interactive report (HTML)", report.encode("utf-8"),
                       file_name=f"{comparison}_Rewiring_Map.html", mime="text/html", key=f"{key_prefix}_dl_html")
    e1, e2, e3, _ = st.columns(4)
    pt = result.pathway_table_display()
    e1.download_button("Pathway-level tests (CSV)", utils.to_download_bytes_csv(labelled(pt, "pathways", (ga, gb)) if len(pt) else pd.DataFrame()),
                       file_name=f"{comparison}_Pathway_Rewiring.csv", mime="text/csv", key=f"{key_prefix}_dl_pw",
                       disabled=not len(pt))
    rt = result.reaction_table_display()
    e2.download_button("Reaction pairs (CSV)", utils.to_download_bytes_csv(labelled(rt) if len(rt) else pd.DataFrame()),
                       file_name=f"{comparison}_Reaction_Pairs.csv", mime="text/csv", key=f"{key_prefix}_dl_rx",
                       disabled=not len(rt))
    idm = result.prep.id_map
    e3.download_button("Metabolite ID / pathway mapping (CSV)",
                       utils.to_download_bytes_csv(idm if idm is not None else pd.DataFrame()),
                       file_name=f"{comparison}_ID_Mapping.csv", mime="text/csv", key=f"{key_prefix}_dl_ids",
                       disabled=idm is None)

    with st.expander("Tables: pathway tests, reaction pairs, rewired pairs, per-metabolite statistics, modules, ID mapping"):
        if len(pt):
            st.markdown(f"**Pathway-level rewiring** ({S['pathway_library']}; sets with ≥ {result.cfg.pathway_min_size} "
                        "measured members; mean z² of the pairs inside the set vs. label permutations; FDR q-value = "
                        "Benjamini–Hochberg across sets)")
            st.dataframe(utils.format_df_for_display(labelled(pt, "pathways", (ga, gb)).set_index("Pathway")), width="stretch", height=300)
        if len(rt):
            st.markdown("**Reaction-linked pairs** (Human-GEM; product/substrate log2 ratio per group and its shift, "
                        "Welch t-test; FDR q-value = Benjamini–Hochberg across linked pairs)")
            st.dataframe(utils.format_df_for_display(labelled(rt)), width="stretch", height=300)
        st.markdown("**Rewired pairs**")
        st.dataframe(utils.format_df_for_display(labelled(result.rewired_pairs())), width="stretch", height=320)
        st.markdown(f"**Per-metabolite level change** (Welch t-test on log2 values, {gb} vs {ga}; FDR q-value = Benjamini–Hochberg)")
        st.dataframe(utils.format_df_for_display(labelled(result.features, "features").set_index("feature")), width="stretch", height=320)
        if result.hypotheses:
            mods = pd.DataFrame([{
                "Rank": h["rank"], "Module": h["pathway_a"] if h["pathway_a"] == h["pathway_b"] else f"{h['pathway_a']} × {h['pathway_b']}",
                "Class": h["class"], "Pairs": h["n_pairs"], "Reaction-linked pairs": len(h.get("reaction_pairs", [])),
                f"Mean r {ga}": h["r_a"], f"Mean r {gb}": h["r_b"],
                "Bootstrap stability": h["stability"], "Method agrees": h["method_check"]["ok"],
                "Outlier-robust": h["influence_check"]["ok"], "Imputation-robust": h["imputation_check"]["ok"],
                "Verdict": h["verdict"],
            } for h in result.hypotheses]).set_index("Rank")
            st.markdown("**All hypothesis modules**")
            st.dataframe(utils.format_df_for_display(mods), width="stretch")
        if idm is not None:
            st.markdown("**Metabolite ID and pathway mapping** (MSEA ID standardization + Human-GEM reaction network)")
            st.dataframe(utils.format_df_for_display(idm), width="stretch", height=300)

    info = rr.reaction_library().info
    with st.expander("How this works"):
        st.markdown(
            f"""
- **Per-group correlation.** For every metabolite pair, correlation is computed separately in {ga} and {gb}
  ({result.cfg.method.title()}).
- **Test.** The difference is scored with Fisher's z, z = (atanh r_B − atanh r_A) / SE, with
  SE = √(1/(n_A−3) + 1/(n_B−3)) (×√1.06 for Spearman). With 3 samples in a group the SE is undefined; the raw
  Fisher-z difference is used and every p-value comes from label permutations.
- **FDR.** Group labels are permuted (on group-centred data) and every pair is recomputed; the
  pooled null gives a permutation FDR that stays calibrated with small groups. When there are fewer distinct label
  splits than permutations requested, all of them are enumerated (exact test; the smallest possible p-value is
  1 / number of splits). The Benjamini–Hochberg FDR q-value of the analytic p-values ("FDR q-value (BH)") is
  reported alongside. Throughout this page, *FDR q-value* means the false-discovery-rate-adjusted p-value.
- **Pairs tested.** *All pairs*, or only pairs linked by a reaction: *direct* (one enzymatic reaction interconverts
  them) or *≤ 2 steps* (A → X → B through one unmeasured, non-hub intermediate, following reaction direction).
  Testing fewer, mechanistically linked pairs gives more power.
- **Reaction library.** {info['n_pairs']:,} substrate→product pairs derived from Human-GEM v{info['version']}
  ({info['source']}): compartments collapsed, currency species removed (NAD(P)(H), ATP/ADP, … only when acting as
  cofactors), substrates matched to products by conserved atoms (maximum common substructure), routes through
  enzyme-bound intermediates contracted (PDH, OGDH). Checked against 89 textbook conversions and 24 negative controls
  (see VALIDATION_REWIRING.md). For every linked pair the product/substrate log2 ratio is compared between groups.
- **Network display.** The ring shows the top 5 / 10 / 15 / 20 (or all) pathways of the pathway-level test that pass
  a p-value or FDR q-value cutoff (default: top 10 by p-value ≤ 0.05); each metabolite is drawn once, under its
  best-ranked selected pathway. The ring grows with the number of metabolites and pathway names sit outside the
  arcs (with a leader line when moved), so no name overlaps a line, a dot or another name. Correlation plots show
  relative abundance on the analysed log2 scale. The buttons under each image download SVG, PNG, JPEG, TIFF or PDF,
  including the legend (edge classes, reactions, log2 fold-change colour scale).
- **Pathway-level test (MSEA-style).** Metabolites are mapped with the same ID standardization and libraries as the
  Pathway Analysis page. For each pathway with ≥ {result.cfg.pathway_min_size} measured members, the mean z² of the
  pairs inside it is compared with the same label permutations (BH across pathways). Pooling pairs gives power
  when groups are small.
- **Classes.** A pair is *rewired* when its {"p-value" if result.cfg.cutoff_on == "p" else "FDR q-value"} ≤ {result.cfg.fdr:g}
  and |Δr| ≥ {result.cfg.min_delta}: *flipped* if the
  sign reverses with |r| ≥ {result.cfg.flip_min_abs_r} in both groups, otherwise *lost* or *gained* by which group is
  stronger. *Intact* pairs are strongly correlated (|r| ≥ {result.cfg.intact_min_abs_r}) with the same sign in both.
- **Hypothesis cards.** Rewired pairs are grouped by pathway pair and class, reaction-linked pairs first. Each module
  is checked for bootstrap stability, agreement with the other correlation method, robustness to removing outlying
  samples and imputed values. Verdict: *robust* (all pass), *mixed* (one fails),
  *fragile* (two or more fail, or stability < 40%), *exploratory* (fewer than 5 samples per group, so resampling
  checks cannot be assessed).
- **Interpretations** are rule-based text matched on pathway names. They suggest what to test next; they are not findings.
"""
        )
