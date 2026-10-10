"""Step 2 — experimental design: groups, treatments, time points, replicates, tracer and correction."""
from __future__ import annotations

import streamlit as st

from core import design as dz, labeling as lb, correction as cr
from ui import common as C

TRACERS = {"[U-¹³C₆]-glucose": 6, "[U-¹³C₅]-glutamine": 5}


def render():
    C.header(2, "Experimental design",
             "Define groups, treatments, labeling time points and replicates, choose the samples to analyse, and "
             "describe the tracer and natural-isotope correction. Groups (group × treatment) and time points drive "
             "every comparison downstream.")
    if not C.require("upload"):
        return
    s = C.ss()
    st.markdown("#### Samples")
    st.caption("Edit Group, Treatment, Time and Replicate if needed (Treatment is optional — when given, each "
               "group × treatment combination becomes its own group). Untick Include to leave a sample out.")
    D = s.design.copy()
    D["TimeUnit"] = D.TimeUnit.fillna("")
    edited = st.data_editor(
        D, key="design_editor", hide_index=True, disabled=["Sample"],
        column_config={
            "Sample": st.column_config.TextColumn("Sample"),
            "SampleType": st.column_config.SelectboxColumn("Type", options=["13C labeled", "12C reference", "Blank", "QC"]),
            "Group": st.column_config.TextColumn("Group / genotype"),
            "Treatment": st.column_config.TextColumn("Treatment (optional)"),
            "TimeValue": st.column_config.NumberColumn("Time", min_value=0.0, step=0.5),
            "TimeUnit": st.column_config.SelectboxColumn("Unit", options=["", "min", "h", "d"]),
            "Replicate": st.column_config.NumberColumn("Biological replicate", min_value=0, step=1),
            "Include": st.column_config.CheckboxColumn("Include"),
        })
    if not edited.equals(D):
        s.design = edited
        if s.design_confirmed:
            C.reset_downstream("design")
            s.design_confirmed = False
    Df = C.design_final()
    An = dz.analysed(Df)
    if An.empty:
        st.error("No labeled samples are included.")
        return
    summary = An.groupby(["Group", "Hours"]).size().unstack(fill_value=0).reindex(dz.groups_of(Df))
    summary.columns = [f"{h:g} h" for h in summary.columns]
    c1, c2 = st.columns([2, 3])
    with c1:
        st.markdown("**Replicates per group and time**")
        C.show_df(summary)
    with c2:
        w = dz.design_warnings(Df)
        for x in w:
            C.note("⚠ " + x)
        if not w:
            C.note("Design looks complete: ≥ 3 replicates per group and time, ≥ 3 time points, unlabeled reference "
                   "present.", ok=True)
        st.caption(f"Groups: {', '.join(dz.groups_of(Df))} · time points: "
                   f"{', '.join(f'{h:g} h' for h in dz.times_of(Df))}.")

    st.markdown("#### Tracer")
    t = s.tracer
    c = st.columns(4)
    names = list(TRACERS)
    t["name"] = c[0].selectbox("Tracer", names, index=names.index(t["name"]) if t["name"] in names else 0,
                               key="tracer_name")
    t["element"] = "13C"
    t["labeled_atoms"] = TRACERS[t["name"]]
    c[1].text_input("Labeled element", "13C", disabled=True)
    c[2].text_input("Labeled atoms in tracer", str(t["labeled_atoms"]), disabled=True)
    t["purity"] = float(c[3].number_input("Isotopic purity", 0.90, 1.0, float(t["purity"]), 0.001, format="%.3f"))
    with st.expander("Tracer enrichment in the medium"):
        t["enrichment"] = float(st.slider("Fraction of the substrate supplied as tracer", 0.0, 1.0,
                                          float(t["enrichment"]), 0.05))
        st.caption("Recorded in the report. Labeling values are not rescaled by this number.")

    st.markdown("#### Natural-isotope correction")
    ref = dz.reference_samples(Df)
    verdict = None
    if ref:
        obs, forms = lb.reference_m1_ct(s.ds.labeling_table, ref[0])
        verdict = cr.check_correction_status(obs, forms, t["element"])
        msg = f"<b>Automatic check using {ref[0]}:</b> {verdict['text']}"
        if verdict["verdict"] == "uncorrected":
            st.warning(f"Automatic check using {ref[0]}: {verdict['text']}")
        else:
            C.note(msg)
    else:
        st.caption("No unlabeled reference sample: correction status cannot be checked automatically.")
    corr = s.correction
    for k_, v_ in dz.CorrectionSettings().to_dict().items():          # settings files from older versions
        corr.setdefault(k_, v_)
    opts = {"already": "Data already corrected by the vendor software (use as supplied)",
            "apply": "Apply natural-abundance correction here (data not yet corrected)",
            "none": "No correction (analyse uncorrected data)"}
    # the user decides; nothing is corrected automatically (Compound Discoverer corrects against the 12C reference)
    corr["mode"] = st.radio("Correction", list(opts), format_func=opts.get,
                            index=list(opts).index(corr["mode"]) if corr["mode"] in opts else 0, key="corr_mode")
    if verdict and verdict["verdict"] == "uncorrected" and corr["mode"] != "apply":
        C.note("The unlabeled reference suggests the data are not yet corrected — choose <b>Apply natural-abundance "
               "correction here</b> and a correction method below.")
    if verdict and verdict["verdict"] == "corrected" and corr["mode"] == "apply":
        st.error("The unlabeled reference indicates the data are already corrected. Correcting again distorts "
                 "isotopologue distributions.")
    methods = list(cr.METHODS)
    corr["method"] = st.selectbox("Correction method", methods, format_func=cr.METHODS.get,
                                  index=methods.index(corr["method"]) if corr["method"] in methods else 0,
                                  key="corr_method", disabled=corr["mode"] != "apply",
                                  help="All methods build the correction matrix from IUPAC natural isotope abundances "
                                       "and the tracer isotopic purity set above.")
    st.caption(("Applied: " if corr["mode"] == "apply" else "Method used when correction is applied: ")
               + cr.METHOD_REFS[corr["method"]])
    with st.expander("Advanced correction settings", expanded=corr["mode"] == "apply" and
                     corr["method"] in ("accucor", "lsq")):
        if corr["method"] == "accucor":
            a = st.columns(3)
            ins = list(cr.INSTRUMENTS)
            corr["instrument"] = a[0].selectbox("Instrument", ins, format_func=cr.INSTRUMENTS.get,
                                                index=ins.index(corr["instrument"]) if corr["instrument"] in ins else 0,
                                                key="corr_instr")
            corr["resolving_power"] = float(a[1].number_input(
                "Resolving power (FWHM)", 5000.0, 2000000.0, float(corr["resolving_power"]), 5000.0, format="%.0f",
                key="corr_rp"))
            corr["ref_mz"] = float(a[2].number_input("…specified at m/z", 50.0, 2000.0, float(corr["ref_mz"]), 10.0,
                                                     key="corr_refmz"))
            st.caption("A natural isotope of a non-tracer element (²H, ¹⁵N, ¹⁷O, ¹⁸O, ³³S, ³⁴S, Si) is corrected only "
                       "when its mass difference from the tracer isotopologue is below 1.66 × FWHM at the "
                       "metabolite mass (AccuCor criterion), so the correction follows the actual resolution.")
        elif corr["method"] == "lsq":
            corr["resolution"] = st.radio(
                "Elements corrected", ["high", "low"], index=["high", "low"].index(corr["resolution"]),
                key="corr_res",
                format_func={"high": "Tracer element only (high-resolution data)",
                             "low": "All elements — C, H, N, O, S, Si (nominal-mass data)"}.get)
            st.caption("Unconstrained least-squares solution of the skewed correction matrix; negative values are set "
                       "to 0 and the distribution is renormalized.")
        else:
            st.caption("IsoCor-type correction: " + ("natural isotopes of the tracer element only (Orbitrap / Q-TOF "
                       "high-resolution data)." if corr["method"] == "isocor_hr" else
                       "natural isotopes of all elements (nominal-mass data)."))
        st.caption("Correction applies to isotopologue fractions and labeling metrics. Intensity analyses use the "
                   "measured channel areas.")

    groups = dz.groups_of(Df)
    # reference / compared group used as defaults later (chosen in each analysis step)
    if s.group_a not in groups:
        s.group_a = groups[0]
    if s.group_b not in groups or s.group_b == s.group_a:
        s.group_b = next((g for g in groups if g != s.group_a), groups[0])
    st.write("")
    if st.button("Confirm design", type="primary"):
        s.design_confirmed = True
        C.reset_downstream("design")
        s.step = "normalize"
        st.rerun()
    if s.design_confirmed:
        st.caption("✓ Design confirmed. Changing it later recalculates downstream steps.")
