"""
core.correction — natural-isotope-abundance correction of mass isotopologue distributions (MIDs).

Correction-matrix method (as in IsoCor / AccuCor):
    measured MID  m = A · c      (c = tracer-derived label distribution)
A[i, j] = probability that a molecule carrying j tracer-derived labels is observed at mass shift +i, given
    * tracer purity:   each of the j tracer positions is truly heavy with probability = purity
    * natural abundance of the (n - j) unlabeled positions of the tracer element
    * (low-resolution mode only) natural heavy isotopes of the other elements (H, N, O, S, Si)
c is solved by non-negative least squares and renormalized to sum 1.

High-resolution mode (default) corrects only for the tracer element, appropriate for Orbitrap / high-res
data where isotopologues of other elements are resolved from the tracer isotopologues.

Correction methods offered (v2.4.2, all with tracer-purity correction):
    isocor_hr   IsoCor-type matrix + NNLS, tracer element only (high resolution)
                Millard et al., Bioinformatics 2012; Millard et al., Bioinformatics 2019 (IsoCor v2)
    isocor_lr   IsoCor / IsoCorrectoR-type matrix + NNLS, natural isotopes of ALL elements (nominal mass)
                Heinrich et al., Sci Rep 2018 (IsoCorrectoR)
    accucor     resolution-dependent matrix + NNLS (AccuCor): a natural isotope of a non-tracer element is corrected
                only when the instrument cannot resolve it from the tracer isotopologue
                (|Δm| < 1.66 × FWHM at the metabolite mass). Su et al., Anal Chem 2017
    lsq         classical skewed correction matrix solved by (unconstrained) least squares, negative values set to 0
                and renormalized. Fernandez et al., J Mass Spectrom 1996; van Winden et al., Biotechnol Bioeng 2002

IMPORTANT: Compound Discoverer can already report corrected isotopologue areas. Correcting twice distorts the
MID. `check_correction_status` uses the unlabeled reference to tell whether data look corrected.
"""
from __future__ import annotations

import re
from math import comb

import numpy as np
from scipy.optimize import nnls

# natural abundances of heavy isotopes (IUPAC representative values): {element: [(mass shift, abundance), ...]}
NATURAL = {
    "C": [(1, 0.0107)],
    "H": [(1, 0.000115)],
    "N": [(1, 0.00364)],
    "O": [(1, 0.00038), (2, 0.00205)],
    "S": [(1, 0.0075), (2, 0.0425), (4, 0.0001)],
    "Si": [(1, 0.0467), (2, 0.0310)],
}
TRACER_ELEMENT = {"13C": "C", "15N": "N", "2H": "H"}
# exact mass difference heavy − light isotope (Da): {element: {mass shift: Δm}}
ISO_DELTA = {
    "C": {1: 1.0033548},
    "H": {1: 1.0062767},
    "N": {1: 0.9970349},
    "O": {1: 1.0042169, 2: 2.0042449},
    "S": {1: 0.9993878, 2: 1.9957959, 4: 3.9950100},
    "Si": {1: 0.9995685, 2: 1.9968437},
}
MONO = {"C": 12.0, "H": 1.00782503, "N": 14.00307401, "O": 15.99491462, "P": 30.97376151, "S": 31.97207069,
        "Si": 27.97692653, "Cl": 34.96885271, "Br": 78.9183376, "F": 18.9984032, "I": 126.904468, "Na": 22.98976966,
        "K": 38.9637069, "Se": 79.9165218, "Fe": 55.9349375, "Mg": 23.9850417, "Co": 58.933195}

METHODS = {
    "isocor_hr": "IsoCor — correction matrix + NNLS, tracer element only (high resolution)",
    "isocor_lr": "IsoCor / IsoCorrectoR — correction matrix + NNLS, all elements (low / nominal resolution)",
    "accucor": "AccuCor — resolution-dependent correction matrix + NNLS (Orbitrap / FT-ICR / TOF)",
    "lsq": "Classical skewed correction matrix — least-squares inversion (Fernandez 1996)",
}
METHOD_REFS = {
    "isocor_hr": "Millard P. et al. Bioinformatics 2012;28:1294 · Millard P. et al. Bioinformatics 2019;35:4484",
    "isocor_lr": "Millard P. et al. Bioinformatics 2012;28:1294 · Heinrich P. et al. Sci Rep 2018;8:17910",
    "accucor": "Su X. et al. Anal Chem 2017;89:5940",
    "lsq": "Fernandez C.A. et al. J Mass Spectrom 1996;31:255 · van Winden W.A. et al. Biotechnol Bioeng 2002;80:477",
}
INSTRUMENTS = {"orbitrap": "Orbitrap (FWHM ∝ m^1.5)", "fticr": "FT-ICR (FWHM ∝ m²)", "tof": "Q-TOF (constant resolving power)"}


def parse_formula(f: str) -> dict:
    d = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", str(f).replace(" ", "")):
        d[el] = d.get(el, 0) + int(n or 1)
    return d


def _element_dist(count: int, isotopes, size: int) -> np.ndarray:
    """Mass-shift distribution for `count` atoms of one element (multinomial via repeated convolution)."""
    single = np.zeros(size)
    single[0] = 1 - sum(a for _, a in isotopes)
    for s, a in isotopes:
        if s < size:
            single[s] += a
    out = np.zeros(size); out[0] = 1.0
    for _ in range(count):
        out = np.convolve(out, single)[:size]
    return out


def _binom(n: int, p: float, size: int) -> np.ndarray:
    v = np.zeros(size)
    for k in range(min(n, size - 1) + 1):
        v[k] = comb(n, k) * p ** k * (1 - p) ** (n - k)
    return v


def monoisotopic_mass(formula: str) -> float:
    comp = parse_formula(formula)
    return float(sum(MONO.get(e, 0.0) * n for e, n in comp.items()))


def fwhm(mz: float, resolving_power: float, ref_mz: float = 200.0, instrument: str = "orbitrap") -> float:
    """Peak width (Da) at m/z for a resolving power specified at ref_mz."""
    R = max(float(resolving_power), 1.0)
    if instrument == "fticr":
        return mz ** 2 / (R * ref_mz)
    if instrument == "tof":
        return mz / R
    return mz ** 1.5 / (R * np.sqrt(ref_mz))          # Orbitrap (AccuCor)


def unresolved_isotopes(formula: str, tracer: str = "13C", resolving_power: float = 1e5, ref_mz: float = 200.0,
                        instrument: str = "orbitrap") -> dict:
    """{element: [(shift, abundance), …]} of non-tracer natural isotopes that the instrument does NOT resolve from the
    tracer isotopologue with the same nominal shift (AccuCor criterion |Δm| < 1.66 × FWHM)."""
    el = TRACER_ELEMENT.get(tracer, "C")
    d_tr = ISO_DELTA[el][1]
    lim = 1.66 * fwhm(monoisotopic_mass(formula) or 200.0, resolving_power, ref_mz, instrument)
    out = {}
    for e in parse_formula(formula):
        if e == el or e not in NATURAL:
            continue
        keep = [(k, a) for k, a in NATURAL[e] if abs(ISO_DELTA[e][k] - k * d_tr) < lim]
        out[e] = keep
    return out


def _element_dist_partial(count: int, isotopes, observed, size: int) -> np.ndarray:
    """Like _element_dist, but only isotopes in `observed` fall into the measured channels; molecules carrying a
    resolved heavy isotope are measured elsewhere (AccuCor)."""
    single = np.zeros(size)
    single[0] = 1 - sum(a for _, a in isotopes)
    for s_, a in observed:
        if s_ < size:
            single[s_] += a
    out = np.zeros(size); out[0] = 1.0
    for _ in range(count):
        out = np.convolve(out, single)[:size]
    return out


def correction_matrix(formula: str, n_iso: int, tracer: str = "13C", purity: float = 0.99,
                      resolution: str = "high", method: str | None = None, resolving_power: float = 1e5,
                      ref_mz: float = 200.0, instrument: str = "orbitrap") -> np.ndarray:
    """A with shape (n_iso, n_iso): rows observed shift 0..n_iso-1, columns tracer labels 0..n_iso-1.
    resolution 'high' = tracer element only, 'low' = all elements; method 'accucor' = resolution-dependent."""
    el = TRACER_ELEMENT.get(tracer, "C")
    comp = parse_formula(formula)
    n_el = comp.get(el, n_iso - 1)
    size = n_iso + 4
    nat_tracer = sum(a for _, a in NATURAL.get(el, [(1, 0.0)]) if _ == 1)
    other = np.zeros(size); other[0] = 1.0
    if method == "isocor_lr":
        resolution = "low"
    elif method == "isocor_hr":
        resolution = "high"
    if method == "accucor":
        unres = unresolved_isotopes(formula, tracer, resolving_power, ref_mz, instrument)
        for e, cnt in comp.items():
            if e != el and e in NATURAL:
                other = np.convolve(other, _element_dist_partial(cnt, NATURAL[e], unres.get(e, []), size))[:size]
    elif resolution == "low":
        for e, cnt in comp.items():
            if e != el and e in NATURAL:
                other = np.convolve(other, _element_dist(cnt, NATURAL[e], size))[:size]
    A = np.zeros((n_iso, n_iso))
    for j in range(n_iso):
        lab = _binom(j, purity, size)                                  # truly heavy among tracer positions
        nat = _binom(max(n_el - j, 0), nat_tracer, size)               # natural heavy among the rest
        col = np.convolve(np.convolve(lab, nat)[:size], other)[:size]
        A[:, j] = col[:n_iso]
    return A


def solve(A: np.ndarray, m: np.ndarray, method: str | None = None) -> np.ndarray:
    """Corrected label distribution (renormalized). NNLS for IsoCor / AccuCor, least squares for 'lsq'."""
    if method == "lsq":
        c = np.linalg.lstsq(A, m, rcond=None)[0]
        c = np.clip(c, 0, None)
    else:
        c, _ = nnls(A, m)
    return c / c.sum() if c.sum() > 0 else np.full_like(m, np.nan)


def correct_mid(measured, formula: str, tracer: str = "13C", purity: float = 0.99,
                resolution: str = "high", method: str | None = None, **kw) -> np.ndarray:
    """Corrected, renormalized MID (fractions summing to 1). `measured` = areas or fractions M+0..M+n."""
    m = np.asarray(measured, float)
    if not np.isfinite(m).all() or m.sum() <= 0:
        return np.full_like(m, np.nan)
    m = m / m.sum()
    A = correction_matrix(formula, len(m), tracer, purity, resolution, method, **kw)
    return solve(A, m, method)


def natural_mid(formula: str, n_iso: int, tracer: str = "13C", resolution: str = "high") -> np.ndarray:
    """Expected MID of an unlabeled molecule (what an uncorrected 12C reference should look like)."""
    return correction_matrix(formula, n_iso, tracer, 1.0, resolution)[:, 0]


def check_correction_status(ref_fractions: dict, formulas: dict, tracer: str = "13C") -> dict:
    """
    ref_fractions: {metabolite: observed M+1 fraction in the unlabeled reference}
    Compares observed M+1 with the natural-abundance expectation.
    Returns verdict 'corrected' / 'uncorrected' / 'unclear' with the supporting numbers.
    """
    ratios = []
    for met, f1 in ref_fractions.items():
        exp1 = natural_mid(formulas[met], 3, tracer)[1] / natural_mid(formulas[met], 3, tracer)[0]
        if exp1 > 0.02 and np.isfinite(f1):
            ratios.append(f1 / exp1)
    if len(ratios) < 5:
        return dict(verdict="unclear", n=len(ratios), median_ratio=np.nan,
                    text="Too few metabolites in the unlabeled reference to judge correction status.")
    med = float(np.median(ratios))
    if med < 0.35:
        v, t = "corrected", ("The unlabeled reference shows far less M+1 than natural abundance predicts "
                             f"(median {med:.0%} of expected): the data look ALREADY corrected. Do not correct again.")
    elif med > 0.7:
        v, t = "uncorrected", ("The unlabeled reference shows M+1 close to natural abundance "
                               f"(median {med:.0%} of expected): the data look UNCORRECTED. Apply correction.")
    else:
        v, t = "unclear", f"M+1 in the unlabeled reference is {med:.0%} of the natural-abundance expectation."
    return dict(verdict=v, n=len(ratios), median_ratio=med, text=t)
