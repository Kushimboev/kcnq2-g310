#!/usr/bin/env python3
"""Read-out of the mixed tetramers (H1/H2): the rule fixed before the runs finished is applied mechanically.
Primary measure: fraction of the G310S effect on the gate water,  f = (W_WT - W_H) / (W_WT - W_GS),  W_WT = 8.73 (n 4), W_GS = 4.73 (n 3).
Additive model: f(H1) = 0.25, f(H2) = 0.50.
Decision:
  the CI contains 0.25 (H1) / 0.50 (H2)                -> "additive: each serine contributes about a quarter"
  f(H1) >= 0.5 AND the CI excludes 0.25                -> "one serine gives most of the effect"
  the CI of f(H1) contains 0 and excludes 0.25         -> "one serine is tolerated"
CI: delta method, in two variants because n = 2: (a) replica SE, (b) 70 ns block means.
The WIDER (more cautious) interval is used for the decision.
Secondary: net filter charge (first 20 ns excluded) -> H/WT ratio, conditional-binomial CI.
usage: hyb_readout.py [--out FILE]
"""
import os, re, sys, glob
import numpy as np
from scipy import stats

P = "${KCNQ2_ROOT}/13_conduction"
M23 = f"{P}/scripts/analysis_2026-09-23/out/mech_0925"      # WT/G310S mechanism npy
MH  = f"{P}/scripts/analysis_2026-09-29_hybrid/mech"          # mixed-tetramer mechanism npy
BH  = f"{P}/scripts/analysis_2026-09-29_hybrid"               # mixed-tetramer books
B23 = f"{P}/scripts/analysis_2026-09-23/out"                  # WT/G310S books
e = 1.602176634e-19; V = 0.565
GW = 6            # column of the mech npy: gate_water (time,zg,d1,d2,car,ogr,GW,cw,ck,bk)

def gw_series(d, run):
    f = f"{d}/{run}.npy"
    if not os.path.exists(f): return None
    R = np.load(f); return R[:, 0] / 1000.0, R[:, GW]     # ns, gate water (1 frame per ns)

def rep_mean(d, runs):
    out = {}
    for r in runs:
        s = gw_series(d, r)
        if s is not None: out[r] = s
    return out

def blocks(t, y, bw=70.0):
    """70 ns blocks (inside the analysis window) -> block means"""
    out = []; t0 = t[0]
    while t0 + bw <= t[-1] + 1e-9:
        m = (t >= t0) & (t < t0 + bw)
        if m.sum() >= 0.8 * bw: out.append(y[m].mean())
        t0 += bw
    return np.array(out)

def mean_se(series, mode):
    """mode 'rep': SE from the replica means; 'blk': SE from all block means"""
    if mode == "rep":
        v = np.array([y.mean() for _, y in series.values()])
    else:
        v = np.concatenate([blocks(t, y) for t, y in series.values()])
    n = len(v)
    return v.mean(), (v.std(ddof=1) / np.sqrt(n) if n > 1 else float("nan")), n

def f_ci(A, sA, H, sH, Gg, sG):
    """delta method: f = (A-H)/(A-G)"""
    den = A - Gg; f = (A - H) / den
    dA = (H - Gg) / den**2; dH = -1.0 / den; dG = (A - H) / den**2
    var = (dA * sA)**2 + (dH * sH)**2 + (dG * sG)**2
    se = np.sqrt(var)
    return f, se, (f - 1.96 * se, f + 1.96 * se)

def books(d, run):
    f = f"{d}/{run}_b20.books.txt"
    if not os.path.exists(f): return None
    t = open(f).read()
    m1 = re.search(r"_b20: ([0-9.]+) ns", t); m2 = re.search(r"NET outward charge .*?= (-?\d+)", t)
    return (int(m2.group(1)), float(m1.group(1))) if m1 and m2 else None

def ratio(k, t, N, T):
    m = k + N; fr = t / (t + T)
    plo = stats.beta.ppf(0.025, k, m - k + 1) if k > 0 else 0.0
    phi = stats.beta.ppf(0.975, k + 1, m - k) if m - k > 0 else 1.0
    rr = lambda p: p / (1 - p) * (T / t) if p < 1 else float("inf")
    return (k / t) / (N / T), rr(plo), rr(phi), stats.binom.cdf(k, m, fr)

WT = ["ANE_r1", "ANE_r2", "ANE_r4", "ANE_r5"]; GS = ["GS_r1", "GS_r2", "GS_r3"]
H1 = ["H1_r1", "H1_r2"]; H2 = ["H2_r1", "H2_r2"]
sWT = rep_mean(M23, WT); sGS = rep_mean(M23, GS); sH1 = rep_mean(MH, H1); sH2 = rep_mean(MH, H2)

L = []
L.append("# MIXED-TETRAMER READ-OUT (the rule written BEFORE the result is applied mechanically)")
L.append(f"# {__import__('datetime').datetime.now():%F %T}  ·  window: first 20 ns discarded, 1 frame/ns")
L.append("")
L.append("## Gate water (primary measure) — replica means")
for nm, s in (("WT", sWT), ("G310S", sGS), ("H1 (1xS310)", sH1), ("H2 (2xS310)", sH2)):
    if not s: L.append(f"   {nm:14s} — not available yet"); continue
    for r, (t, y) in s.items():
        L.append(f"   {nm:14s} {r:8s} {y.mean():5.2f} ± {y.std():4.2f}   ({t[0]:.0f}-{t[-1]:.0f} ns, n={len(y)})")
L.append("")

res = {}
for nm, s, add in (("H1", sH1, 0.25), ("H2", sH2, 0.50)):
    if not s: L.append(f"## {nm}: no result yet"); L.append(""); continue
    rows = []
    for mode, lab in (("rep", "replica SE"), ("blk", "70 ns block means")):
        A, sA, nA = mean_se(sWT, mode); Gg, sG, nG = mean_se(sGS, mode); H, sH, nH = mean_se(s, mode)
        f, se, (lo, hi) = f_ci(A, sA, H, sH, Gg, sG)
        rows.append((lab, f, se, lo, hi, nA, nG, nH, A, Gg, H))
    L.append(f"## {nm} — f = (W_WT - W_{nm}) / (W_WT - W_GS);  additive expectation f = {add:.2f}")
    for lab, f, se, lo, hi, nA, nG, nH, A, Gg, H in rows:
        tag = " (n=1: CI yo`q)" if not np.isfinite(se) else ""
        L.append(f"   {lab+tag:38s} W_WT {A:5.2f} · W_GS {Gg:5.2f} · W_{nm} {H:5.2f}  ->  f = {f:5.2f} ± {se:4.2f}  (95 % {lo:5.2f}-{hi:5.2f})  [n {nA}/{nG}/{nH}]")
    # decision: with the WIDER interval (rows with a nan interval, n = 1, are ignored)
    okrows = [r for r in rows if np.isfinite(r[3]) and np.isfinite(r[4])]
    if not okrows:
        L.append("   >>> DECISION: no CI (too few replicas) — result INCOMPLETE"); L.append(""); continue
    wide = max(okrows, key=lambda r: r[4] - r[3])
    lab, f, se, lo, hi = wide[0], wide[1], wide[2], wide[3], wide[4]
    inc = lambda x: lo <= x <= hi
    if inc(add):
        v = f"ADDITIVE — {'one serine gives about a quarter' if nm == 'H1' else 'two serines give about a half'} of the full effect (the CI contains {add:.2f})"
    elif nm == "H1" and f >= 0.5 and not inc(0.25):
        v = "ONE SERINE GIVES MOST OF THE EFFECT (compatible with a dominant-negative action)"
    elif nm == "H1" and inc(0.0) and not inc(0.25):
        v = "ONE SERINE IS TOLERATED"
    else:
        v = f"none of the three branches of the rule applies exactly (f {f:.2f}, CI {lo:.2f}-{hi:.2f}) — to be reported as such"
    L.append(f"   >>> DECISION ({lab} — the wider CI): {v}")
    L.append("")
    res[nm] = (f, lo, hi)

L.append("## Conduction (secondary measure; not a conclusion on its own)")
pw = [books(B23, r) for r in WT]; pw = [x for x in pw if x]
Nw = sum(n for n, _ in pw); Tw = sum(t for _, t in pw)
L.append(f"   WT   {Nw} / {Tw:.1f} ns  (G {Nw*e/(Tw*1e-9*V)*1e12:.2f} pS)")
for nm, runs in (("G310S", GS), ("H1", H1), ("H2", H2)):
    d = [(r, books(B23 if nm == 'G310S' else BH, r)) for r in runs]
    ok = [(r, x) for r, x in d if x]
    for r, x in d: L.append(f"   {nm:6s} {r:8s} " + (f"{x[0]:3d} / {x[1]:6.1f} ns" if x else "missing"))
    if ok:
        N = sum(x[0] for _, x in ok); T = sum(x[1] for _, x in ok)
        r_, lo_, hi_, p_ = ratio(N, T, Nw, Tw)
        L.append(f"   {nm:6s} POOLED {N} / {T:.1f} ns -> {nm}/WT = {r_:.2f} (95 % {lo_:.2f}-{hi_:.2f}), p = {p_:.3g}"
                 f"   [expected events: ~{Nw/Tw*T:.1f} at the WT rate, ~{2/840*T:.1f} at the G310S rate]")
L.append("")
L.append("NOTE: the blocks may be correlated, so the block-based CI may be optimistic; the decision uses the WIDER CI.")
txt = "\n".join(L) + "\n"
out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
if out: open(out, "w").write(txt)
print(txt)
