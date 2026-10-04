#!/usr/bin/env python3
"""Final summary of the conduction counts (net charge through the selectivity filter, first 20 ns excluded: out/<RUN>_b20.books.txt), read with the decision rules fixed before the runs.
Output: out/FINAL_SUMMARY.txt. Runs that are not available are marked as missing; the file is rewritten on every call."""
import re, os, datetime
import numpy as np
from scipy import stats
O = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out"); e = 1.602176634e-19
def rd(r):
    f = f"{O}/{r}_b20.books.txt"
    if not os.path.exists(f): return None
    t = open(f).read(); m1 = re.search(r"_b20: ([0-9.]+) ns", t); m2 = re.search(r"NET outward charge .*?= (-?\d+)", t)
    return (int(m2.group(1)), float(m1.group(1))) if m1 and m2 else None
def G(n, t, V): return n * e / (t * 1e-9 * V) * 1e12
def ci(n): return (stats.chi2.ppf(0.025, 2 * n) / 2 if n > 0 else 0.0, stats.chi2.ppf(0.975, 2 * n + 2) / 2)
def pool(runs, V):
    d = {r: rd(r) for r in runs}; ok = {r: v for r, v in d.items() if v}
    if not ok: return d, None
    N = sum(v[0] for v in ok.values()); T = sum(v[1] for v in ok.values()); lo, hi = ci(max(N, 0))
    g = [G(n, t, V) for n, t in ok.values()]; se = np.std(g, ddof=1) / np.sqrt(len(g)) if len(g) > 1 else float("nan")
    return d, dict(N=N, T=T, G=G(N, T, V), lo=G(lo, T, V), hi=G(hi, T, V), reps=g, mean=np.mean(g), se=se)
def ratio(k, t, N, T):
    m = k + N; f = t / (t + T)
    plo = stats.beta.ppf(0.025, k, m - k + 1) if k > 0 else 0.0; phi = stats.beta.ppf(0.975, k + 1, m - k) if m - k > 0 else 1.0
    rr = lambda p: p / (1 - p) * (T / t) if p < 1 else float("inf")
    return (k / t) / (N / T), rr(plo), rr(phi), stats.binom.cdf(k, m, f)
L = [f"# FINAL SUMMARY — {datetime.datetime.now():%F %T}; estimator: net charge through the selectivity filter, first 20 ns discarded", ""]
def fmt(name, d, p, V):
    L.append(f"## {name} ({int(V*1000)} mV)")
    for r, v in d.items(): L.append(f"   {r:8s} " + (f"{v[0]:3d} / {v[1]:6.1f} ns  G {G(v[0], v[1], V):5.2f} pS" if v else "missing (not finished)"))
    if p: L.append(f"   POOLED {p['N']} / {p['T']:.1f} ns -> G {p['G']:.2f} pS (95 % {p['lo']:.2f}-{p['hi']:.2f}); replica mean {p['mean']:.2f} +- {p['se']:.2f} SE")
dW, pW = pool(["ANE_r1", "ANE_r2", "ANE_r4", "ANE_r5"], 0.565); fmt("WT (ANE)", dW, pW, 0.565)
dG, pG = pool(["GS_r1", "GS_r2", "GS_r3"], 0.565); fmt("G310S (GS)", dG, pG, 0.565)
if pW and pG:
    r, lo, hi, p = ratio(pG["N"], pG["T"], pW["N"], pW["T"]); n = sum(1 for v in dG.values() if v)
    if r < 0.5 and hi < 1: v = "permeation IMPAIRED (the Mkrtchyan hypothesis)"
    elif lo <= 1 <= hi and r > 0.5: v = "permeation PRESERVED -> loss of function through gating/PIP2 (Mosca 2022)"
    else: v = "neither case -> add replicas"
    L.append(f"   GS/WT = {r:.2f} (95 % {lo:.2f}-{hi:.2f}), one-sided p = {p:.3f}; G310S replicas {n}/3 -> rule: {v}" + ("" if n == 3 else "  [PROVISIONAL: not all G310S replicas have finished]"))
dC, pC = pool(["CR3_r1", "CR3_r2"], 0.565); fmt("CR3 (7CR3 closed pore)", dC, pC, 0.565)
if pC:
    n = sum(1 for v in dC.values() if v)
    L.append(f"   CR3 total net charge {pC['N']} ({n}/2 replicas) -> rule: " + ("specificity CONFIRMED (<= 1)" if pC["N"] <= 1 else "THE PROTOCOL ALSO CONDUCTS THROUGH THE CLOSED PORE -> stop and check") + ("" if n == 2 else "  [PROVISIONAL]"))
dA, pA = pool(["ANC_r1"], 0.565); fmt("ANC (G310 7.6 A)", dA, pA, 0.565)
dV, pV = pool(["V280_r1", "V280_r2", "V280_r3"], 0.282); fmt("V280 (I-V)", dV, pV, 0.282)
if pV and pW:
    # G(282)/G(565) = (rate282/0.282)/(rate565/0.565)
    r, lo, hi, p = ratio(pV["N"], pV["T"], pW["N"], pW["T"]); k = 0.565 / 0.282; n = sum(1 for v in dV.values() if v)
    gr, glo, ghi = r * k, lo * k, hi * k
    L.append(f"   G(282)/G(565) = {gr:.2f} (95 % {glo:.2f}-{ghi:.2f}); V280 replicas {n}/3 -> " + ("the 95 % interval includes 1: no deviation from ohmic behaviour detected" if glo <= 1 <= ghi else
                 "the 95 % interval excludes 1: deviation from ohmic behaviour -> discuss") + ("" if n == 3 else "  [PROVISIONAL]"))
open(f"{O}/FINAL_SUMMARY.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
