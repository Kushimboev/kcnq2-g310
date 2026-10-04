#!/usr/bin/env python3
"""Charge-scaling (ECC) comparison and KCNQ1 benchmark, first 20 ns excluded. Primary estimator: net charge through the selectivity filter (sf_books.py, paperdata.book);
independent check: indep_count.py (net crossings of the lower filter boundary, z = -0.8 nm, and of the gate plane). Output: ecc_b20_summary.txt."""
import sys, re
sys.path.insert(0, "${KCNQ2_ROOT}/10_figures/scripts")
import paperdata as D
ECC = ["AN_r1", "AN_r2", "KN_r1"]; K10 = ["K10_r1"]; FULL = ["ANE_r1", "ANE_r2", "ANE_r3", "KNE_r1", "KNE_r2"]; A9 = ["A9_r1"]
def indep(run):
    t = open(f"{D.BOOKS}/{run}_b20.indep.txt").read()
    g = int(re.search(r"net crossings gate\s*: free [+-]?\d+\s+KSF [+-]?\d+\s+total ([+-]?\d+)", t).group(1))
    lo = int(re.search(r"net crossings SF_low\(-0\.8\)\s*: free [+-]?\d+\s+KSF [+-]?\d+\s+total ([+-]?\d+)", t).group(1))
    return g, lo
L = []
def out(s=""): L.append(s); print(s)
out("run       ns    main(net SF)  indep SF_low  indep gate   (first 20 ns discarded)")
for r in ECC + K10 + A9 + FULL:
    n, ns, _ = D.book(r); g, lo = indep(r); out(f"{r:8s} {ns:6.1f}   {n:+3d}          {lo:+3d}           {g:+3d}")
def block(title, e, f):
    ne, te = D.pool(e); nf, tf = D.pool(f); out(f"\n{title}: ECC {ne}/{te:.1f} ns vs full {nf}/{tf:.1f} ns")
    for lab, get in (("main (net SF charge)", lambda r: D.book(r)[0]), ("indep SF_low", lambda r: indep(r)[1]), ("indep gate", lambda r: indep(r)[0])):
        a = sum(get(r) for r in e); b = sum(get(r) for r in f)
        if a < 0: a = 0
        rr, lo, hi, p = D.rate_ratio(a, te, b, tf)
        out(f"   {lab:22s}: {a}/{te:.1f} vs {b}/{tf:.1f} -> ratio {rr:.3f} ({lo:.3f}-{hi:.3f}), p = {p:.2e}, suppression x{1/rr:.1f}" if rr > 0 else
            f"   {lab:22s}: {a}/{te:.1f} vs {b}/{tf:.1f} -> ratio 0 (0-{hi:.3f}), p = {p:.2e}")
block("MATCHED PAIRS (as in the paper, Fig 6A): AN r1/r2 + KN r1 vs ANE r1/r2/r3 + KNE r1/r2", ECC, FULL)
block("all ECC runs (K10 r1 added)", ECC + K10, FULL)
block("sensitivity: ANE r3 (same starting frame as r1) excluded", ECC, [r for r in FULL if r != "ANE_r3"])
block("KCNQ2 only: AN r1/r2 vs ANE r1/r2/r3", ["AN_r1", "AN_r2"], ["ANE_r1", "ANE_r2", "ANE_r3"])
block("KCNQ1 only: KN r1 vs KNE r1/r2", ["KN_r1"], ["KNE_r1", "KNE_r2"])
nk, tk = D.pool(["KNE_r1", "KNE_r2"]); lo_, hi_ = D.G_ci(nk, tk)
out(f"\nKCNQ1 benchmark (KNE r1+r2, first 20 ns discarded): {nk} net crossings / {tk:.1f} ns -> G = {D.G_pS(nk, tk):.2f} pS (95 % {lo_:.2f}-{hi_:.2f})")
nw, tw = D.pool(["ANE_r1", "ANE_r2", "ANE_r3"]); na, ta = D.pool(A9)
out(f"panel C KCNQ2: A9 {na}/{ta:.1f} ns = {na/ta*1000:.1f}/us vs ANE {nw}/{tw:.1f} ns = {nw/tw*1000:.1f}/us")
open("ecc_b20_summary.txt", "w").write("\n".join(L) + "\n")
