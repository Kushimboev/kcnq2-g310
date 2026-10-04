#!/usr/bin/env python3
"""Independent, visit-based count of gate crossings, from the per-ion records of the coordination analysis
(gate_hydration.py -> *.hyd.npz: for every frame and every on-axis K+: t, atom index, z - z_G310, R, water O, Ser OG, protein O).
  (1) visits to the gate plane (|z| < HW nm; gaps of up to GAP frames are bridged) and their duration;
  (2) outcome of a visit that started from below: continued upwards (lo->hi) or went back (lo->lo);
  (3) duration of K+ - Ser OG contact episodes in G310S (at least one OG within 0.35 nm).
Main definition: HW = 0.25 nm (the plane of gate_hydration.py), GAP = 3 frames (30 ps). Sensitivity: HW 0.15/0.25/0.35 x GAP 1/3/10.
This uses a different axis definition (filter CA atoms), a radial pre-filter (R < 0.5 nm) and different code from gate_commit_verify.py.
usage: python3 gate_commit.py [HYD_NPZ_DIR] > GATE_COMMIT_SUMMARY.txt"""
import os, sys, numpy as np
from scipy import stats
D = (sys.argv[1].rstrip("/") + "/") if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)) + "/../analysis_2026-09-28_paper/out/"
WT, GS = ["ANE_r1", "ANE_r2", "ANE_r4", "ANE_r5"], ["GS_r1", "GS_r2", "GS_r3"]; TNS = dict(WT=1120.0, GS=840.0)

def segs(t, mask, gap):
    out, cur = [], None
    for i, (ti, m) in enumerate(zip(t, mask)):
        if not m: continue
        if cur is None: cur = [i, i]
        elif ti - t[cur[1]] <= gap * 10 + 1e-6: cur[1] = i
        else: out.append(tuple(cur)); cur = [i, i]
    if cur is not None: out.append(tuple(cur))
    return out

def run(r, hw, gap):
    A = np.load(D + r + ".hyd.npz")["rows"]; dw, up, ret, og = [], 0, 0, []
    for k in np.unique(A[:, 1]):
        a = A[A[:, 1] == k]; a = a[np.argsort(a[:, 0])]; t, z = a[:, 0], a[:, 2]
        for s, e in segs(t, np.abs(z) < hw, gap):
            dw.append((t[e] - t[s]) / 1000 + 0.01)
            if s > 0 and t[s] - t[s - 1] <= 100 and z[s - 1] < 0 and e + 1 < len(z) and t[e + 1] - t[e] <= 100:
                if z[e + 1] > 0: up += 1
                else: ret += 1
        for s, e in segs(t, a[:, 6] >= 1, gap): og.append((t[e] - t[s]) / 1000 + 0.01)
    return np.array(dw), up, ret, np.array(og)

print("# gate_commit.py: visits to the gate plane and the fraction that continue upwards (first 20 ns excluded; WT 4 x 280 ns, G310S 3 x 280 ns)\n")
print("## Main definition: |z - z_G310| < 0.25 nm, gaps <= 30 ps")
tot = {}
for lab, rs in (("WT", WT), ("GS", GS)):
    U = Rt = 0; DW = []; per = []
    for r in rs:
        dw, up, ret, og = run(r, 0.25, 3); U += up; Rt += ret; DW += dw.tolist(); per.append(dw.mean() * 1000)
        s = f"{r}: visits {len(dw)} | duration: mean {dw.mean()*1000:.0f} ps, median {np.median(dw)*1000:.0f} ps, longest {dw.max()*1000:.0f} ps | from below: upwards {up}, back {ret}"
        if len(og): s += f" | Ser-OG contacts: {len(og)} episodes, mean {og.mean()*1000:.0f} ps, median {np.median(og)*1000:.0f} ps, 90 % < {np.percentile(og,90)*1000:.0f} ps, longest {og.max()*1000:.0f} ps"
        print(s)
    DW = np.array(DW); tot[lab] = (U, Rt, per)
    print(f"{lab} POOLED: visits {len(DW)} = {len(DW)/TNS[lab]*1000:.0f}/us | mean duration {DW.mean()*1000:.0f} ps | plane occupied {DW.sum()/TNS[lab]*100:.2f} % of the time | "
          f"arrived from below and continued upwards {U}/{U+Rt} = {U/(U+Rt):.3f} ({U/TNS[lab]*1000:.1f}/us)\n")
(u1, r1, p1), (u2, r2, p2) = tot["WT"], tot["GS"]
print(f"Ratio of the fractions, G310S/WT = {(u2/(u2+r2))/(u1/(u1+r1)):.2f}; Fisher p = {stats.fisher_exact([[u1, r1], [u2, r2]])[1]:.4f}; duration (replica means) Welch p = {stats.ttest_ind(p1, p2, equal_var=False).pvalue:.2f}\n")
print("## Sensitivity (HW nm, GAP frames): WT upwards/from below  duration mean/median ps | G310S ... | ratio of fractions, Fisher p | duration Welch p")
for hw in (0.15, 0.25, 0.35):
    for gap in (1, 3, 10):
        o = {}
        for lab, rs in (("WT", WT), ("GS", GS)):
            U = Rt = 0; DW = []; per = []
            for r in rs:
                dw, up, ret, _ = run(r, hw, gap); U += up; Rt += ret; DW += dw.tolist(); per.append(dw.mean() * 1000)
            o[lab] = (U, Rt, np.mean(DW) * 1000, np.median(DW) * 1000, per)
        a, b = o["WT"], o["GS"]
        print(f"HW {hw:.2f} GAP {gap:2d}: WT {a[0]}/{a[0]+a[1]} = {a[0]/(a[0]+a[1]):.3f}  {a[2]:.0f}/{a[3]:.0f} | G310S {b[0]}/{b[0]+b[1]} = {b[0]/(b[0]+b[1]):.3f}  {b[2]:.0f}/{b[3]:.0f} | "
              f"{(b[0]/(b[0]+b[1]))/(a[0]/(a[0]+a[1])):.2f}, p {stats.fisher_exact([[a[0], a[1]], [b[0], b[1]]])[1]:.3f} | p {stats.ttest_ind(a[4], b[4], equal_var=False).pvalue:.2f}")
print("\n## Note\n- The number of upward crossings does not depend on the definition; the fraction of visits that continue upwards does,"
      "\n  because a wider slab counts more visits. The excursion analysis (gate_commit_verify.py) is the one reported in Table S6.")
