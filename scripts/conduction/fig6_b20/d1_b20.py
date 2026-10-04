#!/usr/bin/env python3
"""Ca-Ca diagonals of G310 (KCNQ1: G345): the x axis of Fig. 6C.
In every frame the six 3D distances between the four CA atoms (made whole relative to the first atom) are sorted; the two largest are the diagonals:
d1 = the shorter ("narrow") diagonal, d2 = the longer ("wide") one. The "G310 diagonal" quoted in the paper is the time average of d1.
Reported for the whole trajectory and with the first 20 ns excluded (t >= 20 ns).
Output: d1_b20.json with the number of frames and the sums for each run (groups are pooled as frame-weighted means)."""
import numpy as np, MDAnalysis as mda, warnings, sys, json
warnings.filterwarnings("ignore")
P = "${KCNQ2_ROOT}"; C = f"{P}/13_conduction"; S = f"{C}/min_trajectories_scratch_2026-09-22/s645628ac"
GA = f"{C}/AW/sub_analysis_min_A.gro"; GK = f"{S}/KNE_min_ref.gro"
RUNS = {"AN_r1": (GA, f"{C}/AN/interim/prod_AN_r1_min.xtc", 310), "AN_r2": (GA, f"{C}/AN/interim/prod_AN_r2_min.xtc", 310),
        "A9_r1": (GA, f"{S}/A9_r1_min.xtc", 310), "ANE_r1": (GA, f"{S}/ANE_r1_min.xtc", 310), "ANE_r2": (GA, f"{S}/ANE_r2_min.xtc", 310),
        "ANE_r3": (GA, f"{S}/ANE_r3_min.xtc", 310), "KN_r1": (GK, f"{S}/KN_r1_min.xtc", 345), "K10_r1": (GK, f"{S}/K10_r1_min.xtc", 345),
        "KNE_r1": (GK, f"{S}/KNE_r1_min.xtc", 345), "KNE_r2": (GK, f"{S}/KNE_r2_min.xtc", 345)}
want = sys.argv[1:] or list(RUNS)
def whole(p, L):
    v = p - p[0]; v -= L * np.round(v / L); return p[0] + v
res = {}
for nm in want:
    gro, xtc, gres = RUNS[nm]
    u = mda.Universe(gro, xtc); g = u.select_atoms(f"resid {gres} and name CA"); assert g.n_atoms == 4, nm
    d1, d2, T = [], [], []
    for ts in u.trajectory:
        L = ts.dimensions[:3]; p = whole(g.positions, L)
        d = sorted(np.linalg.norm(p[i] - p[j]) for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)))
        d1.append(d[-2]); d2.append(d[-1]); T.append(ts.time / 1000)
    d1, d2, T = np.array(d1), np.array(d2), np.array(T); T -= T[0]; m = T >= 20.0
    nar, wid = d1, d2          # d1 = shorter (narrow), d2 = longer (wide) diagonal, per frame
    res[nm] = dict(nfull=int(len(T)), nb20=int(m.sum()), nar_full=float(nar.mean()), nar_b20=float(nar[m].mean()), wid_full=float(wid.mean()),
                   wid_b20=float(wid[m].mean()), sd_nar_b20=float(nar[m].std()), sum_nar_b20=float(nar[m].sum()), sum_wid_b20=float(wid[m].sum()),
                   sum_nar_full=float(nar.sum()), sum_wid_full=float(wid.sum()))
    r = res[nm]
    print(f"{nm:7s} frames {r['nfull']:6d}/{r['nb20']:6d}  NARROW full {r['nar_full']:.3f} b20 {r['nar_b20']:.3f} (sd {r['sd_nar_b20']:.2f})  WIDE full {r['wid_full']:.3f} b20 {r['wid_b20']:.3f}", flush=True)
import os
old = json.load(open("d1_b20.json")) if os.path.exists("d1_b20.json") else {}
old.update(res); json.dump(old, open("d1_b20.json", "w"), indent=1)          # merges the results of several calls
