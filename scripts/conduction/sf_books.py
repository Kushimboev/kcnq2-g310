#!/usr/bin/env python3
"""Selectivity-filter bookkeeping (zone-based, radially restricted; no plane counting).
Zones (z relative to SF-O centre, nm): SF = -0.8..+0.45 & R<0.6 ; CAVdeep = gate+0.1 .. -1.2 & R<0.9 ; OUTfar = z > +1.2
For each ion: sequence of 'S' visits. For every SF visit record where the ion CAME FROM (last of CAVdeep/OUTfar before entering)
and where it WENT (first of CAVdeep/OUTfar reached after leaving). Visits still open at t=0 / t=end are flagged.
usage: sf_books.py gro xtc sfres gres ksel nksf label"""
import sys, numpy as np, MDAnalysis as mda, warnings; warnings.filterwarnings("ignore")
gro, xtc, sfres, gres, ksel, nksf, lab = sys.argv[1:8]; nksf = int(nksf)
u = mda.Universe(gro, xtc); K = u.select_atoms(ksel)
sf = u.select_atoms(f"resid {sfres} and name O"); ga = u.select_atoms(f"resid {gres} and name CA")
nK, nfr = K.n_atoms, len(u.trajectory)
Z = np.empty((nfr, nK), np.float32); R = np.empty((nfr, nK), np.float32); ZG = np.empty(nfr); T = np.empty(nfr)
for i, ts in enumerate(u.trajectory):
    L = ts.dimensions[:3]; s = sf.positions; s = s[0] + ((s - s[0]) - L * np.round((s - s[0]) / L)); c = s.mean(0)
    d = K.positions - c; d -= L * np.round(d / L); Z[i] = d[:, 2] / 10; R[i] = np.hypot(d[:, 0], d[:, 1]) / 10
    g = ga.positions - c; g -= L * np.round(g / L); ZG[i] = g[:, 2].mean() / 10; T[i] = ts.time / 1000
T -= T[0]; zg = ZG.mean()
S = (Z >= -0.8) & (Z <= 0.45) & (R < 0.6)
CAV = (Z > zg + 0.1) & (Z < -1.2) & (R < 0.9)
OUT = Z > 1.2
rows = []
for j in range(nK):
    s = S[:, j]
    if not s.any(): continue
    # SF visits = maximal runs of s allowing short gaps (<=20 frames = 0.2 ns) outside SF (flicker at the boundary)
    idx = np.flatnonzero(s); runs = []; st = idx[0]; pv = idx[0]
    for k in idx[1:]:
        if k - pv > 20: runs.append((st, pv)); st = k
        pv = k
    runs.append((st, pv))
    # merge runs separated by excursions that never reach CAV or OUT
    merged = [list(runs[0])]
    for a, b in runs[1:]:
        gap = slice(merged[-1][1] + 1, a)
        if not (CAV[gap, j].any() or OUT[gap, j].any()): merged[-1][1] = b
        else: merged.append([a, b])
    for a, b in merged:
        pre = np.flatnonzero(CAV[:a, j] | OUT[:a, j]); post = np.flatnonzero(CAV[b + 1:, j] | OUT[b + 1:, j])
        came = ("CAV" if CAV[pre[-1], j] else "OUT") if len(pre) else ("t0" if a == 0 else "?")
        went = ("CAV" if CAV[b + 1 + post[0], j] else "OUT") if len(post) else ("end" if b == nfr - 1 else "?")
        rows.append((j, j >= nK - nksf, T[a], T[b], came, went))
def cnt(f): return sum(1 for r in rows if f(r))
print(f"#### {lab}: {T[-1]:.1f} ns; SF visits {len(rows)}")
for tag, f in (("CAV->SF->OUT (conduction)", lambda r: r[4] == "CAV" and r[5] == "OUT"),
               ("OUT->SF->CAV (reverse)", lambda r: r[4] == "OUT" and r[5] == "CAV"),
               ("t0 ->SF->OUT (pre-loaded exits)", lambda r: r[4] in ("t0", "?") and r[5] == "OUT"),
               ("CAV->SF->end (entered from cavity, still inside)", lambda r: r[4] == "CAV" and r[5] in ("end", "?")),
               ("OUT->SF->end (entered from outside, still inside)", lambda r: r[4] == "OUT" and r[5] in ("end", "?")),
               ("OUT->SF->OUT (outside flicker)", lambda r: r[4] == "OUT" and r[5] == "OUT"),
               ("CAV->SF->CAV (bounced back)", lambda r: r[4] == "CAV" and r[5] == "CAV"),
               ("t0 ->SF->end (never left)", lambda r: r[4] in ("t0", "?") and r[5] in ("end", "?")),
               ("t0 ->SF->CAV", lambda r: r[4] in ("t0", "?") and r[5] == "CAV")):
    n_free = cnt(lambda r: f(r) and not r[1]); n_ksf = cnt(lambda r: f(r) and r[1])
    print(f"   {tag:50s} free {n_free:3d}  KSF {n_ksf:2d}")
ent_bot = cnt(lambda r: r[4] == "CAV"); ent_top = cnt(lambda r: r[4] == "OUT")
ex_top = cnt(lambda r: r[5] == "OUT"); ex_bot = cnt(lambda r: r[5] == "CAV")
print(f"   entries from cavity {ent_bot}, from outside {ent_top}; exits to outside {ex_top}, to cavity {ex_bot}")
print(f"   NET outward charge through SF (exits_top - entries_top) = {ex_top - ent_top}   | (entries_bottom - exits_bottom) = {ent_bot - ex_bot}")
occ = S.sum(1); print(f"   SF occupancy mean {occ.mean():.2f} (t<5ns {occ[T<5].mean():.2f}); occupancy histogram " + " ".join(f"{k}:{(occ==k).mean()*100:.0f}%" for k in range(0, 6) if (occ==k).any()))
for r in sorted(rows, key=lambda r: r[2]):
    if r[5] == "OUT" or r[4] == "OUT":
        print(f"      ion {r[0]:4d} {'KSF ' if r[1] else 'free'} in SF {r[2]:7.2f}-{r[3]:7.2f} ns  came {r[4]:4s} went {r[5]}")
