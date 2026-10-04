#!/usr/bin/env python3
"""z(t) traces of the K+ ions along the pore axis (relative to the centre of the filter oxygens) and the gate plane z_g, for Fig. 3.
usage: ion_traces.py NAME GRO XTC OUTDIR [KSEL]   (KSEL: the ion selection; the ions are named POT in the reduced wild-type .gro and K in the G310S systems)"""
import sys, numpy as np, warnings; warnings.filterwarnings("ignore")
import MDAnalysis as mda
name, gro, xtc, od = sys.argv[1:5]
ksel = sys.argv[5] if len(sys.argv) > 5 else "resname K and name K"
u = mda.Universe(gro, xtc)
sf = u.select_atoms("resid 276:281 and name O"); ga = u.select_atoms("resid 310 and name CA")
K = u.select_atoms(ksel)
assert len(K) > 0, f"0 ion: {ksel}"
assert len(sf) == 24 and len(ga) == 4
n = len(u.trajectory); Z = np.empty((n, len(K)), np.float32); R = np.empty((n, len(K)), np.float32)
T = np.empty(n); ZG = np.empty(n)
for i, ts in enumerate(u.trajectory):
    L = ts.dimensions[:3]; s = sf.positions; s = s[0] + ((s - s[0]) - L * np.round((s - s[0]) / L)); c = s.mean(0)
    d = K.positions - c; d -= L * np.round(d / L)
    Z[i] = d[:, 2] / 10; R[i] = np.hypot(d[:, 0], d[:, 1]) / 10
    g = ga.positions - c; g -= L * np.round(g / L); ZG[i] = g[:, 2].mean() / 10; T[i] = ts.time / 1000
np.savez_compressed(f"{od}/{name}.trace.npz", t=T - T[0], z=Z, r=R, zg=ZG, ids=K.indices)
print(name, "frames", n, "ions", len(K), "zg", round(ZG.mean(), 2))
