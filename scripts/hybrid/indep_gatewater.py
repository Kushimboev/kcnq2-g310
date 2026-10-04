#!/usr/bin/env python3
# Independent check of the mixed-tetramer result: the same quantity recomputed along a different code path.
# Differences from mech_gs.py: axis = xy centre of the four G310 CA atoms (not the filter CA atoms); frames at 21, 23, ... ns (mech_gs: 20, 21, 22 ... ns);
# PBC: minimum image of every atom relative to the G310 centre. Same definition: OW, R < 0.5 nm, |z - z(G310 CA)| < 0.3 nm.
# Also reports the residue name at position 310 of every chain (composition of the mixed pore) and the CA distance between the SER chains (adjacent/diagonal).
# usage: indep_gatewater.py NAME GRO XTC
import sys, numpy as np, warnings; warnings.filterwarnings('ignore')
import MDAnalysis as mda
name, gro, xtc = sys.argv[1:4]
u = mda.Universe(gro, xtc)
g = u.select_atoms('resid 310 and name CA'); ow = u.select_atoms('resname SOL and name OW')
og = u.select_atoms('resid 310 and resname SER and name OG')
assert len(g) == 4
rn = [a.resname for a in g]
ser = [i for i, r in enumerate(rn) if r == 'SER']
dt = u.trajectory.dt; t0 = u.trajectory[0].time
i0 = int(round((21000.0) / dt)); st = int(round(2000.0 / dt))
mic = lambda d, L: d - L * np.round(d / L)
W = []; OGR = []; SS = []
for ts in u.trajectory[i0::st]:
    L = ts.dimensions[:3]
    p = g.positions; p = p[0] + mic(p - p[0], L); c = p.mean(0)
    w = mic(ow.positions - c, L)
    W.append(np.sum((np.hypot(w[:, 0], w[:, 1]) < 5.0) & (np.abs(w[:, 2]) < 3.0)))
    if len(og):
        o = mic(og.positions - c, L); OGR.append(np.hypot(o[:, 0], o[:, 1]).mean())
    if len(ser) == 2: SS.append(np.linalg.norm(mic(p[ser[0]] - p[ser[1]], L)))
W = np.array(W)
print(f"{name}\tframes={len(W)}\tt={u.trajectory[i0].time/1000:.0f}..ns\tres310={','.join(rn)}\tgate_water={W.mean():.2f}"
      f"\tOG_R={np.mean(OGR) if OGR else float('nan'):.2f}A\tSER-SER_CA={np.mean(SS) if SS else float('nan'):.2f}A")
