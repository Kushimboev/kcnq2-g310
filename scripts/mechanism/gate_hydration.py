#!/usr/bin/env python3
# Coordination of K+ along the pore axis in the wild-type (ANE) and G310S (GS) production runs.
# Every frame (10 ps, first 20 ns excluded): for each K+ on the axis (R < 0.5 nm): number of water O (0.30 / 0.35 nm), S310 OG (0.35 nm) and other protein O (0.35 nm);
# z = z_ion - z_G310 (mean z of the G310 CA atoms; axis = xy centre of the filter CA atoms, as in mech_gs.py). Bulk K+: more than 3 nm (in z) from the lipid P plane (every 50th frame).
# usage: gate_hydration.py NAME GRO XTC OUTDIR   -> OUTDIR/NAME.hyd.npz + NAME.hyd.txt
import sys, numpy as np, warnings; warnings.filterwarnings('ignore')
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
name, gro, xtc, od = sys.argv[1:5]
u = mda.Universe(gro, xtc)
sf = u.select_atoms('resid 276:281 and name CA'); g = u.select_atoms('resid 310 and name CA')
K = u.select_atoms('resname K and name K'); OW = u.select_atoms('resname SOL and name OW')
OG = u.select_atoms('resid 310 and resname SER and name OG')
nprot = len(u.select_atoms('not (resname SOL K CL PA PC OL)'))
PO = u.atoms[:nprot].select_atoms('name O* and not name OG'); PP = u.select_atoms('resname PC and name P31 P')
assert len(g) == 4 and len(sf) == 24, (len(g), len(sf))
mic = lambda d, L: d - L * np.round(d / L)
rows, bulk = [], []
t0 = u.trajectory[0].time
import os; MAXF = int(os.environ.get('HYD_MAX', '0')); B0 = float(os.environ.get('HYD_B', '20000'))
for k, ts in enumerate(u.trajectory):
    if ts.time - t0 < B0 - 1e-6: continue
    if MAXF and k > MAXF: break
    L = ts.dimensions[:3]; box = ts.dimensions
    p = sf.positions; p = p[0] + mic(p - p[0], L); c = p.mean(0)
    gz = mic(g.positions - c, L)[:, 2].mean()
    kr = mic(K.positions - c, L); R = np.hypot(kr[:, 0], kr[:, 1])
    sel = np.where((R < 5.0) & (kr[:, 2] > gz - 12.0) & (kr[:, 2] < gz + 15.0))[0]
    if len(sel):
        pos = K.positions[sel]
        def cnt(B, cut):
            if len(B) == 0: return np.zeros(len(sel), int)
            pr = capped_distance(pos, B.positions, cut, box=box, return_distances=False)
            return np.bincount(pr[:, 0], minlength=len(sel)) if len(pr) else np.zeros(len(sel), int)
        w30, w35, og, po = cnt(OW, 3.0), cnt(OW, 3.5), cnt(OG, 3.5), cnt(PO, 3.5)
        for j, i in enumerate(sel):
            rows.append((ts.time, K.indices[i], (kr[i, 2] - gz) / 10.0, R[i] / 10.0, w30[j], w35[j], og[j], po[j]))
    if k % 50 == 0 and len(PP):
        zm = PP.positions[:, 2].mean(); dz = mic(K.positions[:, 2] - zm, L[2]); b = np.where(np.abs(dz) > 30.0)[0]
        if len(b):
            pr = capped_distance(K.positions[b], OW.positions, 3.5, box=box, return_distances=False)
            bulk.extend(np.bincount(pr[:, 0], minlength=len(b)).tolist())
A = np.array(rows, float); np.savez_compressed(f"{od}/{name}.hyd.npz", rows=A, bulk=np.array(bulk))
gate = A[np.abs(A[:, 2]) < 0.25] if len(A) else A
f = lambda a, j: f"{a[:, j].mean():.2f}±{a[:, j].std():.2f}" if len(a) else "nan"
line = (f"{name}\tframes_rows={len(A)}\tgate_plane_samples={len(gate)}\tnW35_gate={f(gate,5)}\tnW30_gate={f(gate,4)}\tnOG_gate={f(gate,6)}\tnProtO_gate={f(gate,7)}"
        f"\tbulk_nW35={np.mean(bulk):.2f}±{np.std(bulk):.2f} (n {len(bulk)})")
open(f"{od}/{name}.hyd.txt", "w").write(line + "\n"); print(line)
