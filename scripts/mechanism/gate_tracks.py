#!/usr/bin/env python3
"""Positions of K+ ions near the gate, for the excursion analysis (independent of gate_hydration.py and its *.hyd.npz files).
The pore axis is the centre of the selectivity-filter carbonyl oxygens. Every K+ within 1.2 nm of the axis and within 1.6 nm
(along z) of the G310 CA plane is written for every frame after the first 20 ns, together with its distance to the nearest
Ser310 OG. The number of frames read is printed; 28001 are expected for a 300 ns run saved every 10 ps.
usage: gate_tracks.py NAME GRO XTC OUTDIR   ->  OUTDIR/NAME.tracks.npz (rows: t - t0 [ps], atom index, z - z_G310 [nm], R [nm], d(OG) [nm])"""
import sys, numpy as np, warnings; warnings.filterwarnings("ignore")
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
name, gro, xtc, od = sys.argv[1:5]
u = mda.Universe(gro, xtc)
sf = u.select_atoms("resid 276:281 and name O and not resname SOL"); g = u.select_atoms("resid 310 and name CA")
K = u.select_atoms("resname K and name K"); OG = u.select_atoms("resid 310 and resname SER and name OG")
assert len(sf) == 24 and len(g) == 4, (len(sf), len(g))
mic = lambda d, L: d - L * np.round(d / L)
rows = []; t0 = u.trajectory[0].time; nfr = 0
for ts in u.trajectory:
    if ts.time - t0 < 20000 - 1e-6: continue
    nfr += 1; L = ts.dimensions[:3]
    p = sf.positions; p = p[0] + mic(p - p[0], L); c = p.mean(0)
    gz = mic(g.positions - c, L)[:, 2].mean()
    d = mic(K.positions - c, L); R = np.hypot(d[:, 0], d[:, 1]); dz = d[:, 2] - gz
    sel = np.where((R < 12.0) & (np.abs(dz) < 16.0))[0]
    if not len(sel): continue
    dog = distance_array(K.positions[sel], OG.positions, box=ts.dimensions).min(1) if len(OG) else np.full(len(sel), 99.0)
    for j, i in enumerate(sel): rows.append((ts.time - t0, K.indices[i], dz[i] / 10, R[i] / 10, dog[j] / 10))
np.savez_compressed(f"{od}/{name}.tracks.npz", rows=np.array(rows, np.float32), nframes=nfr)
print(name, "frames", nfr, "rows", len(rows))
