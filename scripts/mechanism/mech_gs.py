#!/usr/bin/env python3
# Mechanism of G310S: G310S vs wild type with identical definitions; first 20 ns excluded, 1 frame per ns.
# Topology: the system .gro (a .tpr renumbers the residues). Axis: xy centre of the filter (276-281) CA atoms, taken with the minimum image.
# z_rel = z - z_SF (centre of the filter CA atoms). Gate plane z_g = mean z_rel of the 310 CA atoms.
#   gate water    : OW, R < 0.5 nm, |z_rel - z_g| < 0.3 nm
#   cavity        : z_g + 0.3 < z_rel < -0.8 nm, R < 0.6 nm; number of OW and K+; fraction of frames with K+ present
#   below the gate: z_g - 0.8 < z_rel < z_g - 0.1, R < 0.8 nm; fraction of frames with at least one K+ (do ions reach the gate?)
#   xy distance of S310 OG (G310S) and of the 310 CA atoms from the axis; d1/d2 = distances between 310 CA atoms of opposite chains (narrow/wide)
# usage: mech_gs.py NAME GRO XTC   -> one tab-separated line per run
import sys, numpy as np, warnings; warnings.filterwarnings('ignore')
import MDAnalysis as mda
name, gro, xtc = sys.argv[1:4]
u = mda.Universe(gro, xtc)
sf = u.select_atoms('resid 276:281 and name CA'); g = u.select_atoms('resid 310 and name CA')
og = u.select_atoms('resid 310 and resname SER and name OG'); ow = u.select_atoms('resname SOL and name OW'); k = u.select_atoms('resname K and name K')
assert len(g) == 4 and len(sf) == 24, (len(g), len(sf))
dt = u.trajectory.dt; i0 = int(round((u.trajectory[0].time + 20000.0 - u.trajectory[0].time) / dt)); st = int(round(1000.0 / dt))
mic = lambda d, L: d - L * np.round(d / L)
R = []
for ts in u.trajectory[i0::st]:
    L = ts.dimensions[:3] / 10.0
    p = sf.positions / 10.0; p = p[0] + mic(p - p[0], L); c = p.mean(0)
    rel = lambda x: mic(x / 10.0 - c, L)
    gr = rel(g.positions); zg = gr[:, 2].mean(); car = np.hypot(gr[:, 0], gr[:, 1]).mean()
    da = np.linalg.norm(gr[0] - gr[2]); db = np.linalg.norm(gr[1] - gr[3]); d1, d2 = min(da, db), max(da, db)
    w = rel(ow.positions); wr = np.hypot(w[:, 0], w[:, 1]); wz = w[:, 2]
    kk = rel(k.positions); kr = np.hypot(kk[:, 0], kk[:, 1]); kz = kk[:, 2]
    gw = np.sum((wr < 0.5) & (np.abs(wz - zg) < 0.3))
    cav = lambda r, z, rr: (r < rr) & (z > zg + 0.3) & (z < -0.8)
    cw = np.sum(cav(wr, wz, 0.6)); ck = np.sum(cav(kr, kz, 0.6))
    bk = np.sum((kr < 0.8) & (kz > zg - 0.8) & (kz < zg - 0.1))
    ogr = np.hypot(*rel(og.positions)[:, :2].T).mean() if len(og) else np.nan
    R.append((ts.time, zg, d1 * 10, d2 * 10, car * 10, ogr * 10, gw, cw, ck, bk))
R = np.array(R)
m = lambda j: R[:, j].mean(); s = lambda j: R[:, j].std()
line = (f"{name}\tn={len(R)}\tt={R[0,0]/1000:.0f}-{R[-1,0]/1000:.0f}ns\td1={m(2):.2f}±{s(2):.2f}\td2={m(3):.2f}\t310CA_R={m(4):.2f}\tOG_R={m(5):.2f}\t"
        f"gate_water={m(6):.2f}±{s(6):.2f}\tcav_water={m(7):.1f}±{s(7):.1f}\tcav_K_mean={m(8):.3f}\tcav_K_frames={100*np.mean(R[:,8]>0):.1f}%\t"
        f"below_gate_K_frames={100*np.mean(R[:,9]>0):.1f}%\tzg={m(1):.2f}nm")
open(f"{sys.argv[4] if len(sys.argv) > 4 else '.'}/{name}.txt", 'w').write(line + '\n'); np.save(f"{sys.argv[4] if len(sys.argv) > 4 else '.'}/{name}.npy", R)
print(line)
