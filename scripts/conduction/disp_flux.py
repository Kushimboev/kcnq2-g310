#!/usr/bin/env python3
"""Independent check of the conduction counts (logic independent of sf_books.py).
(1) Displacement estimator: N_net = sum over ions and frame pairs of [signed length of the ion path inside the segment z in [zlo, zhi]] / (zhi - zlo)
    (R < rmax in both frames; |dz| > 1 nm is a PBC jump and is skipped). Ions already inside the segment contribute fractionally.
    Filter segment [-0.8, +0.45] nm (relative to the centre of the filter oxygens); gate segment [zg-0.3, zg+0.3] nm (R < 0.6).
(2) Complete transits (state machine): +1 if an ion enters the filter segment from below (z < zlo) and leaves through the top (z > zhi) without interruption; -1 for the reverse.
(3) Potential: V = E_z * <Lz> (E_z is an argument, taken from the mdp).
usage: disp_flux.py gro xtc ksel label Ez[V/nm] [skip_ns=20]"""
import sys, numpy as np, MDAnalysis as mda, warnings
warnings.filterwarnings("ignore")
gro, xtc, ksel, lab, Ez = sys.argv[1:6]; Ez = float(Ez)
skip = float(sys.argv[6]) if len(sys.argv) > 6 else 20.0
u = mda.Universe(gro, xtc); K = u.select_atoms(ksel)
sf = u.select_atoms("resid 276:281 and name O"); ga = u.select_atoms("resid 310 and name CA")
assert len(sf) == 24 and len(ga) == 4, (len(sf), len(ga))
nfr, nK = len(u.trajectory), K.n_atoms
Z = np.empty((nfr, nK), np.float32); R = np.empty((nfr, nK), np.float32); ZG = np.empty(nfr); T = np.empty(nfr); LZ = np.empty(nfr)
for i, ts in enumerate(u.trajectory):
    L = ts.dimensions[:3]; s = sf.positions; s = s[0] + ((s - s[0]) - L * np.round((s - s[0]) / L)); c = s.mean(0)
    d = K.positions - c; d -= L * np.round(d / L); Z[i] = d[:, 2] / 10; R[i] = np.hypot(d[:, 0], d[:, 1]) / 10
    g = ga.positions - c; g -= L * np.round(g / L); ZG[i] = g[:, 2].mean() / 10; T[i] = ts.time / 1000; LZ[i] = L[2] / 10
T -= T[0]; w = T >= skip - 1e-6; Zw, Rw, Tw = Z[w], R[w], T[w]; tlen = Tw[-1] - Tw[0]; zg = ZG[w].mean()
def seg_flux(zlo, zhi, rmax):
    z0, z1, r0, r1 = Zw[:-1], Zw[1:], Rw[:-1], Rw[1:]
    ok = (np.abs(z1 - z0) < 1.0) & (r0 < rmax) & (r1 < rmax)
    lo, hi = np.minimum(z0, z1), np.maximum(z0, z1)
    ov = np.clip(hi, zlo, zhi) - np.clip(lo, zlo, zhi)
    return float((np.sign(z1 - z0) * ov * ok).sum() / (zhi - zlo))
def traversals(zlo, zhi, rmax):
    # B = cap of the cavity below the filter [zlo-0.4, zlo), R<1.0 ; T = vestibule above the filter (zhi, zhi+0.6], R<1.5 ; S = segment, R<rmax.
    # B -> (S only) -> T = +1 ; T -> (S only) -> B = -1 ; the state is reset if the ion leaves elsewhere (bulk, PBC jump).
    up = dn = up_ksf = 0
    for j in range(nK):
        z, r = Zw[:, j], Rw[:, j]; state = None
        for k in range(len(z)):
            zk, rk = z[k], r[k]
            if zlo - 0.4 <= zk < zlo and rk < 1.0:
                if state == 'T': dn += 1
                state = 'B'
            elif zhi < zk <= zhi + 0.6 and rk < 1.5:
                if state == 'B': up += 1; up_ksf += (j >= nK - 3)
                state = 'T'
            elif zlo <= zk <= zhi and rk < rmax: pass
            else: state = None
    return up, dn, up_ksf
Nsf = seg_flux(-0.8, 0.45, 0.8); Ngate = seg_flux(zg - 0.3, zg + 0.3, 0.6)
up, dn, upk = traversals(-0.8, 0.45, 0.8)
V = Ez * LZ[w].mean(); G = lambda n: n * 1.602176634e-19 / (tlen * 1e-9 * V) * 1e12
print(f"{lab:9s} window {tlen:6.1f} ns  V {V*1000:5.1f} mV (Lz {LZ[w].mean():.3f})  zg {zg:+.2f}  | SF displacement N {Nsf:+6.2f} (G {G(Nsf):5.2f} pS)"
      f"  gate N {Ngate:+6.2f}  | complete transits up {up} (KSF {upk}) down {dn} -> net {up-dn} (G {G(up-dn):5.2f} pS)")
