#!/usr/bin/env python3
"""Independent recount of ion crossings (does not use the code of the primary estimator).
usage: indep_count.py <gro> <xtc> <sf_resids> <gate_resid> <ksel> <nksf> <efield V/nm> <label>
Computes:
 (1) a continuous z coordinate (PBC unwrapped) for every K+, relative to the centre of the selectivity filter;
 (2) net crossings (Schmitt trigger, hysteresis 0.3 nm) of the gate plane, the lower filter boundary (-0.8), the filter centre (0) and the upper filter boundary (+0.45);
 (3) complete transits: the sequence IN (below the gate) -> CAV -> SF -> OUT (z > +1.2) observed within the window;
 (4) the starting zone of every ion (t = 0): IN/CAV/SF/OUT; (5) filter and cavity occupancy against time; (6) 50 ns blocks.
"""
import sys, numpy as np, MDAnalysis as mda, warnings
warnings.filterwarnings("ignore")
gro, xtc, sfres, gres, ksel, nksf, ef, lab = sys.argv[1:9]
nksf = int(nksf); ef = float(ef)
u = mda.Universe(gro, xtc)
K = u.select_atoms(ksel); sf = u.select_atoms(f"resid {sfres} and name O"); ga = u.select_atoms(f"resid {gres} and name CA")
assert sf.n_atoms == 24 and ga.n_atoms == 4, (sf.n_atoms, ga.n_atoms)
nK, nfr = K.n_atoms, len(u.trajectory)
Zr = np.empty((nfr, nK)); R = np.empty((nfr, nK)); ZG = np.empty(nfr); T = np.empty(nfr); LZ = np.empty(nfr); D1 = np.empty(nfr); D2 = np.empty(nfr)
for i, ts in enumerate(u.trajectory):
    L = ts.dimensions[:3]
    s = sf.positions; s = s[0] + ((s - s[0]) - L * np.round((s - s[0]) / L)); c = s.mean(0)
    d = K.positions - c; d -= L * np.round(d / L)
    Zr[i] = d[:, 2] / 10; R[i] = np.hypot(d[:, 0], d[:, 1]) / 10
    g = ga.positions - c; g -= L * np.round(g / L)
    ZG[i] = g[:, 2].mean() / 10
    D1[i] = np.linalg.norm(g[0, :2] - g[2, :2]) / 10; D2[i] = np.linalg.norm(g[1, :2] - g[3, :2]) / 10
    T[i] = ts.time / 1000; LZ[i] = L[2] / 10
T = T - T[0]
Lz = LZ.mean(); zg = ZG.mean()
dz = np.diff(Zr, axis=0); dz -= LZ[1:, None] * np.round(dz / LZ[1:, None])
Zu = np.vstack([Zr[:1], Zr[:1] + np.cumsum(dz, axis=0)])
free = np.arange(nK) < nK - nksf
def schmitt(zp, h):
    """net upward crossings of plane zp (and periodic images) per ion; returns events list (t, sign, ion)"""
    ev = []
    w = (Zu - zp) / Lz; H = h / Lz
    for j in range(nK):
        c = np.floor(w[0, j]); wj = w[:, j]
        # vectorised enough: loop only over frames where it can change
        for i in np.flatnonzero((wj >= c + 1 + H) | (wj <= c - H)) if False else range(1, nfr):
            p = wj[i]
            while p >= c + 1 + H:
                c += 1
                if R[i, j] < 1.0: ev.append((T[i], +1, j))
            while p <= c - H:
                c -= 1
                if R[i, j] < 1.0: ev.append((T[i], -1, j))
    return ev
planes = {"gate": zg, "SF_low(-0.8)": -0.8, "SF_mid(0)": 0.0, "SF_top(+0.45)": 0.45}
res = {}
for nm, zp in planes.items():
    ev = schmitt(zp, 0.3)
    fe = [e for e in ev if free[e[2]]]; ke = [e for e in ev if not free[e[2]]]
    res[nm] = (sum(e[1] for e in fe), sum(e[1] for e in ke), fe)
# zone at t=0
def zone(z, r, zgate):
    if z > 1.2: return "OUT"
    if -0.8 <= z <= 0.45 and r < 0.6: return "SF"
    if zgate < z < -0.8 and r < 0.9: return "CAV"
    if z <= zgate: return "IN"
    return "other"
# full permeation: last time in IN before reaching OUT, with CAV and SF visited in between
evs = []
for j in range(nK):
    z = Zu[:, j]; r = R[:, j]
    # use winding-aware local coordinate: position relative to nearest SF image
    zl = Zr[:, j]
    inz = (zl <= zg - 0.3) & (zl > -Lz / 2 + 0.01)
    cav = (zl > zg) & (zl < -0.8) & (r < 0.9)
    sfz = (zl >= -0.8) & (zl <= 0.45) & (r < 0.6)
    out = (zl > 1.2)
    st = None; seen_cav = seen_sf = False; start_zone = zone(zl[0], r[0], zg)
    last = "start"
    for i in range(nfr):
        if inz[i]: st = i; seen_cav = seen_sf = False; last = "IN"
        elif cav[i] and last in ("IN", "CAV", "start"): seen_cav = True; last = "CAV"
        elif sfz[i] and last in ("CAV", "SF", "start") : seen_sf = True; last = "SF"
        elif out[i]:
            if last == "SF":
                evs.append(dict(ion=j, t=T[i], full=(st is not None and seen_cav and seen_sf), start=start_zone, ksf=not free[j], t_in=(T[st] if st is not None else None)))
            last = "OUT"; st = None; seen_cav = seen_sf = False
# occupancy
sfocc = ((Zr >= -0.8) & (Zr <= 0.45) & (R < 0.6)).sum(1)
cavocc = ((Zr > zg) & (Zr < -0.8) & (R < 0.9)).sum(1)
V = ef; e = 1.602176634e-19; tt = T[-1]
def G(n): return n * e / (tt * 1e-9 * V) * 1e12
print(f"#### {lab}: {nfr} frames, {tt:.1f} ns, Lz {Lz:.3f} nm, V {V*1000:.0f} mV, gate plane z {zg:+.2f} nm, d1/d2 {np.minimum(D1,D2).mean()*10:.2f}/{np.maximum(D1,D2).mean()*10:.2f} A")
for nm in planes:
    nf, nk, fe = res[nm]
    print(f"  net crossings {nm:14s}: free {nf:+d}  KSF {nk:+d}  total {nf+nk:+d}   (G_free {G(nf):.2f} pS, G_all {G(nf+nk):.2f} pS)")
full = [x for x in evs if x['full'] and not x['ksf']]
part = [x for x in evs if not x['full'] and not x['ksf']]
print(f"  exits to OUT(z>+1.2) after SF: free {len([x for x in evs if not x['ksf']])} (FULL IN->CAV->SF->OUT in window: {len(full)}; pre-positioned/partial: {len(part)}), KSF {len([x for x in evs if x['ksf']])}")
for x in sorted(evs, key=lambda y: y['t']):
    print(f"     t_out {x['t']:7.2f} ns  ion {x['ion']:4d} {'KSF' if x['ksf'] else 'free'}  start={x['start']:5s} full={x['full']}  t_leftIN={x['t_in']}")
print(f"  G(full-only) = {G(len(full)):.2f} pS")
print(f"  SF occupancy: first 5 ns {sfocc[T<5].mean():.2f}, 5-50 ns {sfocc[(T>=5)&(T<50)].mean():.2f}, 50 ns-end {sfocc[T>=50].mean():.2f}; cavity: first 5 ns {cavocc[T<5].mean():.2f}, rest {cavocc[T>=5].mean():.2f}")
z0 = [(int(j), zone(Zr[0, j], R[0, j], zg), round(float(Zr[0, j]), 2)) for j in range(nK) if zone(Zr[0, j], R[0, j], zg) in ("SF", "CAV")]
print(f"  ions in SF/CAV at t=0: {z0}")
fe = res["SF_mid(0)"][2]
blocks = np.arange(0, tt + 50, 50)
cnt = [sum(e[1] for e in fe if b <= e[0] < b + 50) for b in blocks[:-1]]
print(f"  SF_mid free net per 50-ns block: {cnt}")
fe = res["gate"][2]
print(f"  gate  free net per 50-ns block: {[sum(e[1] for e in fe if b <= e[0] < b + 50) for b in blocks[:-1]]}")
