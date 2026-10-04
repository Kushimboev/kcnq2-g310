#!/usr/bin/env python3
# First coordination shell of K+ at the gate plane: wild type (4 replicas) vs G310S (3 replicas).
# Source: gate_hydration.py npz (every frame, every on-axis ion: z_rel, R, nW3.0, nW3.5, nOG3.5, nProtO3.5).
# Gate plane: |z - z_G310| < 0.25 nm. Replica means -> Welch t test.
import numpy as np, glob, os
from scipy import stats
D = os.path.dirname(os.path.abspath(__file__)) + "/out"
WT = ["ANE_r1", "ANE_r2", "ANE_r4", "ANE_r5"]; GS = ["GS_r1", "GS_r2", "GS_r3"]
def load(n):
    z = np.load(f"{D}/{n}.hyd.npz"); A = z["rows"]; b = z["bulk"]
    g = A[np.abs(A[:, 2]) < 0.25]
    return A, g, b
out = []
tab = {}
for n in WT + GS:
    A, g, b = load(n)
    tab[n] = dict(n_gate=len(g), W35=g[:, 5].mean(), W30=g[:, 4].mean(), OG=g[:, 6].mean(), PO=g[:, 7].mean(),
                  tot=(g[:, 5] + g[:, 6] + g[:, 7]).mean(), bulk=b.mean(), nb=len(b))
    out.append(f"{n}\tn_gate={len(g)}\tW3.5={tab[n]['W35']:.2f}\tOG={tab[n]['OG']:.2f}\tprotO={tab[n]['PO']:.2f}\tTOTAL={tab[n]['tot']:.2f}\tbulk_W3.5={b.mean():.2f}")
def cmp(key, lab):
    a = np.array([tab[n][key] for n in WT]); c = np.array([tab[n][key] for n in GS])
    t, p = stats.ttest_ind(a, c, equal_var=False)
    return (f"{lab:22s} WT {a.mean():.2f} ± {a.std(ddof=1)/np.sqrt(len(a)):.2f} (n{len(a)})   "
            f"GS {c.mean():.2f} ± {c.std(ddof=1)/np.sqrt(len(c)):.2f} (n{len(c)})   Welch p = {p:.4f}")
out.append("")
for k, l in (("W35", "water O (< 3.5 A)"), ("W30", "water O (< 3.0 A)"), ("OG", "S310 OG (< 3.5 A)"),
             ("PO", "other protein O"), ("tot", "TOTAL coordination"), ("bulk", "bulk water (control)")):
    out.append(cmp(k, l))
# z profile: coordination and number of samples per bin
zb = np.arange(-1.2, 1.21, 0.1); prof = {}
for grp, runs in (("WT", WT), ("GS", GS)):
    P = []
    for n in runs:
        A, _, _ = load(n); r = []
        for z0 in zb:
            m = (A[:, 2] >= z0 - 0.05) & (A[:, 2] < z0 + 0.05)
            r.append([m.sum(), A[m, 5].mean() if m.sum() else np.nan, A[m, 6].mean() if m.sum() else np.nan,
                      A[m, 7].mean() if m.sum() else np.nan])
        P.append(r)
    prof[grp] = np.array(P, float)
np.savez(f"{D}/hyd_profile.npz", z=zb, WT=prof["WT"], GS=prof["GS"])
out.append("\nz profile (nm, relative to the G310 CA plane; columns: WT water | WT OG | GS water | GS OG | WT samples | GS samples)")
for i, z0 in enumerate(zb):
    w, g = prof["WT"][:, i], prof["GS"][:, i]
    out.append(f"  z {z0:+.1f}: {np.nanmean(w[:,1]):5.2f} | {np.nanmean(w[:,2]):4.2f} | {np.nanmean(g[:,1]):5.2f} | "
               f"{np.nanmean(g[:,2]):4.2f} | {int(np.nansum(w[:,0])):6d} | {int(np.nansum(g[:,0])):6d}")
open(f"{D}/HYD_SUMMARY.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out[:14]))
