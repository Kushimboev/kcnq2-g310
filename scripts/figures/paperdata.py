#!/usr/bin/env python3
"""Single source of the data behind the manuscript figures.
Events are re-read from the sf_books outputs (`*_b20.books.txt`) rather than typed in by hand, so that the figures cannot drift from the analysis.
Coordination: output of gate_hydration.py."""
import os, re, numpy as np
P = "${KCNQ2_ROOT}"
BOOKS = f"{P}/13_conduction/scripts/analysis_2026-09-23/out"
HYD = f"{P}/13_conduction/scripts/analysis_2026-09-28_paper/out"
E_CHARGE = 1.602176634e-19
WT = ["ANE_r1", "ANE_r2", "ANE_r4", "ANE_r5"]; GS = ["GS_r1", "GS_r2", "GS_r3"]
CR3 = ["CR3_r1", "CR3_r2"]; V280 = ["V280_r1", "V280_r2", "V280_r3"]
LABEL = {"ANE": "WT", "GS": "G310S", "CR3": "closed (7CR3)", "ANC": "G310 narrowed", "V280": "WT, 282 mV"}

def events(run):
    """List of (t_ns, +1/-1): 'went OUT' = +1 (outward), 'came OUT' = -1. In the analysis window t0 = 20 ns."""
    ev = []
    for ln in open(f"{BOOKS}/{run}_b20.books.txt"):
        m = re.match(r"\s+ion\s+(\d+)\s+(\S+)\s+in SF\s+([\d.]+)-\s*([\d.]+) ns\s+came (\S+)\s+went (\S+)", ln)
        if not m: continue
        _, kind, t0, t1, came, went = m.groups()
        if went == "OUT": ev.append((float(t1), +1))
        if came == "OUT": ev.append((float(t0), -1))
    return sorted(ev)

def cumulative(run, tmax=280.0):
    ev = events(run); t = [0.0]; y = [0.0]; c = 0
    for te, s in ev:
        t += [te, te]; y += [c, c + s]; c += s
    t.append(tmax); y.append(c)
    return np.array(t), np.array(y)

def net(run): return sum(s for _, s in events(run))

def G_pS(n, t_ns, V_mV=565.0):
    return n * E_CHARGE / (t_ns * 1e-9 * V_mV * 1e-3) * 1e12

def poisson_ci(n, cl=0.95):
    from scipy.stats import chi2
    a = (1 - cl) / 2
    lo = 0.0 if n == 0 else chi2.ppf(a, 2 * n) / 2
    return lo, chi2.ppf(1 - a, 2 * n + 2) / 2

def G_ci(n, t_ns, V_mV=565.0):
    lo, hi = poisson_ci(n); return G_pS(lo, t_ns, V_mV), G_pS(hi, t_ns, V_mV)

def mech(runs):
    """Output of mech_gs.py: gate water, d1, OG_R, cavity."""
    out = {}
    for r in runs:
        f = f"{BOOKS}/mech_0925/{r}.txt"
        if not os.path.exists(f): continue
        d = dict(kv.split("=", 1) for kv in open(f).read().strip().split("\t")[1:] if "=" in kv)
        g = lambda k: float(d[k].split("±")[0])
        out[r] = dict(gate_water=g("gate_water"), d1=g("d1"), OG_R=g("OG_R"), cav_water=g("cav_water"))
    return out

def hyd(runs):
    """gate_hydration.py: coordination at the gate plane (water, OG, other protein O) and the bulk value."""
    out = {}
    for r in runs:
        z = np.load(f"{HYD}/{r}.hyd.npz"); A = z["rows"]; g = A[np.abs(A[:, 2]) < 0.25]
        out[r] = dict(n=len(g), water=g[:, 5].mean(), OG=g[:, 6].mean(), protO=g[:, 7].mean(),
                      total=(g[:, 5] + g[:, 6] + g[:, 7]).mean(), bulk=z["bulk"].mean())
    return out

def hyd_profile():
    z = np.load(f"{HYD}/hyd_profile.npz")
    return z["z"], z["WT"], z["GS"]          # [replica, z, (n, water, OG, protO)]


# ---- charge-scaling comparison (Fig. 6) and KCNQ1 benchmark: first 20 ns excluded, as for all other measurements ----
ECC_PAIRS = [("KCNQ2 (8J01)", ["AN_r1", "AN_r2"], ["ANE_r1", "ANE_r2", "ANE_r3"]),      # matched pair: same structure and gate width, only the charges differ
             ("KCNQ1 (6V01)", ["KN_r1"], ["KNE_r1", "KNE_r2"])]
KCNQ1_FULL = ["KNE_r1", "KNE_r2"]

def book(run):
    """(net, ns, occ{k: %}) from {run}_b20.books.txt: net filter charge (exits_top - entries_top), window length, histogram of the filter occupancy."""
    t = open(f"{BOOKS}/{run}_b20.books.txt").read()
    net = int(re.search(r"\(exits_top - entries_top\) = *(-?\d+)", t).group(1))
    ns = float(re.search(r"#### \S+: *([\d.]+) ns", t).group(1))
    occ = {int(k): float(v) for k, v in re.findall(r"(\d+):(\d+)%", re.search(r"occupancy histogram ([^\n]+)", t).group(1))}
    return net, ns, occ

def pool(runs):
    """(net, ns): pooled net crossings and time, first 20 ns excluded."""
    return sum(book(r)[0] for r in runs), sum(book(r)[1] for r in runs)

def rate_ratio(n1, t1, n2, t2, cl=0.95):
    """Ratio of two Poisson rates (n1/t1)/(n2/t2): exact conditional-binomial 95 % interval and one-sided p (as in fig6_ecc.py)."""
    from scipy.stats import beta, binom
    n, N = n1, n1 + n2; f = t1 / t2
    lo = beta.ppf((1 - cl) / 2, n, N - n + 1) if n > 0 else 0.0
    hi = beta.ppf(1 - (1 - cl) / 2, n + 1, N - n) if n < N else 1.0
    g = lambda p: (p / (1 - p)) / f if p < 1 else float("inf")
    return (n1 / t1) / (n2 / t2), g(lo), g(hi), binom.cdf(n, N, f / (1 + f))

# ---- mixed tetramers (H1 = one S310 subunit, chain A; H2 = two adjacent, A+B), 2 x 300 ns; first 20 ns excluded ----
HYB = {"H1": ["H1_r1", "H1_r2"], "H2": ["H2_r1", "H2_r2"]}
HYBD = f"{P}/13_conduction/scripts/analysis_2026-09-29_hybrid"

def hyb_book(run):
    """(net, ns): net filter charge and window length of a mixed-tetramer run (from the book written by hyb_an.sh)."""
    t = open(f"{HYBD}/{run}_b20.books.txt").read()
    return (int(re.search(r"\(exits_top - entries_top\) = *(-?\d+)", t).group(1)),
            float(re.search(r"#### \S+: *([\d.]+) ns", t).group(1)))

def gate_water(run):
    """Gate-water series (1 frame per ns): mech_gs.py output for WT/G310S and for the mixed tetramers."""
    d = f"{BOOKS}/mech_0925" if run.startswith(("ANE", "GS")) else f"{HYBD}/mech"
    return np.load(f"{d}/{run}.npy")[:, 6]

def gw_stats(runs):
    """(replica means, mean, SE)"""
    v = np.array([gate_water(r).mean() for r in runs])
    return v, v.mean(), (v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else float("nan"))

def f_water(runs):
    """f = (W_WT - W_H) / (W_WT - W_GS) with its 95 % CI (delta method, replica SE), as in hyb_readout.py."""
    _, A, sA = gw_stats(WT); _, Gg, sG = gw_stats(GS); _, H, sH = gw_stats(runs)
    den = A - Gg; f = (A - H) / den
    se = np.sqrt((((H - Gg) / den**2) * sA)**2 + ((1 / den) * sH)**2 + (((A - H) / den**2) * sG)**2)
    return f, se, f - 1.96 * se, f + 1.96 * se

if __name__ == "__main__":
    for r in WT + GS + CR3 + V280: print(f"{r:9s} net={net(r):+d}")
    print("WT total", sum(net(r) for r in WT), "G =", round(G_pS(sum(net(r) for r in WT), 280 * 4), 2), "pS")
    print("GS total", sum(net(r) for r in GS), "G =", round(G_pS(sum(net(r) for r in GS), 280 * 3), 2), "pS")
    ne, te = pool([r for _, e, _ in ECC_PAIRS for r in e]); nf, tf = pool([r for _, _, f in ECC_PAIRS for r in f])
    print(f"ECC matched pairs (first 20 ns excluded): {ne}/{te:.1f} ns vs {nf}/{tf:.1f} ns -> ratio {rate_ratio(ne, te, nf, tf)[0]:.3f}  | KCNQ1 benchmark {pool(KCNQ1_FULL)}")
    for k, rs in HYB.items():
        n = sum(hyb_book(r)[0] for r in rs); t = sum(hyb_book(r)[1] for r in rs); f, se, lo, hi = f_water(rs)
        print(f"{k}: gate water f = {f:.2f} ({lo:.2f}-{hi:.2f}) | crossings {n}/{t:.0f} ns, {k}/WT = {rate_ratio(n, t, 17, 1120.0)[0]:.2f}")
