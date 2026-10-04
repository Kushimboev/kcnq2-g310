#!/usr/bin/env python3
"""Excursions of K+ ions towards the gate (Table S6), computed from the tracks written by gate_tracks.py.
An excursion is defined with absorbing boundaries:
  below (L): z - z_G310 <= -0.6 nm;  cavity (H): z - z_G310 >= +0.5 nm;  in between the ion is in the lumen (R < 0.6 nm).
  An excursion starts when an ion leaves L and ends when it returns to L ("returned") or reaches H ("crossed").
  It "reached the gate plane" if it came within 0.25 nm of the G310 CA plane (the plane used for the coordination analysis).
Output: GATE_COMMIT_VERIFY.txt (readable summary) and gate_commit.json (source of Table S6).   usage: python3 gate_commit_verify.py [TRACKS_DIR]"""
import os, sys, json, numpy as np
from scipy import stats
D = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))      # directory holding *.tracks.npz (output is written there too)
WT, GS = ["ANE_r1", "ANE_r2", "ANE_r4", "ANE_r5"], ["GS_r1", "GS_r2", "GS_r3"]
ZL, ZH, RL, HW, DT = -0.6, 0.5, 0.6, 0.25, 0.01          # nm, nm, nm, nm, ns
DEPTHS = [-0.5, -0.4, -0.3, -0.25, -0.2, -0.1, 0.0, 0.1, 0.25]

def excursions(r):
    z = np.load(f"{D}/{r}.tracks.npz"); A = z["rows"]; nfr = int(z["nframes"]); ex = []; og = []
    for k in np.unique(A[:, 1]):
        a = A[A[:, 1] == k]; a = a[np.argsort(a[:, 0])]
        t, dz, R, dog = a[:, 0], a[:, 2], a[:, 3], a[:, 4]
        st = np.where(dz <= ZL, 0, np.where(dz >= ZH, 2, 1)); st[(st == 1) & (R >= RL)] = np.where(dz[(st == 1) & (R >= RL)] < 0, 0, 2)
        cur = None; last = None
        for i in range(len(t)):
            gap = i > 0 and t[i] - t[i - 1] > 15
            if gap and cur is not None: cur["end"] = "lost"; ex.append(cur); cur = None
            if gap: last = None
            s = st[i]
            if s == 1:
                if cur is None: cur = dict(src=last, zmax=dz[i], zmin=dz[i], n=0, nplane=0, nog=0, t0=float(t[i]))
                cur["zmax"] = max(cur["zmax"], dz[i]); cur["zmin"] = min(cur["zmin"], dz[i]); cur["n"] += 1
                cur["nplane"] += abs(dz[i]) < HW; cur["nog"] += dog[i] < 0.35
            else:
                if cur is not None:
                    cur["end"] = s; cur["zmax"] = max(cur["zmax"], ZH) if s == 2 else cur["zmax"]; ex.append(cur); cur = None   # an ion that crossed has reached ZH (it may jump between two frames)
                last = s
        if cur is not None: cur["end"] = "open"; ex.append(cur)
        # K+ - Ser OG contact episodes (consecutive frames)
        c = dog < 0.35; i = 0
        while i < len(c):
            if c[i]:
                j = i
                while j + 1 < len(c) and c[j + 1] and t[j + 1] - t[j] <= 15: j += 1
                og.append((j - i + 1) * DT); i = j + 1
            else: i += 1
    return ex, np.array(og), nfr * DT

res = {}; lines = ["# gate_commit_verify.py: excursions of K+ towards the gate (first 20 ns of each run excluded; axis = centre of the filter carbonyl oxygens)",
                   f"# below z <= {ZL} nm, cavity z >= +{ZH} nm (relative to the G310 CA plane), lumen R < {RL} nm; gate plane |z| < {HW} nm\n"]
for lab, runs in (("WT", WT), ("GS", GS)):
    tot = dict(ns=0.0, att=0, plane=0, up=0, ret=0, down=0, dw=[], og=[], reach={d: [0, 0] for d in DEPTHS}, per=[])
    for r in runs:
        ex, og, ns = excursions(r)
        fromL = [e for e in ex if e["src"] == 0 and e["end"] in (0, 2)]
        pl = [e for e in fromL if e["nplane"] > 0]
        up = sum(e["end"] == 2 for e in pl); ret = len(pl) - up; upall = sum(e["end"] == 2 for e in fromL)
        down = sum(1 for e in ex if e["src"] == 2 and e["end"] == 0)
        dw = [e["nplane"] * DT for e in pl]
        for d in DEPTHS:
            rr = [e for e in fromL if e["zmax"] >= d]; tot["reach"][d][0] += len(rr); tot["reach"][d][1] += sum(e["end"] == 2 for e in rr)
        tot["ns"] += ns; tot["att"] += len(fromL); tot["plane"] += len(pl); tot["up"] += up; tot["ret"] += ret; tot["down"] += down; tot["dw"] += dw; tot["og"] += og.tolist()
        tot["per"].append(dict(run=r, ns=ns, excursions=len(fromL), plane=len(pl), crossed=up, returned=ret, crossed_all=upall, down=down,
                               dwell_mean_ps=float(np.mean(dw) * 1000) if dw else None, dwell_median_ps=float(np.median(dw) * 1000) if dw else None,
                               og_n=int(len(og)), og_median_ps=float(np.median(og) * 1000) if len(og) else None, og_mean_ps=float(og.mean() * 1000) if len(og) else None,
                               og_max_ps=float(og.max() * 1000) if len(og) else None))
        p = tot["per"][-1]
        lines.append(f"{r}: {ns:.0f} ns | excursions from below {len(fromL)} | reached the plane {len(pl)} | crossed {up}, returned {ret} | from the cavity downwards {down} | "
                     f"time in the plane: mean {p['dwell_mean_ps'] or 0:.0f} ps, median {p['dwell_median_ps'] or 0:.0f} ps"
                     + (f" | Ser-OG contacts: {len(og)} episodes, median {p['og_median_ps']:.0f}, mean {p['og_mean_ps']:.0f}, longest {p['og_max_ps']:.0f} ps" if len(og) else ""))
    dw = np.array(tot["dw"]); og = np.array(tot["og"])
    res[lab] = dict(ns=tot["ns"], excursions=tot["att"], plane=tot["plane"], crossed=tot["up"], returned=tot["ret"], down=tot["down"],
                    plane_per_us=tot["plane"] / tot["ns"] * 1000, crossed_per_us=tot["up"] / tot["ns"] * 1000, frac=tot["up"] / tot["plane"],
                    dwell_mean_ps=float(dw.mean() * 1000), dwell_median_ps=float(np.median(dw) * 1000), plane_time_pct=float(dw.sum() / tot["ns"] * 100),
                    dwell_rep_means=[p["dwell_mean_ps"] for p in tot["per"]],
                    og=dict(n=int(len(og)), median_ps=float(np.median(og) * 1000), mean_ps=float(og.mean() * 1000), p90_ps=float(np.percentile(og, 90) * 1000), max_ps=float(og.max() * 1000)) if len(og) else None,
                    reach={str(d): tot["reach"][d] for d in DEPTHS}, per=tot["per"])
    x = res[lab]
    lines.append(f"{lab} POOLED ({x['ns']:.0f} ns): reached the plane {x['plane']} = {x['plane_per_us']:.0f}/us | crossed {x['crossed']} ({x['crossed_per_us']:.1f}/us), returned {x['returned']} | "
                 f"fraction crossing {x['frac']:.3f} | time in the plane: mean {x['dwell_mean_ps']:.0f} ps, median {x['dwell_median_ps']:.0f} ps | plane occupied {x['plane_time_pct']:.2f} % of the time"
                 + (f" | Ser-OG contacts: median {x['og']['median_ps']:.0f} ps, mean {x['og']['mean_ps']:.0f} ps, 90 % < {x['og']['p90_ps']:.0f} ps, longest {x['og']['max_ps']:.0f} ps" if x["og"] else "") + "\n")
a, b = res["WT"], res["GS"]
odds, pf = stats.fisher_exact([[a["crossed"], a["returned"]], [b["crossed"], b["returned"]]])
pw = stats.ttest_ind(a["dwell_rep_means"], b["dwell_rep_means"], equal_var=False).pvalue
# ratio of the rates of arrival at the plane (conditional binomial, exposure = time)
N = a["plane"] + b["plane"]; p0 = b["ns"] / (a["ns"] + b["ns"]); pr = 2 * min(stats.binom.cdf(b["plane"], N, p0), stats.binom.sf(b["plane"] - 1, N, p0))
res["tests"] = dict(frac_ratio=b["frac"] / a["frac"], fisher_p=float(pf), dwell_welch_p=float(pw), arrival_ratio=b["plane_per_us"] / a["plane_per_us"], arrival_p=float(min(1, pr)))
lines.append(f"Fraction crossing, G310S/WT = {res['tests']['frac_ratio']:.2f} (Fisher p = {pf:.4f}) | rate of arrival at the plane, G310S/WT = {res['tests']['arrival_ratio']:.2f} (conditional binomial p = {min(1, pr):.3f}) | time in the plane, Welch p = {pw:.2f}\n")
lines.append("## Depth profile: excursions from below that reached z* / how many of them crossed")
lines.append("   z* (nm)   WT reached (/us) crossed fraction |  G310S reached (/us) crossed fraction")
for d in DEPTHS:
    wa, wb = a["reach"][str(d)]; ga, gb = b["reach"][str(d)]
    lines.append(f"   {d:+.2f}     {wa:4d} ({wa/a['ns']*1000:5.0f})   {wb:3d}   {wb/max(wa,1):.3f}  |  {ga:4d} ({ga/b['ns']*1000:5.0f})   {gb:3d}   {gb/max(ga,1):.3f}")
open(f"{D}/GATE_COMMIT_VERIFY.txt", "w").write("\n".join(lines) + "\n"); json.dump(res, open(f"{D}/gate_commit.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("\n".join(lines))
