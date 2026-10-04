#!/usr/bin/env python3
"""Fig. S1, independence of the estimators. A: three estimators for every run. B: G310S/WT ratio by estimator.
C: the potential actually realised in each run (V = E*<Lz>). Data: output of the independent estimators (disp_flux.py) and of sf_books.py."""
import sys, os, re, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../style")
import paperdata as D
from palette import C, MM, use_paper_style
use_paper_style()
V = "${KCNQ2_ROOT}/13_conduction/scripts/analysis_2026-09-26_verify/verify_0926.txt"
rows = {}
for ln in open(V):
    m = re.match(r"(\S+)\s+window\s+([\d.]+) ns\s+V\s+([\d.]+) mV.*SF displacement N\s+([+-][\d.]+).*gate N\s+([+-][\d.]+)", ln)
    if m: rows[m.group(1)] = dict(t=float(m.group(2)), V=float(m.group(3)), disp=float(m.group(4)), gate=float(m.group(5)))
ORDER = D.WT + D.GS + D.CR3 + ["ANC_r1"] + D.V280
LAB = {"ANE": "WT", "GS": "G310S", "CR3": "closed", "ANC": "G310 7.6 $\\AA$", "V280": "282 mV"}
fig = plt.figure(figsize=(MM.TWO_COL, 2.8))
gs = fig.add_gridspec(1, 3, width_ratios=[1.55, 0.85, 0.95], wspace=0.42, left=0.065, right=0.985, bottom=0.30, top=0.86)

# --- A: three estimators ---
ax = fig.add_subplot(gs[0])
x = np.arange(len(ORDER)); w = 0.27
sf = [D.net(r) for r in ORDER]; ga = [rows[r]["gate"] for r in ORDER]; dp = [rows[r]["disp"] for r in ORDER]
ax.bar(x - w, sf, w, color="#1B7EBD", label="net filter charge (primary)")
ax.bar(x, ga, w, color="#7AD151", label="gate flux")
ax.bar(x + w, dp, w, color="#5B3A9E", label="displacement")
ax.axhline(0, color="#888888", lw=0.6)
ax.set_xticks(x); ax.set_xticklabels([r.replace("_", " ") for r in ORDER], rotation=60, ha="right", fontsize=6.3)
ax.set_ylabel("Net K$^+$ per 280 ns"); ax.set_ylim(-1.5, 7.5)
ax.legend(frameon=False, loc="upper right", fontsize=6.5, handlelength=0.9, labelspacing=0.2, borderpad=0.1)
ax.set_title("Three independent estimators agree", loc="left", fontsize=8.2)

# --- B: G310S/WT ratio ---
ax = fig.add_subplot(gs[1])
vals = []
for lab, key in (("net filter\ncharge", None), ("gate\nflux", "gate"), ("displace-\nment", "disp")):
    if key is None: a, b = sum(D.net(r) for r in D.GS), sum(D.net(r) for r in D.WT)
    else: a, b = sum(rows[r][key] for r in D.GS), sum(rows[r][key] for r in D.WT)
    vals.append((lab, (a / 840) / (b / 1120)))
ax.bar(range(3), [v for _, v in vals], 0.55, color=["#1B7EBD", "#7AD151", "#5B3A9E"])
for i, (_, v) in enumerate(vals): ax.text(i, v + 0.024, f"{v:.2f}", ha="center", fontsize=7.2, fontweight="bold")
ax.axhline(1.0, color=C.EXPERIMENT, ls="--", lw=0.8)
ax.text(2.4, 1.0, "no effect", fontsize=6.3, color=C.EXPERIMENT, ha="right", va="bottom")
ax.set_xticks(range(3)); ax.set_xticklabels([l for l, _ in vals], fontsize=6.5)
ax.set_ylabel("G310S / WT rate ratio"); ax.set_ylim(0, 1.15)
ax.set_title("G310S effect", loc="left", fontsize=8.2)

# --- C: realised potential ---
ax = fig.add_subplot(gs[2])
for i, r in enumerate(ORDER):
    col = "#7AD151" if r.startswith("V280") else "#1B7EBD"
    ax.plot(i, rows[r]["V"], "o", ms=3.6, color=col)
ax.axhline(565, color=C.EXPERIMENT, ls="--", lw=0.8); ax.axhline(282, color=C.EXPERIMENT, ls="--", lw=0.8)
ax.text(len(ORDER) - 0.5, 578, "nominal 565 mV", fontsize=6.3, color=C.EXPERIMENT, ha="right")
ax.text(len(ORDER) - 0.5, 295, "nominal 282 mV", fontsize=6.3, color=C.EXPERIMENT, ha="right")
ax.set_xticks(x); ax.set_xticklabels([r.replace("_", " ") for r in ORDER], rotation=60, ha="right", fontsize=6.3)
ax.set_ylabel("Realised potential (mV)"); ax.set_ylim(250, 620)
ax.set_title("Potential per run", loc="left", fontsize=8.2)

for ax, L in zip(fig.axes, "ABC"):
    ax.text(-0.13 if L == "A" else -0.26, 1.15, L, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
O = os.path.dirname(os.path.abspath(__file__)) + "/../out"
fig.savefig(f"{O}/FIG_S1_estimators.pdf"); fig.savefig(f"{O}/FIG_S1_estimators.png", dpi=600)
print("FIG_S1 ok |", [(l, round(v, 3)) for l, v in vals])
