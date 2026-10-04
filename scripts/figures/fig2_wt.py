#!/usr/bin/env python3
"""Fig. 2, wild-type conduction. A: cumulative net charge (4 replicas). B: conductance per replica, pooled, and experiment.
C: filter occupancy (two-ion knock-on)."""
import sys, os, re, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../style")
import paperdata as D
from palette import C, MM, use_paper_style
use_paper_style()
REP_C = ["#08306B", "#4292C6", "#1B9E77", "#A6D854"]     # clearly distinguishable colours (dE >= 42)
fig = plt.figure(figsize=(MM.TWO_COL, 2.7))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 0.95], wspace=0.46, left=0.075, right=0.985, bottom=0.16, top=0.86)

# --- A: cumulative net charge ---
ax = fig.add_subplot(gs[0])
for i, r in enumerate(D.WT):
    t, y = D.cumulative(r)
    ax.step(t, y, where="post", color=REP_C[i], lw=1.2, label=f"replica {i+1}")
ax.set_xlabel("Time (ns)"); ax.set_ylabel("Net outward K$^+$ crossings")
ax.set_xlim(0, 280); ax.set_ylim(-0.3, 5.6); ax.set_yticks(range(6))
ax.legend(frameon=False, loc="upper left", handlelength=1.3, borderpad=0.1, labelspacing=0.25)
ax.set_title("Activated pore conducts", loc="left")
ax.text(0.97, 0.06, "17 crossings / 1.12 $\\mu$s", transform=ax.transAxes, ha="right", fontsize=7.5, color=C.ZERO)

# --- B: conductance (replicas + pooled) and experiment ---
ax = fig.add_subplot(gs[1])
ax.axhspan(17.8 - 3.1, 18.4 + 1.5, color=C.EXPERIMENT, alpha=0.13, lw=0)
ax.text(8.6, 17.6, "experiment,\nsymmetric K$^+$", fontsize=6.3, color=C.EXPERIMENT, va="center", ha="right")
ax.axhspan(5.8 - 0.3, 6.2 + 0.3, color=C.EXPERIMENT, alpha=0.20, lw=0)
ax.text(8.6, 7.4, "experiment,\nphysiological K$^+$", fontsize=6.3, color=C.EXPERIMENT, va="bottom", ha="right")
for i, r in enumerate(D.WT):
    n = D.net(r); g = D.G_pS(n, 280); lo, hi = D.G_ci(n, 280)
    ax.errorbar(i + 1, g, yerr=[[g - lo], [hi - g]], fmt="o", ms=4, color=REP_C[i], lw=1.0, capsize=2)
n = sum(D.net(r) for r in D.WT); g = D.G_pS(n, 1120); lo, hi = D.G_ci(n, 1120)
ax.errorbar(5.4, g, yerr=[[g - lo], [hi - g]], fmt="s", ms=5, color=C.ZERO, lw=1.2, capsize=2.5)
ax.annotate(f"{g:.2f} pS", (5.4, g), textcoords="offset points", xytext=(9, -1), fontsize=7.5, fontweight="bold", va="center")
ax.set_xticks([1, 2, 3, 4, 5.4]); ax.set_xticklabels(["1", "2", "3", "4", "pooled"])
ax.set_xlabel("Replica"); ax.set_ylabel("Conductance (pS)")
ax.set_xlim(0.4, 8.7); ax.set_ylim(0, 22); ax.set_yticks([0, 5, 10, 15, 20])
ax.set_title("Conductance", loc="left")

# --- C: filter occupancy ---
ax = fig.add_subplot(gs[2])
occ = {}
for r in D.WT:
    txt = open(f"{D.BOOKS}/{r}_b20.books.txt").read()
    m = re.search(r"occupancy histogram ([^\n]+)", txt)
    occ[r] = {int(k): float(v.rstrip("%")) for k, v in re.findall(r"(\d+):(\d+%)", m.group(1))}
x = np.arange(4)
two = [occ[r].get(2, 0) for r in D.WT]; three = [occ[r].get(3, 0) for r in D.WT]
ax.bar(x, two, 0.62, color=C.K, alpha=0.85)
ax.bar(x, three, 0.62, bottom=two, color=C.KSF, alpha=0.85)
ax.set_xticks(x); ax.set_xticklabels(["1", "2", "3", "4"]); ax.set_xlabel("Replica")
ax.set_ylabel("Filter occupancy (% of frames)"); ax.set_ylim(0, 118)
ax.annotate("2 ions", (x[-1] + 0.32, two[-1] / 2), xytext=(6, 0), textcoords="offset points",
            fontsize=7, color=C.K, va="center", fontweight="bold")
ax.annotate("3 ions", (x[-1] + 0.32, two[-1] + three[-1] / 2), xytext=(6, 0), textcoords="offset points",
            fontsize=7, color=C.KSF, va="center", fontweight="bold")
ax.set_xlim(-0.6, 4.7)
ax.set_title("Filter state", loc="left")
ax.text(0.5, 0.955, "all entries from the cavity", transform=ax.transAxes, fontsize=6.6, color=C.ZERO, ha="center")

for ax, L in zip(fig.axes, "ABC"):
    ax.text(-0.19, 1.16, L, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
O = os.path.dirname(os.path.abspath(__file__)) + "/../out"
fig.savefig(f"{O}/FIG_2_WT.pdf"); fig.savefig(f"{O}/FIG_2_WT.png", dpi=600)
print("FIG_2 ok:", {r: D.net(r) for r in D.WT}, "pooled", round(g, 2), "pS")
