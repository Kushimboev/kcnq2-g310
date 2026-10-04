#!/usr/bin/env python3
"""Fig. S2, CHARMM36m control. A: expected and observed entries from the cavity into the filter. B: cavity occupancy (ions are present but do not enter).
C: gate water and S310 OG in both force fields. Data: C36_FINAL_SUMMARY.txt."""
import sys, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import poisson
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../style")
import paperdata as D
from palette import C, MM, use_paper_style
use_paper_style()
# from C36_FINAL_SUMMARY.txt; per run: (entries from the cavity, gate water, S310 OG, % of frames with K+ in the cavity)
C36 = {"C36W r1": (0, 10.07, np.nan, 21.4), "C36W r2": (0, 9.75, np.nan, 13.2),
       "C36G r1": (0, 4.52, 3.44, 0.0), "C36G r2": (1, 4.17, 3.23, 2.1)}
AMB_W = [D.mech(D.WT)[r]["gate_water"] for r in D.WT]; AMB_G = [D.mech(D.GS)[r]["gate_water"] for r in D.GS]
WTC, GSC = "#1B7EBD", "#D95F02"
rate_amber = 17 / 1120.0                      # wild type, entries per ns
exp560 = rate_amber * 560

fig = plt.figure(figsize=(MM.TWO_COL, 2.7))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.15], wspace=0.44, left=0.07, right=0.985, bottom=0.20, top=0.86)

# --- A: expected vs observed ---
ax = fig.add_subplot(gs[0])
k = np.arange(0, 18)
ax.bar(k, poisson.pmf(k, exp560) * 100, 0.8, color="#CCCCCC", label=f"expected at the\nAmber14sb rate ({exp560:.1f})")
ax.axvline(0, color=WTC, lw=2.0)
ax.annotate("observed: 0", (0, poisson.pmf(0, exp560) * 100 + 4), xytext=(3.2, 14), fontsize=7, color=WTC,
            fontweight="bold", arrowprops=dict(arrowstyle="->", lw=0.8, color=WTC))
ax.set_xlabel("Cavity $\\rightarrow$ filter entries in 560 ns"); ax.set_ylabel("Probability (%)")
ax.set_xlim(-0.8, 17); ax.set_ylim(0, 22)
ax.legend(frameon=False, loc="upper right", fontsize=6.5, handlelength=0.9, borderpad=0.1)
ax.set_title("CHARMM36m WT: no entries", loc="left", fontsize=8.2)
pm, pe = f"{poisson.pmf(0, exp560):.1e}".split("e")
ax.text(0.98, 0.80, f"$P$(0) = {pm}$\\times$10$^{{{int(pe)}}}$\nconditional $p$ = 0.001", transform=ax.transAxes,
        ha="right", va="top", fontsize=6.6, color=C.ZERO)

# --- B: cavity occupancy ---
ax = fig.add_subplot(gs[1])
names = list(C36); occ = [C36[n][3] for n in names]
cols = [WTC if "W" in n else GSC for n in names]
ax.bar(range(4), occ, 0.6, color=cols)
for i, n in enumerate(names):
    ax.text(i, occ[i] + 0.9, f"{C36[n][0]} {'entry' if C36[n][0] == 1 else 'entries'}", ha="center", fontsize=6.3, color=cols[i])
ax.set_xticks(range(4)); ax.set_xticklabels([n.replace("C36W", "WT").replace("C36G", "G310S") for n in names],
                                            rotation=30, ha="right", fontsize=6.6)
ax.set_ylabel("K$^+$ in the cavity (% of frames)"); ax.set_ylim(0, 28)
ax.set_title("Ions reach the cavity", loc="left", fontsize=8.2)
ax.text(0.98, 0.97, "the block is at the filter,\nnot at the gate", transform=ax.transAxes, ha="right", va="top", fontsize=6.3, color=C.ZERO)

# --- C: gate water in both force fields + OG ---
ax = fig.add_subplot(gs[2])
data = [("Amber14sb", AMB_W, AMB_G), ("CHARMM36m", [C36["C36W r1"][1], C36["C36W r2"][1]], [C36["C36G r1"][1], C36["C36G r2"][1]])]
for gi, (nm, w, g) in enumerate(data):
    for xo, vals, col in ((0, w, WTC), (0.42, g, GSC)):
        x = gi + xo
        ax.bar(x, np.mean(vals), 0.34, color=col, alpha=0.30, lw=0)
        ax.plot(np.full(len(vals), x) + np.linspace(-0.06, 0.06, len(vals)), vals, "o", ms=3.4, color=col)
    ax.text(gi + 0.21, -0.155, nm, ha="center", fontsize=7.2, transform=ax.get_xaxis_transform())
    ax.text(gi + 0.21, 11.2, f"ratio {np.mean(g)/np.mean(w):.2f}", ha="center", fontsize=6.6, color=C.ZERO)
ax.set_xticks([0, 0.42, 1, 1.42]); ax.set_xticklabels(["WT", "G310S"] * 2, fontsize=7); ax.tick_params(axis="x", pad=1.5)
ax.set_ylabel("Waters in the G310 plane"); ax.set_ylim(0, 12.6); ax.set_xlim(-0.35, 1.77)
ax.set_title("Mechanism transfers between force fields", loc="left", fontsize=8.2)

for ax, L in zip(fig.axes, "ABC"):
    ax.text(-0.20, 1.15, L, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
O = os.path.dirname(os.path.abspath(__file__)) + "/../out"
fig.savefig(f"{O}/FIG_S2_CHARMM36m.pdf"); fig.savefig(f"{O}/FIG_S2_CHARMM36m.png", dpi=600)
print(f"FIG_S2 ok | expected {exp560:.1f}, P(0) = {poisson.pmf(0, exp560):.2e}")
