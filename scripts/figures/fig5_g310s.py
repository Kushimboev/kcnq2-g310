#!/usr/bin/env python3
"""Fig. 5, G310S. A: cumulative net charge (WT vs G310S). B: gate water (both force fields). C: ion coordination at the gate plane.
D: coordination profile along z."""
import sys, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../style")
import paperdata as D
from palette import C, MM, use_paper_style
use_paper_style()
WTC, GSC = "#1B7EBD", "#D95F02"          # WT = G310 blue, G310S = serine orange
# CHARMM36m gate water (from C36_FINAL_SUMMARY.txt)
C36 = dict(WT=[10.07, 9.75], GS=[4.52, 4.17])

fig = plt.figure(figsize=(MM.TWO_COL, 4.9))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[1.1, 1.0],
                      hspace=0.55, wspace=0.30, left=0.07, right=0.985, bottom=0.09, top=0.92)

# --- A: cumulative net charge ---
ax = fig.add_subplot(gs[0, 0])
for r in D.WT:
    t, y = D.cumulative(r); ax.step(t, y, where="post", color=WTC, lw=1.1, alpha=0.9)
for r in D.GS:
    t, y = D.cumulative(r); ax.step(t, y, where="post", color=GSC, lw=1.4, alpha=0.9)
ax.plot([], [], color=WTC, lw=1.2, label="WT (4 replicas)"); ax.plot([], [], color=GSC, lw=1.4, label="G310S (3)")
ax.set_xlabel("Time (ns)"); ax.set_ylabel("Net outward K$^+$ crossings")
ax.set_xlim(0, 280); ax.set_ylim(-0.3, 9.6); ax.set_yticks(range(6))
ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.79), handlelength=1.1, borderpad=0.1, labelspacing=0.2, fontsize=7)
ax.set_title("Conduction", loc="left")
ax.text(0.02, 0.985, "17 / 1.12 $\\mu$s   vs   2 / 0.84 $\\mu$s\nratio 0.16 (0.02–0.66), $p$ = 0.003",
        transform=ax.transAxes, ha="left", va="top", fontsize=6.8, color=C.ZERO, linespacing=1.35)

# --- B: gate water, both force fields ---
ax = fig.add_subplot(gs[0, 1])
m = D.mech(D.WT + D.GS)
groups = [("Amber14sb", [m[r]["gate_water"] for r in D.WT], [m[r]["gate_water"] for r in D.GS]),
          ("CHARMM36m", C36["WT"], C36["GS"])]
for gi, (name, w, g) in enumerate(groups):
    for xoff, vals, col in ((0, w, WTC), (0.42, g, GSC)):
        x = gi + xoff
        ax.bar(x, np.mean(vals), 0.34, color=col, alpha=0.30, lw=0)
        ax.plot(np.full(len(vals), x) + np.linspace(-0.07, 0.07, len(vals)), vals, "o", ms=3.4, color=col)
    ax.text(gi + 0.21, -0.155, name, ha="center", fontsize=7.2, transform=ax.get_xaxis_transform())
ax.set_xticks([0, 0.42, 1, 1.42]); ax.set_xticklabels(["WT", "G310S"] * 2, fontsize=7)
ax.tick_params(axis="x", pad=1.5)
ax.set_ylabel("Waters in the G310 plane"); ax.set_ylim(0, 12); ax.set_xlim(-0.35, 1.77)
ax.set_title("Gate water", loc="left")
ax.text(0.5, 0.955, "ratio 0.54   /   0.44", transform=ax.transAxes, ha="center", fontsize=6.8, color=C.ZERO)

# --- C: coordination at the gate plane ---
ax = fig.add_subplot(gs[1, 0])
h = D.hyd(D.WT + D.GS)
mk = lambda runs, k: np.mean([h[r][k] for r in runs])
se = lambda runs, k: np.std([h[r][k] for r in runs], ddof=1) / np.sqrt(len(runs))
bulk = np.mean([h[r]["bulk"] for r in D.WT + D.GS])
labs = ["water O", "Ser310 O$\\gamma$", "other protein O"]
cols = [C.WATER, GSC, "#7570B3"]
for xi, runs in ((0, D.WT), (1, D.GS)):
    bot = 0
    for k, lab, col in zip(["water", "OG", "protO"], labs, cols):
        v = mk(runs, k)
        ax.bar(xi, v, 0.52, bottom=bot, color=col, edgecolor="white", lw=0.6, label=lab if xi == 0 else None)
        bot += v
    ax.errorbar(xi, bot, yerr=se(runs, "total"), color=C.ZERO, lw=1.0, capsize=2.5)
    ax.text(xi, bot + 0.55, f"{bot:.2f}", ha="center", fontsize=7.5, fontweight="bold")
ax.axhline(bulk, color=C.EXPERIMENT, ls="--", lw=0.9)
ax.text(1.47, 7.75, f"bulk {bulk:.2f}", fontsize=6.8, color=C.EXPERIMENT, va="bottom", ha="right")
ax.set_xticks([0, 1]); ax.set_xticklabels(["WT", "G310S"]); ax.set_ylim(0, 10.2); ax.set_xlim(-0.6, 1.5)
ax.set_ylabel("O atoms around K$^+$ at G310")
ax.legend(frameon=False, loc="upper left", ncol=1, handlelength=0.9, labelspacing=0.2, fontsize=6.6, borderpad=0.0)
ax.set_title("First shell is conserved", loc="left")

# --- D: coordination profile along z ---
ax = fig.add_subplot(gs[1, 1])
z, W, G = D.hyd_profile()
okW, okG = np.nansum(W[:, :, 0], axis=0) >= 25, np.nansum(G[:, :, 0], axis=0) >= 25     # >= 25 ion-frames per bin, else not drawn
nz = lambda A, j, ok: np.where(ok, np.nanmean(A[:, :, j], axis=0), np.nan)
ax.plot(z, nz(W, 1, okW), color=WTC, lw=1.4, label="WT water")
ax.plot(z, nz(G, 1, okG), color=GSC, lw=1.4, label="G310S water")
ax.plot(z, nz(G, 2, okG), color=GSC, lw=1.4, ls="--", label="G310S Ser O$\\gamma$")
ax.text(0.97, 3.5, "no G310S data\nabove the gate\n(< 25 ion-frames\nper bin)", fontsize=6.0, color=GSC, ha="right", va="top", linespacing=1.25)
ax.axvline(0, color="#BBBBBB", lw=0.8, zorder=0)
ax.text(0.02, 7.6, "G310\nplane", fontsize=6.6, color=C.GLY310, ha="left", va="top")
ax.set_xlabel("z relative to the G310 plane (nm)"); ax.set_ylabel("O atoms around K$^+$")
ax.set_xlim(-1.2, 1.0); ax.set_ylim(0, 8)
ax.legend(frameon=False, loc="center left", bbox_to_anchor=(0.0, 0.36), handlelength=1.4, labelspacing=0.2, fontsize=6.8, borderpad=0.2)
ax.set_title("Where the shell changes", loc="left")
ax.annotate("", xy=(0.97, 7.3), xytext=(0.55, 7.3), arrowprops=dict(arrowstyle="->", lw=0.7, color="#888888"))
ax.text(0.97, 7.55, "extracellular", fontsize=6.2, color="#888888", ha="right")


for ax, L in zip(fig.axes, "ABCD"):
    pos = ax.get_position(); fig.text(pos.x0 - 0.06, pos.y1 + 0.055, L, fontsize=11, fontweight="bold", va="top")
O = os.path.dirname(os.path.abspath(__file__)) + "/../out"
fig.savefig(f"{O}/FIG_5_G310S.pdf"); fig.savefig(f"{O}/FIG_5_G310S.png", dpi=600)
print("FIG_5 ok | WT total", round(mk(D.WT, "total"), 2), "GS", round(mk(D.GS, "total"), 2), "bulk", round(bulk, 2))
