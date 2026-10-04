#!/usr/bin/env python3
"""Fig. 4, controls.
A: cumulative net charge (wild type, closed 7CR3, G310 narrowed); the three colours are clearly distinct (blue / dark purple / black dashed).
B: conductance at 565 mV: replicas (small dots) and pooled value (large marker, 95 % Poisson interval). C: wild-type conductance at 282 and 565 mV; an ohmic pore gives the same value (grey band).
The long bars are 95 % Poisson intervals: only 3 events at 282 mV (2, 0 and 1 in the replicas)."""
import sys, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerTuple
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../style")
import paperdata as D
from palette import C, MM, use_paper_style
use_paper_style()
WTC, GSC, CLC, ANCC = "#1B7EBD", "#D95F02", "#5B3A9E", "#222222"      # WT blue, G310S orange, closed dark purple, narrowed black

fig = plt.figure(figsize=(MM.TWO_COL, 2.8))
gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1.0, 1.0], wspace=0.5, left=0.075, right=0.985, bottom=0.17, top=0.86)

# --- A: cumulative net charge ---
ax = fig.add_subplot(gs[0])
for r in D.WT:
    t, y = D.cumulative(r); ax.step(t, y, where="post", color=WTC, lw=1.0, alpha=0.9)
for r in D.CR3:
    t, y = D.cumulative(r); ax.step(t, y, where="post", color=CLC, lw=2.8, alpha=0.9, zorder=2)
t, y = D.cumulative("ANC_r1"); ax.step(t, y, where="post", color=ANCC, lw=1.5, ls=(0, (3.5, 1.8)), zorder=3)
h = [Line2D([], [], color=WTC, lw=1.4), Line2D([], [], color=CLC, lw=2.6), Line2D([], [], color=ANCC, lw=1.5, ls=(0, (3.5, 1.8)))]
ax.legend(h, ["wild type, G310 free (4 replicas)", "closed pore, 7CR3 (2 replicas)", "G310 restrained at 7.6 $\\AA$ (1 replica)*"],
          frameon=False, loc="upper left", handlelength=1.8, borderpad=0.1, labelspacing=0.25, fontsize=6.6)
ax.set_xlabel("Time (ns)"); ax.set_ylabel("Net outward K$^+$ crossings")
ax.set_xlim(0, 280); ax.set_ylim(-0.35, 6.9); ax.set_yticks(range(6))
ax.set_title("Controls do not conduct", loc="left")

# --- B: conductance at 565 mV: replicas + pooled ---
ax = fig.add_subplot(gs[1])
rows = [("WT", D.WT, WTC), ("G310S", D.GS, GSC), ("closed\n(7CR3)", D.CR3, CLC), ("G310\n7.6 $\\AA$", ["ANC_r1"], ANCC)]
for i, (lab, runs, col) in enumerate(rows):
    n = sum(D.net(r) for r in runs); tt = 280.0 * len(runs)
    g = D.G_pS(n, tt, 565.0); lo, hi = D.G_ci(n, tt, 565.0)
    for j, r in enumerate(runs):
        ax.plot(i - 0.24 + (j - (len(runs) - 1) / 2) * 0.06, D.G_pS(D.net(r), 280.0, 565.0), "o", ms=3, color=col, alpha=0.6, mew=0, clip_on=False)
    ax.errorbar(i + 0.14, g, yerr=[[g - lo], [hi - g]], fmt="o", ms=5, color=col, lw=1.2, capsize=2.5, zorder=4, clip_on=False)
    ax.text(i + 0.14, hi + 0.35, f"{n}{'*' if i == 3 else ''}/{int(tt)} ns", ha="center", fontsize=6.3, color=col)
ax.set_xticks(range(len(rows))); ax.set_xticklabels([r[0] for r in rows], fontsize=6.8)
ax.set_ylabel("Conductance at 565 mV (pS)"); ax.set_xlim(-0.6, len(rows) - 0.35); ax.set_ylim(-0.4, 9.2); ax.set_yticks([0, 2, 4, 6, 8])
ax.set_title("Conductance by system", loc="left")
ax.text(0.99, 0.985, "small dots: replicas\nlarge markers: pooled,\n95 % Poisson interval", transform=ax.transAxes,
        ha="right", va="top", fontsize=5.9, color=C.ZERO, linespacing=1.3)

# --- C: voltage dependence (conductance) ---
ax = fig.add_subplot(gs[2])
n565 = sum(D.net(r) for r in D.WT); g565 = D.G_pS(n565, 1120.0, 565.0); lo565, hi565 = D.G_ci(n565, 1120.0, 565.0)
band = ax.axhspan(lo565, hi565, color="#E4E4E4", zorder=0, lw=0); ln = ax.axhline(g565, color=C.EXPERIMENT, ls="--", lw=0.9, zorder=1)
for V, runs in ((282, D.V280), (565, D.WT)):
    filled = V == 565
    for j, r in enumerate(runs):
        ax.plot(V - 40 + (j - (len(runs) - 1) / 2) * 9, D.G_pS(D.net(r), 280.0, V), "o", ms=3, color=WTC, mfc=WTC if filled else "white", alpha=0.8, mew=0.8, clip_on=False)
    n = sum(D.net(r) for r in runs); tt = 280.0 * len(runs); g = D.G_pS(n, tt, V); lo, hi = D.G_ci(n, tt, V)
    ax.errorbar(V + 14, g, yerr=[[g - lo], [hi - g]], fmt="o", ms=6, color=WTC, mfc=WTC if filled else "white", mew=1.4, lw=1.3, capsize=3, zorder=4, clip_on=False)
    ax.text(V + 14, hi + 0.3, f"{n} crossings\n/ {int(tt)} ns", ha="center", va="bottom", fontsize=6.0, color=WTC, linespacing=1.15)
ax.text(136, g565 + 0.18, "ohmic", fontsize=6.4, color=C.EXPERIMENT, ha="left", va="bottom")
ax.set_xticks([282, 565]); ax.set_xlim(130, 700); ax.set_ylim(-0.4, 9.2)
ax.set_xlabel("Applied potential (mV)"); ax.set_ylabel("Conductance (pS)")
ax.set_title("Voltage dependence", loc="left")
ax.text(0.03, 0.985, "$G$(282)/$G$(565) = 0.47\n(0.09–1.63), $p$ = 0.33", transform=ax.transAxes, ha="left", va="top", fontsize=6.4, color=C.ZERO, linespacing=1.35)

for ax, L in zip(fig.axes, "ABC"):
    ax.text(-0.19, 1.16, L, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
O = os.path.dirname(os.path.abspath(__file__)) + "/../out"
fig.savefig(f"{O}/FIG_4_controls.pdf"); fig.savefig(f"{O}/FIG_4_controls.png", dpi=600)
print("FIG_4 ok")
