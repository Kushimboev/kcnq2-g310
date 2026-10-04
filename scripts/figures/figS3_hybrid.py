#!/usr/bin/env python3
"""Fig. S3, mixed tetramers: H1 = one S310 subunit (chain A), H2 = two adjacent S310 subunits (A+B), 2 x 300 ns each.
A: gate water against the number of variant subunits (additive line). B: f = fraction of the G310S effect (95 % CI; expectation n/4).
C: ratio of net crossing rates (95 % conditional-binomial CI). Numbers come through paperdata.py."""
import sys, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../style")
import paperdata as D
from palette import C, MM, use_paper_style
use_paper_style()
WTC, GSC, HYC = "#1B7EBD", "#D95F02", "#984EA3"   # a separate colour for the mixed pores (CIE76 dE 49 / 97 from WT and G310S)

GROUPS = [(0, "0\nWT", D.WT, WTC), (1, "1\nH1", D.HYB["H1"], HYC), (2, "2\nH2", D.HYB["H2"], HYC), (4, "4\nG310S", D.GS, GSC)]
fig = plt.figure(figsize=(MM.TWO_COL, 2.35))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 0.85, 0.95], wspace=0.42, left=0.065, right=0.99, bottom=0.235, top=0.86)

# --- A: gate water vs number of variant subunits ---
ax = fig.add_subplot(gs[0, 0])
for n, lab, runs, col in GROUPS:
    v, m, se = D.gw_stats(runs)
    ax.plot(np.full(len(v), n) + np.linspace(-0.09, 0.09, len(v)), v, "o", ms=4.0, color=col, zorder=3)
    ax.errorbar(n, m, yerr=(0 if not np.isfinite(se) else se), fmt="_", ms=13, mew=1.8, color=col, ecolor=col, capsize=3, zorder=4)
_, wt, _ = D.gw_stats(D.WT); _, gsm, _ = D.gw_stats(D.GS)
ax.plot([0, 4], [wt, gsm], ls="--", lw=0.9, color="#999999", zorder=1, label="additive expectation")
ax.set_xticks([0, 1, 2, 4]); ax.set_xticklabels([g[1] for g in GROUPS], fontsize=7.0)
ax.set_xlim(-0.45, 4.45); ax.set_ylim(4, 10)
ax.set_ylabel("Waters in the G310 plane"); ax.set_xlabel("Subunits carrying Ser310", labelpad=3.0)
ax.legend(frameon=False, loc="upper right", handlelength=1.6, fontsize=6.6, borderpad=0.0)
ax.set_title("Gate water vs. subunit count", loc="left")

# --- B: f and the additive expectation ---
ax = fig.add_subplot(gs[0, 1])
for i, (k, exp) in enumerate((("H1", 0.25), ("H2", 0.50))):
    f, se, lo, hi = D.f_water(D.HYB[k])
    ax.errorbar(i, f, yerr=[[f - lo], [hi - f]], fmt="o", ms=5, color=HYC, ecolor=HYC, capsize=3.5, lw=1.3, zorder=3)
    ax.plot(i, exp, marker="_", ms=15, mew=1.6, color="#999999", zorder=2)
    ax.text(i + 0.17, f, f"{f:.2f}", fontsize=7, color=HYC, va="center")
ax.plot([], [], marker="_", ms=10, mew=1.6, ls="none", color="#999999", label="additive (n/4)")
ax.axhline(0, color="#CCCCCC", lw=0.7, zorder=0); ax.axhline(1, color="#CCCCCC", lw=0.7, zorder=0)
ax.set_xticks([0, 1]); ax.set_xticklabels(["H1\n(1 Ser)", "H2\n(2 Ser)"], fontsize=7.0)
ax.set_xlim(-0.5, 1.5); ax.set_ylim(-0.1, 1.05)
ax.set_ylabel("Fraction of the G310S\nwater-column effect")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.93), handlelength=1.3, fontsize=6.6, borderpad=0.0)
ax.set_title("Effect is additive", loc="left")

# --- C: ratio of crossing rates ---
ax = fig.add_subplot(gs[0, 2])
NW, TW = 17, 1120.0
pts = [("H1", D.HYB["H1"], HYC), ("H2", D.HYB["H2"], HYC), ("G310S", None, GSC)]
for i, (lab, runs, col) in enumerate(pts):
    if runs is None: n, t = 2, 840.0
    else: n = sum(D.hyb_book(r)[0] for r in runs); t = sum(D.hyb_book(r)[1] for r in runs)
    r, lo, hi, p = D.rate_ratio(n, t, NW, TW)
    hi_ = min(hi, 3.4)
    ax.errorbar(i, r, yerr=[[r - lo], [hi_ - r]], fmt="o", ms=5, color=col, ecolor=col, capsize=3.5, lw=1.3,
                uplims=[hi > 3.4], zorder=3)
    ax.text(i, hi_ + 0.13, f"{n}/{t:.0f} ns", fontsize=6.4, color=col, va="bottom", ha="center")
ax.axhline(1, color=WTC, ls="--", lw=0.9, zorder=1)
ax.text(2.48, 1.08, "wild type", fontsize=6.4, color=WTC, ha="right")
ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["H1", "H2", "G310S"], fontsize=6.8)
ax.set_xlim(-0.55, 2.55); ax.set_ylim(0, 3.5)
ax.set_ylabel("Crossing rate / wild type")
ax.set_title("Conduction: not resolved", loc="left")

for ax, L in zip(fig.axes, "ABC"):
    ax.text(-0.20, 1.20, L, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
O = os.path.dirname(os.path.abspath(__file__)) + "/../out"
fig.savefig(f"{O}/FIG_S3_hybrid.pdf"); fig.savefig(f"{O}/FIG_S3_hybrid.png", dpi=600)
f1 = D.f_water(D.HYB["H1"]); f2 = D.f_water(D.HYB["H2"])
assert 0.13 <= f1[0] <= 0.32 and 0.37 <= f2[0] <= 0.78, (f1, f2)
print(f"FIG_S3 ok | f(H1) {f1[0]:.2f} ({f1[2]:.2f}-{f1[3]:.2f}) · f(H2) {f2[0]:.2f} ({f2[2]:.2f}-{f2[3]:.2f})")
