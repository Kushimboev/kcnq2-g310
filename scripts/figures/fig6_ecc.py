#!/usr/bin/env python3
"""Fig. 6, charge scaling (Amber14sb-ECC, q = 0.78). A: matched pairs (same structure and gate width, only the charges differ).
B: filter occupancy (scaled charges lock the filter with three ions). C: shift of the threshold (rate against the G310 width).
All data use the window with the first 20 ns excluded, as in the other figures.
Source: paperdata.book/pool -> {RUN}_b20.books.txt (event books written by conduction/fig6_b20/b20_ecc.sh);
x axis (narrow diagonal): conduction/fig6_b20/d1_b20.json (d1_b20.py)."""
import sys, os, re, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, HERE + "/../style")
import paperdata as D
from palette import C, MM, use_paper_style
use_paper_style()
ECCC, FULLC = "#B2182B", "#1B7EBD"

book, pool, ratio_ci = D.book, D.pool, D.rate_ratio
D1 = json.load(open(f"{D.P}/13_conduction/scripts/analysis_2026-09-29_fig6_b20/d1_b20.json"))
def d1(runs):
    """Narrow G310 diagonal (first 20 ns excluded), frame-weighted mean."""
    return sum(D1[r]["sum_nar_b20"] for r in runs) / sum(D1[r]["nb20"] for r in runs)

PAIRS = [(lab.replace(" (", "\n("), e, f) for lab, e, f in D.ECC_PAIRS]
fig = plt.figure(figsize=(MM.TWO_COL, 2.75))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.05, 1.15], wspace=0.46, left=0.07, right=0.985, bottom=0.18, top=0.86)

# --- A: matched pairs ---
ax = fig.add_subplot(gs[0])
for i, (lab, ecc, full) in enumerate(PAIRS):
    for xo, runs, col in ((-0.23, ecc, ECCC), (0.23, full, FULLC)):     # 0.23: keeps the "n/T ns" label clear of the neighbouring bar
        n, t = pool(runs); rate = n / t * 1000
        ax.bar(i + xo, rate, 0.34, color=col, alpha=0.85)
        ax.text(i + xo, rate + 1.1, f"{n}/{t:.1f} ns".replace(".0 ns", " ns"), ha="center", fontsize=6.2, color=col)
ne, te = pool([r for _, e, _ in PAIRS for r in e]); nf, tf = pool([r for _, _, f in PAIRS for r in f])
assert (ne, nf) == (2, 23) and abs(te - 840.0) < 0.05 and abs(tf - 1217.8) < 0.05, "Fig 6: counts differ from the reported values (2/840 vs 23/1218, first 20 ns excluded)"
r, lo, hi, p = ratio_ci(ne, te, nf, tf)
ax.set_xticks(range(len(PAIRS))); ax.set_xticklabels([p[0] for p in PAIRS], fontsize=7.2)
ax.set_ylabel("Net K$^+$ transfer (per $\\mu$s)"); ax.set_ylim(0, 36)
ax.set_title("Scaling suppresses transfer", loc="left", fontsize=8.2)
h = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (ECCC, FULLC)]
ax.legend(h, ["scaled charges (ECC)", "full charges"], frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.85),
          handlelength=0.9, handleheight=0.9, labelspacing=0.25, fontsize=6.6, borderpad=0.1)
pm, pe = f"{p:.1e}".split("e")
ax.text(0.5, 0.985, f"pooled ratio {r:.2f}\n({lo:.2f}–{hi:.2f}), $p$ = {pm}$\\times$10$^{{{int(pe)}}}$", transform=ax.transAxes,
        ha="center", va="top", fontsize=6.5, color=C.ZERO)

# --- B: filter occupancy ---
ax = fig.add_subplot(gs[1])
groups = [("scaled\ncharges", ["AN_r1", "AN_r2", "KN_r1"], ECCC), ("full\ncharges", ["ANE_r1", "ANE_r2", "KNE_r1", "KNE_r2"], FULLC)]
for i, (lab, runs, col) in enumerate(groups):
    two = np.mean([book(r)[2].get(2, 0) for r in runs]); three = np.mean([book(r)[2].get(3, 0) for r in runs])
    ax.bar(i, two, 0.5, color=C.K, alpha=0.85); ax.bar(i, three, 0.5, bottom=two, color=C.KSF, alpha=0.85)
    ax.text(i, 103, lab, ha="center", fontsize=7, color=col, fontweight="bold")
ax.set_xticks([]); ax.set_ylim(0, 118); ax.set_xlim(-0.55, 2.15)
ax.set_ylabel("Filter occupancy (% of frames)")
ax.annotate("2 ions", (1.27, 45), xytext=(5, 0), textcoords="offset points", fontsize=7, color=C.K, fontweight="bold", va="center")
ax.annotate("3 ions", (1.27, 94), xytext=(5, 0), textcoords="offset points", fontsize=7, color=C.KSF, fontweight="bold", va="center")
ax.set_title("Scaling locks the filter", loc="left", fontsize=8.2)

# --- C: shift of the threshold ---
ax = fig.add_subplot(gs[2])
PTS = [("AN", ["AN_r1", "AN_r2"], ECCC, "o", "KCNQ2, ECC"),
       ("A9", ["A9_r1"], ECCC, "o", None),
       ("KN", ["KN_r1"], ECCC, "s", "KCNQ1, ECC"),
       ("K10", ["K10_r1"], ECCC, "s", None),
       ("ANE", ["ANE_r1", "ANE_r2", "ANE_r3"], FULLC, "o", "KCNQ2, full"),
       ("KNE", ["KNE_r1", "KNE_r2"], FULLC, "s", "KCNQ1, full")]
XY = {}; DODGE = {"ANE": -0.06, "KNE": +0.06}        # the two full-charge groups have almost the same width (8.39 and 8.37 A) -> offset so that the markers do not coincide
for nm, runs, col, mk, lab in PTS:
    x = d1(runs); n, t = pool(runs); rate = n / t * 1000; XY[nm] = (x, rate, n, t)
    plo, phi = (0.0 if n == 0 else chi2.ppf(0.025, 2 * n) / 2), chi2.ppf(0.975, 2 * n + 2) / 2      # exact 95 % Poisson interval
    ax.errorbar(x + DODGE.get(nm, 0.0), rate, yerr=[[rate - plo / t * 1000], [phi / t * 1000 - rate]], fmt=mk, ms=6, color=col, mfc=col if col == ECCC else "white",
                mew=1.3, elinewidth=0.8, capsize=2, label=lab)
dx = XY["A9"][0] - XY["ANE"][0]
ax.annotate("", xy=XY["A9"][:2], xytext=(XY["ANE"][0] + DODGE["ANE"], XY["ANE"][1]), arrowprops=dict(arrowstyle="->", lw=0.9, color="#888888",
            connectionstyle="arc3,rad=-0.25", shrinkA=7, shrinkB=7))
ax.text(9.62, 45.5, f"KCNQ2: similar rate at\n+{dx:.2f} $\\AA$ with scaling\n(single run)", fontsize=6.3, color="#888888", ha="center", va="top")
ax.set_xlabel("G310 C$\\alpha$–C$\\alpha$ diagonal ($\\AA$)"); ax.set_ylabel("Net K$^+$ transfer (per $\\mu$s)")
ax.set_xlim(8.0, 10.25); ax.set_ylim(0, 60)
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=6.0, handlelength=0.9, labelspacing=0.25, columnspacing=0.8, borderpad=0.1)
ax.set_title("A wider gate offsets scaling", loc="left", fontsize=8.2)

for ax, L in zip(fig.axes, "ABC"):
    ax.text(-0.19, 1.15, L, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
O = os.path.dirname(os.path.abspath(__file__)) + "/../out"
fig.savefig(f"{O}/FIG_6_ECC.pdf"); fig.savefig(f"{O}/FIG_6_ECC.png", dpi=600)
# --- numbers for the text (first 20 ns excluded); also written to DATA_fig6_b20.txt ---
allecc = ["AN_r1", "AN_r2", "KN_r1", "K10_r1"]; na, ta = pool(allecc); ra, loa, hia, pa = ratio_ci(na, ta, nf, tf)
nk, tk = pool(D.KCNQ1_FULL)
occ = lambda runs, k: np.mean([book(r)[2].get(k, 0) for r in runs])
L = [f"FIG_6 (b20) | matched pooled ECC {ne}/{te:.1f} ns vs full {nf}/{tf:.1f} ns -> ratio {r:.3f} ({lo:.3f}-{hi:.3f}), p={p:.2e}, suppression x{1/r:.1f}",
     f"all-ECC incl. K10: {na}/{ta:.1f} ns vs {nf}/{tf:.1f} -> {ra:.3f} ({loa:.3f}-{hia:.3f}), p={pa:.2e}, x{1/ra:.1f}",
     f"KCNQ1 benchmark (KNE r1+r2): {nk}/{tk:.1f} ns -> G = {D.G_pS(nk, tk):.2f} pS",
     f"filter occupancy: ECC 3 ions {occ(groups[0][1], 3):.0f}% (2 ions {occ(groups[0][1], 2):.0f}%); full 2 ions {occ(groups[1][1], 2):.0f}% (3 ions {occ(groups[1][1], 3):.0f}%)"]
for nm, runs, col, mk, lab in PTS:
    x, rate, n, t = XY[nm]; L.append(f"panel C {nm:4s} d1(narrow,b20) {x:.3f} | {n}/{t:.1f} ns = {rate:.1f}/us")
L.append(f"A9 - ANE shift {dx:.3f} A ; K10 - KNE shift {XY['K10'][0] - XY['KNE'][0]:.3f} A")
for r_ in ["AN_r1", "AN_r2", "KN_r1", "K10_r1", "A9_r1", "ANE_r1", "ANE_r2", "ANE_r3", "KNE_r1", "KNE_r2"]:
    L.append(f"run {r_:7s} net {book(r_)[0]:+d} in {book(r_)[1]:.1f} ns")
open(HERE + "/../DATA_fig6_b20.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
