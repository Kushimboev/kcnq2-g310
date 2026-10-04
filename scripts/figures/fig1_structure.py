#!/usr/bin/env python3
"""Fig. 1, structure. A: rendering of the pore (PyMOL; ../panels/fig1A_pore.png is inserted if present).
B: narrow CA-CA diagonal along S6 (activated 8J01 vs closed 7CR3). C: widening on activation (delta d1); G310 is the pivot."""
import sys, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../style")
from palette import C, MM, use_paper_style
use_paper_style()
H = os.path.dirname(os.path.abspath(__file__))
d = np.load(f"{H}/../out/fig1_data.npz")["diag"]          # res, 8j01 d1,d2, 7cr3 d1,d2, delta
m = (d[:, 0] >= 306) & (d[:, 0] <= 319)
res, a1, b1, dd = d[m, 0].astype(int), d[m, 1], d[m, 3], d[m, 5]
OPEN, CLOSED = "#1B7EBD", "#5B3A9E"     # closed state in dark purple (as in Fig. 4; distinct from the blue)
FIELD = {313, 314, 317, 318}                               # residues of the published gate-distance sets (Nappi 2024; Roscioni 2026, d2-d5)

png = f"{H}/../panels/fig1A_pore.png"
fig = plt.figure(figsize=(MM.TWO_COL, 2.75))
gs = fig.add_gridspec(1, 3, width_ratios=[0.72, 1.12, 1.12], wspace=0.52, left=0.035, right=0.985, bottom=0.17, top=0.86)

# --- A: molecular panel ---
ax = fig.add_subplot(gs[0]); ax.axis("off")
if os.path.exists(png):
    im = mpimg.imread(png)
    if im.shape[2] == 4:                                   # crop the transparent background
        nz = np.argwhere(im[:, :, 3] > 0.02)
        (y0, x0), (y1, x1) = nz.min(0), nz.max(0) + 1
        im = im[max(y0 - 8, 0):y1 + 8, max(x0 - 8, 0):x1 + 8]
    ax.imshow(im); ax.set_title("Activated pore", loc="left")
    Hh, Ww = im.shape[0], im.shape[1]
    ax.annotate("selectivity\nfilter", (0.55 * Ww, 0.22 * Hh), xytext=(1.02 * Ww, 0.05 * Hh), fontsize=6.8,
                color="#4A8FBF", ha="right", va="top", arrowprops=dict(arrowstyle="-", lw=0.6, color="#4A8FBF"))
    ax.annotate("G310", (0.52 * Ww, 0.655 * Hh), xytext=(-0.02 * Ww, 0.60 * Hh), fontsize=7.5, fontweight="bold",
                color=C.GLY310, ha="left", va="center", arrowprops=dict(arrowstyle="-", lw=0.6, color=C.GLY310))
    ax.annotate("S314", (0.46 * Ww, 0.755 * Hh), xytext=(-0.02 * Ww, 0.85 * Hh), fontsize=7.5, fontweight="bold",
                color=C.SER314, ha="left", va="center", arrowprops=dict(arrowstyle="-", lw=0.6, color=C.SER314))
    ax.annotate("", xy=(1.0 * Ww, 0.62 * Hh), xytext=(1.0 * Ww, 0.40 * Hh),
                arrowprops=dict(arrowstyle="->", lw=0.8, color="#AAAAAA"))
    ax.text(1.03 * Ww, 0.40 * Hh, "out", fontsize=6.3, color="#999999", ha="left", va="center")
    ax.text(1.03 * Ww, 0.62 * Hh, "in", fontsize=6.3, color="#999999", ha="left", va="center")
else:
    ax.text(0.5, 0.5, "[PyMOL panel:\nfig1A_pore.png]", ha="center", va="center", fontsize=7, color="#AAAAAA")
    ax.set_title("Pore domain", loc="left")

# --- B: diagonal along S6 ---
ax = fig.add_subplot(gs[1])
ax.plot(res, a1, "o-", color=OPEN, ms=3.4, label="activated (8J01)")
ax.plot(res, b1, "s-", color=CLOSED, ms=3.2, label="closed (7CR3)")
ax.fill_between(res, b1, a1, where=a1 > b1, color=OPEN, alpha=0.10, lw=0)
for r, y, col in ((310, a1[res == 310][0], C.GLY310), (314, a1[res == 314][0], C.SER314)):
    dx, dy, ha = ((0, -14, "center") if r == 310 else (10, -9, "left"))
    ax.annotate(f"{'G' if r == 310 else 'S'}{r}", (r, y), xytext=(dx, dy), textcoords="offset points",
                ha=ha, fontsize=7.5, fontweight="bold", color=col)
ax.set_xlabel("S6 residue"); ax.set_ylabel("Narrow C$\\alpha$–C$\\alpha$ diagonal ($\\AA$)")
ax.set_xticks([306, 310, 314, 318]); ax.set_xlim(305.4, 319.6); ax.set_ylim(5, 27.5)
ax.legend(frameon=False, loc="upper left", handlelength=1.4, labelspacing=0.2, fontsize=7)
ax.set_title("G310: narrowest S6 position", loc="left")

# --- C: widening on activation ---
ax = fig.add_subplot(gs[2])
cols = [C.GLY310 if r == 310 else (C.SER314 if r == 314 else ("#9EC9E2" if r in FIELD else "#CCCCCC")) for r in res]
ax.bar(res, dd, 0.72, color=cols)
ax.axhline(0, color="#888888", lw=0.6)
ax.annotate(f"+{dd[res == 310][0]:.2f} $\\AA$", (310, dd[res == 310][0]), xytext=(-4, 7), textcoords="offset points",
            ha="right", fontsize=7, fontweight="bold", color=C.GLY310)
ax.annotate(f"+{dd[res == 314][0]:.2f} $\\AA$", (314, dd[res == 314][0]), xytext=(-4, 7), textcoords="offset points",
            ha="right", fontsize=7, fontweight="bold", color=C.SER314)
ax.set_xlabel("S6 residue"); ax.set_ylabel("Widening on activation ($\\AA$)")
ax.set_xticks([306, 310, 314, 318]); ax.set_xlim(305.4, 319.6); ax.set_ylim(-1.8, 8.6)
ax.set_title("Activation pivots at G310", loc="left")
h = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ("#9EC9E2", "#CCCCCC")]
ax.legend(h, ["used by published gate metrics", "other S6 residues"], frameon=False, loc="upper left",
          handlelength=0.9, handleheight=0.9, labelspacing=0.25, fontsize=6.5, borderpad=0.1)

for ax, L in zip(fig.axes, "ABC"):
    ax.text(-0.16 if L == "A" else -0.19, 1.16, L, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
O = f"{H}/../out"
fig.savefig(f"{O}/FIG_1_structure.pdf"); fig.savefig(f"{O}/FIG_1_structure.png", dpi=600)
print("FIG_1 ok | G310 +%.2f A, S314 +%.2f A" % (dd[res == 310][0], dd[res == 314][0]))
