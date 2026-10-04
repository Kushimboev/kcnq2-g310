#!/usr/bin/env python3
"""Fig. 3, how ions cross. A: z(t) traces of the K+ ions in the wild type (knock-on). B: the same in G310S (ions reach the gate but few pass).
C: exchange of coordination near the gate (water -> Ser310 OG)."""
import sys, os, re, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../style")
import paperdata as D
from palette import C, MM, use_paper_style
use_paper_style()
H = os.path.dirname(os.path.abspath(__file__))
TR = f"{H}/../../13_conduction/scripts/analysis_2026-09-28_paper/out"
WTC, GSC = "#1B7EBD", "#D95F02"

def traces(name):
    z = np.load(f"{TR}/{name}.trace.npz"); return z["t"], z["z"], z["r"], float(z["zg"].mean())

fig = plt.figure(figsize=(MM.TWO_COL, 2.9))
gs = fig.add_gridspec(1, 3, width_ratios=[1.18, 1.18, 1.0], wspace=0.34, left=0.062, right=0.985, bottom=0.155, top=0.86)

def crossers(run):
    """Indices of the ions marked 'went OUT' by sf_books (same order as in the K selection)."""
    out = []
    for ln in open(f"{D.BOOKS}/{run}_b20.books.txt"):
        m = re.match(r"\s+ion\s+(\d+)\s+(\S+)\s+in SF\s+([\d.]+)-\s*([\d.]+) ns\s+came (\S+)\s+went (\S+)", ln)
        if m and m.group(6) == "OUT": out.append((int(m.group(1)), m.group(2) == "KSF", m.group(5)))
    return out

for k, (name, ttl, col) in enumerate([("ANE_r1", "WT: ions cross by knock-on", WTC),
                                      ("GS_r2", "G310S: few ions pass the gate", GSC)]):
    ax = fig.add_subplot(gs[k])
    t, Z, R, zg = traces(name)
    m = t >= 20.0                                          # analysis window: first 20 ns excluded
    t, Z, R = t[m] - 20.0, Z[m], R[m]
    lo, hi = zg - 0.55, 1.15
    ax.axhspan(-0.8, 0.45, color=C.HILITE, alpha=0.30, lw=0)     # filter
    ax.axhline(zg, color=C.GLY310, lw=1.0, ls="--")
    ax.text(3, zg - 0.05, "G310", fontsize=6.6, color=C.GLY310, ha="left", va="top", fontweight="bold")
    ax.text(259, -0.18, "filter", fontsize=6.6, color="#B8860B", ha="right", va="center")
    ax.text(259, (zg + (-0.8)) / 2, "cavity", fontsize=6.6, color="#999999", ha="right", va="center")
    on = R < 0.6
    cross = {i: k for i, k, c in crossers(name)}
    for j in range(Z.shape[1]):
        zz = np.where(on[:, j] & (Z[:, j] > lo) & (Z[:, j] < hi), Z[:, j], np.nan)
        if np.isfinite(zz).sum() < 30: continue
        if j in cross:
            ax.plot(t, zz, lw=1.1, color=col, alpha=0.95, zorder=3)
        else:
            ax.plot(t, zz, lw=0.6, color="#BBBBBB", alpha=0.8, zorder=1)
    ax.set_xlim(0, 262); ax.set_ylim(lo, hi)
    ax.set_xlabel("Time (ns)")
    if k == 0: ax.set_ylabel("z relative to the filter (nm)")
    ax.set_title(ttl, loc="left", fontsize=8.2)
    npre = sum(1 for _, isksf, came in crossers(name) if isksf or came in ("t0", "?"))     # exits of ions already in the filter at t = 0
    n = D.net(name); lab = f"{n} net crossing" + ("s" if n != 1 else "")
    if npre and n > 1: lab += ";\n" + ("both" if npre == n == 2 else "all" if npre == n else str(npre)) + " by filter ions present at $t$ = 0"
    elif npre: lab += "; by a filter ion present at $t$ = 0"
    ax.text(0.015, 0.975, lab, transform=ax.transAxes, fontsize=6.6, color=col, ha="left", va="top", linespacing=1.3,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=1.2))
    if k == 0:
        ax.plot([], [], color=col, lw=1.1, label="ion that exits")
        ax.plot([], [], color="#BBBBBB", lw=0.7, label="other on-axis K$^+$")
        ax.legend(frameon=True, facecolor="white", framealpha=0.85, edgecolor="none", loc="lower right", fontsize=6.6, handlelength=1.3, labelspacing=0.2, borderpad=0.2)

# --- C: exchange of coordination (another view of Fig. 5D: fractions) ---
ax = fig.add_subplot(gs[2])
z, W, G = D.hyd_profile()
nz = lambda A, j: np.nanmean(A[:, :, j], axis=0)
gw, go, gp = nz(G, 1), nz(G, 2), nz(G, 3)
tot = gw + go + gp
ok = (np.nansum(G[:, :, 0], axis=0) >= 30) & (z <= 0.05)     # 30 ion-frames per bin; nothing is drawn where the axis is empty
ax.stackplot(z[ok], gw[ok] / tot[ok] * 100, go[ok] / tot[ok] * 100, gp[ok] / tot[ok] * 100,
             colors=[C.WATER, GSC, "#7570B3"], labels=["water", "Ser310 O$\\gamma$", "other protein"])
ax.axvline(0, color="#666666", lw=0.9, ls="--")
ax.set_xlabel("z relative to the G310 plane (nm)"); ax.set_ylabel("Share of the K$^+$ first shell (%)")
ax.set_xlim(-1.0, 0.05); ax.set_ylim(0, 100)
ax.legend(frameon=True, facecolor="white", framealpha=0.92, edgecolor="none", loc="lower left", fontsize=6.6, handlelength=0.9, labelspacing=0.2, borderpad=0.3)
ax.set_title("G310S: what surrounds the ion", loc="left", fontsize=8.2)

for ax, L in zip(fig.axes, "ABC"):
    ax.text(-0.15, 1.15, L, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
O = f"{H}/../out"
fig.savefig(f"{O}/FIG_3_mechanism.pdf"); fig.savefig(f"{O}/FIG_3_mechanism.png", dpi=600)
print("FIG_3 ok")
