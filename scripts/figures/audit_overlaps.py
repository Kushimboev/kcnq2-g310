#!/usr/bin/env python3
"""Automatic check of the figures for overlapping text.
Runs every figure script with savefig intercepted and, once the figure is drawn, measures all text objects (titles, axis labels, ticks, annotations, legends):
(1) text-text overlap and text outside the figure are reported as errors; (2) a data line or bar under a text is flagged for visual inspection.
usage: python3 audit_overlaps.py [fig*.py ...]   ->  ../out/AUDIT_overlaps.txt"""
import sys, os, runpy, io, contextlib, itertools, numpy as np, matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.text import Text

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/../out/AUDIT_overlaps.txt"
SCRIPTS = sys.argv[1:] or ["fig1_structure.py", "fig2_wt.py", "fig3_mechanism.py", "fig4_controls.py", "fig5_g310s.py", "fig6_ecc.py",
                           "figS1_estimators.py", "figS2_c36m.py", "figS3_hybrid.py"]
_orig, LOG, SEEN = Figure.savefig, [], set()

def to_lab(rgb):
    c = np.array(matplotlib.colors.to_rgb(rgb), float); c = np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    M = np.array([[0.4124564, 0.3575761, 0.1804375], [0.2126729, 0.7151522, 0.0721750], [0.0193339, 0.1191920, 0.9503041]])
    xyz = M @ c / np.array([0.95047, 1.0, 1.08883]); f = np.where(xyz > 0.008856, np.cbrt(xyz), 7.787 * xyz + 16 / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])

def similar_colors(ax, thr=25.0):
    """Pairs of line / bar / marker colours within one axes with CIE76 dE < thr (a reader cannot tell them apart)."""
    cols = {}
    grid = set(ax.xaxis.get_gridlines()) | set(ax.yaxis.get_gridlines())
    for ln in ax.lines:
        if ln in grid or not ln.get_visible() or len(ln.get_xdata()) < 2: continue
        if ln.get_linestyle() in ("None", "", " ") and ln.get_marker() in ("None", "", None): continue
        cols.setdefault(matplotlib.colors.to_hex(ln.get_color()), str(ln.get_label())[:18] if not str(ln.get_label()).startswith("_") else "line")
    for p in ax.patches:
        if isinstance(p, Rectangle) and p is not ax.patch and p.get_visible() and p.get_width() and p.get_height():
            fc = p.get_facecolor()
            if fc[3] > 0.2: cols.setdefault(matplotlib.colors.to_hex(fc[:3]), "bar")
    keys = [k for k in cols if k not in ("#000000", "#ffffff") and abs(to_lab(k)[1]) + abs(to_lab(k)[2]) > 8]      # neutral (grey/black) colours are excluded
    out = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            d = float(np.linalg.norm(to_lab(keys[i]) - to_lab(keys[j])))
            if d < thr: out.append(f'SIMILAR-COLOR  {keys[i]} ("{cols[keys[i]]}") vs {keys[j]} ("{cols[keys[j]]}")  dE={d:.0f}')
    return out

def inter(a, b):
    w = min(a.x1, b.x1) - max(a.x0, b.x0); h = min(a.y1, b.y1) - max(a.y0, b.y0)
    return max(w, 0) * max(h, 0)

def gather(fig, R):
    items = []
    def add(t, kind, ai):
        s = t.get_text()
        if not s or not s.strip() or not t.get_visible(): return
        bb = Text.get_window_extent(t, R)          # the text only (without the annotation arrow)
        if bb.width <= 1 or bb.height <= 1: return
        items.append(dict(t=t, kind=kind, ax=ai, bb=bb, s=s.replace("\n", " / ")[:46], rot=float(t.get_rotation()) % 90))
    for i, ax in enumerate(fig.axes):
        for tt in (ax.title, ax._left_title, ax._right_title): add(tt, "title", i)
        for t in ax.texts: add(t, "text", i)
        lg = ax.get_legend()
        if lg:
            for t in lg.get_texts(): add(t, "legend", i)
        if not ax.axison: continue
        add(ax.xaxis.label, "xlabel", i); add(ax.yaxis.label, "ylabel", i)
        for axis, lim in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            lo, hi = min(lim), max(lim)
            for tk in axis.get_major_ticks():
                if lo - 1e-9 <= tk.get_loc() <= hi + 1e-9: add(tk.label1, "tick", i)
    for t in fig.texts: add(t, "figtext", None)
    return items

def geoms(fig):
    G = []
    for i, ax in enumerate(fig.axes):
        grid = set(ax.xaxis.get_gridlines()) | set(ax.yaxis.get_gridlines())
        for ln in ax.lines:
            if ln in grid or not ln.get_visible() or len(ln.get_xdata()) < 2: continue
            xy = ln.get_transform().transform(np.column_stack([np.asarray(ln.get_xdata(), float), np.asarray(ln.get_ydata(), float)]))
            xy = xy[np.isfinite(xy).all(1)]
            if len(xy) < 2: continue
            pts = np.vstack([xy[:-1] + (xy[1:] - xy[:-1]) * f for f in (0, .25, .5, .75)] + [xy[-1:]])
            G.append(dict(kind="line", ax=i, pts=pts, lab=str(ln.get_label())[:24], ls=ln.get_linestyle(), lw=ln.get_linewidth()))
        for p in ax.patches:
            if isinstance(p, Rectangle) and p.get_visible() and p is not ax.patch and p.get_width() != 0 and p.get_height() != 0:
                G.append(dict(kind="bar", ax=i, bb=p.get_window_extent(), lab="bar"))
    return G

def check(fig, name):
    fig.set_dpi(600)                                # same dpi as when saving: at low dpi the font sizes differ by 5-10 %
    fig.canvas.draw(); R = fig.canvas.get_renderer(); W, Hh = fig.canvas.get_width_height()
    dpi = fig.dpi; pt2 = (dpi / 72.0) ** 2; items = gather(fig, R); G = geoms(fig); err, rev = [], []
    for a, b in itertools.combinations(items, 2):
        if a["t"] is b["t"]: continue
        if a["kind"] == b["kind"] == "tick" and (a["rot"] or b["rot"]): continue      # rotated tick: its bbox gives a false overlap
        ov = inter(a["bb"], b["bb"]) / pt2; small = min(a["bb"].width * a["bb"].height, b["bb"].width * b["bb"].height) / pt2
        if ov > 1.0 and ov / max(small, 1e-9) > 0.04:
            err.append(f'TEXT-TEXT   {a["kind"]}[{a["ax"]}] "{a["s"]}"  x  {b["kind"]}[{b["ax"]}] "{b["s"]}"   ({ov:.1f} pt2, {100 * ov / small:.0f} % of the smaller)')
    for it in items:
        bb = it["bb"]
        if matplotlib.rcParams["savefig.bbox"] != "tight" and (bb.x0 < -1 or bb.y0 < -1 or bb.x1 > W + 1 or bb.y1 > Hh + 1):   # with tight saving, text outside the canvas is not cut off
            err.append(f'OUTSIDE     {it["kind"]}[{it["ax"]}] "{it["s"]}" leaves the canvas (x {bb.x0:.0f}..{bb.x1:.0f} of {W}, y {bb.y0:.0f}..{bb.y1:.0f} of {Hh})')
        if it["kind"] in ("tick", "xlabel", "ylabel"): continue
        for g in G:
            if g["kind"] == "bar":
                ov = inter(bb, g["bb"]) / pt2
                if ov > 2.0: rev.append(f'TEXT-BAR    {it["kind"]}[{it["ax"]}] "{it["s"]}"  over a bar of axes [{g["ax"]}] ({ov:.0f} pt2)')
            else:
                inside = (g["pts"][:, 0] > bb.x0 + 1) & (g["pts"][:, 0] < bb.x1 - 1) & (g["pts"][:, 1] > bb.y0 + 1) & (g["pts"][:, 1] < bb.y1 - 1)
                if inside.sum() >= 2: rev.append(f'TEXT-LINE   {it["kind"]}[{it["ax"]}] "{it["s"]}"  crosses line "{g["lab"]}" of axes [{g["ax"]}] ({inside.sum()} pts)')
    for i, ax in enumerate(fig.axes):                      # clipped data: points/bars outside the axis limits
        if not ax.axison: continue
        x0, x1 = sorted(ax.get_xlim()); y0, y1 = sorted(ax.get_ylim()); dx, dy = (x1 - x0) * 1e-3, (y1 - y0) * 1e-3
        grid = set(ax.xaxis.get_gridlines()) | set(ax.yaxis.get_gridlines())
        for ln in ax.lines:
            if ln in grid or not ln.get_visible(): continue
            xd, yd = np.asarray(ln.get_xdata(), float), np.asarray(ln.get_ydata(), float)
            if xd.size < 2 or xd.size != yd.size or ln.get_transform() != ax.transData: continue
            if xd.size >= 1000: continue                      # time series (traces) are windowed on purpose
            out = ((xd >= x0 - dx) & (xd <= x1 + dx) & ((yd < y0 - dy) | (yd > y1 + dy))) & np.isfinite(yd)      # clipping in y only within the visible x range
            if out.any() and out.sum() < xd.size: rev.append(f"DATA-CLIP   axes[{i}] line \"{str(ln.get_label())[:22]}\": {int(out.sum())} of {xd.size} points outside the axis limits")
        for p in ax.patches:
            if isinstance(p, Rectangle) and p.get_visible() and p is not ax.patch and p.get_transform() == ax.transData:
                if p.get_y() + p.get_height() > y1 + dy or p.get_x() + p.get_width() > x1 + dx: rev.append(f"DATA-CLIP   axes[{i}] bar (x={p.get_x():.2f}, top={p.get_y() + p.get_height():.2f}) beyond the axis limits")
    for i, ax in enumerate(fig.axes):
        if ax.axison: rev += [f'axes[{i}] ' + x for x in similar_colors(ax)]
    LOG.append((name, err, sorted(set(rev)), len(items)))

def patched(self, *a, **k):
    if id(self) not in SEEN:
        SEEN.add(id(self)); check(self, CURRENT[0])
    return _orig(self, *a, **k)
Figure.savefig = patched
CURRENT = [""]
for sc in SCRIPTS:
    CURRENT[0] = sc
    with contextlib.redirect_stdout(io.StringIO()):
        try: runpy.run_path(f"{HERE}/{sc}", run_name="__main__")
        except SystemExit: pass
    import matplotlib.pyplot as plt; plt.close("all")
lines = ["# Text overlaps in the figures (audit_overlaps.py). ERROR = must be fixed; REVIEW = to be checked by eye and accepted.", ""]
for name, err, rev, n in LOG:
    lines.append(f"## {name}  ({n} text objects) — ERROR: {len(err)} · review: {len(rev)}")
    lines += ["  " + e for e in err] + ["  ? " + r for r in rev] + [""]
open(OUT, "w", encoding="utf-8").write("\n".join(lines))
print("\n".join(l for l in lines if l.startswith("##")))
