#!/usr/bin/env python3
"""Final CHARMM36m table (SI). Usage: python3 c36_summary.py <O> > <O>/C36_FINAL_SUMMARY.txt
Input (<O>): disp_flux.txt (independent estimator), <run>_b20.books.txt (primary estimator: net filter charge, first 20 ns excluded),
mech/<run>.txt (mech_c36.py: 1 frame per ns, same window). Amber14sb reference values: WT 17/1120 ns = 4.30 pS, G310S 2/840 ns;
mechanism (WT 4 vs G310S 3 replicas): gate water 8.73 +- 0.24 vs 4.73 +- 0.12, d1 8.39 vs 8.17 A, S310 OG 3.53 A from the axis.
NET top = exits_top - entries_top (as for Amber14sb); NET bottom = entries_bottom - exits_bottom (net charge entering from the cavity);
their difference is a change in filter occupancy (for example the exit of an ion already in the filter)."""
import sys, re, datetime, numpy as np
from scipy import stats
O = sys.argv[1]; QE = 1.602176634e-19
RUNS = [('C36W_r1', 'WT'), ('C36W_r2', 'WT'), ('C36G_r1', 'GS'), ('C36G_r2', 'GS')]
AMB = {'WT': (17, 1120.0), 'GS': (2, 840.0)}
nan = float('nan')

def G(n, t_ns, v_mV): return n * QE / (t_ns * 1e-9) / (v_mV * 1e-3) * 1e12
def ci(k):  # exact (Garwood) 95 % interval
    return (stats.chi2.ppf(0.025, 2 * k) / 2 if k > 0 else 0.0), stats.chi2.ppf(0.975, 2 * k + 2) / 2

rx = re.compile(r"^(\S+)\s+window\s+([\d.]+) ns\s+V\s+([\d.]+) mV.*?SF displacement N\s+([+-]?[\d.]+).*?gate N\s+([+-]?[\d.]+)"
                r".*?transits up (\d+) \(KSF (\d+)\) down (\d+) -> net ([+-]?\d+)")
DF = {}
for l in open(f'{O}/disp_flux.txt'):
    m = rx.match(l)
    if m: DF[m.group(1)] = dict(t=float(m.group(2)), V=float(m.group(3)), disp=float(m.group(4)), gate=float(m.group(5)), full=int(m.group(9)))

def books(p):
    s = open(p).read()
    r = dict(t=float(re.search(r'####\s+\S+:\s+([\d.]+) ns', s).group(1)),
             top=int(re.search(r'\(exits_top - entries_top\) = (-?\d+)', s).group(1)),
             bot=int(re.search(r'\(entries_bottom - exits_bottom\) = (-?\d+)', s).group(1)))
    m = re.search(r'entries from cavity (\d+), from outside (\d+); exits to outside (\d+), to cavity (\d+)', s)
    r.update(ecav=int(m.group(1)), eout=int(m.group(2)), xout=int(m.group(3)), xcav=int(m.group(4)))
    for key, tag in (('cond', 'CAV->SF->OUT (conduction)'), ('pre', 't0 ->SF->OUT (pre-loaded exits)')):
        m = re.search(re.escape(tag) + r'\s+free\s+(\d+)\s+KSF\s+(\d+)', s); r[key] = int(m.group(1)) + int(m.group(2))
    return r

def mech(p):
    d = dict(q.split('=', 1) for q in open(p).read().strip().split('\t')[1:])
    num = lambda k: float(re.match(r'nan|-?[\d.]+', d[k]).group(0))
    return dict(n=int(d['n']), d1=num('d1'), og=num('OG_R'), gw=num('gate_water'), cavK=num('cav_K_frames'), belK=num('below_gate_K_frames'), cavW=num('cav_water'))

B = {r: books(f'{O}/{r}_b20.books.txt') for r, _ in RUNS}; M = {r: mech(f'{O}/mech/{r}.txt') for r, _ in RUNS}
out = []; P = out.append
P(f"# CHARMM36m FINAL — 4 runs x 300 ns, ~565 mV, 0.5 M KCl, POSRES_ENDS (same protocol as the Amber14sb runs)  [{datetime.datetime.now():%Y-%m-%d %H:%M}]")
P("# Estimator: net charge through the selectivity filter (sf_books.py, first 20 ns discarded) + independent disp_flux.py; group Analysis_min_O (filter O). Mechanism: mech_c36.py (1 frame/ns, same window).")
P("\n## 1. Transfer (first 20 ns discarded)")
P("| run | window, ns | V, mV | NET top | NET bottom | entries from the cavity | CAV→SF→OUT | pre-loaded exits | displacement N | gate N | complete transits |")
P("|---|---|---|---|---|---|---|---|---|---|---|")
for r, g in RUNS:
    b, d = B[r], DF.get(r, {})
    P(f"| {r} | {b['t']:.1f} | {d.get('V', nan):.0f} | {b['top']} | {b['bot']} | {b['ecav']} | {b['cond']} | {b['pre']} | {d.get('disp', nan):+.2f} | {d.get('gate', nan):+.2f} | {d.get('full', 0)} |")
P("\n## 2. Pooled, and comparison with Amber14sb (Amber14sb, same window: WT 17/1120 ns, GS 2/840 ns)")
summ = {}
for g in ('WT', 'GS'):
    rs = [r for r, gg in RUNS if gg == g]
    T = sum(B[r]['t'] for r in rs); V = np.mean([DF[r]['V'] for r in rs if r in DF]) if any(r in DF for r in rs) else 565.0
    top, bot, ecav = (sum(B[r][k] for r in rs) for k in ('top', 'bot', 'ecav'))
    disp = sum(DF[r]['disp'] for r in rs if r in DF); gate = sum(DF[r]['gate'] for r in rs if r in DF)
    ka, Ta = AMB[g]; lam = ka / Ta * T
    for name, k in (('NET top', top), ('NET bottom (from the cavity)', bot)):
        kk = max(k, 0); lo, hi = ci(kk)
        P(f"- {g} {name}: {k} / {T:.1f} ns -> G = {G(kk, T, V):.2f} pS (95 % {G(lo, T, V):.2f}-{G(hi, T, V):.2f}); expected at the Amber14sb rate {lam:.1f}"
          f" -> P(<= {kk}) = {stats.poisson.cdf(kk, lam):.1e}; conditional binomial (C36m < Amber) p = {stats.binom.cdf(kk, kk + ka, T / (T + Ta)):.1e}")
    P(f"- {g} displacement N = {disp:+.2f} (G = {G(disp, T, V):.2f} pS) · gate N = {gate:+.2f} · entries from the cavity into the filter {ecav} · mean V {V:.0f} mV")
    summ[g] = dict(T=T, top=top, bot=bot, ecav=ecav, disp=disp, lam=lam)
P("\n## 3. Mechanism (first 20 ns discarded, 1 frame/ns; replica mean +- SE; n = 2 vs 2 -> p is descriptive only; Amber14sb: WT 4 vs GS 3)")
P("| run | frames | d1, Å | S310 OG from the axis, Å | gate water | cavity water | K⁺ in the cavity, % frames | K⁺ below the gate, % frames |")
P("|---|---|---|---|---|---|---|---|")
for r, g in RUNS:
    m = M[r]; P(f"| {r} | {m['n']} | {m['d1']:.2f} | {m['og']:.2f} | {m['gw']:.2f} | {m['cavW']:.1f} | {m['cavK']:.1f} | {m['belK']:.1f} |")
def grp(key, g):
    v = np.array([M[r][key] for r, gg in RUNS if gg == g]); return v.mean(), v.std(ddof=1) / np.sqrt(len(v)), v
for key, lab in (('gw', 'gate water'), ('d1', 'd1, Å'), ('cavW', 'cavity water'), ('cavK', 'K⁺ in the cavity, % frames'), ('belK', 'K⁺ below the gate, % frames')):
    (mw, sw, vw), (mg, sg, vg) = grp(key, 'WT'), grp(key, 'GS')
    P(f"- {lab}: WT {mw:.2f} +- {sw:.2f} vs GS {mg:.2f} +- {sg:.2f} (Welch p = {stats.ttest_ind(vw, vg, equal_var=False).pvalue:.2g})")
og = grp('og', 'GS'); rw = grp('gw', 'GS')[0] / grp('gw', 'WT')[0]
P(f"- S310 OG from the axis: {og[0]:.2f} +- {og[1]:.2f} Å (Amber 3.53) · gate water GS/WT: C36m {rw:.2f} · Amber14sb {4.73 / 8.73:.2f}")
w, s = summ['WT'], summ['GS']
P(f"\nSUMMARY: C36m WT NET top/bottom {w['top']}/{w['bot']} ({w['T']:.0f} ns, entries from the cavity {w['ecav']}, displacement {w['disp']:+.2f}; expected at the Amber14sb rate {w['lam']:.1f})"
  f" · GS {s['top']}/{s['bot']} ({s['T']:.0f} ns, entries {s['ecav']}) · gate water WT {grp('gw', 'WT')[0]:.2f} vs GS {grp('gw', 'GS')[0]:.2f}"
  f" (GS/WT {rw:.2f}; Amber 0.54) · S310 OG {og[0]:.2f} Å")
print('\n'.join(out))
