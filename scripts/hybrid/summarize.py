#!/usr/bin/env python3
# Independent gate water (indep_gatewater.py) -> f(H1), f(H2), for comparison with the primary read-out (hyb_readout.py)
import re, numpy as np
d = {}
for ln in open('verify_out.txt'):
    m = re.match(r'^(\S+)\t.*gate_water=([\d.]+)', ln)
    if m: d[m.group(1)] = float(m.group(2))
G = lambda rs: np.array([d[r] for r in rs])
WT, GS = G(['ANE_r1', 'ANE_r2', 'ANE_r4', 'ANE_r5']), G(['GS_r1', 'GS_r2', 'GS_r3'])
se = lambda v: v.std(ddof=1) / np.sqrt(len(v))
L = ["# Independent check: axis = centre of the G310 CA atoms (primary: filter CA atoms), frames at 21, 23, ... ns (primary: 20, 21, ... ns)",
     f"WT  {WT.mean():.2f} ± {se(WT):.2f} (n4; primary 8.73 ± 0.24)", f"GS  {GS.mean():.2f} ± {se(GS):.2f} (n3; primary 4.73 ± 0.12)"]
for nm, rs, prim in (('H1', ['H1_r1', 'H1_r2'], '0.22 (0.13-0.32)'), ('H2', ['H2_r1', 'H2_r2'], '0.58 (0.37-0.78)')):
    H = G(rs); A, B, h = WT.mean(), GS.mean(), H.mean(); den = A - B; f = (A - h) / den
    var = (((h - B) / den**2) * se(WT))**2 + ((1 / den) * se(H))**2 + (((A - h) / den**2) * se(GS))**2
    L.append(f"{nm}  {h:.2f} ({', '.join(f'{x:.2f}' for x in H)}) -> f = {f:.2f} (95 % {f-1.96*var**.5:.2f}-{f+1.96*var**.5:.2f})   | primary: {prim}")
L.append("Composition checked: H1 = SER,GLY,GLY,GLY (chain A); H2 = SER,SER,GLY,GLY (A+B), SER-SER CA 6.4-6.6 A -> adjacent (the diagonal is ~8.3 A).")
open('VERIFY_SUMMARY.txt', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
