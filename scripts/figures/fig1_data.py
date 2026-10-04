#!/usr/bin/env python3
"""Data for Fig. 1: CA-CA diagonals along S6 and the HOLE radius profile in 8J01 (activated) and 7CR3 (closed).
Diagonal: the CA atoms of the four chains form two opposite pairs, the narrow (d1) and the wide (d2) one. Output: ../out/fig1_data.npz and a text summary."""
import sys, os, numpy as np
P = "${KCNQ2_ROOT}"
G = f"{P}/4_gate_analysis"
def ca(pdb, resids):
    """CA coordinates from a PDB file: {resid: [4 chains x 3]} (the four chains of the pore tetramer)."""
    out = {}
    for ln in open(pdb):
        if not ln.startswith("ATOM") or ln[12:16].strip() != "CA": continue
        try: r = int(ln[22:26])
        except ValueError: continue
        if r not in resids: continue
        out.setdefault(r, []).append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return {r: np.array(v) for r, v in out.items() if len(v) >= 4}
def diagonals(X):
    """4 x 3 -> (narrow, wide) diagonal (the two most distant pairs)."""
    X = X[:4]; c = X.mean(0); ang = np.arctan2(X[:, 1] - c[1], X[:, 0] - c[0]); X = X[np.argsort(ang)]
    d = [np.linalg.norm(X[0] - X[2]), np.linalg.norm(X[1] - X[3])]
    return min(d), max(d)
RES = list(range(300, 320))
A = ca(f"{G}/pdbs/8j01.pdb", set(RES)); B = ca(f"{G}/pdbs/7cr3.pdb", set(RES))
rows = []
for r in RES:
    if r in A and r in B:
        a1, a2 = diagonals(A[r]); b1, b2 = diagonals(B[r])
        rows.append((r, a1, a2, b1, b2, a1 - b1))
rows = np.array(rows)
prof = {}
for k, f in (("8J01", f"{G}/out/refval_8j01_profile.dat"), ("7CR3", f"{G}/out/refval_7cr3_profile.dat")):
    d = np.loadtxt(f); z, rad = d[:, 0], d[:, 1] * 10.0     # nm -> A
    o = np.argsort(z); z, rad = z[o], rad[o]
    zu, idx = np.unique(z, return_index=True)
    prof[k] = np.column_stack([zu, np.minimum.reduceat(rad, idx) if False else np.array([rad[z == t].min() for t in zu])])
np.savez(f"{os.path.dirname(os.path.abspath(__file__))}/../out/fig1_data.npz", diag=rows, p8j01=prof["8J01"], p7cr3=prof["7CR3"])
print(" res   8J01 d1/d2      7CR3 d1/d2     open by (d1)")
for r, a1, a2, b1, b2, dd in rows:
    tag = "  <-- PAG G310" if r == 310 else ("  <-- S314" if r == 314 else "")
    print(f" {int(r)}   {a1:5.2f} {a2:5.2f}   {b1:5.2f} {b2:5.2f}   {dd:+5.2f}{tag}")
