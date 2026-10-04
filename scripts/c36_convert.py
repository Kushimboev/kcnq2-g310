#!/usr/bin/env python3
"""Convert the Amber system (Amber14sb + Lipid21, full charges) to CHARMM36m. Coordinates are kept; only the force field changes.
Uses the official port charmm36-jul2021.ff (MacKerell lab), checked against the CHARMM-GUI toppar files
(POPC: 134 atom types and charges, no differences; POT/CLA and NBFIX identical).
usage: c36_convert.py AMBER_SYS_DIR IN.gro REF.gro OUTDIR NAME
  AMBER_SYS_DIR: topol.top, A_prot_Protein_chain_[A-D].itp, POPC.itp, posre_{bb,ends,s5}.itp
  IN.gro: starting frame (equilibrated Amber system); REF.gro: restraint reference coordinates (-r; the Amber *_system.gro)
Steps: protein heavy atoms -> PDB (HIE->HSE, HID->HSD, HIP->HSP, ILE CD1->CD) -> pdb2gmx -ignh, termini None/None (ACE/NME);
POPC: Lipid21 (PA/PC/OL) -> CHARMM POPC by graph isomorphism (by element, hydrogens included); coordinates are transferred;
SOL unchanged; K -> POT, CL -> CLA; the ion order is kept (the last three POT are the ions placed in the filter).
Output: OUTDIR/{NAME.gro, NAME_ref.gro, topol.top, prot_Protein_chain_X.itp (+POSRES_BB/ENDS/S5), POPC_c36.itp (+POSRES_LIPID), index.ndx}
"""
import sys, os, re, subprocess, shutil
import networkx as nx
from networkx.algorithms import isomorphism as iso
GMX = '/usr/local/gromacs-2025.2/bin/gmx'
FF = '${KCNQ2_ROOT}/13_conduction/build_2026-09-25_c36/ff/charmm36-jul2021.ff'
amb, in_gro, ref_gro, out, name = sys.argv[1:6]
os.makedirs(out, exist_ok=True)
if not os.path.exists(f'{out}/charmm36-jul2021.ff'): os.symlink(FF, f'{out}/charmm36-jul2021.ff')
log = open(f'{out}/convert.log', 'a')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); log.write(s + '\n'); log.flush()

def read_gro(p):
    L = open(p).read().splitlines(); n = int(L[1]); A = []
    for l in L[2:2 + n]:
        A.append(dict(resnr=int(l[0:5]), resn=l[5:10].strip(), name=l[10:15].strip(), x=float(l[20:28]), y=float(l[28:36]), z=float(l[36:44])))
    return A, L[2 + n]
def section(p, sec):
    out, f = [], False
    for l in open(p):
        s = l.split(';')[0].strip()
        if s.startswith('['): f = (s.replace(' ', '') == f'[{sec}]'); continue
        if f and s and not s.startswith('#'): out.append(s.split())
    return out
def itp_atoms(p): return [dict(nr=int(t[0]), type=t[1], resnr=int(t[2]), resn=t[3], name=t[4]) for t in section(p, 'atoms')]
def rtp_block(p, res):
    s = open(p).read(); i = s.index(f'[ {res} ]'); j = s.find('\n[ ', i + 5)
    blk = s[i:j if j > 0 else None]; secs = re.split(r'\n\s*\[\s*(\w+)\s*\]', blk)
    d = {secs[k]: [l.split(';')[0].split() for l in secs[k + 1].splitlines() if l.split(';')[0].strip()] for k in range(1, len(secs) - 1, 2)}
    return d
mols = [(t[0], int(t[1])) for t in section(f'{amb}/topol.top', 'molecules')]
nat = {'POPC': 134, 'SOL': 3, 'K': 1, 'CL': 1}
chains = [m for m, _ in mols if m.startswith('Protein_chain_')]
for c in chains: nat[c] = len(itp_atoms(f'{amb}/A_prot_{c}.itp'))
say(f'=== {name}: molecules {mols}; chain atoms {[nat[c] for c in chains]}')

# ---- POPC mapping (Lipid21 -> CHARMM36) by graph isomorphism
la = itp_atoms(f'{amb}/POPC.itp'); lb = [(int(t[0]), int(t[1])) for t in section(f'{amb}/POPC.itp', 'bonds')]
rc = rtp_block(f'{FF}/lipid.rtp', 'POPC'); cn = [t[0] for t in rc['atoms']]; cb = [(t[0], t[1]) for t in rc['bonds']]
ga = nx.Graph(); [ga.add_node(a['nr'], el=a['name'][0]) for a in la]; ga.add_edges_from(lb)
gc = nx.Graph(); [gc.add_node(n, el=n[0]) for n in cn]; gc.add_edges_from(cb)
assert ga.number_of_nodes() == gc.number_of_nodes() == 134 and ga.number_of_edges() == gc.number_of_edges(), (ga.number_of_nodes(), gc.number_of_nodes(), ga.number_of_edges(), gc.number_of_edges())
gm = iso.GraphMatcher(ga, gc, node_match=iso.categorical_node_match('el', None)); amap = next(gm.isomorphisms_iter())   # amber nr -> charmm name
inv = {v: k for k, v in amap.items()}; lip_order = [inv[n] for n in cn]                                                  # amber nr in CHARMM order
say(f'POPC mapping: 134/134 atoms, {ga.number_of_edges()} bonds; P31 -> {amap[[a["nr"] for a in la if a["name"]=="P31"][0]]}')

RN = {'HIE': 'HSE', 'HID': 'HSD', 'HIP': 'HSP'}
def convert(gro_path, tag):
    A, box = read_gro(gro_path); assert len(A) == sum(nat[m] * n for m, n in mols), (len(A), sum(nat[m] * n for m, n in mols))
    i = 0; blocks = []
    for m, n in mols:
        for _ in range(n): blocks.append((m, A[i:i + nat[m]])); i += nat[m]
    # protein -> PDB (heavy atoms)
    pdb = [];  heavy = {}; k = 0
    for ci, (m, at) in enumerate([b for b in blocks if b[0].startswith('Protein_chain_')]):
        ch = m[-1]
        for a in at:
            if a['name'].startswith('H'): continue
            rn = RN.get(a['resn'], a['resn']); nm = 'CD' if (a['resn'] == 'ILE' and a['name'] == 'CD1') else a['name']
            k += 1; el = nm[0]; nmf = nm if len(nm) == 4 else ' ' + nm
            pdb.append(f"ATOM  {k % 100000:5d} {nmf:<4s} {rn:>3s} {ch}{a['resnr']:4d}    {a['x']*10:8.3f}{a['y']*10:8.3f}{a['z']*10:8.3f}  1.00  0.00          {el:>2s}\n")
            heavy[(ch, a['resnr'], nm)] = (a['x'], a['y'], a['z'])
        pdb.append('TER\n')
    open(f'{out}/prot_{tag}.pdb', 'w').write(''.join(pdb) + 'END\n')
    r = subprocess.run([GMX, 'pdb2gmx', '-f', f'prot_{tag}.pdb', '-o', f'prot_{tag}.gro', '-p', f'prot_{tag}.top', '-i', f'posre_{tag}.itp', '-ff', 'charmm36-jul2021',
                        '-water', 'none', '-ignh', '-ter', '-chainsep', 'id'], input='6\n5\n' * len(chains), text=True, cwd=out, capture_output=True)
    open(f'{out}/pdb2gmx_{tag}.log', 'w').write(r.stdout + r.stderr)
    assert r.returncode == 0, f'pdb2gmx failed ({tag}) — see pdb2gmx_{tag}.log'
    terms = re.findall(r'(Start|End) terminus [A-Z]+-\d+: (\S+)', r.stdout + r.stderr); say(f'pdb2gmx {tag}: termini {sorted(set(t[1] for t in terms))} ({len(terms)} ta)')
    P, _ = read_gro(f'{out}/prot_{tag}.gro')
    # check that the heavy atoms did not move
    dmax = 0.0; ch_idx = 0; prev = None   # a new chain starts where resnr restarts
    for a in P:
        if prev is not None and a['resnr'] < prev: ch_idx += 1
        prev = a['resnr']; a['ch'] = 'ABCD'[ch_idx]
        if not a['name'].startswith('H'):
            key = (a['ch'], a['resnr'], a['name'])
            if key in heavy: x = heavy[key]; dmax = max(dmax, abs(x[0]-a['x']), abs(x[1]-a['y']), abs(x[2]-a['z']))
    say(f'protein {tag}: {len(P)} atoms (Amber {sum(nat[c] for c in chains)}), max. heavy-atom displacement {dmax*10:.4f} A')
    # assemble
    lines = []; resn_c = 0
    def add(resnr, resn, an, x, y, z): lines.append((resnr, resn, an, x, y, z))
    for a in P: add(a['resnr'], a['resn'], a['name'], a['x'], a['y'], a['z'])
    rid = max(a['resnr'] for a in P)
    for m, at in blocks:
        if m.startswith('Protein_chain_'): continue
        rid += 1
        if m == 'POPC':
            for nr, cnm in zip(lip_order, cn): a = at[nr - 1]; add(rid, 'POPC', cnm, a['x'], a['y'], a['z'])
        elif m == 'SOL':
            for a in at: add(rid, 'SOL', a['name'], a['x'], a['y'], a['z'])
        elif m == 'K': a = at[0]; add(rid, 'POT', 'POT', a['x'], a['y'], a['z'])
        elif m == 'CL': a = at[0]; add(rid, 'CLA', 'CLA', a['x'], a['y'], a['z'])
    fn = f'{out}/{name}.gro' if tag == 'in' else f'{out}/{name}_ref.gro'
    with open(fn, 'w') as f:
        f.write(f'{name} CHARMM36m (jul2021), converted from the Amber14sb system ({tag})\n{len(lines)}\n')
        for j, (rr, rn, an, x, y, z) in enumerate(lines, 1): f.write(f'{rr % 100000:5d}{rn:<5s}{an:>5s}{j % 100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n')
        f.write(box + '\n')
    say(f'{fn}: {len(lines)} atom (Amber {len(A)})')
    return P, lines

Pin, Lin = convert(in_gro, 'in'); Pref, Lref = convert(ref_gro, 'ref')
assert len(Lin) == len(Lref) and all(a[2] == b[2] for a, b in zip(Lin, Lref)), 'atom order differs between in and ref'
# ---- topology: chain itp files (+ restraint blocks)
for c in chains:
    shutil.copy(f'{out}/prot_in_{c}.itp', f'{out}/prot_{c}.itp'); shutil.copy(f'{out}/posre_in_{c}.itp', f'{out}/posre_{c}.itp')
    s = open(f'{out}/prot_{c}.itp').read().replace(f'posre_in_{c}.itp', f'posre_{c}.itp'); open(f'{out}/prot_{c}.itp', 'w').write(s)
    Aa = itp_atoms(f'{amb}/A_prot_{c}.itp'); Ca = itp_atoms(f'{out}/prot_{c}.itp'); cidx = {(a['resnr'], a['name']): a['nr'] for a in Ca}
    add = ''
    for pr in ('bb', 'ends', 's5'):
        src = f'{amb}/posre_{pr}.itp'
        if not os.path.exists(src): continue
        rows = [t for t in section(src, 'position_restraints')]; new = []
        for t in rows:
            a = Aa[int(t[0]) - 1]; key = (a['resnr'], 'CD' if (a['resn'] == 'ILE' and a['name'] == 'CD1') else a['name'])
            assert key in cidx, (c, pr, key); new.append(f'{cidx[key]:6d} {" ".join(t[1:])}\n')
        open(f'{out}/posre_{pr}_{c}.itp', 'w').write('[ position_restraints ]\n' + ''.join(new))
        add += f'\n#ifdef POSRES_{pr.upper()}\n#include "posre_{pr}_{c}.itp"\n#endif\n'
        if c == chains[0]: say(f'posre_{pr}: {len(new)} atom (Amber {len(rows)})')
    open(f'{out}/prot_{c}.itp', 'a').write(add)
# ---- POPC itp: pdb2gmx on a single lipid (all atoms)
lp = [];  L0 = [l for l in Lin if l[1] == 'POPC'][:134]
for j, (rr, rn, an, x, y, z) in enumerate(L0, 1):
    nmf = an if len(an) == 4 else ' ' + an
    lp.append(f"ATOM  {j:5d} {nmf:<4s} POPCL{1:4d}    {x*10:8.3f}{y*10:8.3f}{z*10:8.3f}  1.00  0.00          {an[0]:>2s}\n")
open(f'{out}/popc1.pdb', 'w').write(''.join(lp) + 'END\n')
r = subprocess.run([GMX, 'pdb2gmx', '-f', 'popc1.pdb', '-o', 'popc1.gro', '-p', 'popc1.top', '-i', 'posre_popc1.itp', '-ff', 'charmm36-jul2021', '-water', 'none'],
                   text=True, cwd=out, capture_output=True, input='')
open(f'{out}/pdb2gmx_popc.log', 'w').write(r.stdout + r.stderr); assert r.returncode == 0, 'pdb2gmx failed for POPC'
t = open(f'{out}/popc1.top').read(); i = t.index('[ moleculetype ]'); j = t.find('; Include Position restraint', i); j = j if j > 0 else t.find('[ system ]', i)
blk = t[i:j]; blk = re.sub(r'(\[ moleculetype \]\s*\n;[^\n]*\n)\s*\S+', r'\1POPC', blk, count=1)
pidx = [a['nr'] for a in itp_atoms(f'{out}/popc1.top') if a['name'] == 'P'][0]
open(f'{out}/POPC_c36.itp', 'w').write(f'; CHARMM36 POPC (charmm36-jul2021 lipid.rtp, pdb2gmx) — 2026-09-25\n{blk}\n#ifdef POSRES_LIPID\n[ position_restraints ]\n{pidx:6d}     1     0     0  1000\n#endif\n')
say(f'POPC_c36.itp: P indeksi {pidx}, atomlar {len(itp_atoms(out + "/POPC_c36.itp"))}')
# ---- topol.top
molsc = [(('POT' if m == 'K' else 'CLA' if m == 'CL' else m), n) for m, n in mols]
top = ['; CHARMM36m (charmm36-jul2021, official port), converted from the Amber14sb system (same coordinates)', '#include "charmm36-jul2021.ff/forcefield.itp"']
top += [f'#include "prot_{c}.itp"' for c in chains] + ['#include "POPC_c36.itp"', '#include "charmm36-jul2021.ff/tip3p.itp"', '#include "charmm36-jul2021.ff/ions.itp"', '', '[ system ]', name, '', '[ molecules ]']
top += [f'{m:20s} {n}' for m, n in molsc]
open(f'{out}/topol.top', 'w').write('\n'.join(top) + '\n')
# ---- index.ndx
L = Lin; N = len(L); idx = lambda f: [j + 1 for j, a in enumerate(L) if f(j, a)]
nprot = len(Pin); pot = idx(lambda j, a: a[1] == 'POT')
chain_of = []; ch = 0; prev = None
for a in Pin:
    if prev is not None and a['resnr'] < prev: ch += 1
    prev = a['resnr']; chain_of.append('ABCD'[ch])
G = {}
G['System'] = list(range(1, N + 1)); G['Protein'] = list(range(1, nprot + 1)); G['POPC'] = idx(lambda j, a: a[1] == 'POPC')
G['SOL'] = idx(lambda j, a: a[1] == 'SOL'); G['POT'] = pot; G['CLA'] = idx(lambda j, a: a[1] == 'CLA'); G['Filter_K'] = pot[-3:]
G['Protein_POPC'] = G['Protein'] + G['POPC']; G['Water_and_ions'] = G['SOL'] + G['POT'] + G['CLA']
G['SF_O'] = [j + 1 for j, a in enumerate(Pin) if 276 <= a['resnr'] <= 281 and a['name'] == 'O']
for c in 'ABCD': G[f'G310_{c}'] = [j + 1 for j, a in enumerate(Pin) if a['resnr'] == 310 and chain_of[j] == c]
G['Analysis_min'] = [j + 1 for j, a in enumerate(Pin) if ((276 <= a['resnr'] <= 281) or a['resnr'] == 310) and a['name'] == 'CA'] + pot
with open(f'{out}/index.ndx', 'w') as f:
    for k, v in G.items():
        f.write(f'[ {k} ]\n'); [f.write(' '.join(f'{x:6d}' for x in v[i:i + 15]) + '\n') for i in range(0, len(v), 15)]
say('index:', {k: len(v) for k, v in G.items()})
say(f'=== {name} tayyor')
