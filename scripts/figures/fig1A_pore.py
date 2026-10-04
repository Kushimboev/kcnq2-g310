# PyMOL: Fig. 1A, the activated KCNQ2 pore (8J01), two opposite subunits; filter, G310 and S314 marked.
# Channel chains A, B, D, G (C, E, F, H are calmodulin and are removed). Run with: pymol -cq fig1A_pore.py
from pymol import cmd
import itertools, numpy as np
H = "${KCNQ2_ROOT}/10_figures"
G = "${KCNQ2_ROOT}/4_gate_analysis/pdbs"
cmd.reinitialize()
cmd.load(f"{G}/8j01.pdb", "op")
cmd.remove("not polymer")
cmd.remove("op and not chain A+B+D+G")           # remove calmodulin
cmd.remove("op and not resi 258-330")            # pore helix + filter + S6 (S5 removed)
pos = {c: np.array(cmd.get_atom_coords(f"op and chain {c} and resi 310 and name CA")) for c in "ABDG"}
pair = max(itertools.combinations("ABDG", 2), key=lambda p: np.linalg.norm(pos[p[0]] - pos[p[1]]))
cmd.create("two", f"op and chain {'+'.join(pair)}")
cmd.delete("op")
cmd.hide("everything"); cmd.show("cartoon", "two")
cmd.set("cartoon_transparency", 0.05)
cmd.color("0xBBBBBB", "two")
cmd.color("0x7FB3D5", "two and resi 276-281")                  # filter
cmd.show("sticks", "two and resi 310+314 and not hydro")
cmd.set("cartoon_side_chain_helper", 1)
cmd.color("0x1B7EBD", "two and resi 310"); cmd.color("0xD95F02", "two and resi 314")
cmd.show("spheres", "two and resi 310+314 and name CA")
cmd.set("sphere_scale", 0.42, "two and resi 310+314 and name CA")
cmd.bg_color("white"); cmd.set("ray_opaque_background", 0); cmd.set("antialias", 2)
cmd.set("depth_cue", 0); cmd.set("ray_shadows", 0); cmd.set("ray_trace_mode", 0)
# make the pore axis vertical: filter (276-281) at the top, gate (317-319) at the bottom
sf = np.array(cmd.centerofmass("two and resi 276-281 and name CA"))
gt = np.array(cmd.centerofmass("two and resi 317-319 and name CA"))
ax = gt - sf; ax /= np.linalg.norm(ax)                        # extracellular -> intracellular
ref = np.array([0.0, -1.0, 0.0])                              # downwards on the screen
v = np.cross(ax, ref); s = np.linalg.norm(v); c = float(np.dot(ax, ref))
K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
R = np.eye(3) + K + K @ K * ((1 - c) / (s ** 2 + 1e-12))
cmd.orient("two")
view = list(cmd.get_view())
view[:9] = list(R.T.flatten())
cmd.set_view(view)
cmd.zoom("two and resi 265-322", 2.2)
cmd.png(f"{H}/panels/fig1A_pore.png", width=760, height=1450, dpi=600, ray=1)
print("PORE PNG ok, subunits", pair)
