# A permeable G310 constriction in the activated KCNQ2 pore, and its disruption by G310S

Input files, analysis code and processed data for the molecular dynamics study of K<sup>+</sup> permeation
through the activated KCNQ2 (Kv7.2) pore domain and the developmental and epileptic encephalopathy
variant G310S.

**Manuscript:** bioRxiv preprint, [doi:10.64898/2026.10.05.756556](https://doi.org/10.64898/2026.10.05.756556) · **Archived release:** *(Zenodo DOI to be added)*

Everything here is derived from the simulations reported in the manuscript. Full production
trajectories (about 9.3 GB per 300 ns run, roughly 190 GB in total) are not in this repository; they are
available from the corresponding author on request. The reduced trajectories used for the conduction
analysis are deposited with the archived release.

## What the data show

| quantity | result |
|---|---|
| Wild-type pore (PDB 8J01), 4 x 300 ns at ~565 mV | 17 net filter crossings in 1.12 us → **4.30 pS** (95% CI 2.51–6.89) |
| G310S, 3 x 300 ns | 2 in 840 ns → **G310S/WT = 0.16** (95% CI 0.02–0.66), p = 0.003 |
| Closed-state pore (PDB 7CR3), 2 x 300 ns | 0 in 560 ns (p = 0.001) |
| G310 restrained to closed-state width | 0 gate crossings in 280 ns |
| 282 mV, 3 x 300 ns | 3 in 840 ns; G(282)/G(565) = 0.47 (0.09–1.63) |
| KCNQ1 benchmark (PDB 6V01), 2 runs | 11 in 553 ns → 5.6 pS |
| Gate-plane water, WT → G310S | 8.73 ± 0.24 → 4.73 ± 0.12 (Amber14sb); 9.91 → 4.34 (CHARMM36m) |
| First coordination shell at the constriction | 6.30 ± 0.05 (WT) vs 6.52 ± 0.07 (G310S) oxygens — conserved |
| Excursions of ions towards the gate (Table S6) | reached the gate plane: 88 in 1.12 us (WT) vs 38 in 0.84 us (G310S), rate ratio 0.58; of those, crossed: 20 (23%) vs 2 (5%); time in the plane 141 vs 117 ps; K–serine Oγ contacts median 30 ps — the ion is turned back at the serine ring, not retained |
| Charge scaling (ECC, q = 0.78), structure- and width-matched | 2 in 840 ns vs 23 in 1.22 us → 0.13 (0.01–0.51), p = 3.2e-4 |
| **Mixed tetramers**: Ser310 in 1 and in 2 adjacent subunits, 2 x 300 ns each | gate water 7.83 and 6.43 → **f = 0.22** (0.13–0.32) and **0.58** (0.37–0.78) of the G310S effect, against the additive expectation n/4 → **additive**; conduction 7 net crossings in 560 ns each (not resolved) |

All counts use the same analysis window: the first 20 ns of every run is discarded.

## Layout

```
inputs/          one directory per simulated system (Table S1): coordinates (.gro.gz), topologies
                 (.top/.itp), index files (index.ndx.gz); inputs/mdp/ holds the GROMACS parameter files
forcefield/      Amber14sb (protein) + Lipid21 + TIP3P + Joung-Cheatham ions, and the CHARMM36m set
scripts/
  conduction/    sf_books.py (primary estimator: net charge through the selectivity filter),
                 indep_count.py, disp_flux.py (independent displacement estimator), drivers;
                 fig6_b20/ holds the charge-scaling comparison and the KCNQ1 benchmark
  mechanism/     gate water and gate geometry (mech_gs.py), ion coordination (gate_hydration.py,
                 hyd_stats.py), ion traces for Fig. 3 (ion_traces.py), excursions of ions towards the
                 gate (gate_tracks.py, gate_commit_verify.py) and an independent visit-based count
                 (gate_commit.py), CHARMM36m summary
  hybrid/        mixed tetramers (one and two S310 subunits): hyb_an.sh, hyb_readout.py (applies the
                 decision rule fixed before those runs finished), indep_gatewater.py
  figures/       one script per manuscript figure; all numbers come through paperdata.py
  style/         colour dictionary and matplotlib style used by the figure scripts
processed/
  conduction/    per-run event books (*_b20.books.txt), independent estimators, summaries
  mechanism/     gate water and geometry, coordination (per-ion records and summary), excursions
  figures/       the manuscript figures (PDF) and the numbers behind Fig. 6
MANIFEST.txt     every file with its size
```

Run names: ANE = wild type (full charges), GS = G310S, CR3 = closed-state pore (7CR3), ANC = narrowed gate,
V280 = wild type at 282 mV, H1 / H2 = mixed tetramers, KNE = KCNQ1 (full charges), AN / KN / A9 / K10 =
scaled-charge (ECC) runs, C36W / C36G = CHARMM36m wild type / G310S.

**Paths.** The scripts are the ones that were run for the study. They locate their input through the environment
variable `KCNQ2_ROOT` and the directory layout of the authors' working tree; the corresponding files are under
`processed/` and `inputs/` here (for example `13_conduction/scripts/analysis_2026-09-23/out/` corresponds to
`processed/conduction/`).

**Rule labels in the summaries.** The text that follows "-> rule:" in `processed/conduction/FINAL_SUMMARY.txt`, the "DECISION" lines of
`HYBRID_READOUT.txt` apply reading rules that were written down before the
corresponding runs finished. The 282 mV comparison is reported
as the ratio of conductances with its 95 % interval.

## Reproducing the numbers

1. Conduction counts: `scripts/conduction/sf_books.py <system.gro> <trajectory.xtc> 276:281 310 "resname K and name K" 3 <label>`.
   The last argument before the label is the number of ions placed in the filter at build time; these are
   tracked separately from free ions, and counting them as permeation events is a known source of error.
2. Conductance and confidence intervals: `scripts/conduction/final23_summary.py` (Poisson intervals for
   pooled counts, conditional binomial tests for ratios between systems).
3. Mechanism: `scripts/mechanism/mech_gs.py <label> <system.gro> <trajectory.xtc> <outdir>`. Use the
   system `.gro` as topology, not the `.tpr`, which renumbers residues.
4. Figures: `scripts/figures/fig*.py` (they import `scripts/style/palette.py`); the numbers are read through `paperdata.py`.
5. Excursions towards the gate (Table S6): `scripts/mechanism/gate_tracks.py <label> <system.gro> <trajectory.xtc> <outdir>` writes the
   positions of every K+ within 1.2 nm of the pore axis for each frame of the analysis window (it prints the number of frames
   read; 28,001 are expected for a 300 ns run). `scripts/mechanism/gate_commit_verify.py processed/mechanism/gate_excursions`
   then reproduces `GATE_COMMIT_VERIFY.txt` and `gate_commit.json` from the deposited `*.tracks.npz`. An excursion starts when an
   ion rises above 0.6 nm below the G310 plane and ends with its return or its arrival in the cavity, 0.5 nm above the plane.
   `scripts/mechanism/gate_commit.py processed/mechanism/coordination` is the independent visit-based count.
6. Mixed tetramers: `scripts/hybrid/hyb_readout.py` reproduces `processed/conduction/HYBRID_READOUT.txt`, which applies
   the rule fixed before those runs finished (additive expectation f = n/4 for n variant subunits; gate water is the
   primary measure, conduction the secondary one). `scripts/hybrid/indep_gatewater.py` recomputes the gate water with a
   different pore-axis definition and a different set of frames; its result is `processed/conduction/VERIFY_SUMMARY.txt`.

Analyses were run with GROMACS 2025.2, MDAnalysis 2.9.0, NumPy, SciPy and Matplotlib.

## Citation

Qoshimboyev S, Marimuthu P, Makhkamov M, Razzokov J. A permeable G310 constriction limits K<sup>+</sup> flux in the activated
KCNQ2 pore and is disrupted by the encephalopathy variant G310S. *bioRxiv* (2026). [https://doi.org/10.64898/2026.10.05.756556](https://doi.org/10.64898/2026.10.05.756556)

## Licence

Code: MIT License (see `LICENSE`) · Data and documentation: Creative Commons Attribution 4.0 International (CC BY 4.0).
