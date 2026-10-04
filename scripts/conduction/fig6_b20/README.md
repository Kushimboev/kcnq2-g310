# Charge-scaling (ECC) comparison and KCNQ1 benchmark (Fig. 6)

All counts use the same analysis window as the rest of the study: the first 20 ns of every run is excluded. The protocol is that
of the main analysis (`../an23.sh`): `gmx trjconv -b t0+20 ns`, then `sf_books.py` (primary estimator: net charge through the
selectivity filter) and `indep_count.py` (independent estimator).

| file | content |
|---|---|
| `b20_ecc.sh` | per-run event books for the seven runs of the comparison (AN r1/r2, A9 r1, KN r1, KNE r1/r2, K10 r1); the results are in `processed/conduction/{RUN}_b20.books.txt` and `{RUN}_b20.indep.txt` |
| `d1_b20.py` → `d1_b20.json`, `d1_b20_part1.txt` | x axis of Fig. 6C: the narrow G310 (KCNQ1: G345) Cα–Cα diagonal, taken per frame and averaged over the window |
| `ecc_b20_summary.py` → `ecc_b20_summary.txt` | all numbers: matched pairs, all scaled-charge runs, sensitivity to excluding ANE r3, KCNQ2 and KCNQ1 separately, KCNQ1 benchmark; the primary and both independent estimators |
| `../../figures/fig6_ecc.py` | the figure |

Run names: AN = KCNQ2 with scaled charges; ANE = KCNQ2 with full charges; KN / KNE = KCNQ1 with scaled / full charges;
A9 and K10 = scaled-charge runs with the G310 (G345) diagonals restrained at 0.90 and 1.00 nm.

## Result

| comparison | first 20 ns excluded |
|---|---|
| matched pairs (AN r1/r2 + KN r1 vs ANE r1/r2/r3 + KNE r1/r2) | 2 / 840 ns vs 23 / 1217.8 ns → ratio 0.13 (0.01–0.51), p = 3.2×10⁻⁴ |
| all scaled-charge runs (K10 r1 added) | 4 / 1120 ns vs 23 / 1217.8 ns → 0.19 (0.05–0.55), p = 3.4×10⁻⁴ |
| KCNQ1 benchmark (KNE r1 + r2) | 11 / 553.3 ns = 5.64 pS (2.81–10.09) |
| filter occupancy | scaled charges: three ions in 97 % of frames; full charges: two ions in 88 % |
| Fig. 6C | A9 7 / 280 ns vs ANE 12 / 664.5 ns (+0.55 Å); K10 2 / 280 ns vs KNE 11 / 553.3 ns (+1.5 Å) |

**Why the first 20 ns matter here.** The full-charge runs show six net crossings within their first 20 ns, all of them exits of
ions placed in the filter at build time (a relaxation transient); the scaled-charge runs show none. Excluding that window lowers
the full-charge rate from 22 to 19 per µs.

**Independent estimators.** Crossings of the lower filter boundary: 2 / 840 vs 24 / 1218 → 0.12 (0.01–0.49); crossings of the
gate plane: 4 / 840 vs 27 / 1218 → 0.22 (0.06–0.62). Without ANE r3, which shares its starting frame with r1: 0.13 (0.01–0.52).
Pairs taken separately: KCNQ2 1 / 560 vs 12 / 664.5 → 0.10 (0.002–0.67), p = 0.004; KCNQ1 1 / 280 vs 11 / 553.3 → 0.18
(0.004–1.24), p = 0.052.

**To recompute:** `bash b20_ecc.sh`, `python3 d1_b20.py RUN…`, `python3 ecc_b20_summary.py`, `python3 ../../figures/fig6_ecc.py`.
Paths inside the scripts refer to the authors' working tree through `${KCNQ2_ROOT}`.
