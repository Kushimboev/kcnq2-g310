#!/bin/bash
# Conduction analysis: sf_books.py (primary estimator: net charge through the selectivity filter) + indep_count.py,
# (a) over the whole trajectory and (b) with the first 20 ns excluded (the reduced xtc is cut with gmx trjconv -b t0+20 ns; the scripts are unchanged).
# usage: an23.sh NAME [NAME...]   NAME: ANE_r1 ANE_r2 ANE_r3 ANE_r4 ANE_r5 ANC_r1 GS_r1 GS_r2 GS_r3 CR3_r1 CR3_r2 V280_r1 V280_r2
export LC_ALL=C GMX_MAXBACKUP=-1
P=${KCNQ2_ROOT}; C=$P/13_conduction; A=$C/scripts/analysis_2026-09-22_audit; O=$C/min_trajectories_2026-09-23
S=$C/min_trajectories_scratch_2026-09-22/s645628ac; OUT=$C/scripts/analysis_2026-09-23/out; mkdir -p $OUT
G=/usr/local/gromacs-2025.2/bin/gmx; GA=$C/AW/sub_analysis_min_A.gro; KA="resname POT and name POT"; KK="resname K and name K"
log(){ echo "$(date '+%F %T') $*" >> $OUT/an23.log; }
one(){ NM=$1 GRO=$2 XTC=$3 KS=$4 V=$5
  [ -s "$XTC" ] || { log "$NM: no xtc ($XTC)"; return; }
  python3 $A/sf_books.py $GRO $XTC 276:281 310 "$KS" 3 $NM > $OUT/$NM.books.txt 2>&1
  python3 $A/indep_count.py $GRO $XTC 276:281 310 "$KS" 3 $V $NM > $OUT/$NM.indep.txt 2>&1
  t0=$(python3 -c "import MDAnalysis as m,warnings; warnings.filterwarnings('ignore'); u=m.Universe('$GRO','$XTC'); print(u.trajectory[0].time)" 2>/dev/null)
  B=$OUT/${NM}_min_b20.xtc
  echo 0 | $G trjconv -s $GRO -f $XTC -b $(python3 -c "print($t0+20000)") -o $B > $OUT/${NM}_b20.trjconv.log 2>&1
  python3 $A/sf_books.py $GRO $B 276:281 310 "$KS" 3 ${NM}_b20 > $OUT/${NM}_b20.books.txt 2>&1
  python3 $A/indep_count.py $GRO $B 276:281 310 "$KS" 3 $V ${NM}_b20 > $OUT/${NM}_b20.indep.txt 2>&1
  rm -f $B $OUT/.${NM}_min_b20.xtc_offsets.*
  log "$NM (t0 $t0 ps): $(grep -m1 NET $OUT/$NM.books.txt | sed 's/ *| .*//; s/.*= //') full / $(grep -m1 NET $OUT/${NM}_b20.books.txt | sed 's/ *| .*//; s/.*= //') b20"; }
for NM in "$@"; do
  case $NM in
    ANE_r1|ANE_r2|ANE_r3) one $NM $GA $S/${NM}_min.xtc "$KA" 0.565 ;;
    ANE_r4|ANE_r5|ANC_r1) one $NM $GA $O/${NM}_min.xtc "$KA" 0.565 ;;
    V280_r*)              one $NM $GA $O/${NM}_min.xtc "$KA" 0.282 ;;
    GS_r*|CR3_r*)         one $NM $O/ref_min_$NM.gro $O/${NM}_min.xtc "$KK" 0.565 ;;
  esac
done
log "ALLDONE $*"
