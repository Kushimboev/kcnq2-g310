#!/bin/bash
# Charge-scaling (ECC) comparison (Fig. 6) and KCNQ1 benchmark in the window with the first 20 ns excluded,
# with exactly the protocol of the main analysis (an23.sh): gmx trjconv -b t0+20 ns -> sf_books.py (primary) + indep_count.py (independent).
# The scripts sf_books.py and indep_count.py are unchanged. Output: {RUN}_b20.books.txt and {RUN}_b20.indep.txt
# (read by paperdata.py). The reduced trajectory of KN r1 is regenerated here from the raw trajectory,
# and its full-window book is compared with the stored one as a consistency check.
export LC_ALL=C GMX_MAXBACKUP=-1
P=${KCNQ2_ROOT}; C=$P/13_conduction; A=$C/scripts/analysis_2026-09-22_audit
S=$C/min_trajectories_scratch_2026-09-22/s645628ac; OUT=$C/scripts/analysis_2026-09-23/out; W=$C/scripts/analysis_2026-09-29_fig6_b20
K=$P/14_kcnq1_control; G=/usr/local/gromacs-2025.2/bin/gmx
GA=$C/AW/sub_analysis_min_A.gro; GK=$S/KNE_min_ref.gro; KA="resname POT and name POT"; KK="resname K and name K"
log(){ echo "$(date '+%F %T') $*" >> $W/b20_ecc.log; }
mkdir -p $W/work

# --- 0) regenerate KN_r1_min.xtc from the raw trajectory (as for KNE r1: group Analysis_min, 193 atoms)
if [ ! -s $S/KN_r1_min.xtc ]; then
  echo Analysis_min | $G trjconv -f $K/prod_KN_r1.xtc -s $K/prod_KN_r1.tpr -n $K/sys/index.ndx -o $S/KN_r1_min.xtc > $W/work/KN_r1_trjconv.log 2>&1
  log "KN_r1_min.xtc regenerated: $(ls -la $S/KN_r1_min.xtc | awk '{print $5}') bytes"
fi
python3 $A/sf_books.py $GK $S/KN_r1_min.xtc 311:316 345 "$KK" 3 KN_r1 > $W/work/KN_r1.full.books.txt 2>&1
a=$(grep -m1 -E 'NET|####' $W/work/KN_r1.full.books.txt | head -n1); 
sn1=$(grep -m1 'NET outward' $W/work/KN_r1.full.books.txt | sed 's/.*= *//;s/ .*//'); sn0=$(grep -m1 'NET outward' $A/KN_r1.books.txt | sed 's/.*= *//;s/ .*//')
sv1=$(grep -m1 '####' $W/work/KN_r1.full.books.txt | sed 's/.*visits //'); sv0=$(grep -m1 '####' $A/KN_r1.books.txt | sed 's/.*visits //')
if diff <(grep -v '^####' $W/work/KN_r1.full.books.txt | sed 's/KN_r1//') <(grep -v '^####' $A/KN_r1.books.txt | sed 's/KN_r1//') > /dev/null; then v=IDENTICAL; else v="DIFFERENT"; fi
log "KN_r1 check (regenerated vs stored book): NET $sn1 vs $sn0, SF visits $sv1 vs $sv0 -> $v"

# --- 1) event books with the first 20 ns excluded
one(){ NM=$1 GRO=$2 XTC=$3 SFR=$4 GR=$5 KS=$6 V=$7
  [ -s $OUT/${NM}_b20.books.txt ] && grep -q NET $OUT/${NM}_b20.books.txt && { log "$NM: book already present — skipped"; return; }
  t0=$(python3 -c "import MDAnalysis as m,warnings; warnings.filterwarnings('ignore'); u=m.Universe('$GRO','$XTC'); print(u.trajectory[0].time)" 2>/dev/null)
  B=$W/work/${NM}_min_b20.xtc
  echo 0 | $G trjconv -s $GRO -f $XTC -b $(python3 -c "print($t0+20000)") -o $B > $W/work/${NM}_b20.trjconv.log 2>&1
  python3 $A/sf_books.py $GRO $B $SFR $GR "$KS" 3 ${NM}_b20 > $OUT/${NM}_b20.books.txt 2>&1
  python3 $A/indep_count.py $GRO $B $SFR $GR "$KS" 3 $V ${NM}_b20 > $OUT/${NM}_b20.indep.txt 2>&1
  rm -f $B $W/work/.${NM}_min_b20.xtc_offsets.*
  log "$NM (t0 $t0 ps): $(grep -m1 '####' $OUT/${NM}_b20.books.txt | sed 's/.*_b20: //') | $(grep -m1 'NET outward' $OUT/${NM}_b20.books.txt | sed 's/.*= *//')"; }
one AN_r1  $GA $C/AN/interim/prod_AN_r1_min.xtc 276:281 310 "$KA" 0.565
one AN_r2  $GA $C/AN/interim/prod_AN_r2_min.xtc 276:281 310 "$KA" 0.565
one A9_r1  $GA $S/A9_r1_min.xtc                 276:281 310 "$KA" 0.565
one KN_r1  $GK $S/KN_r1_min.xtc                 311:316 345 "$KK" 0.565
one KNE_r1 $GK $S/KNE_r1_min.xtc                311:316 345 "$KK" 0.565
one KNE_r2 $GK $S/KNE_r2_min.xtc                311:316 345 "$KK" 0.565
one K10_r1 $GK $S/K10_r1_min.xtc                311:316 345 "$KK" 0.565
echo ALLDONE > $W/b20_ecc.done; log "ALLDONE"
