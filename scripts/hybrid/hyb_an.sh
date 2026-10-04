#!/bin/bash
# Analysis of the mixed tetramers (H1/H2). The estimators are exactly those used for G310S and the wild type (the read-out rule was fixed before the runs finished):
#   (1) conduction: sf_books.py (net filter charge) + indep_count.py, on the reduced (Analysis_min) trajectory
#   (2) mechanism (primary measure, gate water): mech_gs.py; needs the FULL trajectory and the system .gro (the reduced trajectory has no water),
#       so it runs only where the full trajectory is available locally.
# usage: hyb_an.sh [H1_r1 H1_r2 H2_r1 H2_r2 ...]     flags: NOCOND=1 (mechanism only), NOMECH=1 (conduction only)
set -u; export LC_ALL=C GMX_MAXBACKUP=-1
P=${KCNQ2_ROOT}; C=$P/13_conduction; A=$C/scripts/analysis_2026-09-22_audit
OUT=$C/scripts/analysis_2026-09-29_hybrid; W=$OUT/work; mkdir -p $W $OUT/mech
G=/usr/local/gromacs-2025.2/bin/gmx
SSHC="ssh -i $HOME/.ssh/ssh_key -o BatchMode=yes -o ConnectTimeout=10"
NOCOND=${NOCOND:-0}; NOMECH=${NOMECH:-0}
log(){ echo "$(date '+%F %T') $*" >> $OUT/hyb.log; }
box(){ case $1 in H1_r1) echo local;; H1_r2) echo user@host2;; H2_r1) echo user@host3;;
       H2_r2) sed 's/^L$/local/; s/^T1$/user@host2/; s/^T2$/user@host3/' $C/queue_2026-09-23/q28_H2_r2.box 2>/dev/null;; esac; }
bdir(){ case $1 in H1_r2) echo 4070Ti_1;; H2_r1) echo 4070Ti_2;; H2_r2) sed 's/^T1$/4070Ti_1/; s/^T2$/4070Ti_2/; s/^L$/local/' $C/queue_2026-09-23/q28_H2_r2.box 2>/dev/null;; esac; }

fullxtc(){ # NM -> path of the local full trajectory (empty if absent)
  local NM=$1 N=${1%_*} R=${1#*_} b
  [ -s $C/$N/run_$R/prod_$NM.xtc ] && { echo $C/$N/run_$R/prod_$NM.xtc; return; }
  b=$(bdir $NM); [ -n "${b:-}" ] && [ -s $P/backup_2026-09-23/fleet/$b/8J01_$NM/prod_$NM.xtc ] && echo $P/backup_2026-09-23/fleet/$b/8J01_$NM/prod_$NM.xtc; }

prep(){ # NM -> $W/${NM}_min.xtc + $W/ref_min_$NM.gro
  local NM=$1 N=${1%_*} R=${1#*_} H S; H=$(box $NM); S=$C/$N/run_$R/sys
  [ -n "${H:-}" ] || { log "$NM: host not found"; return 1; }
  [ -s $W/ref_min_$NM.gro ] || echo Analysis_min | $G trjconv -s $S/${N}_system.gro -f $S/${N}_system.gro -n $S/index.ndx -o $W/ref_min_$NM.gro > $W/ref_$NM.log 2>&1
  local F; F=$(fullxtc $NM)
  # The .tpr is taken from next to the trajectory and a stale reduced trajectory is deleted first: otherwise a failed trjconv
  # would leave an older, shorter reduced trajectory in place and the analysis would silently use it.
  rm -f $W/${NM}_min.xtc $W/.${NM}_min.xtc_offsets.*
  if [ -n "${F:-}" ]; then
    local T=$C/$N/run_$R/prod_$NM.tpr; [ -s $T ] || T=$(dirname $F)/prod_$NM.tpr
    echo Analysis_min | $G trjconv -f $F -s $T -n $S/index.ndx -o $W/${NM}_min.xtc > $W/trjconv_$NM.log 2>&1
  elif [ "$H" != local ]; then
    timeout 2400 $SSHC -n $H "cd ~/8J01_$NM && export GMX_MAXBACKUP=-1 LC_ALL=C && echo Analysis_min | $G trjconv -f prod_${NM}.xtc -s prod_${NM}.tpr -n sys/index.ndx -o ${NM}_min.xtc > trjconv_min.log 2>&1; ls -l ${NM}_min.xtc" >> $W/trjconv_$NM.log 2>&1
    rsync -a -e "$SSHC" $H:8J01_$NM/${NM}_min.xtc $W/ >> $W/trjconv_$NM.log 2>&1
  fi
  [ -s $W/${NM}_min.xtc ] || { log "$NM: min.xtc was not created"; return 1; }
  log "$NM min.xtc last frame: $(grep -a 'Last written' $W/trjconv_$NM.log | tail -1)"
  log "$NM min.xtc $(stat -c %s $W/${NM}_min.xtc) bytes"; }

cond(){ local NM=$1 GRO=$W/ref_min_$1.gro XT=$W/$1_min.xtc t0 B
  python3 $A/sf_books.py   $GRO $XT 276:281 310 "resname K and name K" 3 $NM        > $OUT/$NM.books.txt 2>&1
  python3 $A/indep_count.py $GRO $XT 276:281 310 "resname K and name K" 3 0.565 $NM > $OUT/$NM.indep.txt 2>&1
  t0=$(python3 -c "import MDAnalysis as m,warnings; warnings.filterwarnings('ignore'); u=m.Universe('$GRO','$XT'); print(u.trajectory[0].time)" 2>/dev/null)
  B=$W/${NM}_min_b20.xtc
  echo 0 | $G trjconv -s $GRO -f $XT -b $(python3 -c "print($t0+20000)") -o $B > $W/${NM}_b20.trjconv.log 2>&1
  python3 $A/sf_books.py   $GRO $B 276:281 310 "resname K and name K" 3 ${NM}_b20        > $OUT/${NM}_b20.books.txt 2>&1
  python3 $A/indep_count.py $GRO $B 276:281 310 "resname K and name K" 3 0.565 ${NM}_b20 > $OUT/${NM}_b20.indep.txt 2>&1
  rm -f $B $W/.${NM}_min_b20.xtc_offsets.*
  log "$NM COND: NET $(grep -m1 NET $OUT/$NM.books.txt | sed 's/ *| .*//; s/.*= //') full / $(grep -m1 NET $OUT/${NM}_b20.books.txt | sed 's/ *| .*//; s/.*= //') b20"; }

mech(){ local NM=$1 N=${1%_*} R=${1#*_} F; F=$(fullxtc $NM)
  [ -n "${F:-}" ] || { log "$NM MECH: full xtc not available locally yet — skipped"; return 0; }
  python3 $C/scripts/analysis_2026-09-23/mech_gs.py $NM $C/$N/run_$R/sys/${N}_system.gro $F $OUT/mech > $OUT/mech/$NM.log 2>&1
  log "$NM MECH: $(tr '\t' ' ' < $OUT/mech/$NM.txt 2>/dev/null | grep -oE 't=[^ ]+|d1=[^ ]+|OG_R=[^ ]+|gate_water=[^ ]+' | tr '\n' ' ')"; }

SEL="${*:-H1_r1 H1_r2 H2_r1}"
log "=== start: $SEL (NOCOND=$NOCOND NOMECH=$NOMECH)"
for NM in $SEL; do
  [ $NOCOND = 1 ] || { prep $NM && cond $NM; }
  [ $NOMECH = 1 ] || mech $NM
done
log "ALLDONE $SEL"
