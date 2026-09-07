#!/bin/bash
# Refit everything that had HallA_y21 in the fit, after the xB correction.
# C0 first: it is the production vector and the generator is waiting on it.
set -e
cd /Users/vpk/rc_iter1
PY=~/.venv/bin/python3
ALL11="bsa_clas12,bsa_demasi_phi,bsa_zhao,clas6_eta,clas6_pi0,eg1,halla_n,halla_y11,halla_y16,halla_y21,compass"
ALL10="bsa_clas12,bsa_demasi_phi,bsa_zhao,clas6_eta,clas6_pi0,eg1,halla_n,halla_y11,halla_y16,halla_y21"
SEED=runs/C0_with_compass/fitpar.npy
run () { echo "--- $1 $(date +%H:%M)"; $PY run_phase_fix.py "$1" "$2" $SEED 2>&1 | grep -E '^===|^    norm'; }

run C1_with_compass "$ALL11"          # the production fit
SEED=runs/C1_with_compass/fitpar.npy
run D00_base "$ALL10"                 # series D: leave-one-out, y21 in, corrected xB
i=1
for drop in clas6_pi0 clas6_eta halla_n halla_y16 halla_y11 halla_y21 bsa_demasi_phi bsa_zhao bsa_clas12 eg1; do
  KEEP=$($PY -c "print(','.join(k for k in '$ALL10'.split(',') if k!='$drop'))")
  run "$(printf 'D%02d_drop_%s' $i $drop)" "$KEEP"
  i=$((i+1))
done
TIE_BET=1 $PY run_phase_fix.py T2_tie_bET "$ALL11" runs/C1_with_compass/fitpar.npy 2>&1 | grep -E '^==='
echo "ПЕРЕФИТ ЗАКОНЧЕН $(date +%H:%M)"
