#!/bin/bash
# Series B: leave-one-out with Hall-A 2021 inside the fit.  y21 is the only set
# above Q2 = 4.3, so this is the first series in which the Q2 lever arm is not
# an extrapolation.
set -e
cd /Users/vpk/rc_iter1
ALL="bsa_clas12,bsa_demasi_phi,bsa_zhao,clas6_eta,clas6_pi0,eg1,halla_n,halla_y11,halla_y16,halla_y21"
SEED=runs/A0_with_y21/fitpar.npy
PY=~/.venv/bin/python3
run () { echo "--- $1"; $PY run_phase_fix.py "$1" "$2" $SEED 2>&1 | grep -E '^===|^    (phases|norm)'; }
run B00_base "$ALL"
i=1
for drop in clas6_pi0 clas6_eta halla_n halla_y16 halla_y11 halla_y21 bsa_demasi_phi bsa_zhao bsa_clas12 eg1; do
  KEEP=$($PY -c "print(','.join(k for k in '$ALL'.split(',') if k!='$drop'))")
  run "$(printf 'B%02d_drop_%s' $i $drop)" "$KEEP"
  i=$((i+1))
done
echo "СЕРИЯ B ГОТОВА"
