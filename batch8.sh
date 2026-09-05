#!/bin/bash
# Series 8: leave-one-out on the corrected data (De Masi and Hall-A systematics
# from the papers, normalisation nuisances, phase bounds repaired).
set -e
cd /Users/vpk/rc_iter1
ALL="bsa_clas12,bsa_demasi_phi,bsa_zhao,clas6_eta,clas6_pi0,eg1,halla_n,halla_y11,halla_y16"
SEED=runs/83_systtrue/fitpar.npy
PY=~/.venv/bin/python3
run () { echo "--- $1"; $PY run_phase_fix.py "$1" "$2" $SEED 2>&1 | grep -E '^===|^    (phases|norm)'; }
run 90_base "$ALL"
i=91
for drop in clas6_pi0 clas6_eta halla_n halla_y16 halla_y11 bsa_demasi_phi bsa_zhao bsa_clas12 eg1; do
  KEEP=$($PY -c "print(','.join(k for k in '$ALL'.split(',') if k!='$drop'))")
  run "${i}_drop_${drop}" "$KEEP"
  i=$((i+1))
done
echo "СЕРИЯ 8 ГОТОВА"
