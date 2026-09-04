#!/bin/bash
# Fourth series.  De Masi enters at the phi level, 703 published points instead of
# 56 digitised amplitudes whose errors were 1.5 times too small.  b(H_T^d) stays
# tied to b(H_T^u).  Base excludes the preliminary CLAS12 cross sections, COMPASS
# and HallA y21.
set -u
PY=~/.venv/bin/python3
cd ~/rc_iter1
export TIE_BD=1
BASE="clas6_pi0,clas6_eta,halla_n,halla_y16,halla_y11,bsa_demasi_phi,bsa_zhao,bsa_clas12,eg1"
run () {
  echo "### $(date +%H:%M:%S)  $1"
  $PY fitrun.py "$1" "$2" > "runs_$1.log" 2>&1 || { echo "   FIT FAILED"; return; }
  $PY plots.py  "$1"      >> "runs_$1.log" 2>&1 || echo "   PLOTS FAILED"
  mv "runs_$1.log" "runs/$1/run.log"
  $PY -c "
import json;r=json.load(open('runs/$1/summary.json'))
print('   chi2/ndf = %.4f  (%.1f/%d)'%(r['chi2_ndf'],r['chi2_fitted'],r['ndf']))"
}
drop () {
  local keys=$(echo "$BASE" | tr ',' '\n' | grep -v "^${1}$" | paste -sd, -)
  run "$2" "$keys"
}
run 70_base "$BASE"
drop clas6_pi0       71_drop_clas6_pi0
drop clas6_eta       72_drop_clas6_eta
drop halla_n         73_drop_halla_n
drop halla_y16       74_drop_halla_y16
drop halla_y11       75_drop_halla_y11
drop bsa_demasi_phi  76_drop_bsa_demasi
drop bsa_zhao        77_drop_bsa_zhao
drop bsa_clas12      78_drop_bsa_clas12
drop eg1             79_drop_eg1
$PY compare_runs.py > runs/COMPARISON.txt 2>&1
echo "### $(date +%H:%M:%S)  SERIES 5 DONE"
