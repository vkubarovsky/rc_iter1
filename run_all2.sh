#!/bin/bash
# Second series: leave-one-out from a base that excludes the preliminary CLAS12
# cross sections AND COMPASS.  COMPASS sits at xB ~ 0.1 and W ~ 4.4 GeV, far
# outside anything the model was built on, and in the first series it dragged
# every fit to the same distorted place.
set -u
PY=~/.venv/bin/python3
cd ~/rc_iter1
BASE="clas6_pi0,clas6_eta,halla_n,halla_y16,halla_y11,halla_y21,bsa_demasi,bsa_zhao,bsa_clas12,eg1"
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
run 20_base "$BASE"
drop clas6_pi0   21_drop_clas6_pi0
drop clas6_eta   22_drop_clas6_eta
drop halla_n     23_drop_halla_n
drop halla_y16   24_drop_halla_y16
drop halla_y11   25_drop_halla_y11
drop halla_y21   26_drop_halla_y21
drop bsa_demasi  27_drop_bsa_demasi
drop bsa_zhao    28_drop_bsa_zhao
drop bsa_clas12  29_drop_bsa_clas12
drop eg1         30_drop_eg1
$PY compare_runs.py > runs/COMPARISON.txt 2>&1
echo "### $(date +%H:%M:%S)  SERIES 2 DONE"
