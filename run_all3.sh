#!/bin/bash
# Third series.  b(H_T^d) tied to b(H_T^u).  Base excludes the preliminary CLAS12
# cross sections, COMPASS and HallA y21.  Then leave one out.
set -u
PY=~/.venv/bin/python3
cd ~/rc_iter1
export TIE_BD=1
BASE="clas6_pi0,clas6_eta,halla_n,halla_y16,halla_y11,bsa_demasi,bsa_zhao,bsa_clas12,eg1"
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
run 50_base "$BASE"
drop clas6_pi0   51_drop_clas6_pi0
drop clas6_eta   52_drop_clas6_eta
drop halla_n     53_drop_halla_n
drop halla_y16   54_drop_halla_y16
drop halla_y11   55_drop_halla_y11
drop bsa_demasi  56_drop_bsa_demasi
drop bsa_zhao    57_drop_bsa_zhao
drop bsa_clas12  58_drop_bsa_clas12
drop eg1         59_drop_eg1
$PY compare_runs.py > runs/COMPARISON.txt 2>&1
echo "### $(date +%H:%M:%S)  SERIES 3 DONE"
