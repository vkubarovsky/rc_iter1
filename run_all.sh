#!/bin/bash
# Overnight series: fit, then plot, one directory per run.
set -u
PY=~/.venv/bin/python3
cd ~/rc_iter1
ALL="clas6_pi0,clas6_eta,halla_n,halla_y16,halla_y11,halla_y21,clas12_xs,compass,bsa_demasi,bsa_zhao,bsa_clas12,eg1"
NOC12="clas6_pi0,clas6_eta,halla_n,halla_y16,halla_y11,halla_y21,compass,bsa_demasi,bsa_zhao,bsa_clas12,eg1"

run () {   # run <tag> <comma list>
  echo "### $(date +%H:%M:%S)  $1"
  $PY fitrun.py "$1" "$2" > "runs_$1.log" 2>&1 || { echo "   FIT FAILED"; return; }
  $PY plots.py  "$1"      >> "runs_$1.log" 2>&1 || echo "   PLOTS FAILED"
  mv "runs_$1.log" "runs/$1/run.log"
  $PY -c "
import json;r=json.load(open('runs/$1/summary.json'))
print('   chi2/ndf = %.4f  (%.1f/%d)'%(r['chi2_ndf'],r['chi2_fitted'],r['ndf']))"
}

run 01_all              "$ALL"
run 02_no_clas12xs      "$NOC12"

# leave one out, starting from 02
drop () {
  local d=$1
  local keys=$(echo "$NOC12" | tr ',' '\n' | grep -v "^${d}$" | paste -sd, -)
  run "$2" "$keys"
}
drop clas6_pi0   03_drop_clas6_pi0
drop clas6_eta   04_drop_clas6_eta
drop halla_n     05_drop_halla_n
drop halla_y16   06_drop_halla_y16
drop halla_y11   07_drop_halla_y11
drop halla_y21   08_drop_halla_y21
drop compass     09_drop_compass
drop bsa_demasi  10_drop_bsa_demasi
drop bsa_zhao    11_drop_bsa_zhao
drop bsa_clas12  12_drop_bsa_clas12
drop eg1         13_drop_eg1

$PY compare_runs.py > runs/COMPARISON.txt 2>&1
echo "### $(date +%H:%M:%S)  ALL DONE"
