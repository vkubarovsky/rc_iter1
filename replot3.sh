#!/bin/bash
cd ~/rc_iter1
for t in 50_base 51_drop_clas6_pi0 52_drop_clas6_eta 53_drop_halla_n 54_drop_halla_y16 \
         55_drop_halla_y11 56_drop_bsa_demasi 57_drop_bsa_zhao 58_drop_bsa_clas12 59_drop_eg1; do
  ~/.venv/bin/python3 plots.py "$t" >/dev/null 2>&1 && echo "ok $t" || echo "FAIL $t"
done
~/.venv/bin/python3 compare_runs.py > runs/COMPARISON.txt 2>&1
echo REPLOT3 DONE
