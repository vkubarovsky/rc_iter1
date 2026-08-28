#!/bin/bash
# One full RC iteration: given a parameter file, compute RC, recorrect the
# published cross sections, refit phi, refit amplitudes.
#   ./iterate.sh <iter_label> <seed_par.npy>
set -e
IT=$1; SEED=$2
PY=~/.venv/bin/python3
cd ~/rc_iter1
# install the seed parameters as the exclurad_py amp2021 parameters.
# exclurad_py/models/_amplitude_fit.py still uses the OLD slope convention
# (b at xB = 0.15), so convert - a raw cp would silently change the model.
$PY -c "import numpy as np, sys; sys.path.insert(0,'$HOME/rc_iter1'); import reparam; \
np.save('$HOME/exclurad_py/exclurad_py/models/amp2021_par.npy', reparam.to_old(np.load('$SEED')))"
for ch in pi0 eta; do
  $PY compute_rc.py $ch 12 > rc_${ch}_${IT}.log 2>&1
  mv rc_${ch}.txt rc_${ch}_${IT}.txt
  $PY refit_phi.py $ch rc_${ch}_${IT}.txt sf_${ch}_${IT}.txt
  $PY to_strfun.py $ch sf_${ch}_${IT}.txt data/strfun_${ch}.data
done
SEED=$SEED OUTP=fitpar_${IT}.npy $PY fit_clas12.py > fit_${IT}.log 2>&1
echo "ITER $IT DONE"
