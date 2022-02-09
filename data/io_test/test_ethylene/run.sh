#!/bin/bash

module purge
module load state

ln -nsf $STATE_POT/C_pbe3/#vnew.data pot.C_pbe3
ln -nsf $STATE_POT/H_lda3/#vnew.data pot.H_lda3
mpirun -np 1 $STATE_COMMAND < nfinp_gdiis > nfout_gdiis


