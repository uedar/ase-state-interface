#!/bin/bash

module purge
module load state

ln -nsf $STATE_POT/Al_pbe1/#vnew.data pot.Al_pbe1
mpirun -np 1 $STATE_COMMAND < nfinp_scf > nfout_scf


