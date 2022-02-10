#!/bin/bash

module purge
module load state

ln -nsf $STATE_POT/Al_pbe1/#vnew.data pot.Al_pbe1
ln -nsf $STATE_POT/Cl_pbe1/#vnew.data pot.Cl_pbe1
mpirun -np 1 $STATE_COMMAND < nfinp_gdiis_pbc > nfout_gdiis_pbc


