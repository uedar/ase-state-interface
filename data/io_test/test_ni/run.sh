#!/bin/bash

module purge
module load state

ln -nsf $STATE_POT/Ni_pbe4/#vnew.data pot.Ni_pbe4
mpirun -np 1 $STATE_COMMAND < nfinp_scf > nfout_scf


