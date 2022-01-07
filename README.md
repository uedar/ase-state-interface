# ase-state-interface

### Install
```
pip --no-cache-dir install --upgrade --index-url https://test.pypi.org/simple/ state-interface
```

### Examples
```
from state_interface.state import STATE

import os
from ase.build import bulk

super_cell = bulk('Si')
# ---------------- set up executable ----------------
label = 'Si'
input_file = label+'.in'
output_file = label+'.out'
no_cpus = 32
pw_loc = './STATE'

## parallel qe using srun (for slurm system)
os.environ['ASE_STATE_COMMAND'] = f'mpirun -np {no_cpus} {pw_loc} < {input_file} > {output_file}'
# -------------- set up input parameters --------------
input_data = {'WF_OPT'       :    'DAV' ,
              'TYPE'         :       0  ,
              'NSPG'         :       1  ,
              'GMAX'         :     4.00 ,
              'GMAXP'        :     8.00 ,
              'KPOINT_MESH'  : [8, 8, 8],
              'KPOINT_SHIFT' : ["OFF", "OFF", "OFF"],
              'WIDTH'        : 0.0002,
              'EDELTA'       : 0.5e-9,
              'NEG'          : 8,
              'PSEUDOS'      : [['Si', 28.00, 'pot.Si_pbe1']]
              }
# ----------------  pseudo-potentials -----------------
ion_pseudo = {'Si': 'pot.Si_pbe1'}
# -------------- create ASE calculator ----------------
dft_calc = STATE(label=label, input_data=input_data)
super_cell.calc = dft_calc
super_cell.get_forces()
super_cell.get_potential_energy()
```
