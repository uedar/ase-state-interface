from state_interface.state import STATE

import os
from ase.build import fcc100, add_adsorbate

slab = fcc100('Al', size=(1,2,3))
add_adsorbate(slab, 'Cl', 1.5, 'hollow')
slab.cell = [
  [2.86378246381, 0.0, 0.0],
  [0.0, 2.86378246381, 0.0],
  [0.0, 0.0, 5.72756492761]
]

# ---------------- set up executable ----------------
label = 'AlCl'
input_file = label+'.in'
output_file = label+'.out'
no_cpus = 12
pw_loc = './STATE'


## parallel qe using srun (for slurm system)
os.environ['ASE_STATE_COMMAND'] = f'mpirun -np {no_cpus} {pw_loc} < {input_file} > {output_file}'
# -------------- set up input parameters --------------


# -------------- set up input parameters --------------
input_data = {'WF_OPT'       :    'DAV' ,
              'GEO_OPT'      :    'GDIIS',
              'NSPG'         :       1  ,
              'GMAX'         :     4.00 ,
              'GMAXP'        :     10.00 ,
              'SMEARING'     : 'MP',
              'KPOINT_MESH'  : [4, 4, 1],
              'KPOINT_SHIFT' : ["OFF", "OFF", "OFF"],
              'WIDTH'        : 0.002,
              'EDELTA'       : "1.000D-09",
              'NEG'          : 16,
              'MIX'          : 'BROYDEN2',
              'MIX_ALPHA'    : 0.80,
              'FMAX'         : "1.000D-03",
              'DTIO'         : 600.00,
              'PSEUDOS'      : [['Al', 26.9815, 'pot.Al_pbe1'],
                                ['Cl', 35.4527, 'pot.Cl_pbe1']],
              'SUBSTRATE_IDX': ['1:4']
              }

dft_calc = STATE(label=label, input_data=input_data)
# dft_calc.write_input(slab)
slab.calc = dft_calc
slab.get_potential_energy()
# => it will not coverge