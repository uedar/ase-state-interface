import numpy as np

from ase import Atoms
from ase.io.trajectory import Trajectory
from state_interface.state import STATE

a = 4.0  # approximate lattice constant
b = a / 2

input_data = {}

dft_calc = STATE(label=label,input_data=input_data_

ag = Atoms('Ag',
           cell=[(0, b, b), (b, 0, b), (b, b, 0)],
           pbc=1,
           calculator=dft_calc)  # use EMT potential
cell = ag.get_cell()
traj = Trajectory('Ag.traj', 'w')
for x in np.linspace(0.95, 1.05, 5):
    ag.set_cell(cell * x, scale_atoms=True)
    ag.get_potential_energy()
    traj.write(ag)
