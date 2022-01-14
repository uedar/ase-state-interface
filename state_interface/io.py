import re
from ase import Atoms

def read_state_input(input_file):
  file = open(input_file, 'r').read()

  atomic_species_str = re.search(r'(?<=&ATOMIC_SPECIES)(.+?\n)(?=&END)', file, flags=re.DOTALL).group().strip()
  atomic_species = [c.split() for c in atomic_species_str.split('\n')]

  cell_str = re.search(r'(?<=&CELL)(.+?\n)(?=&END)', file, flags=re.DOTALL).group().strip()
  cell = [c.split() for c in cell_str.split('\n')]

  coord_str = re.search(r'(?<=&ATOMIC_COORDINATES CARTESIAN)(.+?\n)(?=&END)', file, flags=re.DOTALL).group().strip()
  coord = [c.split() for c in coord_str.split('\n')]
  positions = [c[:3] for c in coord]

  species_idx = {lst[0]: i+1 for i, lst in enumerate(atomic_species)}
  symbols = [k for k,v in species_idx.items() for c in coord if int(c[5]) == v]
  structure = Atoms(symbols=symbols, positions=positions, cell=cell)
  return structure