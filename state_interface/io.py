import re
from ase import Atoms
import warnings
import numpy as np
from ase.units import Hartree, Bohr

force_unit = Hartree/Bohr

def read_state_input(input_file):
  file = open(input_file, 'r').read()
  svar, sblk = varblock(input_file)

  atomic_species_str = re.search(r'(?<=&ATOMIC_SPECIES)(.+?\n)(?=&END)', file, flags=re.DOTALL).group().strip()
  atomic_species = [c.split() for c in atomic_species_str.split('\n')]

  cell_str = re.search(r'(?<=&CELL)(.+?\n)(?=&END)', file, flags=re.DOTALL).group().strip()
  cell = Bohr*np.array([c.split() for c in cell_str.split('\n')], dtype=float)
  coord_str = re.search(r'(?<=&ATOMIC_COORDINATES CARTESIAN)(.+?\n)(?=&END)', file, flags=re.DOTALL).group().strip()
  coord = [c.split() for c in coord_str.split('\n')]
  positions = Bohr*np.array([c[:3] for c in coord], dtype=float)

  species_idx = {lst[0]: i+1 for i, lst in enumerate(atomic_species)}
  symbols = [k for k,v in species_idx.items() for c in coord if int(c[5]) == v]
  structure = Atoms(symbols=symbols, positions=positions, cell=cell)
  return structure

def read_state_output(output_file):  
  data = []
  extract = False
  with open (output_file) as fd:
      lines = fd.readlines()
      if 'The calculation has converged' not in str(lines):
          warnings.warn('The calculation has not converged')
      for line in lines:
          if re.search('CONVERGED', line):
              extract = True
          if extract:
              data.append (line.split())
              if re.search('EXIT', line):
                  extract = False
  energy, f_max, f_rms = Hartree*float(data[2][1]), force_unit*float(data[2][2]), force_unit*float(data[2][3])
  force_data = []
  positions = []
  species = []
  for line in data[6:]:
      if (line == []):
          break
      species.append (line[2])
      positions.append (line[3:6])
      force_data.append (line[6:9])
  positions = Bohr*np.array(positions, dtype=float)
  structure = Atoms(symbols=species, positions = positions)
  forces = force_unit*np.array(force_data, dtype=float)
  return structure, energy, forces


def varblock(file):
  """
  Remove comments, captures input variable and block inputs separately
  Input
    STATE input file object 
  Returns
    svar = dictionary of input variables
    sblk = list containing lines of block inputs
  """
  with open(file, 'r') as f:
      lines = f.read().splitlines()

  # Clean comments and trailing spaces
  clean = []
  for line in lines:
      li = line.strip()
      if not li.startswith('#'): clean.append(line.strip())
  lines = clean.copy()
  
  # Seprate variables and blocks
  for idx in range(len(lines)):
      if lines[idx].startswith('&'):
          blkstart_idx = idx
          break
  varlines = lines[:blkstart_idx]
  sblk = lines[blkstart_idx:]

  # Process varlines as dictionary
  svar = dict()
  for line in varlines:
      parts = line.split()
      varname = parts[0].upper()
      svar[varname] = [x.upper() for x in parts[1:]]
      
  return svar, sblk

