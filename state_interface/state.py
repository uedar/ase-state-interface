"""STATE Calculator
export ASE_STATE_COMMAND="/path/to/STATE -in PREFIX.in > PREFIX.out"
Run STATE jobs.
"""
import warnings
from ase import io
from ase.calculators.calculator import FileIOCalculator, PropertyNotPresent
from ase.atoms import Atoms
from ase.calculators.singlepoint import (SinglePointDFTCalculator)
import re
import numpy as np
from ase.data import chemical_symbols

error_template = 'Property "%s" not available. Please try running STATE\n' \
                 'first by calling Atoms.get_potential_energy().'
warn_template = 'Property "%s" is None. Typically, this is because the ' \
                'required information has not been printed by STATE ' \
                'at a "low" verbosity level (the default). ' \
                'Please try running STATE with "high" verbosity.'
Ha2eV = 27.2114
Bohr2Ang = 0.529177
force_unit = Ha2eV/Bohr2Ang

class STATE(FileIOCalculator):
    """
    """
    implemented_properties = ['energy', 'forces']
    command = 'STATE < PREFIX.in > PREFIX.out'
    discard_results_on_any_change = True
    def __init__(self, label='state', atoms=None, **kwargs):
        FileIOCalculator.__init__(self,
                                  label, atoms, **kwargs)
        self.calc = None
        self.label= label
    
    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        NATM      = atoms.get_global_number_of_atoms()
        COORDS    = atoms.get_positions() * 1/Bohr2Ang
        SPECIES   = list(set(atoms.get_chemical_symbols()))
        ATOMIC_NUMS = atoms.get_atomic_numbers()
        NTYP      = len(SPECIES)
        CELL      = atoms.get_cell() * 1/Bohr2Ang
        
        with open (self.label + '.in', 'w') as fd:
            input_data = self.parameters['input_data']
            PSEUDO_IDX = {lst[0]: i+1 for i, lst in enumerate(input_data['PSEUDOS'])}
            print ('#', file = fd)
            print ('#', file = fd)
            print ('#', file = fd)
            print("WF_OPT", input_data['WF_OPT'], file = fd)
            if 'GEO_OPT' in input_data:
                print("GEO_OPT",   input_data['GEO_OPT'], file = fd)
            print("NTYP",   NTYP, file = fd)
            print("NATM",   NATM, file = fd)
            if 'TYPE' in input_data:
                print("TYPE",   input_data['TYPE'], file = fd)
            print("NSPG",   input_data['NSPG'], file = fd)
            if 'VERBOSITY' in input_data:
                print("VERBOSITY",   input_data['VERBOSITY'], file = fd)
            print("GMAX",   input_data['GMAX'], file = fd)
            print("GMAXP",  input_data['GMAXP'], file = fd)
            if 'NSCF' in input_data:
                print("NSCF",   input_data['NSCF'], file = fd)
            print("KPOINT_MESH",  *input_data['KPOINT_MESH'], file = fd)
            print("KPOINT_SHIFT", *input_data['KPOINT_SHIFT'], file = fd)
            if 'SMEARING' in input_data:
                print("SMEARING",   input_data['SMEARING'], file = fd)
            print("WIDTH",  input_data['WIDTH'], file = fd)
            if 'MIX' in input_data:
                print("MIX",   input_data['MIX'], file = fd)
            if 'MIX_ALPHA' in input_data:
                print("MIX_ALPHA",   input_data['MIX_ALPHA'], file = fd)
            if 'DTIO' in input_data:
                print("DTIO",   input_data['DTIO'], file = fd)
            if 'FMAX' in input_data:
                print("FMAX",   input_data['FMAX'], file = fd)
            print("EDELTA", input_data['EDELTA'], file = fd)
            print("NEG",    input_data['NEG'], file = fd)
            if 'XTYPE' in input_data:
                print("XTYPE",   input_data['XTYPE'], file = fd)
            print ("&CELL", file = fd)
            for cell_vector in CELL:
                print (*cell_vector, file = fd)
            print ("&END", file = fd)
            print("&ATOMIC_SPECIES", file = fd)
            for line in input_data['PSEUDOS']:
                print(*line, file = fd)
            print("&END", file = fd)
            print("&ATOMIC_COORDINATES", "CARTESIAN", file = fd)
            for i,line in enumerate(COORDS):
                cs = chemical_symbols[ATOMIC_NUMS[i]]
                print(*line, 1,1,PSEUDO_IDX[cs], file=fd)
            print ("&END", file = fd)
            if 'VDW-DF' in input_data:
                print('&VDW-DF', file = fd)
                print("QCUT",  input_data['VDW-DF']['QCUT'], file = fd)
                print("NQ",  input_data['VDW-DF']['NQ'], file = fd)
                print('&END', file = fd)
            
    def read_results(self):
        
        data = []
        extract = False
        with open (self.label + '.out') as fd:
            for line in fd:
                if re.search('CONVERGED', line):
                    extract = True
                if extract:
                    data.append (line.split())
                    if re.search('EXIT', line):
                        extract = False
        energy, f_max, f_rms = Ha2eV*float(data[2][1]), force_unit*float(data[2][2]), force_unit*float(data[2][3])
        force_data = []
        positions = []
        species = []
        for line in data[6:]:
            if (line == []):
                break
            species.append (line[2])
            positions.append (line[3:6])
            force_data.append (line[6:9])
        positions = Bohr2Ang*np.array(positions, dtype=float)
        structure = Atoms(symbols=species, positions = positions)
        forces = force_unit*np.array(force_data, dtype=float)
        calc = SinglePointDFTCalculator(structure, energy=energy,
                                        forces=forces)
        
        self.calc = calc
        self.results = calc.results
    
    def get_stress(self, atoms):
        return np.zeros((3,3))
