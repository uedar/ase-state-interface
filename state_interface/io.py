import re
from ase import Atoms
import warnings
import numpy as np
from ase.units import Hartree, Bohr

force_unit = Hartree/Bohr


def read_state_input(input_file):
    file = open(input_file, 'r').read()
    # USE PARSER
    # input = InputParser(input_file)
    # input.InitParse()
    # GET POSITION
    # positions = input.blocks[ATOMIC_COORDINATES]['input']['cps']
    # cell = input.blocks['CELL']['input']

    atomic_species_str = re.search(
        r'(?<=&ATOMIC_SPECIES)(.+?\n)(?=&END)', file, flags=re.DOTALL).group().strip()
    atomic_species = [c.split() for c in atomic_species_str.split('\n')]

    cell_str = re.search(r'(?<=&CELL)(.+?\n)(?=&END)',
                         file, flags=re.DOTALL).group().strip()
    cell = Bohr*np.array([c.split()
                         for c in cell_str.split('\n')], dtype=float)
    print("CELL \n", cell, "\n\n")
    
    
    coord_str = re.search(r'(?<=&ATOMIC_COORDINATES CARTESIAN)(.+?\n)(?=&END)',
                          file, flags=re.DOTALL).group().strip()
    coord = [c.split() for c in coord_str.split('\n')]
    positions = Bohr*np.array([c[:3] for c in coord], dtype=float)
    print("POSITION \n", positions, "\n\n")

    species_idx = {lst[0]: i+1 for i, lst in enumerate(atomic_species)}
    print('species_idx', species_idx)
    symbols = [k for k, v in species_idx.items()
               for c in coord if int(c[5]) == v]
    print('symboles', symbols)
    structure = Atoms(symbols=symbols, positions=positions, cell=cell)
    print(structure)
    return structure


def read_state_output(output_file):
    data = []
    extract = False
    with open(output_file) as fd:
        lines = fd.readlines()
        if 'The calculation has converged' not in str(lines):
            warnings.warn('The calculation has not converged')
        for line in lines:
            if re.search('CONVERGED', line):
                extract = True
            if extract:
                data.append(line.split())
                if re.search('EXIT', line):
                    extract = False
    energy, f_max, f_rms = Hartree * \
        float(data[2][1]), force_unit * \
        float(data[2][2]), force_unit*float(data[2][3])
    force_data = []
    positions = []
    species = []
    for line in data[6:]:
        if (line == []):
            break
        species.append(line[2])
        positions.append(line[3:6])
        force_data.append(line[6:9])
    positions = Bohr*np.array(positions, dtype=float)
    structure = Atoms(symbols=species, positions=positions)
    forces = force_unit*np.array(force_data, dtype=float)
    return structure, energy, forces


class InputParser:
    '''
    Conducts Full parsing of all inputs in file.
    Used to output a dictionary containing all variables and block inputs.
    '''

    def __init__(self, FileObject):
        self.FileObject = FileObject
        self.input = dict()

    def InitParse(self):
        self.VarBlock()
        self.CatchVariable()
        self.CatchBlock()
        self.SubroutineManager()
        self.ReadBlock()
        self.ErrorCheck()
        self.ASECall()

    def VarBlock(self):
        '''
        Remove comments, captures input variable and block inputs separately
        '''
        with open(self.FileObject, 'r') as f:
            lines = f.read().splitlines()

        # Clean comments and trailing spaces
        clean = []
        for line in lines:
            li = line.strip()
            if not li.startswith('#'):
                clean.append(line.strip())
        lines = clean.copy()

        # Separate variables and blocks
        for idx in range(len(lines)):
            if lines[idx].startswith('&'):
                blkstart_idx = idx
                break
        self.varlines = lines[:blkstart_idx]
        self.blklines = lines[blkstart_idx:]

    def CatchVariable(self):
        # Process the variable lines to a dictionary
        self.variables = dict()
        for line in self.varlines:
            parts = line.split()
            varname = parts[0].upper()
            self.variables[varname] = [x.upper() for x in parts[1:]]

    def CatchBlock(self):
        # Process block lines to a dictionary
        # Parsing only separates blocks into dictionary
        # Needs specific block reader (provided by SubroutineManager)
        self.blocks = dict()
        for line in self.blklines:
            linestr = line.upper()

            if linestr.startswith('&END'):
                # Finalize block dictionary entry
                self.blocks[blkname] = {'opt': blkopt.upper() if blkopt is not None else blkopt, 'input': blkinp}
                print(self.blocks[blkname])

            if linestr.startswith('&'):
                # Process first line of a block for name and option(if available)
                parts = line.split()
                blkname = parts[0].lstrip('&')
                blkinp, blkopt = [], None

                # For Blocks with additional options (ex. &COORDINATE + CRYSTAL)
                if len(parts) == 2:
                    blkopt = parts[1]
                print(parts, '...', blkname, '...', blkopt)

            else:
                blkinp.append(line)

    def ReadBlock(self):
        '''
        Applies the appopriate read subroutine to each block
        '''
 
        for key in self.blocks.keys():
            entry = self.blocks[key]
            if key+'_subroutine' in self.subroutine.keys():
                apply = self.subroutine[f"{key.upper()}_subroutine"]
                entry = apply(entry)
            else: 
                print(f"{key} block has no support yet. Skipping ...")

    def SubroutineManager(self):
        '''
        Scalable reader for specific input blocks
        '''

        def ATOMIC_COORDINATES_subroutine(entry):
            # Takes the dictionary entry for ATOMIC_COORDINATES
            if entry['opt'] in [None, 'CARTESIAN', 'CART']:
                
                lines = entry['input'].copy()
                entry['input'] = dict()

                # The following are the EXACT syntax naming in 
                _cps = [line.split()[:3] for line in lines]
                entry['input']['cps'] = Bohr*np.array(_cps, dtype=float)
                _iwei = [line.split()[3] for line in lines] 
                entry['input']['iwei'] = _iwei
                _imdtyp = [line.split()[4] for line in lines] 
                entry['input']['imdtyp'] = _imdtyp
                _ityp = [line.split()[5] for line in lines] 
                entry['input']['ityp'] = _ityp
                
                return entry
            
            else: 
                raise ValueError(f"BLOCK OPTION ERROR: No support found for ATOMIC_COORDINATE option = '{entry['opt']}'")
                
            

        def ATOMIC_SPECIES_subroutine(entry):
            # take the list of lines
            print("ATOMIC SPECIES WAS RUN")

        def CELL_subroutine(entry):
            _vector = [line.split() for line in entry['input'] ]
            entry['input'] = Bohr*np.array(_vector, dtype=float)
            return entry
            
            
        self.subroutine = {
            'ATOMIC_COORDINATES_subroutine': ATOMIC_COORDINATES_subroutine,
            'ATOMIC_SPECIES_subroutine': ATOMIC_SPECIES_subroutine,
            'CELL_subroutine': CELL_subroutine,
        }

    def ASECall(self):
        self.positions = self.blocks['ATOMIC_COORDINATES']['input']['cps']
        self.cell = self.blocks['CELL']['input']
        
        
    # def ErrorCheck(self):