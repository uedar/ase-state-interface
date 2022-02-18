import re
from ase import Atoms
import warnings
import numpy as np
from ase.units import Hartree, Bohr
import sys

force_unit = Hartree/Bohr


def read_state_input(input_file):
    # file = open(input_file, 'r').read()
    # USE PARSER to get information
    input = InputParser(input_file)
    input.InitParse()
    
    # Initialize ASE atoms input argments
    kwargs = {}
    
    # Load atomic coordinate (cartesian/crystal) to kwargs
    pos_opt = input.blocks['ATOMIC_COORDINATES']['opt']
    if  pos_opt in [None, 'CART', 'CARTESIAN']:
        pos = Bohr * input.blocks['ATOMIC_COORDINATES']['input']['cps'] ## converted
        kwargs['positions'] = pos
    elif pos_opt in ['CRYSTAL','CRYS']:
        pos = input.blocks['ATOMIC_COORDINATES']['input']['cps']
        kwargs['scaled_positions'] = pos
    else:
        raise ValueError(f"Invalid atomic coordinate option = {pos_opt}.")
        
    # Load cell information 
    if 'CELL' in list(input.variables.keys()):
        clist = list(map(float, input.variables['CELL'])) 
        lena, lenb, lenc = Bohr*clist[0], Bohr*clist[1], Bohr*clist[2] ## converted
        angbc, angac, angab = clist[3], clist[4], clist[5]
        kwargs['cell'] = [lena, lenb, lenc, angbc, angac, angab]
    elif 'CELL' in list(input.blocks.keys()):
        cell_matrix = input.blocks['CELL']['input']
        kwargs['cell'] = Bohr * cell_matrix ## converted
    else:
        raise ValueError(f"Cell information not found.")
        
    # Load Atomic species
    species = input.blocks['ATOMIC_SPECIES']['input']['symbol']
    ityp = input.blocks['ATOMIC_COORDINATES']['input']['ityp']
    symlist = [species[int(i)-1] for i in ityp]
    kwargs['symbols'] = symlist
    
    # Initiate ASE Atoms module
    structure = Atoms(**kwargs)
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
    Interpretation of inputs are not done here.
    '''
    def __init__(self, FileObject):
        self.FileObject = FileObject
        self.input = dict()

    def InitParse(self):
        self.VarBlock()
        self.CatchVariable()
        self.CatchBlock()
        self.SubroutineManager()
        self.ErrorCheck()
        self.ReadBlock()
        
    def Show(self):
        # Shows the main component of the parsed stuff.
        print("VARIABLE \n", self.variables)
        print("BLOCK \n", self.blocks)

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
            if linestr.startswith('&'):
                # Process first line of a block for name and option(if available)
                parts = line.split()
                blkname = parts[0].lstrip('&')
                blkinp, blkopt = [], None
                # For Blocks with additional options (ex. &COORDINATE + CRYSTAL)
                if len(parts) == 2:
                    blkopt = parts[1]
            else:
                blkinp.append(line)

    def ReadBlock(self):
        # Applies the appopriate read subroutine to each block
        for key in self.blocks.keys():
            if key in self.subroutine.keys():
                self.subroutine[key.upper()]()


    def SubroutineManager(self):
        '''
        Scalable reader for specific input blocks
        No interpretation of input done. 
        '''
        def ATOMIC_COORDINATES_read():           
            entry = self.blocks['ATOMIC_COORDINATES']
            lines = entry['input'].copy()
            entry['input'] = dict()
            # The following are the EXACT syntax naming in 
            _cps = [line.split()[:3] for line in lines]
            entry['input']['cps'] = np.array(_cps, dtype=float)
            _iwei = [line.split()[3] for line in lines] 
            entry['input']['iwei'] = _iwei
            _imdtyp = [line.split()[4] for line in lines] 
            entry['input']['imdtyp'] = _imdtyp
            _ityp = [line.split()[5] for line in lines] 
            entry['input']['ityp'] = _ityp

        def ATOMIC_SPECIES_read():   
            entry = self.blocks['ATOMIC_SPECIES']
            lines = entry['input'].copy()
            entry['input'] = dict()
            _symbol, _mass, _file = [], [], []
            for line in lines:
                # For the moment atomic number is not supported 
                _symbol.append(line.split()[0])
                _mass.append(line.split()[1])
                _file.append(line.split()[2])
                
            entry['input']['symbol'] = _symbol
            entry['input']['mass'] = _mass
            entry['input']['file'] = _file
       

        def CELL_read():
            entry = self.blocks['CELL']
            _vector = [line.split() for line in entry['input'] ]
            entry['input'] = np.array(_vector, dtype=float)
            
        self.subroutine = {
            'ATOMIC_COORDINATES': ATOMIC_COORDINATES_read,
            'ATOMIC_SPECIES': ATOMIC_SPECIES_read,
            'CELL': CELL_read,
        }
        
    def ErrorCheck(self):
        
        #CELL Block & CELL Variable incomptibility
        if 'CELL' in self.variables.keys() and 'CELL' in self.blocks.keys():
            raise ValidationError("Simultaneous declaration of CELL in Block and variable is not allowed. ")
        
        # ATOMIC_COORDINATES not declared
        if 'ATOMIC_COORDINATES' not in list(self.blocks.keys()):
            raise ValidationError("ATOMIC_COORDINATES is not declared.")
        