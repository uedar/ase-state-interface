import re
from ase import Atoms
import warnings
from nbformat import ValidationError
import numpy as np
from ase.units import Hartree, Bohr
from ase.geometry.cell import cellpar_to_cell
import sys

force_unit = Hartree/Bohr


def read_state_input(input_file):
    file = open(input_file, 'r').read()
    # USE PARSER
    # input = InputParser(input_file)
    # input.InitParse()
    # GET POSITION
    # positions = input.blocks[ATOMIC_COORDINATES]['input']['cps']
    # cell = input.blocks['CELL']['input']

    input = InputParser(input_file)
    input.InitParse()
    
    kwargs = {}
    
    # Load atomic coordinate (cartesian/crystal) to kwargs
    pos_opt = input.blocks['ATOMIC_COORDINATES']['opt']
    if  pos_opt in [None, 'CART', 'CARTESIAN']:
        pos = Bohr * input.blocks['ATOMIC_COORDINATES']['input']['cps']
        kwargs['positions'] = pos
    elif pos_opt in ['CRYSTAL','CRYS']:
        print(input.blocks['ATOMIC_COORDINATES']['input']['cps'])
        pos = input.blocks['ATOMIC_COORDINATES']['input']['cps']
        kwargs['scaled_positions'] = pos
    else:
        raise ValueError(f"Invalid atomic coordinate option = {pos_opt}.")
        
    # Load cell information 
    if 'CELL' in list(input.variables.keys()):
        clist = list(map(float, input.variables['CELL'])) 
        lena, lenb, lenc = Bohr*clist[0], Bohr*clist[1], Bohr*clist[2]
        angbc, angac, angab = clist[3], clist[4], clist[5]
        kwargs['cell'] = [lena, lenb, lenc, angbc, angac, angab]
    elif 'CELL' in list(input.blocks.keys()):
        cell_matrix = input.blocks['CELL']['input']
        kwargs['cell'] = Bohr * cell_matrix
    else:
        raise ValueError(f"Cell information not found.")
        
    # Load Atomic species
    species = input.blocks['ATOMIC_SPECIES']['input']['symbol']
    ityp = input.blocks['ATOMIC_COORDINATES']['input']['ityp']
    symlist = [species[int(i)-1] for i in ityp]
    kwargs['symbols'] = symlist
    
    
    print(kwargs)
    # Initiate ASE Atoms module
    structure = Atoms(**kwargs)
    
    

    # atomic_species_str = re.search(
    #     r'(?<=&ATOMIC_SPECIES)(.+?\n)(?=&END)', file, flags=re.DOTALL).group().strip()
    # atomic_species = [c.split() for c in atomic_species_str.split('\n')]

    # cell_str = re.search(r'(?<=&CELL)(.+?\n)(?=&END)', file, flags=re.DOTALL).group().strip()
    # cell = Bohr*np.array([c.split() for c in cell_str.split('\n')], dtype=float)
    # print("CELL \n", cell, "\n\n")
    
    
    # coord_str = re.search(r'(?<=&ATOMIC_COORDINATES CARTESIAN)(.+?\n)(?=&END)',
    #                       file, flags=re.DOTALL).group().strip()
    # coord = [c.split() for c in coord_str.split('\n')]
    # positions = Bohr*np.array([c[:3] for c in coord], dtype=float)
    # print("POSITION \n", positions, "\n\n")

    # species_idx = {lst[0]: i+1 for i, lst in enumerate(atomic_species)}
    # print('species_idx', species_idx)
    # symbols = [k for k, v in species_idx.items()
    #            for c in coord if int(c[5]) == v]
    # print('symbols', symbols)
    # structure = Atoms(symbols=symbols, positions=positions, cell=cell)
    # print(structure)
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
        self.ErrorCheck()
        self.ReadBlock()
        # self.ASECall()
        
    def Show(self):
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
                # print(self.blocks[blkname])

            if linestr.startswith('&'):
                # Process first line of a block for name and option(if available)
                parts = line.split()
                blkname = parts[0].lstrip('&')
                blkinp, blkopt = [], None

                # For Blocks with additional options (ex. &COORDINATE + CRYSTAL)
                if len(parts) == 2:
                    blkopt = parts[1]
                # print(parts, '...', blkname, '...', blkopt)

            else:
                blkinp.append(line)

    def ReadBlock(self):
        # Applies the appopriate read subroutine to each block
        for key in self.blocks.keys():
            if key+'_read' in self.subroutine.keys():
                self.subroutine[f"{key.upper()}_read"]()
        
        # self.subroutine['ATOMIC_COORDINATES_read']()
        # self.subroutine['CELL_read']()
        # self.subroutine['ATOMIC_SPECIES_read']()
        
        # for key in self.blocks.keys():
        #     if key+'_read' in self.subroutine.keys():
        #         self.subroutine[f"{key.upper()}_read"]
        #         # apply = self.subroutine[f"{key.upper()}_read"]
        #         # entry = apply(entry)
        #     else: 
        #         print(f"{key} block has no support yet. Skipping ...")

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
            'ATOMIC_COORDINATES_read': ATOMIC_COORDINATES_read,
            'ATOMIC_SPECIES_read': ATOMIC_SPECIES_read,
            'CELL_read': CELL_read,
        }

    # def ASECall(self):
    #     # Simplify calling variables for ASE usage
        
        
        
        
    #     self.positions = self.blocks['ATOMIC_COORDINATES']['input']['cps']
        
    #     if 'CELL' in self.variables.keys():
    #         parts = self.variables['CELL']
    #         abc = (Bohr*np.array(parts[:3], dtype=float)).tolist()
    #         angles = np.array(parts[3:], dtype=float).tolist()
    #         cellpar = abc + angles
    #         self.cell = cellpar_to_cell(cellpar)
    #     else:
    #         self.cell = self.blocks['CELL']['input']
        
        
        
    def ErrorCheck(self):
        
        #CELL Block & CELL Variable incomptibility
        if 'CELL' in self.variables.keys() and 'CELL' in self.blocks.keys():
            raise ValidationError("Simultaneous declaration of CELL in Block and variable is not allowed. ")
        
        # HARD REQUIREMENTS
        # ATOMIC_COORDINATES not declared
        if 'ATOMIC_COORDINATES' not in list(self.blocks.keys()):
            raise ValidationError("ATOMIC_COORDINATES is not declared.")
        