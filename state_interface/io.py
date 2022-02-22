import re
import warnings
import numpy as np
from ase import Atoms
from ase.units import Hartree, Bohr

force_unit = Hartree/Bohr


def read_state_input(input_file):

    # USE PARSER to get information
    info = InputParser(input_file)
    info.parse()

    # Initialize ASE atoms info argments
    kwargs = {}

    # Load atomic coordinate (cartesian/crystal) to kwargs
    pos_opt = info.blocks['ATOMIC_COORDINATES']['opt']
    if pos_opt in [None, 'CART', 'CARTESIAN']:
        # Converted position
        pos = Bohr * info.blocks['ATOMIC_COORDINATES']['input']['cps']
        kwargs['positions'] = pos
    elif pos_opt in ['CRYSTAL', 'CRYS']:
        pos = info.blocks['ATOMIC_COORDINATES']['input']['cps']
        kwargs['scaled_positions'] = pos
    else:
        raise ValueError(f"Invalid atomic coordinate option = {pos_opt}.")

    # Load cell information
    if 'CELL' in list(info.variables.keys()):
        clist = [float(cnum) for cnum in info.variables['CELL']]
        # Converted lattice lengths
        lena, lenb, lenc = Bohr*clist[0], Bohr*clist[1], Bohr*clist[2]
        angbc, angac, angab = clist[3], clist[4], clist[5]
        kwargs['cell'] = [lena, lenb, lenc, angbc, angac, angab]
    elif 'CELL' in list(info.blocks.keys()):
        cell_matrix = info.blocks['CELL']['input']
        # Converted cell matrix
        kwargs['cell'] = Bohr * cell_matrix
    else:
        raise ValueError('Cell information not found.')

    # Load Atomic species
    species = info.blocks['ATOMIC_SPECIES']['input']['symbol']
    ityp = info.blocks['ATOMIC_COORDINATES']['input']['ityp']
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
        if line == []:
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
    def __init__(self, file_object):
        self.file_object = file_object
        self.info = dict()
        # Initialize class attributes
        self.varlines = None
        self.blklines = None
        self.variables = None
        self.blocks = None
        self.subroutine = None

    def parse(self):
        self.varblock()
        self.catchvariable()
        self.catchblock()
        self.manager()
        self.errorcheck()
        self.readblock()

    def show(self):
        # Shows the main component of the parsed stuff.
        print("VARIABLE \n", self.variables)
        print("BLOCK \n", self.blocks)

    def varblock(self):
        '''
        Remove comments, captures info variable and block inputs separately
        '''
        with open(self.file_object, 'r') as f:
            lines = f.read().splitlines()

        # Clean comments and trailing spaces
        clean = []
        for line in lines:
            _line = line.strip()
            if not _line.startswith('#'):
                clean.append(line.strip())
        lines = clean.copy()

        # Separate variables and blocks
        for idx, line in enumerate(lines):
            if line.startswith('&'):
                blkstart_idx = idx
                break
        self.varlines = lines[:blkstart_idx]
        self.blklines = lines[blkstart_idx:]

    def catchvariable(self):
        # Process the variable lines to a dictionary
        self.variables = dict()
        for line in self.varlines:
            parts = line.split()
            varname = parts[0].upper()
            self.variables[varname] = [x.upper() for x in parts[1:]]

    def catchblock(self):
        # Process block lines to a dictionary
        # Parsing only separates blocks into dictionary
        # Needs specific block reader (provided by manager)
        self.blocks = dict()
        for line in self.blklines:
            linestr = line.upper()

            if linestr.startswith('&') and not linestr.startswith('&END'):
                # Process first line of a block for name
                # and option(if available)
                parts = line.split()
                blkname = parts[0].lstrip('&')
                blkinp, blkopt = [], None
                # For Blocks with additional options
                # (ex. &COORDINATE + CRYSTAL)
                if len(parts) == 2:
                    blkopt = parts[1]

            elif linestr.startswith('&END'):
                # Finalize block dictionary entry
                if blkopt:
                    _opt = blkopt.upper()
                self.blocks[blkname] = {'opt': _opt, 'input': blkinp}

            else:
                blkinp.append(line)

    def readblock(self):
        # Applies the appopriate read subroutine to each block
        for key in self.blocks.keys():
            if key in self.subroutine.keys():
                self.subroutine[key.upper()]()


    def manager(self):
        '''
        Scalable reader for specific info blocks
        No interpretation of info done.
        '''
        def atomic_coordinates_read():
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

        def atomic_species_read():
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

        def cell_read():
            entry = self.blocks['CELL']
            _vector = [line.split() for line in entry['input']]
            entry['input'] = np.array(_vector, dtype=float)

        self.subroutine = {
            'ATOMIC_COORDINATES': atomic_coordinates_read,
            'ATOMIC_SPECIES': atomic_species_read,
            'CELL': cell_read,
        }

    def errorcheck(self):

        # CELL Block & CELL Variable incomptibility
        if 'CELL' in self.variables.keys() and 'CELL' in self.blocks.keys():
            raise ValueError("Double CELL declaration in block and variable.")

        # ATOMIC_COORDINATES not declared
        if 'ATOMIC_COORDINATES' not in list(self.blocks.keys()):
            raise ValueError("ATOMIC_COORDINATES is not declared.")
