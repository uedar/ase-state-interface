"""
Read strucutures from STATE input files
and read results from STATE output files.
"""
import re
import warnings
import numpy as np
from ase import Atoms
from ase.units import Hartree, Bohr

force_unit = Hartree / Bohr


def read_state_input(input_file):
    """ASE module to read STATE input files.

    Args:
        input_file (file-like object or str): Points to the STATE input file

    Returns:
        object: Object based on ASE's atoms module using the read input.
    """
    # USE PARSER to get information
    info = InputParser(input_file)
    info.parse()

    # Initialize ASE atoms info argments
    kwargs = {}

    # Load atomic coordinate (cartesian/crystal) to kwargs
    pos_opt = info.blocks["ATOMIC_COORDINATES"]["opt"]
    if pos_opt in [None, "CART", "CARTESIAN"]:
        # Converted position
        pos = Bohr * info.blocks["ATOMIC_COORDINATES"]["input"]["cps"]
        kwargs["positions"] = pos
    elif pos_opt in ["CRYSTAL", "CRYS"]:
        pos = info.blocks["ATOMIC_COORDINATES"]["input"]["cps"]
        kwargs["scaled_positions"] = pos
    else:
        raise ValueError(f"Invalid atomic coordinate option = {pos_opt}.")

    # Load cell information
    if "CELL" in list(info.variables):
        clist = [float(cnum) for cnum in info.variables["CELL"]]
        # Converted lattice lengths
        lena, lenb, lenc = Bohr * clist[0], Bohr * clist[1], Bohr * clist[2]
        kwargs["cell"] = [lena, lenb, lenc, clist[3], clist[4], clist[5]]
    elif "CELL" in list(info.blocks):
        cell_matrix = info.blocks["CELL"]["input"]
        # Converted cell matrix
        kwargs["cell"] = Bohr * cell_matrix
    else:
        raise ValueError("Cell information not found.")

    # Load Atomic species
    species = info.blocks["ATOMIC_SPECIES"]["input"]["symbol"]
    ityp = info.blocks["ATOMIC_COORDINATES"]["input"]["ityp"]
    symlist = [species[int(i) - 1] for i in ityp]
    kwargs["symbols"] = symlist

    # Initiate ASE Atoms module
    structure = Atoms(**kwargs)
    return structure


def read_state_output(output_file):
    """ASE module to read STATE output files."""

    data = []
    extract = False
    with open(output_file, encoding="utf-8") as file:
        lines = file.readlines()
        if "The calculation has converged" not in str(lines):
            warnings.warn("The calculation has not converged")
        for line in lines:
            if re.search("CONVERGED", line):
                extract = True
                count_empty_line = 0
            if extract:
                data.append(line.split())
                if len(line.split()) == 0:
                    count_empty_line += 1
                if count_empty_line == 2:
                    extract = False
            if extract is False and len(data) != 0:
                energy = Hartree * float(data[2][1])
                #                f_max      = force_unit * float(data[2][2])
                #                f_rms      = force_unit * float(data[2][3])
                force_data = []
                positions = []
                species = []
                for vip_line in data[6:]:
                    if vip_line == []:
                        break
                    if vip_line[0] != "MD:":
                        break
                    species.append(vip_line[2])
                    positions.append(vip_line[3:6])
                    force_data.append(vip_line[6:9])
                positions = Bohr * np.array(positions, dtype=float)
                structure = Atoms(symbols=species, positions=positions)
                forces = force_unit * np.array(force_data, dtype=float)
                data = []
    return structure, energy, forces


class InputParser:
    """Fully parse all input parameters in file."""

    def __init__(self, file_object):
        """Intialize the class."""
        self.file_object = file_object
        self.info = dict()
        # Initialize class attributes
        self.varlines = None
        self.blklines = None
        self.variables = None
        self.blocks = None
        self.subroutine = None

    def parse(self):
        """Launch parsing procedure."""
        self.varblock()
        self.catchvariable()
        self.catchblock()
        self.manager()
        self.errorcheck()
        self.readblock()

    def show(self):
        """Show the variable and block component of the parsed file."""
        print("VARIABLE \n", self.variables)
        print("BLOCK \n", self.blocks)

    def varblock(self):
        """Read and clean input file, and separate variable and block lines."""
        with open(self.file_object, "r", encoding="utf-8") as file:
            lines = file.read().splitlines()

        # Clean comments and trailing spaces
        clean = []
        for line in lines:
            _line = line.strip()
            if not _line.startswith("#"):
                clean.append(line.strip())
        lines = clean.copy()

        # Separate variables and blocks
        for idx, line in enumerate(lines):
            if line.startswith("&"):
                blkstart_idx = idx
                break
        self.varlines = lines[:blkstart_idx]
        self.blklines = lines[blkstart_idx:]

    def catchvariable(self):
        """Create a dictionary from variable lines."""
        self.variables = dict()
        for line in self.varlines:
            parts = line.split()
            varname = parts[0].upper()
            self.variables[varname] = [x.upper() for x in parts[1:]]

    def catchblock(self):
        """Create a dictionary from block lines (not specific)."""
        self.blocks = dict()
        for line in self.blklines:
            linestr = line.upper()

            if linestr.startswith("&") and not linestr.startswith("&END"):
                # Process first line of a block for name
                # and option(if available)
                parts = line.split()
                blkname = parts[0].lstrip("&")
                blkinp, blkopt = [], None
                # For Blocks with additional options
                # (ex. &COORDINATE + CRYSTAL)
                if len(parts) == 2:
                    blkopt = parts[1]

            elif linestr.startswith("&END"):
                # Finalize block dictionary entry
                if blkopt:
                    blkopt = blkopt.upper()
                self.blocks[blkname] = {"opt": blkopt, "input": blkinp}

            else:
                blkinp.append(line)

    def readblock(self):
        """Further process the block dictionary using specific subroutines."""
        # Applies the appopriate read subroutine to each block
        for key in self.blocks:
            if key.upper() in self.subroutine:
                self.subroutine[key.upper()]()

    def manager(self):
        """Define the subroutine that further processes the block dictionary."""

        def atomic_coordinates_read():
            entry = self.blocks["ATOMIC_COORDINATES"]
            lines = entry["input"].copy()
            entry["input"] = {}
            # The following are the EXACT syntax naming in
            _cps = [line.split()[:3] for line in lines]
            entry["input"]["cps"] = np.array(_cps, dtype=float)
            _iwei = [line.split()[3] for line in lines]
            entry["input"]["iwei"] = _iwei
            _imdtyp = [line.split()[4] for line in lines]
            entry["input"]["imdtyp"] = _imdtyp
            _ityp = [line.split()[5] for line in lines]
            entry["input"]["ityp"] = _ityp

        def atomic_species_read():
            entry = self.blocks["ATOMIC_SPECIES"]
            lines = entry["input"].copy()
            entry["input"] = {}
            _symbol, _mass, _file = [], [], []
            for line in lines:
                # For the moment atomic number is not supported
                _symbol.append(line.split()[0])
                _mass.append(line.split()[1])
                _file.append(line.split()[2])

            entry["input"]["symbol"] = _symbol
            entry["input"]["mass"] = _mass
            entry["input"]["file"] = _file

        def cell_read():
            entry = self.blocks["CELL"]
            _vector = [line.split() for line in entry["input"]]
            entry["input"] = np.array(_vector, dtype=float)

        self.subroutine = {
            "ATOMIC_COORDINATES": atomic_coordinates_read,
            "ATOMIC_SPECIES": atomic_species_read,
            "CELL": cell_read,
        }

    def errorcheck(self):
        """Error definition and check of code breaking circumstances."""
        # CELL Block & CELL Variable incomptibility
        if "CELL" in self.variables and "CELL" in self.blocks:
            raise ValueError("Double CELL declaration in block and variable.")

        # ATOMIC_COORDINATES not declared
        if "ATOMIC_COORDINATES" not in list(self.blocks):
            raise ValueError("ATOMIC_COORDINATES is not declared.")
