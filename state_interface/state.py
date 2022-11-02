"""STATE Calculator
export ASE_STATE_COMMAND="/path/to/STATE < PREFIX.in > PREFIX.out"
Run STATE jobs.
"""

from ase.units import Bohr, Hartree
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.data import chemical_symbols
import numpy as np
from .io import read_state_output

ERROR_TEMPLATE = (
    'Property "%s" not available. Please try running STATE\n'
    "first by calling Atoms.get_potential_energy()."
)

WARN_TEMPLATE = (
    'Property "%s" is None. Typically, this is because the '
    "required information has not been printed by STATE "
    'at a "low" verbosity level (the default). '
    'Please try running STATE with "high" verbosity.'
)

force_unit = Hartree / Bohr


class STATE(FileIOCalculator):
    """
    STATE calculator interface class
    """

    implemented_properties = ["energy", "forces"]
    command = "STATE < PREFIX.in > PREFIX.out"
    discard_results_on_any_change = True

    def __init__(self, label="state", atoms=None, **kwargs):
        FileIOCalculator.__init__(self, label, atoms, **kwargs)
        self.calc = None
        self.label = label

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        natom = atoms.get_global_number_of_atoms()
        coords = atoms.get_positions() * 1 / Bohr
        species = list(set(atoms.get_chemical_symbols()))
        atomic_nums = atoms.get_atomic_numbers()
        ntyp = len(species)
        cell = atoms.get_cell() * 1 / Bohr
        constraints = atoms.constraints
        with open(self.label + ".in", "w", encoding="utf-8") as file:
            input_data = self.parameters["input_data"]
            print("#\n#\n#", file=file)
            print("NTYP", ntyp, file=file)
            print("NATM", natom, file=file)
            for parameter in input_data.keys():
                write_parameter(parameter, input_json=input_data, file=file)
            print("&CELL", file=file)
            for cell_vector in cell:
                print(
                    f"{cell_vector[0]: >18.12f}",
                    f"{cell_vector[1]: >18.12f}",
                    f"{cell_vector[2]: >18.12f}",
                    file=file,
                )
            print("&END", file=file)
            print("&ATOMIC_SPECIES", file=file)
            for line in input_data["PSEUDOS"]:
                print(*line, file=file)
            print("&END", file=file)
            print("&ATOMIC_COORDINATES", "CARTESIAN", file=file)
            pseudo_idx = {
                lst[0]: i + 1 for i, lst in enumerate(input_data["PSEUDOS"])
            }

            constraints_idx = (
                constraints[0].get_indices() if len(constraints) > 0 else []
            )
            for i, line in enumerate(coords):
                chemical_symbol = chemical_symbols[atomic_nums[i]]
                imdtyp = 0 if i in constraints_idx else 1
                print(
                    f"{line[0]: >18.12f}",
                    f"{line[1]: >18.12f}",
                    f"{line[2]: >18.12f}",
                    1,
                    imdtyp,
                    pseudo_idx[chemical_symbol],
                    file=file,
                )
            print("&END", file=file)
            if "VDW-DF" in input_data:
                print("&VDW-DF", file=file)
                print("QCUT", input_data["VDW-DF"]["QCUT"], file=file)
                print("NQ", input_data["VDW-DF"]["NQ"], file=file)
                print("&END", file=file)

    def read_results(self):
        structure, energy, forces = read_state_output(self.label + ".out")
        calc = SinglePointDFTCalculator(
            structure, energy=energy, forces=forces
        )

        self.calc = calc
        self.results = calc.results
        self.results["updt_pos"] = structure.get_positions()

    def get_stress(self, atoms):
        """Dummy class method to get stress"""
        return np.zeros((3, 3))

    def get_updt_pos(self):
        """Return final atoms positions"""
        return self.results["updt_pos"]


def write_parameter(parameter, input_json, file):
    """Helper method to write input parameters"""
    if parameter == "KPOINT_MESH":
        print("KPOINT_MESH", *input_json["KPOINT_MESH"], file=file)
    elif parameter == "KPOINT_SHIFT":
        print("KPOINT_SHIFT", *input_json["KPOINT_SHIFT"], file=file)
    elif parameter in ["PSEUDOS", "SUBSTRATE_IDX", "VDW-DF"]:
        pass
    else:
        print(parameter, input_json[parameter], file=file)
