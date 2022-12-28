import os
from ase import Atoms
from ase.db import connect
import numpy as np
from state_interface.state import STATE


def test_scf():
    """
    test SCF calculation
    """
    atoms_obj = Atoms("He")
    atoms_obj.set_cell(10.0 * np.identity(3))
    label = "He"
    input_file = label + ".in"
    output_file = label + ".out"
    pw_loc = "./STATE"
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    os.environ["ASE_STATE_COMMAND"] = f"{pw_loc} < {input_file} > {output_file}"
    pw_origin = "/state-5.6.10/src/STATE"
    pw_link = os.path.join(os.path.dirname(os.path.abspath(__file__)), "STATE")
    if os.path.exists(pw_link):
        os.remove(pw_link)
    os.symlink(pw_origin, pw_link)
    input_data = {
        "GMAX": 5,
        "GMAXP": 20,
        "KPOINT_MESH": [1, 1, 1],
        "KPOINT_SHIFT": ["OFF", "OFF", "OFF"],
        "XCTYPE": "ggapbe",
        "NEG": 8,
        "PSEUDOS": [["He", 4.0, "pot.He_lda1TM"]],
    }

    dft_calc = STATE(label=label, input_data=input_data)
    atoms_obj.calc = dft_calc
    atoms_obj.get_potential_energy()

    db_file = connect("result.json")
    db_file.write(atoms=atoms_obj)
