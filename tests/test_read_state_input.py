"""This is a test of state-interface/io.read_state_input"""
import os
from state_interface.io import read_state_input

base = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.normpath(os.path.join(base, "../data/io_test"))


def test_read_al_input():
    """test Al input file"""
    data_path = os.path.join(TEST_DIR, "test_al/nfinp_scf")
    atoms = read_state_input(data_path)
    assert len(atoms.cell) == 3
    assert atoms.symbols == ["Al"]
    assert len(atoms.positions) == 1
