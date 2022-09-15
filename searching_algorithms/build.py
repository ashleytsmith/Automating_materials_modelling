import ase
from ase.io import *
from ase import Atoms
import sys

import math
import numpy as np


def replace_atom(atoms,index,new_symbol):

    '''
    Replaces the atom at a given index by the new symbol keeping the original indices intact.
    '''
    
    pos = atoms.get_positions()
    sym = atoms.get_chemical_symbols()
    #cell = atoms.get_cell()

    sym.pop(index)
    sym.insert(index, new_symbol)

    atoms = Atoms(sym, positions=pos)
    #atoms.set_cell(cell)
    #atoms.set_pbc([True, True, True])

    return atoms


def add_atom(atoms,new_atom):

    '''
    Add atom to atoms object.
    '''

    atoms.extend(new_atom)

    return atoms

    