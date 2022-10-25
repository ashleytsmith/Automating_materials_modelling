from ase import Atoms


def repeat_structure(atoms, rep):
    
    atoms_repeated = atoms.repeat(rep)

    return atoms_repeated


def repeat_structures(atoms, reps):

    structures = []

    for rep in reps:
    
        atoms_repeated = repeat_structure(atoms, rep)
        structures.append(atoms_repeated)

    return structures


def replace_atom(atoms,index,new_symbol):

    '''
    Replaces the atom at a given index by the new symbol keeping the original indices intact.
    '''
    
    pos = atoms.get_positions()
    sym = atoms.get_chemical_symbols()
    cell = atoms.get_cell()

    sym.pop(index)
    sym.insert(index, new_symbol)

    atoms = Atoms(sym, positions=pos)
    atoms.set_cell(cell)
    atoms.set_pbc([True, True, True])

    return atoms


def pack_atom_objects(atoms_list, cell):

    '''
    Packs a list of atom objects back into one atoms object.
    '''

    pos = [a.position for a in atoms_list]
    sym = [a.symbol for a in atoms_list]
    atoms = Atoms(sym, positions=pos)
    atoms.set_cell(cell)
    atoms.set_pbc([True, True, True])

    return atoms







