from ase.io import read


def read_trajfile(infile):

    atoms = read(infile)

    return atoms


def get_atom_properties(atoms):

    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    cell = atoms.get_cell()

    return positions, symbols, cell


def save_file(atoms, name=None):

    if name == None:

        atoms.write('in.traj')

    else:

        atoms.write(name)
