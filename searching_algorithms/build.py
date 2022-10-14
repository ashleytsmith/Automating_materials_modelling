
def repeat_structure(atoms, reps):

    structures = []

    for rep in reps:
    
        atoms_repeated = atoms.repeat(rep)
        structures.append(atoms_repeated)

    return structures


