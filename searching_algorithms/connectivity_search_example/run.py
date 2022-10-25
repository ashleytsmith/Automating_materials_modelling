import pickle

from searching_algorithms.connectivity_search import breadth_first_search


def run(center, skip = True):

    '''
    Walks through the structure starting from a given atom returning the neighbour distances of all other atoms relative to it.
    '''

    infile = open('Results/bonding_info.pkl','rb')
    bond_dict = pickle.load(infile)
    infile.close()

    if skip:

        bond_dict = skip_bridges(bond_dict)

    neighbour_dict = breadth_first_search(bond_dict, center)

    return neighbour_dict

    

def skip_bridges(bond_dict):

    '''
    Removes bridging atoms (double cooridinated) from the bond dictinary.
    '''

    reduced_bond_dict = {}

    number_of_atoms = len(bond_dict)

    for i in range(0,number_of_atoms):

        current_neighbours = bond_dict[i]
        number_of_current_neighbours = len(bond_dict[i])

        if number_of_current_neighbours > 2:

            next_neighbours = []

            for j in current_neighbours:

                next_neighbour = [n for n in bond_dict[j] if n!=i][0]
                next_neighbours.append(next_neighbour)

            reduced_bond_dict[i] = next_neighbours
           

    return reduced_bond_dict