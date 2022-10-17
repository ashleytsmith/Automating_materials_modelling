import pickle


def run(center):

    '''
    Walks through the structure starting from a given atom returning the neighbour distances of all other atoms relative to it.
    '''

    infile = open('bonding_info.pkl','rb')
    bond_dict = pickle.load(infile)
    infile.close()

    reduced_bond_dict = skip_bridges(bond_dict)

    neighbour_dict = breadth_first_search(reduced_bond_dict, center)

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
                


def breadth_first_search(bond_dict, search_center):  

    '''
    Walks breadth first through the structure identifying which shells the other atoms are in relative to the starting atom.
    '''

    neighbour_dict = {}

    number_of_atoms = len(bond_dict)

    found_already = [search_center]
    old_neighbours = [search_center]

    shell = 1

    while len(found_already) < number_of_atoms - 1:

        neighbour_dict[shell] = []
        new_neighbours = []

        for n in old_neighbours:
            
            neighbours = bond_dict[n]

            for m in neighbours:

                if m not in found_already:

                    neighbour_dict[shell].append(m)
                    found_already.append(m)
                    new_neighbours.append(m)


        old_neighbours = new_neighbours
        shell += 1

    return neighbour_dict

        




    # found_already = list(old_neighbours)
# while len(found_already) < number_of_atoms_in_structure - 1:
# current paths
   

    



