import numpy as np
import itertools


def neighbour_search(bins,bins_shape,number_of_atoms,cutoff):

    '''
    Takes the atoms sorted into bins with one extra layer of bins as the input and
    loops through the bins calling the distance calculator.
    '''

    original_shape = (bins_shape[0] -2,bins_shape[1] - 2,bins_shape[2] -2)

    bond_dict = {}

    for i in range(0,number_of_atoms):

        bond_dict[i] = []

    # loop over all unique pairs only once.  

    for i,j,k in np.ndindex(original_shape):

        i = i + 1
        j = j + 1
        k = k + 1

        current_bin = bins[i][j][k]

        if current_bin:

            bin_neighbours = []

            # slicing routine using slices along the z axis
            
            bin_neighbours.extend(bins[i+1][j-1][k-1:k+2])
            bin_neighbours.extend(bins[i+1][j][k-1:k+2])
            bin_neighbours.extend(bins[i+1][j+1][k-1:k+2])
            bin_neighbours.extend(bins[i][j+1][k-1:k+2])
            bin_neighbours.extend([bins[i][j][k+1]])

            # local term

            bonds = get_distance(current_bin,current_bin,cutoff)
            bond_dict = add_to_bond_dict(bonds,bond_dict, bothways = False)

            # neighbour terms

            for neighbouring_bin in bin_neighbours:

                if neighbouring_bin:

                    bonds = get_distance(current_bin,neighbouring_bin,cutoff)
                    bond_dict = add_to_bond_dict(bonds,bond_dict, bothways = True)
                             
    return bond_dict

 
def get_distance(bin_info_1,bin_info_2,cutoff):

    '''
    Calculates the distance matrix for all combinations of positions pos1 and pos2.
    The distance cutoff controls which atoms are considered to be bonded.
    '''

    pos1 = bin_info_1.positions
    pos2 = bin_info_2.positions

    indices_1 = bin_info_1.indices
    indices_2 = bin_info_2.indices

    D = (np.linalg.norm(pos2[np.newaxis, :, :] - pos1[:, np.newaxis, :], axis=2)).reshape(-1)

    Indices = list(itertools.product(indices_1, indices_2))

    bonded = np.where((0.01 < D) & (D < cutoff))[0]

    bonds = []

    for i in bonded:

        bonds.append(Indices[i])

    return bonds



def add_to_bond_dict(bonds,bond_dict, bothways = True):

    '''
    Add both the forwards and backwards pairs to the bond dictionary if bothways is true.
    '''

    if bothways:

        for bond in bonds:

            bond_dict[bond[0]].append(bond[1])
            bond_dict[bond[1]].append(bond[0])

    else:

         for bond in bonds:

            bond_dict[bond[0]].append(bond[1])

    return bond_dict










    

   



    