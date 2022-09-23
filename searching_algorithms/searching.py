import numpy as np
import itertools


def neighbour_search(bins,bins_shape):

    '''
    Takes the atoms sorted into bins with one extra layer of bins as the input and
    loops through the bins calling the distance calculator.
    '''

    original_shape = (bins_shape[0] -2,bins_shape[1] - 2,bins_shape[2] -2)

    bond_dict = {}

    # loop over the unique pairs excluding edge cases for 3 of the lower faces

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

            for neighbouring_bin in bin_neighbours:

                if neighbouring_bin:

                    bonds = get_distance(current_bin,neighbouring_bin)
                    #bond_dict[index_i] = index_j
                    #bond_dict[index_j] = index_x
                     

    # fill in the missed edge cases by looping over the 3 lower edges which were missed by the first loop

    bin_neighbours = []

    for i,j in np.ndindex(original_shape[0],original_shape[1]):

        i = i + 1
        j = j + 1

        bin_neighbours.extend([bins[i][j][0]])

    for i,k in np.ndindex(original_shape[0],original_shape[2]):

        i = i + 1
        k = k + 1

        bin_neighbours.extend(bins[i][0][k-1:k+2])

    for j,k in np.ndindex(original_shape[1],original_shape[2]):

        j = j + 1
        k = k + 1

        bin_neighbours.extend(bins[0][j-1][k-1:k+2])
        bin_neighbours.extend(bins[0][j][k-1:k+2])
        bin_neighbours.extend(bins[0][j+1][k-1:k+2])

    # still need to call get distance for the edge cases

    #for atoms_info in bin_neighbours:

        #if atoms_info:

            #print('  ')
        

 
def get_distance(bin_info_1,bin_info_2):

    '''
    Calculates the distance matrix for all combinations of positions pos1 and pos2.
    The distance cutoff controls which atoms are considered to be bonded.
    '''

    cutoff = 2.4

    pos1 = bin_info_1.positions
    pos2 = bin_info_2.positions

    indices_1 = bin_info_1.indices
    indices_2 = bin_info_2.indices

    D = (np.linalg.norm(pos2[np.newaxis, :, :] - pos1[:, np.newaxis, :], axis=2)).reshape(-1)
 
    Indices = list(itertools.product(indices_1, indices_2))

    bonded = np.where(D< cutoff)[0]
    
    bonds = []

    for i in bonded:

        bonds.append(Indices[i])

    return bonds

    


    

   



    