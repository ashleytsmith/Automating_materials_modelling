import copy
import numpy as np
import datetime

from searching_algorithms import geometry


def sort_into_bins(positions, symbols, indices, cell, bins_shape, repeat):
    
    '''
    Splits atoms into 3D bins by binning along each of the specified basis vectors.
    '''

    a, b, c = geometry.get_cell_vectors(cell)

    slicing_order = ('a', 'b', 'c')

    scaled_positions = geometry.convert_to_scaled_positions(cell, positions)

    atoms_info = create_atom_info_array(
        positions, scaled_positions, symbols, indices)

    single_sliced_positions = slice_along_axis(
        slicing_order[0], atoms_info, bins_shape[0])

    bins = []

    for stuff in single_sliced_positions:

        sub_bin = slice_along_axis(slicing_order[1], stuff, bins_shape[1])

        triple_sliced_positions = []

        for things in sub_bin:

            sub_sub_bin = slice_along_axis(
                slicing_order[2], things, bins_shape[2])
            sub_sub_bin = remove_empty_lists(sub_sub_bin)

            triple_sliced_positions.append(sub_sub_bin)

        bins.append(triple_sliced_positions)

    if repeat:

        bins, bins_shape = add_neighbouring_bins(bins, bins_shape, a, b, c)

    bins = vectorise_bin_data(bins, bins_shape)

    return bins, bins_shape


def create_atom_info_array(positions, scaled_positions, symbols, indices):
   
    '''
    Creates a structured array (very similar to a struct in C) with all required variables for each atom accesible from the same element.
    '''

    data_types = np.dtype([('x', 'f'), ('y', 'f'), ('z', 'f'),
                          ('a', 'f'), ('b', 'f'), ('c', 'f'), ('i', 'i'), ('s', 'U10')])
    array_length = len(positions)
    atoms_info = np.empty(array_length, dtype=data_types)

    x, y, z = zip(*[(pos[0], pos[1], pos[2]) for pos in positions])
    a, b, c = zip(*[(pos[0], pos[1], pos[2]) for pos in scaled_positions])
    i = [i for i in indices]
    s = [s for s in symbols]

    atoms_info['x'] = x
    atoms_info['y'] = y
    atoms_info['z'] = z
    atoms_info['a'] = a
    atoms_info['b'] = b
    atoms_info['c'] = c
    atoms_info['i'] = i
    atoms_info['s'] = s

    return atoms_info


def slice_along_axis(axis, atoms_info, number_of_bins):
   
    '''
    Seperates the inputs into different bins based on their positions.
    '''

    bins = [[] for stuff in range(number_of_bins)]

    for stuff in atoms_info:

        scaled_position = stuff[axis]

        bin_index = assign_to_bin(scaled_position, number_of_bins)

        bins[bin_index].append(stuff)

    return bins


def assign_to_bin(projection, number_of_bins):
    
    '''
    Splits the values along an axis into bins based on the size of the projection.
    '''

    if projection < 0:

        bin_index = 0

    elif projection > 1:

        bin_index = number_of_bins - 1

    else:

        bin_index = np.floor(projection * number_of_bins).astype(int)

    return bin_index


def remove_empty_lists(list_of_lists):
    
    '''
    Replaces empty lists with None.
    '''

    for i, contents in enumerate(list_of_lists):

        if not contents:

            list_of_lists[i] = None

    return list_of_lists


def add_neighbouring_bins(bins, bins_shape, a, b, c):
    
    '''
    Ensures periodic boundary conditions by adding one layer of bins around the outside.
    Works for any complete rectanglar array of bins with shape bins_shape.
    '''

    original_array_size = bins_shape
    bins = np.asarray(bins,dtype='object')
    
    # pad with shifting

    # x padding
   
    x_bottom_slice = apply_shift(bins[0, :, :], a)
    x_bottom_slice = x_bottom_slice[np.newaxis,:,:]

    x_top_slice = apply_shift(bins[-1, :, :], -a)
    x_top_slice = x_top_slice[np.newaxis,:,:]

    x_padded = np.concatenate([x_top_slice, bins, x_bottom_slice], axis = 0)

    # y padding

    y_bottom_slice = apply_shift(x_padded[:, 0, :], b)
    y_bottom_slice = y_bottom_slice[:,np.newaxis,:]

    
    y_top_slice = apply_shift(x_padded[:, 0, :], -b)
    y_top_slice = y_top_slice[:,np.newaxis,:]

    x_y_padded = np.concatenate([y_top_slice, x_padded, y_bottom_slice], axis = 1)

    # z padding

    z_bottom_slice = apply_shift(x_y_padded[:, :, 0], c)
    z_bottom_slice = z_bottom_slice[:,:,np.newaxis]

    z_top_slice = apply_shift(x_y_padded[:, :, -1], -c)
    z_top_slice = z_top_slice[:,:,np.newaxis]

    bins = np.concatenate([z_top_slice, x_y_padded, z_bottom_slice], axis = 2)
    
    bins_shape = (
        original_array_size[0] + 2, original_array_size[0] + 2, original_array_size[0] + 2)

    return bins, bins_shape



def apply_shift(bins_to_copy, shift):

    '''
    Shift all menbers of each bin by the shift vector.
    '''

    shape = np.shape(bins_to_copy)
    bins = np.empty(shape, dtype='object')
    data_types = np.dtype([('x', 'f'), ('y', 'f'), ('z', 'f'), ('i', 'i'), ('s', 'U10')])

    x_shift = shift[0]
    y_shift = shift[1]
    z_shift = shift[2]

    for i,j in np.ndindex(shape):

        bin_contents_to_copy = bins_to_copy[i][j]
    
        if bin_contents_to_copy:

            array_length = len(bin_contents_to_copy)
            bins[i][j] = []

            for k in range(0,array_length):

                atom_info = np.empty(1, dtype = data_types)
                atom_info['x'] = bin_contents_to_copy[k]['x'] + x_shift
                atom_info['y'] = bin_contents_to_copy[k]['y'] + y_shift
                atom_info['z'] = bin_contents_to_copy[k]['z'] + z_shift
                atom_info['i'] = bin_contents_to_copy[k]['i']
                atom_info['s'] = bin_contents_to_copy[k]['s']
                bins[i][j].append(atom_info[0])

        else:

            bins[i][j] = None

    return bins


def vectorise_bin_data(bins, shape):
   
    '''
    Vectorises the atoms_info data in each 3D bin.
    Returns numpy arrays of positions, indices amd symbols.
    '''

    for i, j, k in np.ndindex(shape):

        bin_contents = bins[i][j][k]

        if bin_contents:

            positions = []
            indices = []
            symbols = []

            for stuff in bin_contents:

                position = (stuff['x'], stuff['y'], stuff['z'])
                index = stuff['i']
                symbol = stuff['s']

                positions.append(position)
                indices.append(index)
                symbols.append(symbol)

            positions = np.asarray(positions)
            indices = np.asarray(indices)
            symbols = np.asarray(symbols)

            atoms_info = AtomsInfo(positions, indices, symbols)

            bins[i][j][k] = atoms_info

    return bins


class AtomsInfo:

    __slots__ = ["positions", "indices", "symbols"]

    def __init__(self, positions, indices, symbols):

        self.positions = positions
        self.indices = indices
        self.symbols = symbols
