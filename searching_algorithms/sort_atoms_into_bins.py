import copy
import numpy as np

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
    Seperates the input positions into different bins.
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

    inner_layer_lower = 1
    outer_layer_lower = 0

    inner_layer_higher_x = original_array_size[0]
    inner_layer_higher_y = original_array_size[1]
    inner_layer_higher_z = original_array_size[2]

    outer_layer_higher_x = original_array_size[0] + 1
    outer_layer_higher_y = original_array_size[1] + 1
    outer_layer_higher_z = original_array_size[2] + 1

    bins = pad_array(bins, inner_layer_higher_x,
                     inner_layer_higher_y, inner_layer_higher_z)

    # faces

    # x reflections
    for j, k in np.ndindex((inner_layer_higher_x, inner_layer_higher_z)):

        bins[outer_layer_higher_x][j+1][k +
                                        1] = apply_shift(copy.deepcopy(bins[inner_layer_lower][j+1][k+1]), a)
        bins[outer_layer_lower][j+1][k +
                                     1] = apply_shift(copy.deepcopy(bins[inner_layer_higher_x][j+1][k+1]), -a)

    # y reflections
    for i, k in np.ndindex((inner_layer_higher_x, inner_layer_higher_z)):

        bins[i+1][outer_layer_higher_y][k +
                                        1] = apply_shift(copy.deepcopy(bins[i+1][inner_layer_lower][k+1]), b)
        bins[i+1][outer_layer_lower][k +
                                     1] = apply_shift(copy.deepcopy(bins[i+1][inner_layer_higher_y][k+1]), -b)

    # z reflections
    for i, j in np.ndindex((inner_layer_higher_x, inner_layer_higher_y)):

        bins[i+1][j+1][outer_layer_higher_z] = apply_shift(
            copy.deepcopy(bins[i+1][j+1][inner_layer_lower]), c)
        bins[i+1][j+1][outer_layer_lower] = apply_shift(
            copy.deepcopy(bins[i+1][j+1][inner_layer_higher_z]), -c)

    # edges

    # x edges
    for i in range(0, inner_layer_higher_x):

        bins[i][outer_layer_higher_y][outer_layer_higher_z] = apply_shift(
            copy.deepcopy(bins[i][inner_layer_lower][inner_layer_lower]), b + c)
        bins[i][outer_layer_lower][outer_layer_higher_z] = apply_shift(
            copy.deepcopy(bins[i][inner_layer_higher_y][inner_layer_lower]), -b + c)
        bins[i][outer_layer_higher_y][outer_layer_lower] = apply_shift(
            copy.deepcopy(bins[i][inner_layer_lower][inner_layer_higher_z]), b - c)
        bins[i][outer_layer_lower][outer_layer_lower] = apply_shift(
            copy.deepcopy(bins[i][inner_layer_higher_y][inner_layer_higher_z]), -b - c)

    # y edges
    for j in range(0, inner_layer_higher_y):

        bins[outer_layer_higher_x][j][outer_layer_higher_z] = apply_shift(
            copy.deepcopy(bins[inner_layer_lower][j][inner_layer_lower]), a + c)
        bins[outer_layer_lower][j][outer_layer_higher_z] = apply_shift(
            copy.deepcopy(bins[inner_layer_higher_x][j][inner_layer_lower]), -a + c)
        bins[outer_layer_higher_x][j][outer_layer_lower] = apply_shift(
            copy.deepcopy(bins[inner_layer_lower][j][inner_layer_higher_z]), a - c)
        bins[outer_layer_lower][j][outer_layer_lower] = apply_shift(
            copy.deepcopy(bins[inner_layer_higher_x][j][inner_layer_higher_z]), -a - c)

    # z edges
    for k in range(0, inner_layer_higher_z):

        bins[outer_layer_higher_x][outer_layer_higher_y][k] = apply_shift(
            copy.deepcopy(bins[inner_layer_lower][inner_layer_lower][k]), a + b)
        bins[outer_layer_lower][outer_layer_higher_y][k] = apply_shift(
            copy.deepcopy(bins[inner_layer_higher_x][inner_layer_lower][k]), -a + b)
        bins[outer_layer_higher_x][outer_layer_lower][k] = apply_shift(
            copy.deepcopy(bins[inner_layer_lower][inner_layer_higher_y][k]), a - b)
        bins[outer_layer_lower][outer_layer_lower][k] = apply_shift(
            copy.deepcopy(bins[inner_layer_higher_x][inner_layer_higher_y][k]), -a - b)

    # corners

    bins[outer_layer_higher_x][outer_layer_higher_y][outer_layer_higher_z] = apply_shift(
        copy.deepcopy(bins[inner_layer_lower][inner_layer_lower][inner_layer_lower]), a + b + c)
    bins[outer_layer_lower][outer_layer_higher_y][outer_layer_higher_z] = apply_shift(
        copy.deepcopy(bins[inner_layer_higher_x][inner_layer_lower][inner_layer_lower]), -a + b + c)
    bins[outer_layer_higher_x][outer_layer_lower][outer_layer_higher_z] = apply_shift(
        copy.deepcopy(bins[inner_layer_lower][inner_layer_higher_y][inner_layer_lower]), a - b + c)
    bins[outer_layer_higher_x][outer_layer_higher_y][outer_layer_lower] = apply_shift(
        copy.deepcopy(bins[inner_layer_lower][inner_layer_lower][inner_layer_higher_z]), a + b - c)
    bins[outer_layer_lower][outer_layer_lower][outer_layer_higher_z] = apply_shift(
        copy.deepcopy(bins[inner_layer_higher_x][inner_layer_higher_y][inner_layer_lower]), -a - b + c)
    bins[outer_layer_lower][outer_layer_higher_y][outer_layer_lower] = apply_shift(
        copy.deepcopy(bins[inner_layer_higher_x][inner_layer_lower][inner_layer_higher_z]), -a + b - c)
    bins[outer_layer_higher_x][outer_layer_lower][outer_layer_lower] = apply_shift(
        copy.deepcopy(bins[inner_layer_lower][inner_layer_higher_y][inner_layer_higher_z]), a - b - c)
    bins[outer_layer_lower][outer_layer_lower][outer_layer_lower] = apply_shift(copy.deepcopy(
        bins[inner_layer_higher_x][inner_layer_higher_y][inner_layer_higher_z]), -a - b - c)

    bins_shape = (
        original_array_size[0] + 2, original_array_size[0] + 2, original_array_size[0] + 2)

    return bins, bins_shape


def pad_array(bins, x_max, y_max, z_max):
    
    '''
    Pad the bins array with None values.
    '''

    # pad z values

    for i, j in np.ndindex((x_max, y_max)):

        bins[i][j].insert(0, None)
        bins[i][j].insert(z_max + 1, None)

    # pad y values

    pad = [None for stuff in range(0, y_max+2)]

    for i in range(0, x_max):

        bins[i].insert(0, copy.deepcopy(pad))
        bins[i].insert(y_max + 1, copy.deepcopy(pad))

    pad = [copy.deepcopy((pad)) for stuff in range(0, x_max+2)]

    # pad x values

    bins.insert(0, copy.deepcopy(pad))
    bins.insert(x_max + 1, copy.deepcopy(pad))

    return bins


def apply_shift(bin_contents, shift):

    if bin_contents:

        for stuff in bin_contents:

            stuff['x'] += shift[0]
            stuff['y'] += shift[1]
            stuff['z'] += shift[2]

    return bin_contents


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
