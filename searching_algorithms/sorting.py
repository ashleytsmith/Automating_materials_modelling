import copy
import numpy as np

from searching_algorithms import geometry


def sort_into_bins(positions,symbols,indices,cell,bins_shape):

    '''
    Splits atoms into 3D bins by binning along each of the specified basis vectors.
    '''

    a,b,c = geometry.get_cell_vectors(cell)

    slicing_order = (a,b,c)

    atoms_info = create_atom_info_array(positions,symbols,indices)

    single_sliced_positions = slice_along_axis(slicing_order[0],atoms_info,bins_shape[0])
    
    bins = []
    
    for stuff in single_sliced_positions:

        sub_bin = slice_along_axis(slicing_order[1],stuff,bins_shape[1])

        triple_sliced_positions = []

        for things in sub_bin:

            sub_sub_bin = slice_along_axis(slicing_order[2],things,bins_shape[2])
            sub_sub_bin = remove_empty_lists(sub_sub_bin)

            triple_sliced_positions.append(sub_sub_bin)

        bins.append(triple_sliced_positions)

    new_bins = copy.deepcopy(bins)

    add_neighbouring_bins(new_bins,4,a,b,c)

    bins = vectorise_bin_data(bins,bins_shape)
    
    return bins

  
def create_atom_info_array(positions,symbols,indices):

    '''
    Creates a structured array (very similar to a struct in C) with all required variables for each atom accesible from the same element.
    '''

    data_types = np.dtype([('x','f'),('y','f'),('z','f'),('i','i'),('s','U10')])
    array_length = len(positions)
    atoms_info = np.empty(array_length,dtype=data_types)

    x,y,z = zip(*[(pos[0],pos[1],pos[2]) for pos in positions])
    i = [i for i in indices]
    s = [s for s in symbols]

    atoms_info['x']=x
    atoms_info['y']=y
    atoms_info['z']=z
    atoms_info['i']=i
    atoms_info['s']=s

    return atoms_info



def slice_along_axis(axis,atoms_info,number_of_bins):

    '''
    Seperates the input positions into different bins.
    '''

    bins = [ [] for stuff in range(number_of_bins) ]

    for stuff in atoms_info:

        projection = geometry.project_along_axis(stuff['x'],stuff['y'],stuff['z'],axis)

        bin_index = assign_to_bin(projection,number_of_bins)

        bins[bin_index].append(stuff)

    return bins



def assign_to_bin(projection,number_of_bins):

    '''
    Splits the values along an axis into bins based on the size of the projection.
    '''

    if projection < 0:

        bin_index = 0

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


def add_neighbouring_bins(bins,bins_per_dimension,a,b,c):

    '''
    Ensures periodic boundary conditions by adding one layer of bins around the outside. 
    Currently implemented assuming the same number of bins per dimension. 
    '''

    original_array_size = bins_per_dimension

    bins = pad_array(bins,original_array_size)
    
    #bins = np.asarray(bins, dtype = 'object')  # could use if have to but likely a bit slower for elementwise accessing

    for i,j in np.ndindex((original_array_size,original_array_size)):

        # faces

        # x reflections
        bins[original_array_size + 1][i+1][j+1] = apply_shift(copy.deepcopy(bins[1][i+1][j+1]),a)
        bins[0][i+1][j+1] = apply_shift(copy.deepcopy(bins[original_array_size][i+1][j+1]),-a)
        # y reflections
        bins[i+1][original_array_size + 1][j+1] = apply_shift(copy.deepcopy(bins[i+1][1][j+1]),b)
        bins[i+1][0][j+1] = apply_shift(copy.deepcopy(bins[i+1][original_array_size][j+1]),-b)
        # z reflections
        bins[i+1][j+1][original_array_size + 1] = apply_shift(copy.deepcopy(bins[i+1][j+1][1]),c)
        bins[i+1][j+1][0] = apply_shift(copy.deepcopy(bins[i+1][j+1][original_array_size]),-c)

        # edges

        bins[0][0][i] = apply_shift(...)
        bins[0][original_array_size + 1][i] = apply_shift(...) *2 '++' 

    #corners

    bins[0][0][0] = apply_shift(...)*3

    
        

    # printing

    print(bins)
    #print(bins[1])
    #print('   ')
    #print(bins[0])
    #print(bins[5])

 

def pad_array(bins,original_array_size):

    '''
    Pad the bins array with None values.
    '''
   
    # pad z values

    for i,j in np.ndindex((original_array_size,original_array_size)):

        bins[i][j].insert(0, None)
        bins[i][j].insert(original_array_size + 1, None)

    pad = [None for stuff in range(0,original_array_size+2)]

    # pad y values

    for i in range(0,original_array_size):

        bins[i].insert(0,copy.deepcopy(pad))
        bins[i].insert(original_array_size +1, copy.deepcopy(pad))
    
    pad = [copy.deepcopy((pad)) for stuff in range(0,original_array_size+2)]

    # pad x values

    bins.insert(0,copy.deepcopy(pad))
    bins.insert(original_array_size + 1, copy.deepcopy(pad))

    return bins




def apply_shift(bin_contents,shift):

    if bin_contents:

        for stuff in bin_contents:

            stuff['x'] += shift[0]
            stuff['y'] += shift[1]
            stuff['z'] += shift[2]

    return bin_contents





def vectorise_bin_data(bins,shape):

    '''
    Vectorises the atoms_info data in each 3D bin.
    Returns numpy arrays of positions, indices amd symbols.
    '''

    for i,j,k in np.ndindex(shape):
 
        bin_contents = bins[i][j][k]
        
        if bin_contents:

            positions = []
            indices = []
            symbols = []

            for stuff in bin_contents:

                position = (stuff['x'],stuff['y'],stuff['z'])
                index = stuff['i']
                symbol = stuff['s']

                positions.append(position)
                indices.append(index)
                symbols.append(symbol)

            positions = np.asarray(positions)
            indices = np.asarray(indices)
            symbols = np.asarray(symbols)

            atoms_info = AtomsInfo(positions,indices,symbols)

            bins[i][j][k] = atoms_info

    return bins


class AtomsInfo:

  __slots__ = ["positions","indices","symbols"]

  def __init__(self, positions, indices, symbols):

    self.positions = positions
    self.indices = indices
    self.symbols = symbols
        

