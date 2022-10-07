import collections
import numpy as np

from ase import Atoms
from ase.io.trajectory import Trajectory

from searching_algorithms import input_and_ouput as io
from searching_algorithms import geometry

def view_the_grid(input_structure):

    '''
    Displays the grid by putting a hydrogen atom at each of the grid positions.
    '''

    points_per_dimension = 4
    atoms,pos,sym, cell = io.read_trajfile(input_structure)
    bin_positions,_ = geometry.generate_grid(cell,points_per_dimension)
    
    n_dim = 3
    number_of_points = pow(points_per_dimension,n_dim) 

    positions = np.reshape(bin_positions,(number_of_points,n_dim))
    positions = np.ndarray.tolist(positions)

    symbols = np.repeat("H", len(positions))

    atoms = Atoms(symbols, positions= positions, cell=cell)

    io.save_file(atoms,'test.traj')


def count_occupied_bins(bins,bins_shape):

    '''
    Test number of occupied bins and if the sum over the bins yields the correct atom count.
    '''

    total_bins = 0
    bin_count = 0
    atom_count = 0
    atom_count_per_bin = []
    bin_ids = []

    for i,j,k in np.ndindex(bins_shape):

        bin_contents = bins[i][j][k]
        total_bins += 1
        
        if bin_contents:

            bin_count +=1
            atoms_in_bin = len(bin_contents.indices)
            atom_count += atoms_in_bin
            atom_count_per_bin.append(atoms_in_bin)
            bin_ids.append(id(bin_contents))

   
    
    unique = [item for item, count in collections.Counter(bin_ids).items() if count == 1]
    duplicates = [item for item, count in collections.Counter(bin_ids).items() if count > 1]
   
    print('atom count per ocuupied bin is ' + str(atom_count_per_bin))
    print('total number of bins is ' + str(total_bins))
    print('number of occipied bins is ' + str(bin_count))
    print('atom count is ' + str(atom_count))
    print(len(bin_ids), 'bin ids', bin_ids )
    print(len(unique),  ' unique id(s)', unique)
    print(len(duplicates) , 'duplicate id(s)', duplicates)


def test_bin_assignment_along_first_axis(bins,bin_positions,cell):

    '''
    Test bin assignment along the first axis of the bins object. 
    Returns a movie which shows which bin each atom was assigned to.
    '''
    
    all_positions = []
    all_symbols = []
    
    traj = Trajectory('test.traj','w') 

    for stuff in bins:

        single_sliced = stuff[:][:]

        positions = []
        symbols = []

        # add atoms positions

        for j,k in np.ndindex(np.shape(single_sliced)):

            atoms_info = single_sliced[j][k]

            if atoms_info:

                pos_list,sym_list = get_bin_info(atoms_info)
                positions.extend(pos_list)
                symbols.extend(sym_list)
                all_positions.extend(pos_list)
                all_symbols.extend(sym_list)

                
        #add bin positions

        for l,m,n in np.ndindex(np.shape(bins)):

            positions.append(bin_positions[l][m][n])
            symbols.append("H")

            all_positions.append(bin_positions[l][m][n])
            all_symbols.append("H")

        atoms = Atoms(symbols, positions= positions, cell=cell)
        traj.write(atoms=atoms,mode='a')

    atoms = Atoms(all_symbols, positions= all_positions, cell=cell)
    traj.write(atoms=atoms,mode='a')

    traj.close()


def get_bin_info(atoms_info):

    positions = []
    symbols = []

    pos_list = atoms_info.positions

    for i,symbol in enumerate(atoms_info.symbols):

        symbols.append(symbol)
        positions.append(pos_list[i])
       
    return positions, symbols


def show_atoms_in_each_bin(bins,cell,show_empty_bins = True):

    '''
    Shows a movie of the atoms in all the different bins. 
    '''

    all_positions = []
    all_symbols = []
    
    traj = Trajectory('test.traj','w') 

    shape = np.shape(bins)

    for i,j,k in np.ndindex(shape):

        positions = []
        symbols = []

        atoms_info = bins[i][j][k]

        if atoms_info:

            pos_list,sym_list = get_bin_info(atoms_info)
            positions.extend(pos_list)
            symbols.extend(sym_list)
            all_positions.extend(pos_list)
            all_symbols.extend(sym_list)

            if not show_empty_bins:

                atoms = Atoms(symbols, positions= positions, cell=cell)
                traj.write(atoms=atoms,mode='a')

        if show_empty_bins:

            atoms = Atoms(symbols, positions= positions, cell=cell)
            traj.write(atoms=atoms,mode='a')

    atoms = Atoms(all_symbols, positions= all_positions, cell=cell)
    traj.write(atoms=atoms,mode='a')

    traj.close()



def show_outer_layer_of_bins(bins,cell):

    '''
    Shows a movie of the atoms in all the different bins. 
    '''

    all_positions = []
    all_symbols = []
    
    traj = Trajectory('test.traj','w') 

    shape = np.shape(bins)

    for i,j,k in np.ndindex(shape):

        edge_member = i == 0 or i == shape[0] - 1 or j == 0 or j == shape[1] - 1 or k == 0 or k == shape[2] - 1

    
        if edge_member:

            positions = []
            symbols = []

            atoms_info = bins[i][j][k]

            if atoms_info:

                pos_list,sym_list = get_bin_info(atoms_info)
                positions.extend(pos_list)
                symbols.extend(sym_list)
                all_positions.extend(pos_list)
                all_symbols.extend(sym_list)
    
            atoms = Atoms(symbols, positions= positions, cell=cell)
            traj.write(atoms=atoms,mode='a')

    atoms = Atoms(all_symbols, positions= all_positions, cell=cell)
    traj.write(atoms=atoms,mode='a')

    traj.close()




def show_neighbours_of_each_bin(bins,bins_shape,cell):

    '''
    Shows a movie of the unique neighbour pairs for every occupied bin.
    '''

    traj = Trajectory('test.traj','w') 

    all_positions = []
    all_symbols = []

    original_shape = (bins_shape[0] -2,bins_shape[1] - 2,bins_shape[2] -2)

    # loop over the unique pairs excluding edge cases for 3 of the lower faces

    for i,j,k in np.ndindex(original_shape):

        i = i + 1
        j = j + 1
        k = k + 1

        all_neighbours_positions = []
        all_neighbours_symbols = []

        bin = bins[i][j][k]
        bin_neighbours = []
        bin_neighbours.append(bin)

        # slicing routine using slices along the z axis
        bin_neighbours.extend(bins[i+1][j-1][k-1:k+2])
        bin_neighbours.extend(bins[i+1][j][k-1:k+2])
        bin_neighbours.extend(bins[i+1][j+1][k-1:k+2])
        bin_neighbours.extend(bins[i][j+1][k-1:k+2])
        bin_neighbours.extend([bins[i][j][k+1]])

        indices = []

        for atoms_info in bin_neighbours:

            positions = []
            symbols = []
            
            if atoms_info:

                indices.append(atoms_info.indices)

                pos_list,sym_list = get_bin_info(atoms_info)
                positions.extend(pos_list)
                symbols.extend(sym_list)
                all_neighbours_positions.extend(pos_list)
                all_neighbours_symbols.extend(sym_list)
                all_positions.extend(pos_list)
                all_symbols.extend(sym_list)

                atoms = Atoms(symbols, positions= positions, cell=cell)
                traj.write(atoms=atoms,mode='a')

        atoms = Atoms(all_neighbours_symbols, positions= all_neighbours_positions, cell=cell)
        traj.write(atoms=atoms,mode='a')

        print('   ')
        print('current bin ', i,j,k)
        print('neighbouring bins ')
        print(i + 1, j - 1, k - 1)
        print(i + 1, j - 1, k)
        print(i + 1, j - 1, k + 1)
        print('  ')
        print(i + 1, j, k - 1)
        print(i + 1, j, k)
        print(i + 1, j, k + 1)
        print('  ')
        print(i + 1, j + 1, k - 1)
        print(i + 1, j + 1, k)
        print(i + 1, j + 1, k + 1)
        print('  ')
        print(i, j + 1, k - 1)
        print(i, j + 1, k)
        print(i, j + 1, k + 1)
        print('  ')
        print(i, j, k + 1)
        
        if bin:

            print('bin indices', bin.indices)
            print('neighbour indices',indices)

        else:
            
            print('bin empty')
        
        
    

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

    indices = []

    for atoms_info in bin_neighbours:

        if atoms_info:

            indices.append(atoms_info.indices)

            pos_list,sym_list = get_bin_info(atoms_info)
            all_positions.extend(pos_list)
            all_symbols.extend(sym_list)

    atoms = Atoms(all_symbols, positions= all_positions, cell=cell)
    traj.write(atoms=atoms,mode='a')

    print('        ')
    print('neighbour indices for edge cases',indices)

    traj.close()



def check_for_repeat_bonds(bond_dict,indices):

    fail_count = 0
    pass_count = 0

    for i in indices:

        neighbours = bond_dict[i]
        unique_neighbours = [item for item, count in collections.Counter(neighbours).items()]
    
        number_of_neighbours = len(neighbours)
        number_of_unique_neighbours = len(unique_neighbours)
        
        if number_of_neighbours > number_of_unique_neighbours:

            print( 'fail', i, number_of_neighbours, ' neighbours ', number_of_unique_neighbours, ' unique neighbours ' )
            fail_count +=1  

        else: 

            print('pass',  i, number_of_neighbours, ' neighbours ', number_of_unique_neighbours, ' unique neighbours ' )
            pass_count +=1

    print(fail_count, ' total fails')
    print(pass_count, ' total passes')





    




  
         

   
        

    
  

    