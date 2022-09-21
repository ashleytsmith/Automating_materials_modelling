from multiprocessing import shared_memory
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

    for i,j,k in np.ndindex(bins_shape):

        bin_contents = bins[i][j][k]
        total_bins += 1
        
        if bin_contents:

            bin_count +=1
            atoms_in_bin = len(bin_contents.indices)
            atom_count += atoms_in_bin
            atom_count_per_bin.append(atoms_in_bin)

    print('atom count per ocuupied bin is ' + str(atom_count_per_bin))
    print('total number of bins is ' + str(total_bins))
    print('number of occipied bins is ' + str(bin_count))
    print('atom count is ' + str(atom_count))


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


def show_atoms_in_each_bin(bins,cell,showall = True):

    '''
    Shows a movie of the atoms in all the different bins. 
    '''

    all_positions = []
    all_symbols = []
    
    traj = Trajectory('test.traj','w') 

    for i,j,k in np.ndindex(np.shape(bins)):

        positions = []
        symbols = []

        atoms_info = bins[i][j][k]

        if atoms_info:

            pos_list,sym_list = get_bin_info(atoms_info)
            positions.extend(pos_list)
            symbols.extend(sym_list)
            all_positions.extend(pos_list)
            all_symbols.extend(sym_list)

            if not showall:

                atoms = Atoms(symbols, positions= positions, cell=cell)
                traj.write(atoms=atoms,mode='a')

        if showall:

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

    for i,j,k in np.ndindex(original_shape):

        i = i + 1
        j = j + 1
        k = k + 1

        all_neighbours_positions = []
        all_neighbours_symbols = []

        bin = bins[i][j][k]
        bin_neighbours = []
        bin_neighbours.append(bin)

        # slicing routine along z axis
        bin_neighbours.extend(bins[i+1][j-1][k-1:k+2])
        bin_neighbours.extend(bins[i+1][j][k-1:k+2])
        bin_neighbours.extend(bins[i+1][j+1][k-1:k+2])
        bin_neighbours.extend(bins[i][j+1][k-1:k+2])
        bin_neighbours.extend([bins[i+1][j][k]])

        # works fine for all bins two layers inside but need to do edge cases for 3 faces one layer in
       
        for atoms_info in bin_neighbours:

            positions = []
            symbols = []

            if atoms_info:

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

    atoms = Atoms(all_symbols, positions= all_positions, cell=cell)
    traj.write(atoms=atoms,mode='a')

    traj.close()



        

    
  

    