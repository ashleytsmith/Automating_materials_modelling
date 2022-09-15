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

    bin_count = 0
    atom_count = 0
    atom_count_per_bin = []

    for i,j,k in np.ndindex(bins_shape):

        bin_contents = bins[i][j][k]
        
        if bin_contents:

            bin_count +=1
            atoms_in_bin = len(bin_contents.indices)
            atom_count += atoms_in_bin
            atom_count_per_bin.append(atoms_in_bin)

    print('atom count per bin is ' + str(atom_count_per_bin))
    print('bin count is ' + str(bin_count))
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

        print(single_sliced)

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


def test_periodic_boundary_conditions(bins,cell):

    '''
    Show all the atoms with one layer of bins around the whole cell.
    '''

    '''
    Repeat the bins for one layer around the cell. 
    '''
    a,b,c = geometry.get_cell_vectors(cell)
    min_index = 0
    max_index = np.shape(bins)[0] - 1 

    # faces 

    slice_x_upper = bins[min_index][:][:]
    slice_x_lower = bins[max_index][:][:]
    slice_y_upper = bins[:][min_index][:]
    slice_y_lower = bins[:][max_index][:]
    slice_z_upper = bins[:][:][min_index]
    slice_z_lower = bins[:][:][max_index]

    new_slice = np.copy(slice)
    
    #print(np.shares_memory(slice,bins))


    for i,j in np.ndindex(np.shape(new_slice)):

        atoms_info = slice[i][j]

        #print(np.shares_memory(slice,bins))
    
        if atoms_info:

            for atom_info in atoms_info:

                atom_info['x'] += c[0]
                atom_info['y'] += c[1]
                atom_info['z'] += c[2]
            

    for i,j in np.ndindex(np.shape(new_slice)):

        atoms_info = slice[i][j]
    
        if atoms_info:

            for atom_info in atoms_info:

                print(atom_info)

        

    



def show_neighbours_of_each_bin(bins,cell):

    '''
    Shows a movie of all 26 neighbours for every occupied bin.
    '''

    original_shape = np.shape(bins)
    bins = np.pad(bins, ((1,1), (1,1), (1, 1)), mode='wrap')
    print(np.shape(bins))

    bounding_indices = np.ndindex(np.shape(bins))
    
   # traj = Trajectory('test.traj','w') 

    for i,j,k in np.ndindex(np.shape(bins)):

        positions = []
        symbols = []

        #neighbouring_bins = bins[i-1:i+1][j-1:j+1][k-1:k+1]
        #print(i,j,k)


        if 0 < j < 5 and k < 5: 
        
           neighbouring_bins = bins[i][j-1:j+1][k+1]
           print(np.shape(neighbouring_bins))

        #print(bins[i-1:i+1][j-1:j+1][k-1:k+1])

        #for atom_info in bins[i][j][k]:

                #pos = [atom_info[0],atom_info[1],atom_info[2]]
                #positions.append(pos)
                #symbols.append(atom_info[4])

            
                #print(' sub bins')
                #print(nearest_neighbours)



   # atoms = Atoms(symbols, positions= positions, cell=cell)
   # traj.write(atoms=atoms,mode='a')

   # traj.close()

        

    
  

    