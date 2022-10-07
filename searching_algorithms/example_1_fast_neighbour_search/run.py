import datetime

import os

from ase import Atoms
import numpy as np

from searching_algorithms import input_and_ouput as io
from searching_algorithms import sorting
from searching_algorithms import geometry
from searching_algorithms import searching

from searching_algorithms.example_1_fast_neighbour_search import test
from searching_algorithms.example_1_fast_neighbour_search import parameters


'''
Searches fo the nearest neighbours of each atom in a 3D cell.
'''


def run():

    # input params

    bins_per_dimension = parameters.bins_per_dimension
    cut_off_distance = parameters.cut_off_distance
    input_structure = parameters.get_file_path()
    
    # read infile
    
    atoms,positions,symbols, cell = io.read_trajfile(input_structure)

    
    # generate some bin positions. Just for visualsing. Not needed for actual computation.
    
    #test.view_the_grid(input_structure)
    bin_positions, edges = geometry.generate_grid(cell,bins_per_dimension)


    # binning procedure
 
    number_of_atoms = len(positions)
    indices = np.arange(number_of_atoms)
    bins_shape_before_sorting = (bins_per_dimension,bins_per_dimension,bins_per_dimension)
    bins, bins_shape = sorting.sort_into_bins(positions,symbols,indices,cell,bins_shape_before_sorting, True)
    


    # test the binning procedure

    #test.count_occupied_bins(bins,bins_shape)

    #bins, bins_shape = sorting.sort_into_bins(positions,symbols,indices,cell,bins_shape_before_sorting, False)
    #test.test_bin_assignment_along_first_axis(bins, bin_positions,cell) # only run without neighbouring bins being added

    #test.show_atoms_in_each_bin(bins,cell,show_empty_bins = False)
    #test.show_outer_layer_of_bins(bins,cell)

    # perform neighbour search

    bond_dict = searching.neighbour_search(bins,bins_shape,number_of_atoms,cut_off_distance)

    
    
    # test neighbour search

    #test.show_neighbours_of_each_bin(bins,bins_shape,cell)
    
    
    
    
    
    
    
    
    #test.check_ids_of_bin_neighbours(bins,bins_shape)

    # test search

    #test.check_number_of_bonds_found(bond_dict,symbols,indices)

    #test.check_for_repeat_bonds(bond_dict,indices)

    #save file

    #atoms.extend(bin_atoms)

    #io.save_file(atoms,'test.traj')

    
    

    start = datetime.datetime.now()

    finish = datetime.datetime.now()

    print(finish-start)


    





    





