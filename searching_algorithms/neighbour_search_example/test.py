import collections
import numpy as np

from searching_algorithms import test_neighbour_search
from searching_algorithms import input_and_ouput as io
from searching_algorithms import sort_atoms_into_bins
from searching_algorithms import geometry
from searching_algorithms import neighbour_search

from searching_algorithms.neighbour_search_example import parameters
from searching_algorithms.neighbour_search_example import benchmarking
from searching_algorithms import build
from searching_algorithms.neighbour_search_example.run import run as run_neighbour_search


def show_the_grid(cell, bin_positions):

    test_neighbour_search.show_the_grid(cell, bin_positions)


def count_occupied_bins(bins, bins_shape):

    test_neighbour_search.count_occupied_bins(bins, bins_shape)


def test_bin_assignment_along_first_axis(bins, bin_positions, cell):

    test_neighbour_search.test_bin_assignment_along_first_axis(bins, bin_positions, cell)


def show_atoms_in_each_bin(bins, cell, bin_positions, show_empty_bins, only_show_edges):

    test_neighbour_search.show_atoms_in_each_bin(
        bins, cell, bin_positions, show_empty_bins, only_show_edges)


def show_neighbours_of_each_bin(bins, bins_shape, cell, bin_positions):

    test_neighbour_search.show_neighbours_of_each_bin(bins, bins_shape, cell, bin_positions)


def check_for_repeat_bonds(bond_dict, indices):

    test_neighbour_search.check_for_repeat_bonds(bond_dict, indices)


def check_number_of_bonds_found(bond_dict, symbols, indices):

    fail_count = 0
    pass_count = 0

    for i in indices:

        unique_neighbours = [item for item,
                             count in collections.Counter(bond_dict[i]).items()]
        symbol = symbols[i]

        expected_neighbours = parameters.expected_neighbours(symbol)

        if len(unique_neighbours) == expected_neighbours:

            print('pass', i, 'unique neighbours found:', unique_neighbours)
            pass_count += 1

        else:

            print('fail', i, 'unique neighbours found:', unique_neighbours)
            fail_count += 1

    print(pass_count, ' total passes')
    print(fail_count, ' total fails')



def test_scaling_behavour(reps):

    # base case parameters

    base_bins_shape = parameters.bins_shape
    cut_off_distance = parameters.cut_off_distance
    input_structure = parameters.get_file_path()
    base_atoms = io.read_trajfile(input_structure)
        
    # generate inputs for larger and larger supercells

    all_reps = list(range(1, reps + 1))
    repetitions = len(all_reps) 
    structures = build.repeat_structure(base_atoms,all_reps)
    bins_shapes = benchmarking.genertate_bin_shapes(base_bins_shape,repetitions)

  
    for rep in range(0,reps):

        atoms = structures[rep]
        symbols = atoms.get_chemical_symbols()
        indices = np.arange(len(symbols))
        bond_dict, run_time = run_neighbour_search(bins_shapes[rep], cut_off_distance, atoms)
        check_number_of_bonds_found(bond_dict,symbols,indices)


def run_tests(bins_shape, cut_off_distance, atoms):


   # get atom properties

    positions, symbols, cell = io.get_atom_properties(atoms)
    number_of_atoms = len(positions)
    indices = np.arange(number_of_atoms)
    

    # generate some bin positions. Just for visualsing. Not needed for actual computation.

    bin_positions = geometry.generate_grid(cell, bins_shape)
    test_neighbour_search.show_the_grid(cell, bin_positions)


    # binning procedure

    bins, bins_shape = sort_atoms_into_bins.sort_into_bins(
    positions, symbols, indices, cell, bins_shape, repeat=True)

    
    # test the binning procedure

    test_neighbour_search.count_occupied_bins(bins,bins_shape)

    test_neighbour_search.test_bin_assignment_along_first_axis(bins, bin_positions,cell)

    test_neighbour_search.show_atoms_in_each_bin(bins,cell, bin_positions, show_empty_bins = False, only_show_edges = False)


    # perform neighbour search

    bond_dict = neighbour_search.neighbour_search(
        bins, bins_shape, number_of_atoms, cut_off_distance)


    # test neighbour search

    test_neighbour_search.show_neighbours_of_each_bin(bins,bins_shape,cell, bin_positions)

    check_number_of_bonds_found(bond_dict,symbols,indices)

    test_neighbour_search.check_for_repeat_bonds(bond_dict,indices)