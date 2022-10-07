import collections

from searching_algorithms import test
from searching_algorithms.example_1_fast_neighbour_search import parameters

def view_the_grid(input_structure):

    test.view_the_grid(input_structure)

def count_occupied_bins(bins,bins_shape):

    test.count_occupied_bins(bins,bins_shape)


def test_bin_assignment_along_first_axis(bins, bin_positions,cell): 

    test.test_bin_assignment_along_first_axis(bins, bin_positions,cell)


def show_atoms_in_each_bin(bins,cell,show_empty_bins):

    test.show_atoms_in_each_bin(bins,cell,show_empty_bins)


def show_outer_layer_of_bins(bins,cell):

    test.show_outer_layer_of_bins(bins,cell)


def show_neighbours_of_each_bin(bins,bins_shape,cell):

    test.show_neighbours_of_each_bin(bins,bins_shape,cell)

def check_ids_of_bin_neighbours(bins,bins_shape):

    test.check_ids_of_bin_neighbours(bins,bins_shape)


def check_for_repeat_bonds(bond_dict,indices):

    test.check_for_repeat_bonds(bond_dict,indices)


def check_number_of_bonds_found(bond_dict,symbols,indices):

    fail_count = 0

    for i in indices:

        unique_neighbours = [item for item, count in collections.Counter(bond_dict[i]).items()]
        symbol = symbols[i]

        expected_neighbours = parameters.expected_neighbours(symbol)

        if len(unique_neighbours) == expected_neighbours:

            continue

        else:

            print('fail',i, unique_neighbours)
            fail_count +=1

    print(fail_count, ' total fails')



