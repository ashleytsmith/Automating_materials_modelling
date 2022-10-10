import collections

from searching_algorithms import test
from searching_algorithms.example_1_fast_neighbour_search import parameters


def show_the_grid(cell, bin_positions):

    test.show_the_grid(cell, bin_positions)


def count_occupied_bins(bins, bins_shape):

    test.count_occupied_bins(bins, bins_shape)


def test_bin_assignment_along_first_axis(bins, bin_positions, cell):

    test.test_bin_assignment_along_first_axis(bins, bin_positions, cell)


def show_atoms_in_each_bin(bins, cell, bin_positions, show_empty_bins, only_show_edges):

    test.show_atoms_in_each_bin(
        bins, cell, bin_positions, show_empty_bins, only_show_edges)


def show_neighbours_of_each_bin(bins, bins_shape, cell, bin_positions):

    test.show_neighbours_of_each_bin(bins, bins_shape, cell, bin_positions)


def check_for_repeat_bonds(bond_dict, indices):

    test.check_for_repeat_bonds(bond_dict, indices)


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
