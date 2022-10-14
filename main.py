#! /usr/bin/env python3

from searching_algorithms.example_1_fast_neighbour_search.run import run as run_fast_neighbour_search
from searching_algorithms.example_1_fast_neighbour_search.test import run_tests as run_fast_neighbour_search_tests
from searching_algorithms.example_1_fast_neighbour_search.benchmarking import run_benchmarking as run_fast_neighbour_search_benchmarking
#from searching_algorithms.example_1_fast_neighbour_search.benchmarking import run_ase 
from searching_algorithms.example_1_fast_neighbour_search.plot import plot as plot_benchmarking_results

from searching_algorithms.example_1_fast_neighbour_search import parameters
from searching_algorithms import input_and_ouput as io
from searching_algorithms.example_1_fast_neighbour_search.test import test_scaling_behavour


# reccommended input params

bins_shape = parameters.bins_shape
cut_off_distance = parameters.cut_off_distance
input_structure = parameters.get_file_path()

# read infile

atoms = io.read_trajfile(input_structure)

# run

#bond_dict, run_time = run_fast_neighbour_search(bins_shape, cut_off_distance, atoms)

# testing 

#run_fast_neighbour_search_tests(bins_shape, cut_off_distance, atoms)

# benchmarking

#run_fast_neighbour_search_benchmarking(run_fast_neighbour_search,reps = 8,runs = 5)

# plot 

plot_benchmarking_results()




















