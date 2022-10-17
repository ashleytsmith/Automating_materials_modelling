#! /usr/bin/env python3

import cProfile
import pickle

from searching_algorithms.neighbour_search_example.run import run as run_neighbour_search
from searching_algorithms.neighbour_search_example.test import run_tests as run_neighbour_search_tests
from searching_algorithms.neighbour_search_example.benchmarking import run_benchmarking as run_neighbour_search_benchmarking
from searching_algorithms.neighbour_search_example.benchmarking import run_ase_neighbour_search 
from searching_algorithms.neighbour_search_example.benchmarking import run_KDTree
from searching_algorithms.neighbour_search_example.plot import plot as plot_benchmarking_results
from searching_algorithms.neighbour_search_example import parameters
from searching_algorithms.neighbour_search_example.test import test_scaling_behavour
from searching_algorithms import input_and_ouput as io


from searching_algorithms.connectivity_search import run as run_connectivtiy_search



# reccommended input params

bins_shape = parameters.bins_shape
cut_off_distance = parameters.cut_off_distance
input_structure = parameters.get_file_path()

# read infile

#atoms = io.read_trajfile(input_structure)

# run neighbour search

#bond_dict, run_time = run_neighbour_search(bins_shape, cut_off_distance, atoms)

# save

#outfile = open('bonding_info.pkl','wb')
#pickle.dump(bond_dict,outfile)

# # testing 

#run_neighbour_search_tests(bins_shape, cut_off_distance, atoms)

# profiling

#cProfile.run('run_neighbour_search(bins_shape, cut_off_distance, atoms)')

# benchmarking

#run_neighbour_search_benchmarking(run_neighbour_search,reps = 8,runs = 5)
#run_neighbour_search_benchmarking(run_ase_neighbour_search,reps = 8,runs = 5)
#run_neighbour_search_benchmarking(run_KDTree,reps = 8,runs = 5)





#run connectivity search

center_atom = 7
run_connectivtiy_search(center_atom)






























