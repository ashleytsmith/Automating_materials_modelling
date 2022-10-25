#! /usr/bin/env python3

import cProfile
import pickle

from searching_algorithms.neighbour_search_example.run import run as run_neighbour_search
from searching_algorithms.neighbour_search_example.test import run_tests as run_neighbour_search_tests
from searching_algorithms.neighbour_search_example.benchmarking import run_benchmarking as run_neighbour_search_benchmarking
from searching_algorithms.neighbour_search_example.benchmarking import run_ase_neighbour_search 
from searching_algorithms.neighbour_search_example.benchmarking import run_KDTree
from searching_algorithms.neighbour_search_example.plot import plot as plot_benchmarking_results
from searching_algorithms.neighbour_search_example import parameters as neighbour_search_params
from searching_algorithms.neighbour_search_example.test import test_scaling_behavour
from searching_algorithms import input_and_ouput as io


from searching_algorithms.connectivity_search_example.run import run as run_connectivtiy_search
from searching_algorithms.connectivity_search_example import parameters as connectivity_search_params
from searching_algorithms.connectivity_search_example.plot import plot as plot_connectivity
from searching_algorithms.connectivity_search_example.run_rendering import run as run_rendering


# reccommended input params

bins_shape = neighbour_search_params.bins_shape
cut_off_distance = neighbour_search_params.cut_off_distance
input_structure = neighbour_search_params.get_file_path()

# read infile

atoms = io.read_trajfile(input_structure)

# run neighbour search

bond_dict, run_time = run_neighbour_search(bins_shape, cut_off_distance, atoms)

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

center_atom = connectivity_search_params.center_atom
neighbour_info = run_connectivtiy_search(center_atom, skip = False)

# plot connectivity results

plot_connectivity(neighbour_info, atoms)
run_rendering(1,10)
































