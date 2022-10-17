import datetime
import pickle
import os

from ase.neighborlist import *

from searching_algorithms import input_and_ouput as io
from searching_algorithms import build

from searching_algorithms.neighbour_search_example import parameters
from searching_algorithms.neighbour_search_example.run import run as run_neighbour_search

from scipy.spatial import cKDTree



def run_benchmarking(method,reps,runs):

    '''
    Run a given method for larger and larger super cells for a certain number of repetitions to test the scaling.
    '''

    # assign filename

    if method == run_neighbour_search:

        file_name = 'results.pkl'

    if method == run_ase_neighbour_search:

        file_name = 'results_ase.pkl'

    if method == run_KDTree:

        file_name = 'results_kdtree.pkl'

    # get previously computed data if available
    
    if os.path.isfile(file_name):

        system_sizes, bond_dicts, run_times, runs_list, current_rep, current_run = get_results(file_name,runs)

    else:

        system_sizes = []
        run_times = []
        runs_list = []
        bond_dicts = []

        current_rep = 0
        current_run = 0

    # base case parameters

    base_bins_shape = parameters.bins_shape
    cut_off_distance = parameters.cut_off_distance
    input_structure = parameters.get_file_path()
    base_atoms = io.read_trajfile(input_structure)
        
    # generate inputs for larger and larger supercells

    all_reps = list(range(1, reps + 1))
    repetitions = len(all_reps) 
    structures = build.repeat_structure(base_atoms,all_reps)
    bins_shapes = generate_bin_shapes(base_bins_shape,repetitions)

    # run for different system sizes

    for rep in range(current_rep,reps):

        run_times.append(runs_list)

        for run in range(current_run,runs):

            if method == run_neighbour_search:

                bond_dict, run_time = method(bins_shapes[rep], cut_off_distance, structures[rep])

            if method == run_ase_neighbour_search:

                bond_dict, run_time = method(structures[rep])

            if method == run_KDTree:

                bond_dict, run_time = method(structures[rep])

            run_times[rep].append(run_time)
            save_results(system_sizes,run_times,bond_dicts, file_name) 
            print('super cell size ' + str(rep) + ' run number ' + str(run) + ' at ' + str(datetime.datetime.now()) )


        bond_dicts.append(bond_dict)
        system_sizes.append(len(bond_dict))
        runs_list = []
 
    

def get_results(file_name,runs):

    infile = open(file_name,'rb')
    results_dict = pickle.load(infile)
    infile.close()

    system_sizes = results_dict['system_sizes']
    run_times = results_dict['run_times']
    runs_list = results_dict['run_times'][-1]
    bond_dicts = results_dict['bond_dicts']

    current_rep = len(run_times) - 1
    current_run = len(runs_list) 

    if current_run == runs:

        current_rep += 1 
        current_run = 0
        runs_list = []

    return system_sizes, bond_dicts, run_times, runs_list, current_rep, current_run



def save_results(system_sizes,run_times,bond_dicts,filename):

    outfile = open(filename,'wb')

    collect_all = {}

    collect_all['system_sizes'] = system_sizes
    collect_all['run_times'] = run_times
    collect_all['bond_dicts'] = bond_dicts

    pickle.dump(collect_all,outfile)

    outfile.close()



def generate_bin_shapes(base_bins_shape,reps):

    bins_shapes = []

    current_bins_shape = base_bins_shape

    bins_shapes.append(current_bins_shape)

    reps = reps - 1 

    for i in range(0,reps):
    
        current_bins_shape = tuple([2*x for x in current_bins_shape])
        bins_shapes.append(current_bins_shape)

    return bins_shapes



def run_ase_neighbour_search(atoms):

    '''
    Run Neighbourlist from ase with the latest linear scaling version.
    '''

    start = datetime.datetime.now()

    number_of_atoms = len(atoms)
    cutoffs = natural_cutoffs(atoms,mult = 1.25)

    nl = build_neighbor_list(atoms,cutoffs,self_interaction= False,skin = 0,bothways = True,primitive=NewPrimitiveNeighborList) 

    bond_dict = {}

    for i in range(0,number_of_atoms):
        
        bond_dict[i] = []

    for i in range(0,number_of_atoms):
        
        indices, offsets = nl.get_neighbors(i)
        
        for j in indices:

            bond_dict[i].append(j)

    finish = datetime.datetime.now()

    run_time = finish-start

    return bond_dict, run_time


def run_KDTree(atoms):

    number_of_atoms = len(atoms)
    positions, symbols, cell = io.get_atom_properties(atoms)

    start = datetime.datetime.now()

    kdtree=cKDTree(positions)

    bond_dict = {}

    for i in range(0,number_of_atoms):
        
        bond_dict[i] = []

    for i,pos in enumerate(positions):

        dist, neighbours = kdtree.query(pos,parameters.expected_neighbours(symbols[i]) + 1)
       
        for j in neighbours:

            bond_dict[i].append(j)

    finish = datetime.datetime.now()

    run_time = finish-start

    return bond_dict, run_time














    







