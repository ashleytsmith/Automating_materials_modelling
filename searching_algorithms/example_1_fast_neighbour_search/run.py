import numpy as np
import datetime

from searching_algorithms import sort_atoms_into_bins
from searching_algorithms import neighbour_search
from searching_algorithms import input_and_ouput as io


def run(bins_shape, cut_off_distance, atoms):

    '''
    Searches fo the nearest neighbours of each atom in a 3D cell.
    '''

    # get atom properties

    positions, symbols, cell = io.get_atom_properties(atoms)
    number_of_atoms = len(positions)
    indices = np.arange(number_of_atoms)


    start = datetime.datetime.now()

    # binning procedure

    bins, bins_shape = sort_atoms_into_bins.sort_into_bins(
        positions, symbols, indices, cell, bins_shape, repeat=True)


    # perform neighbour search

    bond_dict = neighbour_search.neighbour_search(
        bins, bins_shape, number_of_atoms, cut_off_distance)


    finish = datetime.datetime.now()

    run_time = finish-start

    return bond_dict, run_time


  
