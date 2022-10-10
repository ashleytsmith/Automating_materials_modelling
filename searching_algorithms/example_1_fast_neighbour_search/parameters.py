import os


def expected_neighbours(symbol):
    
    '''
    Returns expected number of neighbours for a given atom type.
    '''

    if symbol == "Si" or symbol == "Al":

        x = silicon_expected_neighbours

    elif symbol == "O":

        x = oxygen_expected_neighbours

    elif symbol == "H":

        x = hydrogen_expected_neighbours

    return x


def get_file_path():

    input_structure = os.getcwd() + "/Input_structures/CHA.traj"

    return input_structure


silicon_expected_neighbours = 4
oxygen_expected_neighbours = 2
hydrogen_expected_neighbours = 1


bins_shape = (3, 3, 3)
cut_off_distance = 2.4
