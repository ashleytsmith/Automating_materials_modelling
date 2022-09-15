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



number_of_atoms_in_structure = 108
oxygen_of_the_first_acid_site = 107  
Al_of_the_first_acid_site = 7
hydrogen_acid_site = 108  


cut_off_distance = 2
silicon_expected_neighbours = 4
oxygen_expected_neighbours = 2
hydrogen_expected_neighbours = 1



input_structure = os.getcwd() + "../Input_structures/CHA.traj"
output_folder = "Generated_structures"