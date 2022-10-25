def breadth_first_search(bond_dict, search_center):  

    '''
    Walks breadth first through the structure identifying which shells the other atoms are in relative to the starting atom.
    '''

    neighbour_dict = {}

    number_of_atoms = len(bond_dict)

    neighbour_dict[0] = [search_center]
    found_already = [search_center]
    old_neighbours = [search_center]
  
    shell = 1

    while len(found_already) < number_of_atoms - 1:

        neighbour_dict[shell] = []
        new_neighbours = []

        for n in old_neighbours:
            
            neighbours = bond_dict[n]

            for m in neighbours:

                if m not in found_already:

                    neighbour_dict[shell].append(m)
                    found_already.append(m)
                    new_neighbours.append(m)


        old_neighbours = new_neighbours
        shell += 1

    return neighbour_dict

        

   

    



