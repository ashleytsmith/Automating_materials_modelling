import os
from fnmatch import fnmatch


def run(view, largest_frame, height, ratio):

    '''
    Render a movie.
    '''

    paths = get_filepaths( '*view_' + str(view) + '.pov')
    base_pov_file = 'Rendering/' + 'view_' + str(view) + '/frame_' + str(largest_frame) + '_view_' + str(view) + '.pov'

    head, all_atom_positions, num_lines_in_largest_file = read_base_pov_file(base_pov_file)

    home = os.getcwd()
    povray_path = '/Users/ashley/povray/include'     
    os.chdir(povray_path)   # hacky but most pragmatic way to get POVRAY to run on mac without official support.

    for path in paths:

        ini_file = path.split('.')[0] + '.ini'
        edit_ini_file(ini_file, height, ratio)
        edit_pov_file(path, head, all_atom_positions, num_lines_in_largest_file)
        os.system('povray ' + ini_file + ' ' + path)
       
    os.chdir(home)

    # repeat frames

    new_num_frames = 2*(largest_frame + 1)

    for path in paths:

        frame = path.split('frame_')[-1].split('_')[0]
        repeat_and_reverse = str(new_num_frames - 1 - int(frame))
        forwards_frame_path = path.replace('.pov','.png')
        backwards_frame_path = path.split('frame_')[0] + 'frame_' + repeat_and_reverse + '_view_' + str(view) + '.png'

        os.system('cp -r ' + forwards_frame_path + ' ' + backwards_frame_path)


def get_filepaths(pattern):

    target_directory = os.getcwd() + '/Rendering'

    paths = []
    
    for path, subdirs, files in os.walk(target_directory):

        for name in files:
    
            if fnmatch(name, pattern):

                paths.append(os.path.join(path, name))

    return paths



def edit_ini_file(path, height, ratio):

    '''
    Edit the ini file produced by ase.io.pov. 
    By default ase only allows either width or height to be set leading to different image sizes.
    '''

    locs = [5,7]

    with open(path, 'r') as file:
   
        data = file.readlines()

    data[locs[0]] = '; Width / Height ='+ ratio + '\n'
    data[locs[1]] = 'Height=' + height + '\n'

    with open(path, 'w') as file:

        file.writelines(data)



def read_base_pov_file(path):

    '''
    Read the pov file produced by ase.io.pov and store the lines.
    '''

    tail_length = 2
    cut_position = 49

    with open(path, 'r') as file:
   
        data = file.readlines()

    number_of_atoms = len(data) - cut_position - tail_length

    head = data [:cut_position]
    atom_positions = data[cut_position:cut_position + number_of_atoms]

    num_lines = len(head) + len(atom_positions) + tail_length

    return head, atom_positions, num_lines

    

def edit_pov_file(path, head, all_atom_positions, num_lines_largest_file):

    '''
    Edit the pov file produced by ase.io.pov. 
    By default ase sets the view based upon the atoms supplied causing the cell to move around in the movie.
    '''

    locs = [7]
    z_location = '100.00'

    data = list(head) 

    num_lines = sum(1 for line in open(path))
    length_difference = num_lines_largest_file - num_lines

    data[locs[0]] = '  direction ' + z_location + '*z' + '\n'

    if not length_difference == 0:

        atom_positions = all_atom_positions[:-length_difference]
        data.extend(atom_positions)

    else:

        data.extend(all_atom_positions)

    with open(path, 'w') as file:

        file.writelines(data)




        

    



    











