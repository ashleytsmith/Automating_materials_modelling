import copy
import numpy as np
import os

from searching_algorithms import build
from searching_algorithms import input_and_ouput as io
from ase.io.trajectory import Trajectory
from ase import Atoms

from ase.io import write
from ase.utils import hsv


def write_rendering_input_files(neighbour_info, atoms):

    connectivity_search_movie(neighbour_info, atoms)



def connectivity_search_movie(neighbour_info, atoms):

    '''
    Movie of the shells being filled for a 2 by 2 super cell from various different views.
    '''

    views = [0,1,2]

    shells = []
    number_of_atoms = len(atoms)

    # repeat the atoms and the neighbour info dictionary 

    atoms = build.repeat_structure(atoms, 2)
    reps = range(1,8)

    original_neighbour_info = copy.deepcopy(neighbour_info)

    for i in range(0,len(neighbour_info)):

        for rep in reps:

            rep_inds = [x+number_of_atoms*rep for x in original_neighbour_info[i]]
            neighbour_info[i].extend(rep_inds)


    # get atom shells in a list based on the neighbour info

    shells = []

    for i in range(0,len(neighbour_info)):

        current_shell = [atoms[j] for j in neighbour_info[i]]

        shells.append(current_shell)


    # progressivly fill the shells

    filling_shells = copy.deepcopy([shells[0]])
    filled_shells = copy.deepcopy(shells[0])

    for i in range(1,len(shells)):

        filled_shells.extend(shells[i])
        filling_shells.append(copy.deepcopy(filled_shells))

    # pack nested lists of atom objects into list of atoms objects

    cell = atoms.get_cell()

    packed_shells = []

    for i,shell in enumerate(filling_shells):

        current_shell = build.pack_atom_objects(shell, cell)
        packed_shells.append(current_shell)


    # make movie

    main_folder = 'Rendering/'

    for view in views:

        sub_folder = 'view_' + str(view) + '/'
            
        if not os.path.isdir(main_folder + sub_folder):
                
            os.mkdir(main_folder + sub_folder)

        for i,shell in enumerate(packed_shells):

            file_path = sub_folder + 'frame_' + str(i) + '_view_' + str(view) 
            generate_render_files(shell, file_path, view)


    



def check_output(shells):

    for i,stuff in enumerate(shells):
    
        print('shell ' + str(i))
        print(stuff.symbols)
        print('   ')




def generate_render_files(atoms, file_path, view, preview = False):

    '''
    Prepare inputs files for rendering using POVRAY.
    '''

    # prepare generic keywords

    number_of_atoms = len(atoms)
    symbols = atoms.get_chemical_symbols()

    views = ['0x, 0y, 0z','-90x, 0y, 0z','0x, 90y, 90z']
    zooms = [100, 110, 160]

    rotation = views[view] 
    zoom = zooms[view]
    
    red = [229,23,23]
    blue = [11,50,255]
    cream = [239,199,159]

    colors = [(red if s == 'O' else (cream if s =='Si' else (blue if s == 'Al' else None))) for s in symbols]
    colors = [ [r/ 255, g/255, b/255] for r,g,b in colors]

    # prepare povray kwargs

    width = 800 # option to set width in pixels
    texture = ['glass',] * number_of_atoms
    cell_line_thickness = 0.25
    background_color = 'White'
    transparent_background = False 
    camera_dist = zoom # distance from camera to front atom
    camera_type = 'perspective' # perspective, ultra_wide_angle
    point_lights = [[(-1.,-2.,-3.),'Red'],[(-1,-4,-5), 'Blue']] # [[loc1, color1], [loc2, color2],...]
    area_light = [(2., 3., 40.) ,'White', .7, .7, 3, 3]  # loc, color width, height, Nlamps_x, Nlamps_y
    
    # set keywords

    generic_projection_settings = { 
    'rotation': rotation,
    'colors': colors,
    }

    povray_settings = { 

    'canvas_width' : width,  
    'camera_dist'  : camera_dist,                     
    'camera_type'  : camera_type, 
    'point_lights' : point_lights,
    'area_light'   : area_light,
    'transparent'  : transparent_background, 
    'background'   : background_color,       
    'textures'     : texture, 
    'celllinewidth': cell_line_thickness, 
    }

    if preview:

        # make prerender

        write('Rendering/' + file_path + '.png', atoms, **generic_projection_settings)

    else:

         # prepare povray input

       write('Rendering/' + file_path + '.pov', atoms,
             **generic_projection_settings,
                povray_settings=povray_settings)



