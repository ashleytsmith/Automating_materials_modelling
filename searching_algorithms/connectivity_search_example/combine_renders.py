import os 
import itertools

import cv2
import numpy as np
import imageio


from searching_algorithms.connectivity_search_example.run_rendering import get_filepaths 



def combine_images(neighbour_info):

    '''
    Comine the individual scenes into one larger PNG file.
    '''

    main_source_folder = 'Rendering/' 
    main_target_folder = 'Rendering/gif/'
    file_path = 'Rendering/gif/connectivity_search_movie.gif'
    number_of_views = 3

    crop_list = [None, (10,0), (180,150)] # (top,bottom)
    shift_list = [None, ('left',40), ('right',10)] # (left,right)
    
    legend = cv2.imread(  'Rendering/gif/legend.png')
    legend_offset = (45,630)
    annotations_offset = (400,670) 
    
    blank_canvas = np.zeros([800,800,3],dtype=np.uint8)
    blank_canvas[:] = 255

    legend = scale_down_image(legend, 0.4)

    paths = []

    # get paths and sort
    
    for view in range(0,number_of_views):

        path_list = get_filepaths( main_source_folder + 'view_' + str(view), '*view_' + str(view) + '.png')
        frames = [int(path.split('frame_')[-1].split('_')[0]) for path in path_list]
        path_list = [x for _, x in sorted(zip(frames, path_list))]
        paths.append(path_list)

    #output folder

    if not os.path.isdir(main_target_folder):
        
        os.mkdir(main_target_folder)

    # get frames

    num_frames = len(paths[0])
    movie = []

    for view in range(0,number_of_views):

        act = []

        for i in range(0,num_frames):

            img = cv2.imread(paths[view][i])
            act.append(img)

        movie.append(act)

    # cropping and chopping

    for i, act in enumerate(movie):

        if crop_list[i]:

            for j, frame in enumerate(act):

                movie[i][j] = crop_image(movie[i][j], crop_list[i][0], crop_list[i][1])

        if shift_list[i]:

            for j, frame in enumerate(act):
                
                movie[i][j] = shift_image(movie[i][j], shift_list[i][0], shift_list[i][1])

    # adjust color scheme
    
    legend = adjust_color_scheme(legend)

    for i, act in enumerate(movie):

        for j, frame in enumerate(act):

            movie[i][j] = adjust_color_scheme(movie[i][j])

    # add images to canvas

    for i, act in enumerate(movie):

        for j, frame in enumerate(act):

            canvas = np.copy(blank_canvas)
            movie[i][j] = add_image_to_canvas(canvas,[movie[i][j],legend], [(0,0),legend_offset])

    # add annotations

    shell_numbers, atom_count_by_shell, current_face = get_annotations(neighbour_info, number_of_views)

    for i, act in enumerate(movie):

        for j, frame in enumerate(act):

            movie[i][j] = annotate_image(movie[i][j], 'View:', annotations_offset[0], annotations_offset[1] + 100)
            movie[i][j] = annotate_image(movie[i][j], 'Shell number:', annotations_offset[0], annotations_offset[1])
            movie[i][j] = annotate_image(movie[i][j], 'Atom count by shell:', annotations_offset[0], annotations_offset[1] + 25)

            movie[i][j] = annotate_image(movie[i][j], shell_numbers[i][j], annotations_offset[0] + 135, annotations_offset[1])
            movie[i][j] = annotate_image(movie[i][j], atom_count_by_shell[i][j], annotations_offset[0] + 200, annotations_offset[1] + 25)
            movie[i][j] = annotate_image(movie[i][j], current_face[i][j], annotations_offset[0] + 50, annotations_offset[1] + 100)
            
    # write each frame to a file

    for i, act in enumerate(movie):

        for j, frame in enumerate(act):

            cv2.imwrite(main_target_folder + 'frame_' + str(j) + '_view_' + str(i) + '.png', movie[i][j]) 
    
    # flatten movie into list of frames

    frames = list(itertools.chain(*movie))

    # combine into gif

    make_gif(frames, file_path)



def crop_image(img,top,bottom):

    height = np.shape(img)[0]
    img = img[top:height - bottom,:,:]

    return img


def shift_image(img, side, cut_size):

    '''
    Cut whitespace and add to the other side to shift the image on the canvas.
    '''

    if side == 'left':

        left_slice = img[:,0:cut_size,:]
        main_slice = img[:,cut_size:,:]
    
        img = cv2.hconcat([main_slice, left_slice])

    if side == 'right':

        width = np.shape(img)[1]
        right_slice = img[:,width - cut_size:,:]
        main_slice = img[:,:width - cut_size,:]
    
        img = cv2.hconcat([right_slice, main_slice])

    return img


def scale_down_image(img, scale_factor):

    scaled_width = int(img.shape[1] * scale_factor)
    scaled_height = int(img.shape[0] * scale_factor)
    dim = (scaled_width, scaled_height)
    
    img = cv2.resize(img, dim, interpolation = cv2.INTER_AREA)

    return img


def adjust_color_scheme(img):

    '''
    Adjust to convert from cv2 to imageio color scheme.
    '''

    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

    return img


def add_image_to_canvas(canvas, image_list, offset_list):

    for i, img in enumerate(image_list):
    
        canvas[offset_list[i][1]:offset_list[i][1]+img.shape[0], offset_list[i][0]:offset_list[i][0]+img.shape[1]] = img

    return canvas


def get_annotations(neighbour_info, number_of_views):

    number_of_neighbours = len(neighbour_info)
    separator = ','

    # get shell info

    shell_numbers = []
    atom_count_by_shell = []

    for i in range(0, number_of_neighbours):

        shell_numbers.append(str(i))  
        atom_count_by_shell.append(str(len(neighbour_info[i])))

        if i < number_of_neighbours - 1:
            
            shell_numbers.append(separator)
            atom_count_by_shell.append(separator)

    # progressivly build string

    number_of_frames = len(shell_numbers)

    shell_number = ''
    atom_count = ''
    progressivly_filling_shell_numbers = []
    progressivly_filling_atoms_count_by_shell = []

    for i in range(0, number_of_frames):

        shell_number = shell_number + shell_numbers[i]
        atom_count = atom_count + atom_count_by_shell[i]
        progressivly_filling_shell_numbers.append(shell_number)
        progressivly_filling_atoms_count_by_shell.append(atom_count)

    # extend shell info to full cycle

    reversed_shell_numbers = list(reversed(progressivly_filling_shell_numbers))
    reversed_atom_count_by_shell = list(reversed(progressivly_filling_atoms_count_by_shell))

    full_cycle_shell_numbers = progressivly_filling_shell_numbers
    full_cycle_shell_numbers.extend(reversed_shell_numbers)
    full_cycle_atom_count_by_shell = progressivly_filling_atoms_count_by_shell
    full_cycle_atom_count_by_shell.extend(reversed_atom_count_by_shell)

   # repeat for all views shown in movie

    repeated_shell_numbers = []
    repeated_atom_count_by_shell = []

    for i in range(0,number_of_views):

        repeated_shell_numbers.append(full_cycle_shell_numbers)
        repeated_atom_count_by_shell.append(full_cycle_atom_count_by_shell)

    # face information

    faces = ['Top','Front','Left']
    current_face = []

    for i in range(0,number_of_views):

        faces_list = []

        for j in range(0,2*number_of_frames):

            faces_list.append(faces[i])

        current_face.append(faces_list)


    shell_numbers = repeated_shell_numbers
    atom_count_by_shell = repeated_atom_count_by_shell

    return shell_numbers, atom_count_by_shell, current_face
  


def annotate_image(img,text,x,y):

    ft = cv2.freetype.createFreeType2()
    ft.loadFontData(fontFileName='Rendering/gif/calibri-regular.ttf',id=0) 
    ft.putText(img, text, (x,y), fontHeight = 20, color = (0,0,0), thickness = -1, line_type=cv2.LINE_AA, bottomLeftOrigin=True)
    
    return img


def make_gif(frames, file_path):

    with imageio.get_writer(file_path , mode="I", fps = 2) as writer:

        for frame in frames:

            writer.append_data(frame)

    writer.close()