import os 
import itertools

import cv2
import numpy as np
import imageio


from searching_algorithms.connectivity_search_example.run_rendering import get_filepaths 



def combine_images():

    '''
    Comine the individual scenes into one larger PNG file.
    '''

    main_source_folder = 'Rendering/' 
    main_target_folder = 'Rendering/gif/'
    crop_list = [None, (10,0), (180,150)] # (top,bottom)
    shift_list = [None, ('left',40), ('right',10)] # (left,right)
    file_path = 'Rendering/gif/connectivity_search_movie.gif'


    blank_canvas = np.zeros([800,800,3],dtype=np.uint8)
    blank_canvas[:] = 255

    paths = []

    # get paths and sort
    
    for view in range(0,3):

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

    for view in range(0,3):

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


    # add images to canvas

    for i, act in enumerate(movie):

        for j, frame in enumerate(act):

            canvas = np.copy(blank_canvas)
            movie[i][j] = add_image_to_canvas(canvas,movie[i][j], 0, 0)
   
    # write each frame to a file

    for i, act in enumerate(movie):

        for j, frame in enumerate(act):

            cv2.imwrite(main_target_folder + 'frame_' + str(i) + '_view_' + str(j) + '.png', movie[i][j]) 
    
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



def add_image_to_canvas(canvas,img,x_offset,y_offset):
    
    canvas[y_offset:y_offset+img.shape[0], x_offset:x_offset+img.shape[1]] = img

    return canvas


def annotate_image(img,text,x,y):

    font = cv2.FONT_HERSHEY_SIMPLEX
    fontScale = 1
    color = (255,0,255)
    thickness = 2
    img = cv2.putText(img, text, (x,y) , font, fontScale, color, thickness, cv2.LINE_AA)
    
    return img


def make_gif(frames, file_path):

    with imageio.get_writer(file_path , mode="I", fps = 2) as writer:

        for frame in frames:

            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB) 
            writer.append_data(frame)

    writer.close()