import os 

import cv2
import numpy as np
import imageio


from searching_algorithms.connectivity_search_example.run_rendering import get_filepaths 



def combine_images():

    '''
    Comine the individual scenes into one larger PNG file.
    '''

    main_folder = 'Rendering/gif/'
    img0_loc = (0,0)
    img1_loc = (0,250)
    img2_loc = (400,0)
    file_path = 'Rendering/gif/connectivity_search_movie.gif'


    blank_canvas = np.zeros([800,800,3],dtype=np.uint8)
    blank_canvas[:] = 255

    paths = []


    # get paths and sort
    
    for view in range(0,3):

        path_list = get_filepaths( '*view_' + str(view) + '.png')
        frames = [int(path.split('frame_')[-1].split('_')[0]) for path in path_list]
        path_list = [x for _, x in sorted(zip(frames, path_list))]
        paths.append(path_list)


    # add images to canvas

    num_frames = len(paths[0])
    frames = []

    if not os.path.isdir(main_folder):
        
        os.mkdir(main_folder)

    for i in range(0,num_frames):

        canvas = np.copy(blank_canvas)

        img0 = cv2.imread(paths[0][i])
        img1 = cv2.imread(paths[1][i])
        img2 = cv2.imread(paths[2][i])

        canvas = add_image_to_canvas(canvas,img0,img0_loc[0],img0_loc[1])
        canvas = add_image_to_canvas(canvas,img1,img1_loc[0],img1_loc[1])
        canvas = add_image_to_canvas(canvas,img2,img2_loc[0],img2_loc[1])
        
        cv2.imwrite(main_folder + 'frame_' + str(i) + '.png', canvas) 
        frames.append(canvas)

    # combine into gif

    make_gif(frames, file_path)



def crop_image(img,top,bottom):

    height = np.shape(img)[0]
    img = img[top:height - bottom,:,:]

    return img


def crop_images(paths,top,bottom):

    images = []

    for path in paths:

        img = crop_image(img,top,bottom)
        images.aapend(img)

    return images


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

    with imageio.get_writer(file_path , mode="I") as writer:

        for frame in frames:

            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB) 
            writer.append_data(frame)

    writer.close()