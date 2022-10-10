import numpy as np


def generate_grid(cell, shape):
    
    '''
    Returns a 3D grid with scaled positions as the elements. 
    '''

    x_points = np.linspace(0, 1, num=shape[0], endpoint=False)
    y_points = np.linspace(0, 1, num=shape[1], endpoint=False)
    z_points = np.linspace(0, 1, num=shape[2], endpoint=False)

    a, b, c = get_cell_vectors(cell)

    grid_vectors = [[[a*i + b*j + c*k for k in z_points]
                     for j in y_points] for i in x_points]
    grid_vectors = np.asarray(grid_vectors)

    return grid_vectors


def get_cell_vectors(cell):
    
    '''
    Get the cell vectors from the cell matrix.
    '''

    a = cell[0]
    b = cell[1]
    c = cell[2]

    return a, b, c


def convert_to_scaled_positions(cell, positions):
   
    '''
    Projects the input vector along the desired axis.
    '''

    scaled_positions = np.linalg.solve(cell.T, np.transpose(positions)).T

    return scaled_positions
