import numpy as np

def generate_grid(cell,points_per_dimension):

    '''
    Returns a 3D grid with scaled positions as the elements. 
    '''

    points = np.linspace(0, 1, num=points_per_dimension, endpoint= False)
    
    a,b,c = get_cell_vectors(cell)

    grid_vectors = [[[a*i + b*j + c*k for k in points] for j in points] for i in points]
    grid_vectors = np.asarray(grid_vectors)
    
    return grid_vectors, points



def get_cell_vectors(cell):

    '''
    Get the cell vectors from the cell matrix.
    '''

    a = cell[0]
    b = cell[1]
    c = cell[2]
 
    return a,b,c


def project_along_axis(x_in,y_in,z_in,axis):

    '''
    Projects the input vector along the desired axis.
    '''

    projection = np.dot([x_in,y_in,z_in],axis) / np.dot(axis,axis)

    return projection