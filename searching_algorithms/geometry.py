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


def convert_to_scaled_positions(cell,positions):

    '''
    Projects the input vector along the desired axis.
    '''

    scaled_positions = np.linalg.solve(cell.T, np.transpose(positions)).T
   

    return scaled_positions

