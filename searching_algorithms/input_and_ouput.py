import os
from ase.io import read


def read_trajfile(infile):

  atoms=read(infile)

  positions=atoms.get_positions()
  symbols=atoms.get_chemical_symbols()
  cell = atoms.get_cell()
  
  return atoms, positions, symbols, cell
  
  
def save_file(atoms, name = None):
    
    if name == None:

        atoms.write('in.traj')

    else:

        atoms.write(name)


def change_folder(
    folder_name,
    ):  

    '''
    Change the folder and create a folder if the path doesnt exist already.
    '''

    if not os.path.isdir(folder_name):

        os.mkdir(folder_name)

    os.chdir(folder_name)



