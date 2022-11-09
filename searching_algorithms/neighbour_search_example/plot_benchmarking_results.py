import pickle
from matplotlib import pyplot as plt
import matplotlib.lines as mlines


def plot():


    X, Y = get_results('results.pkl')
    X_ase, Y_ase = get_results('results_ase.pkl')
    X_kdtree, Y_kdtree = get_results('results_kdtree.pkl')

    # convert 

    Y = [[t.total_seconds() for t in list] for list in Y]
    Y_ase = [[t.total_seconds() for t in list] for list in Y_ase]
    Y_kdtree = [[t.total_seconds() for t in list] for list in Y_kdtree]
    
    # plotting params

    reps = range(0,6)
    runs_per_structure = range(0,5)

    markers = ["o","x","+"]
    colors = ['b', 'k', 'r']

    file_name = 'scaling_plot.png'

    # plotting
    
    fig = plt.figure()
    fig.set_size_inches(5,4.2)

    for rep in reps:

        for run in runs_per_structure:
        
            plt.plot( X[rep], Y[rep][run], marker = markers[0], markersize = 6, color = colors[0], fillstyle = 'none')
            plt.plot( X_ase[rep], Y_ase[rep][run], marker = markers[1], markersize = 6, color = colors[1], fillstyle = 'none')
            plt.plot( X_kdtree[rep], Y_kdtree[rep][run], marker = markers[2], markersize = 6, color = colors[2], fillstyle = 'none')


    #axes and titles
     
    plt.xlabel('#Atoms')
    plt.ylabel('Run time (s)')

    # legend

    this_work = mlines.Line2D([], [], color = colors[0], marker = markers[0],linestyle='none' ,fillstyle='none', label = 'Neighbour search')
    ase = mlines.Line2D([], [], color = colors[1], marker = markers[1],linestyle='none' ,fillstyle='none', label = 'ASE Neighborlist')
    kdtree = mlines.Line2D([], [], color = colors[2], marker = markers[2],linestyle='none' ,fillstyle='none', label = 'KDTree')

    handles_list=[this_work,ase,kdtree]
    plt.legend(handles=handles_list, loc = 'upper left')

    # save
    
    plt.savefig(file_name)
   

def get_results(file_name):

    infile = open(file_name,'rb')
    results_dict = pickle.load(infile)
    infile.close()

    system_sizes = results_dict['system_sizes']
    run_times = results_dict['run_times']

    return system_sizes, run_times









