import pickle
from matplotlib import pyplot as plt
import matplotlib.lines as mlines


def plot():


    X, Y = get_results('results.pkl')

    print(X,Y)

    # convert 

    Y = [[t.total_seconds() for t in list] for list in Y]
    
    # plotting params

    reps = range(0,5)
    runs_per_structure = range(0,5)

    markers = ["o","x"]
    colors = ['b', 'k', 'r']

    file_name = 'scaling_plot.png'

    # plotting
    
    fig = plt.figure()
    fig.set_size_inches(5,4.2)

    for rep in reps:

        for run in runs_per_structure:
        
            plt.plot( X[rep], Y[rep][run], marker = markers[0], markersize = 6, color = colors[0], fillstyle = 'none')


    #axes and titles
     
    plt.xlabel('#Atoms')
    plt.ylabel('Run time (s)')

    # legend

    this_work = mlines.Line2D([], [], color = colors[0], marker = markers[0],linestyle='none' ,fillstyle='none', label = 'Fast neighbour search')
    #ase = mlines.Line2D([], [], color = colors[1], marker = markers[1],linestyle='none' ,fillstyle='none', label = 'ASE neighbour search')
    #KDTree = 

    handles_list=[this_work]
    plt.legend(handles=handles_list, loc = 'lower right')

    # save
    
    plt.savefig(file_name)
   

def get_results(file_name):

    infile = open(file_name,'rb')
    results_dict = pickle.load(infile)
    infile.close()

    system_sizes = results_dict['system_sizes']
    run_times = results_dict['run_times']

    return system_sizes, run_times









