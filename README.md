## Summary

This repo contains some examples of searching algorithms for materials simulations including,
* a neighbour searching algorithm for finding the nearest neighbour information for atoms in 3D.
* A breadth-first-search based connectivity search algorithm for finding the number of bonded neighbours between each atom and the starting atom.

## Requirements

* atomic simulation environment
* numpy

This project is built with atomic simulation environment which is a GUI built on Tkinter which supports preprocessing input for a broad variety of different simulation methods. There is a good guide on how to do this in the ase documentation. If you prefer using anaconda like me it is possible to install ase via anaconda which also has great documentation.

conda install -c conda-forge ase

## Find the nearest neighbour information for atoms in 3D

**Background:**

* Nearest neighbour search algorithms are a class of algorithms for finding points in a set which are closest or in a more general sense, most similar, to a given point. Many imaginative and open ended solutions such as graph based methods, clustering and space partitioning based methods have been proposed to solve problems of this type across a broad variety of disciplines.

* A brute force approach which searches through all other atoms for every atom scales like O(N<sup>2</sup>). Computing the full distance matrix in this manner is extremely quick on a personal computer and if we only wish to perform this analysis a few times or for a small number of atoms this straightforward implementation is more than adequate. However, a more careful computational approach becomes extremely useful for large scale simulations with 1000+ atoms where coordination information may be useful at every simulation step. 

* Sorting/partitioning the atoms into bins can bring the scaling down to order O(N) because for a fixed cutoff, the neighbourhood searched for each atom is constant. 




**Implementation and design choices:**

* The code was written as far as possible with a test-driven approach with many sections existing as a test first before being moved into the main code. In fact, to understand the code as quickly as possible I would recommend simply running the tests one by one viewing the output. 

* NumPy was used wherever appropriate to try to speed up the computations. Python itself is a highly-dynamic and flexible language and is considered by most to be at least a factor of 10 slower in most applications ( Ball park figure. In fact, the Python builtins are often written in C with some codes nearly matching pure C in terms of speed. ). The idea in numPy is to use vectorisation, broadcasting and fast indexing routines written in pure C instead of Python for-loops and iterators. We want to cram as much of the logic/computation into operations on a numPy arrays with much larger chunks of memory being visited at a time compared to slow one at a time python iterators. 

* However, in this project there are many instances where we can’t really benefit from NumPy’s vectorised operations either because we are pushed to access the elements individually e.g. when finding all the nearest neighbours of a bin, or we wish to use objects. When accessing individual elements of an array it is well known that NumPy tends to be slower than Python. We are asking Python to access elements in the NumPy array stored in C memory scope one by one and then allocating a new Python object in memory and perhaps creating a pointer to this object in a list, this overhead is more expensive than just using vanilla Python where we can just access a list or a list of lists directly. Similarly, whenever you see dtype=object in some numPy code you should be suspicious. Most of the speed benefits are lost because the array now is filled with pointers to Python objects which are being stored elsewhere in memory much like a Python list.

* Here the atoms information was bundled into a 3D nested list with each element containing a bin object which contains information on all atoms in that bin i.e. positions, atomic symbols and indices. Slots were used in order that each bin uses less memory and the padding routine, periodic boundary condition implementation and indexing turned out to be similarly straightforward to if it had been implemented in numPy which is also known for its easy indexing routines. The attributes of each bin object were represented by numPy arrays providing a nice speed boost to the computations e.g. finding the distance matrix in the neighbour search.

**Code summary:**

<p align="center">
<img src="https://github.com/ashleytsmith/Useful_algorithms_for_materials_modelling/blob/main/Images_for_GitHub/neighbour_search_algo_overview.png" width="400" alt="see images folder if image doesn't show"> 
</p>


## Connectivity search with ray tracing visualisation

*	First the neighbour searching algorithm was applied to the zeolite crystal SSZ-13 which CHA framework type.
*	Bin sizes of above 4 angstroms and atom counts per bin of between 1 and 8 tended to perform best for this system.
*	The result is a bond dictionary which, in this application, can be thought of as tree structure, where atoms are the vertices and bonds are the edges. The oxygen bridges in the bond dictionary were skipped and used as the input for a breadth-first-search (BFS) searching algorithm.
*	The BFS identifies the shells an atom belongs to relative to the starting atom and is visualized below for a 2 by 2 by 2 supercell. 
*	The movie was made using POV-Ray with the choice of camera perspectives and use of light making it much easier to view the pores and cavities of the structure.

<p align="center">
<img src="https://github.com/ashleytsmith/Useful_algorithms_for_materials_modelling/blob/main/Images_for_GitHub/connectivity_search_movie.gif" width="400" alt="see images folder if image doesn't show"> 
</p>


## Applications and extension ideas: ##

* The neighbour search should work for any periodic system regardless of cell shape and could easily be adapted to work for non-periodic systems. Simply trying different systems and input parameters would be a quick way to get your hands dirty with the code. A good way to extend the work would be to use the neighbour searching algorithm in more applications. The ase.neighborlist docs could be a good source of inspiration and has examples of applications involving, coordination counting, pair distribution functions, pair potentials and dynamical matrices.

* A nice and potentially short extension to the breadth first search algorithm would be to extend it to do topological searches. [One way of doing this](https://www.sciencedirect.com/science/article/abs/pii/S1387181105003756) is via a series of breath first searches and utilising cyclicality conditions.

* Lastly, another interesting direction would be build your own GUI for viewing chemical structures and exporting them to POV-Ray.  A skeleton all the features you may want is present in other open source GUIs, so a good start would be to try and reproduce these key features. Then one could take it further by really focusing on making it easy to make ray traced pictures and movies perhaps using POV-Ray commands directly which already provide an intuitive way to visualise the scene. 

