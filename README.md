## Summary

This repo contains some examples of searching algorithms for materials simulations including,
* a neighbour searching algorithm for finding the nearest neighbour information for atoms in 3D.
* A breadth-first-search based connectivity search algorithm for finding the number of bonded neighbours between each atom and the starting atom.

## Requirements

* atomic simulation environment
* NumPy
* POV-Ray

This project utilises atomic simulation environment (ase) which is a GUI built on Tkinter which supports preprocessing input for a broad variety of different simulation methods. There is a good guide on how to do this in the ase documentation. If you prefer using anaconda like me it is possible to install ase via anaconda which also has great documentation.

conda install -c conda-forge ase

## Find the nearest neighbour information for atoms in 3D

**Background:**

* Nearest neighbour search algorithms are a class of algorithms for finding points in a set which are closest or in a more general sense, most similar, to a given point. Many imaginative and open ended solutions such as graph based methods, clustering and space partitioning based methods have been proposed to solve problems of this type across a broad variety of disciplines.

* A brute force approach which searches through all other atoms for every atom scales like O(N<sup>2</sup>). Computing the full distance matrix in this manner is extremely quick on a personal computer and if we only wish to perform this analysis a few times or for a small number of atoms this straightforward implementation is more than adequate. However, a more careful computational approach becomes extremely useful for large scale simulations with 1000+ atoms where coordination information may be useful at every simulation step. 

* Sorting/partitioning the atoms into bins can bring the scaling down to order O(N) because for a fixed cutoff, the neighbourhood searched for each atom is constant. 




**Design considerations:**

* The code was written as far as possible with a test-driven approach with many sections existing as a test first before being moved into the main code. In fact, to understand the code as quickly as possible I would recommend simply running the tests one by one viewing the output. 

* NumPy was used wherever appropriate to try to speed up the computations. The idea in numPy is to use vectorisation, broadcasting and fast indexing routines written in pure C instead of Python for-loops and iterators. We want to cram as much of the logic/computation into operations on a numPy arrays with much larger chunks of memory being visited at a time compared to slow one at a time python iterators. 

* However, in this project there are many instances where we can’t really benefit from NumPy’s vectorised operations e.g., we are pushed to access the elements individually when finding all the nearest neighbours of a bin. When accessing individual elements of an array it is well known that NumPy tends to be slower than Python. We are asking Python to access elements in the NumPy array stored in C memory scope one by one and then allocating a new Python object in memory and perhaps creating a pointer to this object in a list, this overhead is more expensive than just using vanilla Python where we can just access a list or a list of lists directly.

**Code summary:**

<p align="center">
<img src="https://github.com/ashleytsmith/Useful_algorithms_for_materials_modelling/blob/main/Images_for_GitHub/neighbour_search_algo_overview.png" width="400" alt="see images folder if image doesn't show"> 
</p>

**Performance:**

Needs work. This was the first attempt to get a clean, easy to follow and working solution to a surprisingly tricky searching problem. Currently it works well for 100s of atoms but does not scale well for larger system sizes and should not be expected to after just one iteration. It is a good starting template and there is plenty of low hanging fruit (see the "Applications and extension ideas" section for some suggestions) which could yield speed boosts. 


## Connectivity search with ray tracing visualisation

*	First the neighbour searching algorithm was applied to the zeolite crystal SSZ-13 which has CHA framework type.
*	Bin sizes of around 4 angstroms and atom counts per bin of between 1 and 8 tended to perform best for this system.
*	The result is a bond dictionary which, in this application, can be thought of as tree structure where the atoms the vertices and the bonds are the edges. The oxygen bridges in the bond dictionary were skipped and used as the input for a breadth-first-search (BFS) searching algorithm.
*	The BFS identifies the shells an atom belongs to relative to the starting atom and is visualized below for a 2 by 2 by 2 supercell. 
*	The movie was made using POV-Ray with the choice of camera perspectives and use of light making it much easier to view the pores and cavities of the structure.

<p align="center">
<img src="https://github.com/ashleytsmith/Useful_algorithms_for_materials_modelling/blob/main/Images_for_GitHub/connectivity_search_movie.gif" width="400" alt="see images folder if image doesn't show"> 
</p>


## Applications and extension ideas ##

* The neighbour search algorithm is a pet project which will likely require a fair amount of work to get it to scale well. To improve the scaling the most pressing things to focus on, in order, would be: Improving the memory layout e.g., by just considering the positions alone and only allowing float data types instead of putting all information in one bin. Implementing a much better indexing routine for the neighbour pairs, likely flattening it out before looping through in the distance calc. And lastly, incorporating the periodic boundary condition routine into the distance calculation rather than using explicit padding.


* A nice and potentially short extension to the breadth first search algorithm would be to extend it to do topological searches. [One way](https://www.sciencedirect.com/science/article/abs/pii/S1387181105003756) of doing this is via a series of BFSs and utilising cyclicality conditions.

* Lastly, another interesting direction would be build your own GUI for viewing chemical structures and exporting them to POV-Ray.  A skeleton all the features you may want is present in other open source GUIs, so a good start would be to try and reproduce these key features. Then one could take it further by really focusing on making it easy to make ray traced pictures and movies perhaps using POV-Ray commands directly which already provide an intuitive way to visualise the scene. 

