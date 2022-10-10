## Summary

* Collection of useful algorithms for materials simulations.
* Code cleanup and write up still in progress.


## Requirements

* atomic simulation environment
* numpy

This project is built with atomic simulation environment which is a GUI built on Tkinter which supports preprocessing input for a broad variety of different simulation methods. There is a good guide on how to do this in the ase documentation. If you prefer using anaconda like me it is possible to install ase via anaconda which also has great documentation.

```
conda install -c conda-forge ase
```

## Example 1: Find the nearest neighbour information for atoms in 3D

**Background:**

* Nearest neighbour search algorithms are a class of algorithms for finding points in a set which are closest or in a more general sense, most similar, to a given point. Many imaginative and open ended solutions such as graph based methods, clustering and space partitioning based methods have been proposed to solve problems of this type across a broad variety of disciplines.

* A brute force approach, for example, using a KDtree and searching through all other atoms for every atom scales like O(N<sup>2</sup>). Computing the full distance matrix in this manner is extremely quick on a personal computer and if we only wish to perform this analysis a few times or for a small number of atoms this straightforward implementation is more than adequate. However, a more careful computational approach becomes extremely useful for large scale simulations with 1000+ atoms where coordination information may be useful at every simulation step. 

* Sorting/partitioning the atoms into bins can bring the scaling down to order O(N) because for a fixed cutoff, the neighbourhood searched for each atom is constant. Bin sorting based approaches to the 3D neighbour search problem were originally pioneered in molecular dynamics simulations but also can be generally quite useful in any applications involving coordination counting, pair distribution functions, pair potentials or dynamical matrices. The speed increase for a 10000 atom periodic system can be around 100X and can be as high as 1000X. This comes at the cost of using more memory to store the bin information and eventually for a very large number of cuts per axis (1000+) the number of bins needs to be limited leading to less effective scaling.



**Implementation and design choices**


* The code was written as far as possible with a test-driven approach with many sections existing as a test first before being moved into the main code. Putting effort into writing intuitive and stable tests saved a lot of time when debugging. In fact to understand the code as quickly as possible I would recommend simply running the tests one by one viewing the output. 

* NumPy was used wherever appropriate to try to speed up the computations. Python itself is a highly-dynamic and flexible language and is considered by most to be at leaset a factor of 10 slower in most applications ( Ball park figure. In fact, the Python builtins are often written in C with some codes nearly matching pure C in terms of speed. ). The idea in numpy is to use vectorisation, broadcasting and fast indexing routines written in pure C instead of Python for-loops and iterators. We want to cram as much of the logic/computation into operations on a numpy arrays with much larger chunks of memory being visited at a time compared to slow one at a time python iterators. 

* However, in this project there are many instances where we can’t really benefit from NumPy’s vectorised operations either because we are pushed to access the elements individually e.g. find all the nearest neighbours of a bin, or we wish to use objects. When accessing individual elements of an array it is well known that NumPy tends to be slower than Python. We are asking Python to access elements in the NumPy array stored in C memory scope one by one, and then allocating a new Python object in memory and perhaps creating a pointer to this object in a list, this overhead is more expensive than just using vanilla Python where we can just access a list or a list of lists directly. A good rule of thumb is if you are trying to do absolutely everything with NumPy you’re likely over optimising and it may even turn out to cause much bigger problems for you later on. Similarly, whenever you see dtype=object in some numPy code you should be suspicious. Most of the speed benefits are lost because the array now is filled with pointers to Python objects which are being stored elsewhere in memory much like a Python list.

* Here the atoms information was bundled into a 3D nested list with each element containing a bin object which contains information on all atoms in that bin i.e. positions, atomic symbols and indices. Slots were used in order that each bin uses less memory and the padding routine, periodic boundary condition implementation and indexing turned out to be similarly straight forward to numPy which is also known for its easy indexing routines. The attributes of the bins object were represented by numPy arrays providing a nice speed boost to the computations e.g. finding the distance matrix in the neighbour search.