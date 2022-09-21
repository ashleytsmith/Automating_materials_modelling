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


* Sorting/partitioning the atoms into bins can bring the scaling down to order O(N) because for a fixed cutoff, the neighbourhood searched for each atom is constant. Bin sorting based approaches to the 3D neighbour search problem were originally pioneered in molecular dynamics simulations but can be generally quite useful in any applications involving coordination counting, pair distribution functions, pair potentials or dynamical matrices. The speed increase for a 10000 atom periodic system can be around a factor of 100 and can be as high as 1000. This comes at the cost of using more memory to store the bin information and eventually for a very large number of cuts per axis (1000+) the number of bins needs to be limited leading to less effective scaling.

* Sorting/partitioning the atoms into bins can bring the scaling down to order O(N) because for a fixed cutoff, the neighbourhood searched for each atom is constant. Bin sorting based approaches to the 3D neighbour search problem were originally pioneered in molecular dynamics simulations but also can be generally quite useful in any applications involving coordination counting, pair distribution functions, pair potentials or dynamical matrices. The speed increase for a 10000 atom periodic system can be around 100X and can be as high as 1000X. This comes at the cost of using more memory to store the bin information and eventually for a very large number of cuts per axis (1000+) the number of bins needs to be limited leading to less effective scaling.
