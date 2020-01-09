# Hard spheres packings within cylinders: sequential linear programming

A python implementation.

### References:
- J. Norman, E. L. Sorrell, Y. Hu, V. Siripurapu, J. Garcia, J. Bagwell, P. Charbonneau, S. R. Lubkin, and M. Bagnat, Philosophical Transactions of the Royal Society B: Biological Sciences, 373, 1759, 20170320. doi:10.1098/rstb.2017.0320
- L. Fu, W. Steinhardt, H. Zhao, J. E. Socolar and P. Charbonneau, Soft matter, 2016, 12, 2505-2514. doi:10.1039/C5SM02875B
- S. Torquato and Y. Jiao, Physical Review E, 2010, 82, 061302. doi:10.1103/PhysRevE.82.061302
- Torquato-Jiao Sequential Linear Programming Sphere Packing Algorithm
https://github.com/sdatkinson/TJ

## Summary
Finding the dense packing structure of hard spheres in cylindrical pores by sequantial linear programming  (SLP) mthod. This python implementation applied scipy.optimize.linprog

## File Description

### Folders

- CircularCylinder/: circular cylinder solver files
- EllipticalCylinder/: elliptical cylinder solver files

### File specification
- batchrun.py : sample batch running script
- config.json : configuration file
- iomodule.py : module related to input and outputs
- lp.py : sequantial linear programming solver
- main.py : file to begin with the calculation
- mcmove.py : conducting random moves
- neilist.py : neighbor list module
- randomconfig.py : generating random configurartion
- tools.py : useful helper functions
- visualize.py : visualizing the particle coordinates

### Prerequisite:
- Python (2 or 3)
- numpy, scipy
- mayavi (for visualize): mayavi may need wxpython

## Usage
- python main.py : conduct one simulation, where input parameters is from config.json
- python main.py 10 : conduct 10 repeat simulations
- python visualize.py resultConfig_0.dat : visualize the final configuration, read from resultConfig_0.dat