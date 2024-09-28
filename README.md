# Getting started

Report: https://www.overleaf.com/9276872947cwwjgyvrmhjp#c8cfa7

# What is in this repository
1. Module `geometry_utils` with functions to create, add charges to, scale, and visualize lammps geometry files, You can import this module to use these functions in your simulation as you want. All functions have been documented with docstring. Any function that deals with geometry files should be written here. To user this module in `simulations` folder, add these lines
```
import sys
import os

# Add the directory containing geometry_utils to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# example how to import
from geometry_utils import convert_lammps_to_cif
```
2. `literature` (incomplete) is an archive of the literature professor Abhishek has sent to us and is relevant to the project.
3. `simulations` contain folders of simulations in the following structure that should be kept. The idea is that we can run the simulation in that directory right away without having to worry about paths and each of these simulations can be zipped and run by others who have access to the module `geometry_utils`. This pattern should be strictly kept unless we agree on adjusting them:
   1. `run-lammps.lmp` lammps input file
   2. `run.py` python file that calls the input file
   4. geometry files
   5. `ffield.reax` force field parameter file
   6. `job_script.sh` shell script to submit job in the cluster (optional) If you have this please make sue that your python code's essential modules are included.
   7. `generate_geometry.py` python script called in order to generate the geometry files (optional)
   8. `process_result.py` python script called in order to process the simulation results (optional)
   9. dump files (optional)
   10. plot of the simulation results (optional)

# What to do next?
1. [x] Get surface energy correct
2. [x] Discuss: Avoid hardcoding numbers even though it is easier and more convenient.
3. [ ] Idea: check [MDTraj](https://mdtraj.org/1.9.3/index.html) python library that can analyze each frame of the simulation
4. [x] Organize generate_geometry.py in the energy volume scan so that it produces the geometry files present in the energy volume scan simulation (Andres) -> between 13.08.2024 - 15.08.2024
5. [x] Write a function that can convert the position in lammpstrj. file in a specific time step to cif file (Use atom type mapping from lammps geometry file)
6. [ ] Set up an NPT and possibly NVT simulation and run them local on Andres's and Yiyang's computer and get a more or less OK function of lattice constant or temperature against time.
7. [ ] Handle the HPC slurm job submission. Try running NPT or NVT on cluster and plot the same function without certainty in how correct the simulation works: Warisa
8. [ ] Set up HPC account and connect to a node in HPC cluster (no need to go further)
9. [ ] Function that randomly generates vacancy
10. [ ] Run NPT and NVT on cluster
11. [ ] Observe heating behavior with and without vacancy
