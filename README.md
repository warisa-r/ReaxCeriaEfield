# Getting started

1. Install [lammps11Aug1017](https://download.lammps.org/tars/lammps-11Aug2017.tar.gz) (version with both reax/fortran and reax/c) as shared library on your operating system.
2. Click [here](https://help.itc.rwth-aachen.de/service/rhr4fjjutttf/article/598d0f7f78cb4ab8b81af1b3f68ba831/) to read instruction on how to Connect to hpc node.
3. Follow the instruction in install_lammps_mpi to install it to your cluster node.
4. Additionally install packages according to your simulations
   1. MISC for applying electric field and using `fix efield` command
   2. reax/c for running `energy_volume_scan_sim` simulation
<br>

Note: If you decide to install lammps in your os system, the instruction might not be applicable and you might have to find your own way around tackling with dependencies.

# What is in this repository
1. Module `geometry_utils` with functions to create, add charges to, scale, and visualizie lammps geometry files, You can import this module to use these functions in your simulation as you want. All functions have been documented with docstring. Any function that deals with geometry files should be written here.
2. `literature` (incomplete) is an archive of the literatures professor Abihshek has sent to us and are relevent to the project.
3. `simulations` contain folders of simulations in the following structure that should be kept. The idea is that we can run the simulation in that directory right away without having to worry about paths and each of these simulation can be zipped and run by others who have access to the module `geometry_utils`. This pattern should be strictly kept unless we agree on adjusting them:
   1. `run-lammps.sh` input file
   2. `run.py` python file that calls the input file
   4. geometry files
   5. `ffield.reax` force field parameter file
   6. `generate_geometry.py`python script called in order to generate the geometry files (optional)
   7. `process_result.py` python script called in order to process the simulation results (optional)
   8. dump files (optional)
   9. plot of the simulation results (optional)
