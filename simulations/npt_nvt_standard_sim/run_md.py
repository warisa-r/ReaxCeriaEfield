
from lammps import lammps
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt

# Initialize the MPI environment
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

lmp = lammps(comm=comm)

lmp.file("run-lammps-md.lmp")
#lmp.file("run-lammps-md-efield.lmp")

#Number of steps tbd?
lmp.command(f"run {500000}")


from extract_data import temperatures_log
from extract_data import steps_log
temperatures = temperatures_log("log.lammps")
steps = steps_log("log.lammps")

plt.figure(figsize=(10, 6))
plt.plot(steps, temperatures, marker='o', linestyle='-', color='b')
plt.xlabel('Step')
plt.ylabel('Temperature (K)')
plt.title('Temperature vs. Step')
plt.grid(True)

# Save the plot as a PNG file
plt.savefig('temperature_vs_step.png')
