from lammps import lammps
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt


# Initialize the MPI environment
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

'''
lmp = lammps(comm=comm)

# Run lammps npt simulation
lmp.file("run-lammps-npt.lmp")

# Some method to calculate mean lattice constant
from extract_lattice import lattice_constant_log
timesteps, lx_values, ly_values, lz_values = lattice_constant_log("log.lammps")

lx_mean = np.mean(lx_values)
ly_mean = np.mean(ly_values)
lz_mean = np.mean(lz_values)

# Re-initialize lammps instance for NVT (optional, for continuation)
lmp = lammps(comm=comm)
# apply mean lattice constant to NVT ( TEST )
_lx = lx_mean
_ly = ly_mean
_lz = lz_mean
# Create simulation box
lmp.file("run-lammps-nvt.lmp")
lmp.command(f"change_box all x final 0 {_lx} y final 0 {_ly} z final 0 {_lz} remap")
lmp.command(f"run {500000}")  # Run in chunks of 500000 steps
lmp.command(f"write_restart   restart.NVT")
'''
lmp = lammps(comm=comm)
lmp.file("run-lammps-efield.lmp")

# Function to extract temperature from log file
def extract_temperature(log_file):
    timesteps = []
    temperatures = []
    start_parsing = False
    with open(log_file, 'r') as file:
        for line in file:
            if "Step" in line and "Temp" in line:
                start_parsing = True
                continue
            if start_parsing:
                if line.strip() == "":
                    break
                data = line.split()
                try:
                    timesteps.append(int(data[0]))
                    temperatures.append(float(data[1]))
                except (IndexError, ValueError):
                    continue
    return np.array(timesteps), np.array(temperatures)

# Extract temperature data from log file
timesteps, temperatures = extract_temperature("log.lammps")

# Plot temperature vs. time
plt.figure()
plt.plot(timesteps, temperatures, label='Temperature')
plt.xlabel('Timestep')
plt.ylabel('Temperature (K)')
plt.title('Temperature vs. Timestep')
plt.legend()
plt.savefig('temperature_plot.png')

