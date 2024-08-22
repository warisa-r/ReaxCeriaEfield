from mpi4py import MPI
from lammps import lammps
import matplotlib.pyplot as plt

# Initialize the MPI environment
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Initialize LAMMPS with the MPI communicator
lmp = lammps(comm=comm)

# Define your input script (load it from a file)
lmp.file("run-lammps-npt.lmp")  # Ensure this script doesn't include a final "run" command

# Total number of timesteps to run for NPT
total_steps = 500000
steps_per_run = 500

# Lists to store time (steps) and lattice constants
time_steps = []
lx_values = []
ly_values = []
lz_values = []

# Initialize old lattice constants with None
lx_old, ly_old, lz_old = None, None, None

# Run the simulation in chunks and collect lattice constants
for i in range(0, total_steps, steps_per_run):
    lmp.command(f"run {steps_per_run} pre no post no")  # Run in chunks of 500 steps
    
    # Get the box dimensions after each chunk of 500 steps
    lx_new = lmp.get_thermo("lx")
    ly_new = lmp.get_thermo("ly")
    lz_new = lmp.get_thermo("lz")

    # Store the data
    time_steps.append(i + steps_per_run)
    lx_values.append(lx_new)
    ly_values.append(ly_new)
    lz_values.append(lz_new)

    # Update old lattice constants
    lx_old, ly_old, lz_old = lx_new, ly_new, lz_new

# Only the root process (rank 0) should handle plotting and saving files
if rank == 0:
    # Plotting the results
    plt.figure(figsize=(10, 6))
    plt.plot(time_steps, lx_values, label='lx', marker='o')
    plt.plot(time_steps, ly_values, label='ly', marker='s')
    plt.plot(time_steps, lz_values, label='lz', marker='^')

    plt.xlabel('Time Steps')
    plt.ylabel('Lattice Constant (Angstroms)')
    plt.title('Lattice Constant as a Function of Time')
    plt.legend()
    plt.grid(True)

    # Save the plot as a picture file
    plt.savefig('plot_of_lattice_constant.png')

    # Save the simulation results
    lmp.command("write_restart restart.NPT")  # Save the restart file

# Re-initialize lammps instance for NVT (optional, for continuation)
lmp = lammps(comm=comm)
lmp.file("run-lammps-nvt.lmp")

