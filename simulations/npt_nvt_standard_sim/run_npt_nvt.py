from lammps import lammps
from mpi4py import MPI
import numpy as np

# Initialize the MPI environment
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

lmp = lammps(comm=comm)

# Run lammps npt simulation
lmp.file("run-lammps-npt.lmp")

# Some method to calculate mean lattice constant
from extract_data import lattice_constant_log
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