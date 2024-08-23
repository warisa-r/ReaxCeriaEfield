from mpi4py import MPI
from lammps import lammps
import matplotlib.pyplot as plt

# Initialize the MPI environment
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Initialize LAMMPS with the MPI communicator
lmp = lammps(comm=comm)

# Run lammps npt simulation
lmp.file("run-lammps-npt.lmp")

# Re-initialize lammps instance for NVT (optional, for continuation)
lmp = lammps(comm=comm)
lmp.file("run-lammps-nvt.lmp")
