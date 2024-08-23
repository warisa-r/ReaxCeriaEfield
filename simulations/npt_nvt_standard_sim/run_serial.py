from lammps import lammps

lmp = lammps()

# Run lammps npt simulation
lmp.file("run-lammps-npt.lmp")

# Some method to calculate mean lattice constant

# Re-initialize lammps instance for NVT (optional, for continuation)
lmp = lammps()
# apply mean lattice constant to NVT ( TEST )
_lx = 10.8
_ly = 10.8
_lz = 10.8
# Create simulation box
lmp.file("run-lammps-nvt.lmp")
lmp.command(f"change_box all x final 0 {_lx} y final 0 {_ly} z final 0 {_lz} remap")
lmp.command(f"run {500000}")  # Run in chunks of 500000 steps
lmp.command(f"write_restart   restart.NVT")