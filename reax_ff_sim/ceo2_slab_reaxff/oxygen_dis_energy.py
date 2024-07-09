import matplotlib.pyplot as plt
from lammps import lammps
import numpy as np

# Even with varying electric field intensity, the dissociation energy of O2 remains constant
import os

file_path = "run-lammps_O2.sh"
# Read the content of the file
with open(file_path, 'r') as file:
    lines = file.readlines()

# Modify the specific line
with open(file_path, 'w') as file:
    for line in lines:
        if line.strip().startswith('read_data'):
            file.write(f'read_data geometry_files/O2.lmp\n')
        else:
            file.write(line)

conv = 0.0433641

lmp = lammps()
lmp.file("run-lammps_O2.sh")

E_normal = lmp.get_thermo("pe") * conv

# change the magnitude of the electric field
# Read the content of the file
with open(file_path, 'r') as file:
    lines = file.readlines()

# Modify the specific line
with open(file_path, 'w') as file:
    for line in lines:
        if line.strip().startswith('read_data'):
            file.write('read_data geometry_files/O2_dis.lmp\n')
        else:
            file.write(line)

# Second LAMMPS instance
lmp2 = lammps()

lmp2.file("run-lammps_O2.sh")

E_dis = lmp2.get_thermo("pe") * conv - E_normal

print(f"Disassociation energy of O2 = {E_dis} eV")