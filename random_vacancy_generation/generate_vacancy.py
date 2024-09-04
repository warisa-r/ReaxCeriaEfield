from ase.io import read, write
from ase.build import surface, make_supercell
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.structure import Structure

import os
import numpy as np
import random

# Use the slab in npt_nvt_simulation as a base structure before vacancy generation

structure = read("bulk.cif")
slab = surface(structure, (1,1,1), layers=14, vacuum=0.0)
supercell = make_supercell(slab,[[2,0,0],[0,2,0],[0,0,1]])
write('111slab_ase.cif', supercell, format='cif')

# Get cell length for slab size in slab generation by pymatgen
def get_cell_length_c(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('_cell_length_c'):
                # Extract the value after the key
                _, value = line.split()
                return float(value)
    raise ValueError("_cell_length_c not found in the file.")

file_path = '111slab_ase.cif'
min_slab_size = get_cell_length_c(file_path)

# Delete the file '111slab_ase.cif'
os.remove(file_path)

# Load the generated slab structure in ase 
atoms = read('111slab.cif')

# List to store the indices of oxygen atoms
oxygen_indices = []

# Iterate through the atoms and store the indices of oxygen atoms
for i, atom in enumerate(atoms):
    if atom.symbol == 'O':
        oxygen_indices.append(i)

# Pick a random oxygen atom from the list
random_oxygen_index = random.choice(oxygen_indices)
print("Randomly selected oxygen atom index:", random_oxygen_index)

# Get the position of the randomly selected oxygen atom
oxygen_position = atoms[random_oxygen_index].position

# List to store the distances and indices of cerium atoms
cerium_distances = []

# Iterate through the atoms and calculate distances to cerium atoms
for i, atom in enumerate(atoms):
    if atom.symbol == 'Ce':
        distance = np.linalg.norm(atom.position - oxygen_position)
        cerium_distances.append((distance, i))

# Sort the cerium atoms by distance and select the 4 nearest ones
cerium_distances.sort()
nearest_cerium_indices = [index for _, index in cerium_distances[:4]]

# Print the indices of the 4 nearest cerium atoms
print("Indices of the 4 nearest cerium atoms:", nearest_cerium_indices)

# Write the indices to a file
with open(f'vacancy_indices_{random_oxygen_index}.txt', 'w') as f:
    f.write(f"Index of the oxygen vacancy: {random_oxygen_index}\n")
    f.write(f"Indices of the 4 nearest cerium atoms: {nearest_cerium_indices}\n")
    f.write(f"Index of the random oxygen atom: {random_oxygen_index}\n")

# Use min_slab_size or the surface thickness without vacancy by generating a new supercell without vacuum

# Extract the positions of all atoms
positions = atoms.positions

# Calculate the minimum and maximum values for x and y coordinates
x_min, y_min = np.min(positions[:, 0]), np.min(positions[:, 1])
x_max, y_max = np.max(positions[:, 0]), np.max(positions[:, 1])

# Compute the middle points for x and y coordinates
x_middle = (x_min + x_max) / 2
y_middle = (y_min + y_max) / 2

# Print the middle coordinates
print("Middle of the x and y plane:", (x_middle, y_middle))

# The vacancy is 120.0

with open('111slab.cif','r') as file:
    lines = file.readlines()

cell = atoms.get_cell()
a, b ,c = cell.lengths()

keyword = "O" + str(random_oxygen_index)
_z_new = round(60.0 / c, 8)

for i, line in enumerate(lines):
    if keyword in line:
        lines[i] = "  O2-  " + keyword + "  1  " + "0.5  0.5  " + str(_z_new) + "  1.0 \n"
        break

with open(f'vacancy_slab_{random_oxygen_index}.cif','w') as file:
    file.writelines(lines)

#TODO: Turn the cif files to lammps geometry files
# !! the index of the atoms in lammps file should start from 1 when cif index starts from 0 !!
vacancy_slab = read(f'vacancy_slab_{random_oxygen_index}.cif')
write(f'vacancy_slab_{random_oxygen_index}.lmp', vacancy_slab, format='lammps-data')
# No need to change index right? If we wanna track the atom, just modify the index by +1?