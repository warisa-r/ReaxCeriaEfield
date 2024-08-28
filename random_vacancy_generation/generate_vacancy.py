from ase.io import read, write
from ase.build import surface, make_supercell
import numpy as np
import random

# Read the CIF file of the slab
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
