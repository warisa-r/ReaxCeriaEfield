# Prepare a slab geometry from the bulk geometry
from ase.io import read, write
from ase.build import surface, make_supercell
import numpy as np

#### 1. Generate a (111) slab geometry from the bulk geometry #### 
structure = read("bulk.cif")
slab = surface(structure, (1,1,1), layers=7, vacuum=120.0)
supercell = make_supercell(slab,[[2,0,0],[0,2,0],[0,0,1]])
write('111slab.cif', supercell, format='cif')

#### 2. Locate the vacancy ####

# Read the slab
supercell = read('111slab.cif')

# Extract the positions of all atoms
positions = supercell.positions

# Find the middle points in x, y, and z directions
x_middle = (np.min(positions[:, 0]) + np.max(positions[:, 0])) / 2
y_middle = (np.min(positions[:, 1]) + np.max(positions[:, 1])) / 2
z_middle = (np.min(positions[:, 2]) + np.max(positions[:, 2])) / 2

# Define a tolerance to consider atoms within the middle region
tolerance = (np.max(positions[:, 2]) - np.min(positions[:, 2])) / 14  # Adjust if needed

# Initialize variables to find the closest oxygen atom to the middle point
min_distance = float('inf')
middle_oxygen_atom = None

# Identify the oxygen atom closest to the middle point in all three dimensions
for atom in supercell:
    if atom.symbol == 'O':
        distance = np.sqrt(
            (atom.position[0] - x_middle)**2 +
            (atom.position[1] - y_middle)**2 +
            (atom.position[2] - z_middle)**2
        )
        if distance < min_distance:
            min_distance = distance
            middle_oxygen_atom = atom

# Print the coordinates of the closest oxygen atom to the middle point
if middle_oxygen_atom:
    # Get fractional coordinates
    fractional_positions = supercell.get_scaled_positions()[middle_oxygen_atom.index]
    print("The oxygen atom closest to the middle point is:")
    print(f"Index: {middle_oxygen_atom.index}, Position: {middle_oxygen_atom.position}, Fractional Position: {fractional_positions}")
else:
    print("No oxygen atoms found.")