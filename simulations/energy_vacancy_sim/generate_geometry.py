# Prepare a slab geometry from the bulk geometry
from ase.io import read, write
from ase.build import surface, make_supercell
import numpy as np
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.structure import Structure

#### 1. Generate a (111) slab geometry from the bulk geometry ####

# First generate a slab using ASE to get the size of the slab with 7 layers and 0 vacuum
structure = read("bulk.cif")
slab = surface(structure, (1,1,1), layers=7, vacuum=0.0)
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

bulk_structure = Structure.from_file("bulk.cif")

#  Create SlabGenerator instance
slabgen = SlabGenerator(
    bulk_structure,
    miller_index=(1, 1, 1),  # Miller index of the surface plane
    min_slab_size=min_slab_size,      # Minimum slab thickness in Angstroms
    min_vacuum_size=120.0,    # Minimum vacuum size in Angstroms
    lll_reduce=False         # Whether to use the LLL algorithm to reduce the cell
)

# Generate slabs with different termination types
slabs = slabgen.get_slabs(
    symmetrize=True,     # Whether to symmetrize the slab
    tol=0.1,              # Symmetrization tolerance
    bonds=None,           # Custom bonds to pass to the slab generator
    max_broken_bonds=0   # Maximum number of broken bonds allowed in the slab
)

# Create a 4x4 supercell in the xy-axis for each slab and save them
for i, slab in enumerate(slabs):
    # Create the 2x2 supercell for the slab
    supercell_scaling_matrix = [[4, 0, 0], [0, 4, 0], [0, 0, 1]]
    supercell = slab.make_supercell(supercell_scaling_matrix)
    
    # Check if the total number of atoms in the supercell is equal to 336
    if len(supercell) == 336:
        # Optionally, write the supercell slab to a file named 111slab.cif only if the condition is met
        supercell.to(fmt="cif", filename="111slab.cif")
        print(f"Supercell Slab {i+1} with 336 atoms saved as 111slab.cif")
    else:
        print(f"Supercell Slab {i+1} does not have 336 atoms and was not saved.")

#### 2. Locate the vacancy and cerium atoms surrounding it ####

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
    print("The oxygen atom closest to the middle point is:") # +1 index because lammps' atom index starts at 1 when python array index starts at 0
    print(f"Index: {middle_oxygen_atom.index + 1}, Position: {middle_oxygen_atom.position}, Fractional Position: {fractional_positions}")

    cerium_distances = []

    for atom in supercell:
        if atom.symbol == 'Ce':
            distance = np.linalg.norm(atom.position - middle_oxygen_atom.position)
            cerium_distances.append((distance, atom.index, atom.position))

    # Sort the list by distance (the first element of the tuple)
    cerium_distances.sort(key=lambda x: x[0])

    # Select the top 4 closest cerium atoms
    closest_cerium_atoms = cerium_distances[:4]

    print("Indices and positions of the 4 closest cerium atoms to the middle oxygen atom:") # This is to locate the cerium atoms for structure deformation
    for distance, index, position in closest_cerium_atoms:
        print(f"Index: {index + 1}, Position: {position}, Distance: {distance}")

    print(f"z_max is {np.max(positions[:, 2])}") # To locate where the vacuum is
    print(f"z_min is {np.min(positions[:, 2])}")
else:
    print("No oxygen atoms found.")

#### 3. After manually adding vacancy in lammps file, visualize to check if everything is correct ####

import sys
import os

# Add the directory containing geometry_utils to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from geometry_utils import convert_lammps_to_cif, add_charges_to_atoms

convert_lammps_to_cif('data_vacancy.CeO2_111slab', 'charge', '111slab_vacancy.cif')
