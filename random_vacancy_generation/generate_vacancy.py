from ase.io import read, write
from ase.build import surface, make_supercell
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.structure import Structure

import os
import numpy as np
import random

# Define the number of vacancies
N_vac = 2       # default is 1 vacancy

# First generate a slab using ASE to get the size of the slab with 7 layers and 0 vacuum
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

bulk_structure = Structure.from_file("bulk.cif")

#  Create SlabGenerator instance
slabgen = SlabGenerator(
    bulk_structure,
    miller_index=(1, 1, 1),  # Miller index of the surface plane
    min_slab_size=min_slab_size,      # Minimum slab thickness in Angstroms
    min_vacuum_size=N_vac*120.0,    # Minimum vacuum size in Angstroms
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
    if len(supercell) == 672:
        # Optionally, write the supercell slab to a file named 111slab.cif only if the condition is met
        supercell.to(fmt="cif", filename="111slab.cif")
        print(f"Supercell Slab {i+1} with 672 atoms saved as 111slab.cif")
    else:
        print(f"Supercell Slab {i+1} does not have 672 atoms and was not saved.")

# Load the generated slab structure in ase 
atoms = read('111slab.cif')

# List to store the indices of oxygen atoms
oxygen_indices = []

# Iterate through the atoms and store the indices of oxygen atoms
for i, atom in enumerate(atoms):
    if atom.symbol == 'O':
        oxygen_indices.append(i)

# Pick N_vac random oxygen atom from the list
random_oxygen_index = []
for i in range(N_vac):
    random_oxygen_index.append(random.choice(oxygen_indices))
    print(f"{i+1}. Randomly selected oxygen atom index:", random_oxygen_index[i])

# Get the position of the randomly selected oxygen atoms
oxygen_position = []
for i in range(N_vac):
    oxygen_position.append(atoms[random_oxygen_index[i]].position)

nearest_cerium_indices = []
for j in range(N_vac):
    # List to store the distances and indices of cerium atoms
    cerium_distances = []

    # Iterate through the atoms and calculate distances to cerium atoms
    for i, atom in enumerate(atoms):
        if atom.symbol == 'Ce':
            distance = np.linalg.norm(atom.position - oxygen_position[j])
            cerium_distances.append((distance, i))

    # Sort the cerium atoms by distance and select the 4 nearest ones
    cerium_distances.sort()
    nearest_cerium_indices.append([index for _, index in cerium_distances[:4]])

    # Print the indices of the 4 nearest cerium atoms
    print(f"Indices of the 4 nearest cerium atoms around {random_oxygen_index[j]}. oxygen atom:", nearest_cerium_indices[j])

# Write the indices to a file
for i in range(N_vac):
    with open(f'vacancy_indices_{random_oxygen_index[i]}.txt', 'w') as f:
        f.write(f"Index of the oxygen vacancy: {random_oxygen_index[i]}\n")
        f.write(f"Indices of the 4 nearest cerium atoms: {nearest_cerium_indices[i]}\n")
        f.write(f"Index of the random oxygen atom: {random_oxygen_index[i]}\n")

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

# The vacancy is 120.0 * N_vac

with open('111slab.cif','r') as file:
    lines = file.readlines()

cell = atoms.get_cell()
a, b ,c = cell.lengths()

for j in range(N_vac):
    keyword = "O" + str(random_oxygen_index[j])
    _z_new = round((60.0 + 120. * j) / c, 8)

    for i, line in enumerate(lines):
        if keyword in line:
            lines[i] = "  O2-  " + keyword + "  1  " + "0.5  0.5  " + str(_z_new) + "  1.0 \n"
            break

with open(f'vacancy_slab_{random_oxygen_index}.cif','w') as file:
    file.writelines(lines)

# !! the index of the atoms in lammps file should start from 1 when cif index starts from 0 !!
vacancy_slab = read(f'vacancy_slab_{random_oxygen_index}.cif')
write(f'vacancy_slab_{random_oxygen_index}.lmp', vacancy_slab, format='lammps-data')

# Add the charges to the atoms in the LAMMPS file
import sys
import os

# Add the directory containing geometry_utils to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from geometry_utils import add_charges_to_atoms
from geometry_utils import add_masses_to_lmp

add_charges_to_atoms(f'vacancy_slab_{random_oxygen_index}.lmp', f'vacancy_slab_{random_oxygen_index}.lmp')
add_masses_to_lmp(f'vacancy_slab_{random_oxygen_index}.lmp', f'vacancy_slab_{random_oxygen_index}.lmp')
