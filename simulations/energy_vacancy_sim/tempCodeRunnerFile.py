
import sys
import os

# Add the directory containing geometry_utils to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from geometry_utils import convert_lammps_to_cif, add_charges_to_atoms

convert_lammps_to_cif('data_vacancy.CeO2_111slab', 'atomic', '111slab_vacancy.cif')

#### 4. Add charges to both geometry files to prepare them for the simulation ####

add_charges_to_atoms('data_vacancy.CeO2_111slab', 'data_vacancy.CeO2_111slab')
add_charges_to_atoms('data.CeO2_111slab', 'data.CeO2_111slab')