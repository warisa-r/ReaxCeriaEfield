import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from data_utils.vacancy_track import read_nearest_cerium_indices, parse_lammpstrj_of_atoms, highlight_vacancy_site

vacant_oxygen_index = 316
cerium_indices = read_nearest_cerium_indices(vacant_oxygen_index)
coordinates = parse_lammpstrj_of_atoms('dump_heat_efield.lammpstrj', cerium_indices)

# Print the coordinates for verification
for frame in coordinates:
    print(f"Timestep: {frame['timestep']}")
    for atom_id, atom_data in frame["atoms"].items():
        print(f"Atom ID: {atom_id}, Coordinates: ({atom_data['x']}, {atom_data['y']}, {atom_data['z']})")

highlight_vacancy_site('dump_heat_efield.lammpstrj', cerium_indices, 'dump_heat_efield_highlighted.lammpstrj')