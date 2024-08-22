from lammps import lammps

conv = 0.0433641 # kcal/mol to eV # Reference: http://wild.life.nctu.edu.tw/class/common/energy-unit-conv-table-detail.html

def calculate_O2_dis_energy(oxygen_input_file):
    """
    Calculates the dissociation energy of an oxygen molecule using LAMMPS.

    This function modifies the input LAMMPS script to use a specific data file for a dissociated oxygen molecule,
    runs the simulation, and calculates the dissociation energy based on the potential energy difference
    between the dissociated state and a previously defined ground state.

    Parameters:
    - oxygen_input_file (str): The path to the LAMMPS input script. This script must contain a line starting with
      'read_data' to specify the data file for the simulation. This function modifies this line to use the
      'data_dis.O2' file, which should represent the dissociated state of the oxygen molecule.

    Returns:
    - D_O2 (float): The calculated dissociation energy of the oxygen molecule in eV.

    Note:
    - The 'data_dis.O2' file must be present in the working directory and represent the dissociated state of
      the oxygen molecule.
    """
    
    # Read the content of the file
    with open(oxygen_input_file, 'r') as file:
        lines = file.readlines()

    # Modify the specific line
    with open(oxygen_input_file, 'w') as file:
        for line in lines:
            if line.strip().startswith('read_data'):
                file.write(f'read_data data.O2\n')
            else:
                file.write(line)

    lmp = lammps()
    lmp.file(oxygen_input_file)
    energy_oxygen_ground = lmp.get_thermo("pe")

    # Read the content of the file
    with open(oxygen_input_file, 'r') as file:
        lines = file.readlines()

    # Modify the specific line
    with open(oxygen_input_file, 'w') as file:
        for line in lines:
            if line.strip().startswith('read_data'):
                file.write(f'read_data data_dis.O2\n')
            else:
                file.write(line)

    lmp = lammps()
    lmp.file(oxygen_input_file)

    energy_oxygen_dis = lmp.get_thermo("pe")

    D_O2 = (energy_oxygen_dis - energy_oxygen_ground) * conv

    return D_O2

def calculate_CeO2_slab_energy(slab_input_file):
    """
    Calculates the energy of a CeO2 slab from a LAMMPS input file.

    This function modifies the provided LAMMPS input file to ensure it reads the correct data file
    for the CeO2 slab and outputs the trajectory in a specified format. It then runs the simulation
    using the modified input file and calculates the slab's energy by extracting the potential energy
    from the LAMMPS thermodynamic output. The calculated energy is adjusted by a predefined conversion
    factor ('conv') to convert it to the desired units.

    Parameters:
    - slab_input_file (str): The path to the LAMMPS input file for the CeO2 slab simulation.

    Returns:
    - float: The calculated energy of the CeO2 slab, adjusted by the conversion factor.
    """

    # Read the content of the file
    with open(slab_input_file, 'r') as file:
        lines = file.readlines()

    # Modify the specific line
    with open(slab_input_file, 'w') as file:
        for line in lines:
            if line.strip().startswith('read_data'):
                file.write(f'read_data data.CeO2_111slab\n')
            elif line.strip().startswith('dump            1'):
                file.write(f'dump            1 all custom 10000 steps-ut.lammpstrj id xu yu zu  \n')
            else:
                file.write(line)
    
    lmp = lammps()
    lmp.file(slab_input_file)
    energy = lmp.get_thermo("pe")

    E_slab = energy * conv

    return E_slab

def calculate_CeO2_slab_with_vacancy_energy(slab_input_file):
    """
    Calculates the potential energy of a CeO2 slab with a vacancy.

    This function modifies the input LAMMPS file to use a specific data file
    for a CeO2 slab with a vacancy. It then runs a LAMMPS simulation using the
    modified input file and calculates the potential energy of the system.

    Parameters:
    - slab_input_file (str): The path to the LAMMPS input file for the simulation.

    Returns:
    - float: The calculated potential energy of the CeO2 slab with a vacancy.
    """

    # Read the content of the file
    with open(slab_input_file, 'r') as file:
        lines = file.readlines()

    # Modify the specific line
    with open(slab_input_file, 'w') as file:
        for line in lines:
            if line.strip().startswith('read_data'):
                file.write(f'read_data data_vacancy.CeO2_111slab\n')
            elif line.strip().startswith('dump            1'):
                file.write(f'dump            1 all custom 10000 steps-ut-vacancy.lammpstrj id xu yu zu  \n')
            else:
                file.write(line)
    
    lmp = lammps()
    lmp.file(slab_input_file)
    energy_vacancy = lmp.get_thermo("pe")

    E_slab_with_vacancy = energy_vacancy * conv

    return E_slab_with_vacancy

oxygen_input_file = "run-lammps-oxygen.lmp"
slab_input_file = "run-lammps-slab.lmp"

D_O2 = calculate_O2_dis_energy(oxygen_input_file)

E_slab = calculate_CeO2_slab_energy(slab_input_file)

E_slab_with_vacancy = calculate_CeO2_slab_with_vacancy_energy(slab_input_file)
print(E_slab_with_vacancy)

E_vacancy = E_slab_with_vacancy - E_slab - D_O2/2

print("Energy of CeO2 slab with vacancy:", E_slab_with_vacancy)
print("Energy of CeO2 slab:", E_slab)
print("Dissociation energy of O2:", D_O2)

print("Vacancy formation energy:", E_vacancy)

print("\n")

### Calculate the structure deformation ###

# Read the content of the file
import math
import itertools

def read_last_timestep_data(filename, atom_ids):
    """
    Reads the last timestep data for specified atoms from a LAMMPS trajectory file.

    This function opens a LAMMPS trajectory file, reads it in reverse to find the last timestep,
    and extracts the data for specified atoms based on their IDs. It assumes the file format
    includes an "ITEM: ATOMS" line followed by atom data lines, where each atom data line starts
    with the atom ID.

    Parameters:
    - filename (str): The path to the LAMMPS trajectory file.
    - atom_ids (list of int): A list of atom IDs for which data is to be extracted.

    Returns:
    - dict: A dictionary where each key is an atom ID (from the specified atom_ids list) and
      the value is another dictionary with 'x', 'y', and 'z' keys representing the atom's
      coordinates from the last timestep.
    """

    with open(filename, 'r') as file:
        lines = file.readlines()

    # Reverse the file lines to find the last timestep from the end
    lines.reverse()

    # Variables to hold the data for the atoms of interest
    atom_data = {}

    # Flags to identify when we are reading the desired atoms
    # When reversed lines are read, the first line read is atom index and coordinate
    reading_atoms = True

    for line in lines:
        if line.startswith("ITEM: ATOMS"):
            reading_atoms = False
        elif reading_atoms:
            parts = line.split()
            atom_id = int(parts[0])
            if atom_id in atom_ids:
                # Store the atom data
                atom_data[atom_id] = {
                    'x': float(parts[1]),
                    'y': float(parts[2]),
                    'z': float(parts[3])
                }
                # Check if all desired atoms have been found
                if all(id in atom_data for id in atom_ids):
                    break  # Stop reading since we found all the needed data

    return atom_data

# Function to calculate Euclidean distance between two points
def calculate_distance(point1, point2):
    return math.sqrt((point1['x'] - point2['x'])**2 + (point1['y'] - point2['y'])**2 + (point1['z'] - point2['z'])**2)

# Specify the file and the atom IDs of interest
vacancy_trajectory = 'steps-ut-vacancy.lammpstrj'
normal_trajectory = 'steps-ut.lammpstrj'

# Atom IDs of Cerium surrounding the vacancy. The indices are one of the output from 'generate_geometry.py'
atom_ids = [54, 55, 58, 74]

# Read the data
atom_data_slab = read_last_timestep_data(normal_trajectory, atom_ids)
atom_data_vacancy = read_last_timestep_data(vacancy_trajectory, atom_ids)

# Ensure all atoms' data were found before calculating distances
if all(id in atom_data_slab for id in atom_ids):
    # Iterate through all unique pairs of atoms
    for atom_id1, atom_id2 in itertools.combinations(atom_ids, 2):
        distance_slab = calculate_distance(atom_data_slab[atom_id1], atom_data_slab[atom_id2])
        print(f"Distance between Cerium atom {atom_id1} and {atom_id2}: {distance_slab}")
else:
    print("Data for one or more atoms not found.")

print("\n After deformation: \n")
# Ensure all atoms' data were found before calculating distances
if all(id in atom_data_vacancy for id in atom_ids):
    # Iterate through all unique pairs of atoms
    for atom_id1, atom_id2 in itertools.combinations(atom_ids, 2):
        distance_vacancy = calculate_distance(atom_data_vacancy[atom_id1], atom_data_vacancy[atom_id2])
        print(f"Distance between Cerium atom {atom_id1} and {atom_id2}: {distance_vacancy}")
else:
    print("Data for one or more atoms not found.")