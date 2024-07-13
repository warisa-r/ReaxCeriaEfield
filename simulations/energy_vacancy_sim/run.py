from lammps import lammps

conv = 0.0433641 # kcal/mol to eV

oxygen_input_file = "run-lammps-oxygen.sh"

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

D_O2 = calculate_O2_dis_energy(oxygen_input_file)
print(D_O2)