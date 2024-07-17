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

oxygen_input_file = "run-lammps-oxygen.sh"
slab_input_file = "run-lammps-slab.sh"

D_O2 = calculate_O2_dis_energy(oxygen_input_file)

E_slab = calculate_CeO2_slab_energy(slab_input_file)

E_slab_with_vacancy = calculate_CeO2_slab_with_vacancy_energy(slab_input_file)
print(E_slab_with_vacancy)

E_vacancy = E_slab_with_vacancy - E_slab - D_O2/2

print("Energy of CeO2 slab with vacancy:", E_slab_with_vacancy)
print("Energy of CeO2 slab:", E_slab)
print("Dissociation energy of O2:", D_O2)

print("Vacancy formation energy:", E_vacancy)