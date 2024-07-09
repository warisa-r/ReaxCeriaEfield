import matplotlib.pyplot as plt
from lammps import lammps
import numpy as np

import os

# Directory where the plot will be saved
plot_directory = '/mnt/c/Users/user/Documents/School/FFparamOpt/NewAttempt_supercell/energy_volume_with_efield/'

# Check if the directory exists, and create it if it does not
if not os.path.exists(plot_directory):
    os.makedirs(plot_directory)

# Define the path to the shell script
file_path = 'run-lammps.sh'
file_path2 = 'run-lammps_vacancy.sh'

conv = 0.0433641 # kcal/mol to eV
D_e_O2 = 4.876615876761253 # in eV Dissociation energy of O2 molecule fromm running oxygen_dis_energy.py

data_per_e_field = {}

e_intensities = np.arange(30, 40, 1)
E_slab_values = []
E_vacancy_values = []

for e_intensity in e_intensities:

    # change the magnitude of the electric field
    # Read the content of the file
    with open(file_path, 'r') as file1:
        lines1 = file1.readlines()

    # Modify the specific line
    with open(file_path, 'w') as file1:
        for line in lines1:
            if line.strip().startswith('variable fz'):
                file1.write(f'variable fz atom q*{e_intensity:.1f}\n')
            elif line.strip().startswith('variable efz'):
                file1.write(f'variable efz atom q*z*{e_intensity:.1f}\n')
            else:
                file1.write(line)
    
    try:
        lmp1 = lammps()
        lmp1.file("run-lammps.sh")
        E_slab = lmp1.get_thermo("pe") * conv
        E_slab_values.append(E_slab)

        with open(file_path2, 'r') as file2:
            lines2 = file2.readlines()

        # Modify the specific line
        with open(file_path2, 'w') as file2:
            for line in lines2:
                if line.strip().startswith('variable fz'):
                    file2.write(f'variable fz atom q*{e_intensity:.1f}\n')
                elif line.strip().startswith('variable efz'):
                    file2.write(f'variable efz atom q*z*{e_intensity:.1f}\n')
                else:
                    file2.write(line)

        lmp2 = lammps()
        lmp2.file("run-lammps_vacancy.sh")
        E_slab_vacancy = lmp2.get_thermo("pe") * conv

        E_vacancy = E_slab_vacancy - E_slab - D_e_O2 / 2
        E_vacancy_values.append(E_vacancy)

    except Exception as e:
        print(f'Error with electric field {e_intensity}: error {e}')
        continue

print(f'E_slab_values = {E_slab_values}')
print(f'E_vacancy_values = {E_vacancy_values}')
# Plotting
plt.plot(e_intensities, E_vacancy_values, 's-', label='E-vacancy')
plt.xlabel('Electric Field Intensity (eV/Ang)')
plt.ylabel('Energy (eV)')
plt.legend()
plt.savefig('energy_vacancy_all_efields.png')