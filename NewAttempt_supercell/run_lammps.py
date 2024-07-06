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

conv = 0.0433641 # eV to kcal/mol
nCeO2 = 32 # Number of CeO2 in the simulation

data_per_e_field = {}

for e in np.arange(30, 100, 10):
    volumes_per_CeO2 = []
    energies_per_CeO2 = []

    # change the magnitude of the electric field
    # Read the content of the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Modify the specific line
    with open(file_path, 'w') as file:
        for line in lines:
            if line.strip().startswith('variable fz'):
                file.write(f'variable fz atom q*{e:.1f}\n')
            elif line.strip().startswith('variable efz'):
                file.write(f'variable efz atom q*z*{e:.1f}\n')
            else:
                file.write(line)
            

    for i in np.arange(0.95, 1.05, 0.01):
        # Read the content of the file
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Modify the specific line
        with open(file_path, 'w') as file:
            for line in lines:
                if line.strip().startswith('read_data'):
                    file.write(f'read_data data_{i:.2f}.CeO2\n')
                else:
                    file.write(line)

        lmp = lammps()
        try:
            lmp.file("run-lammps.sh")

            volume = lmp.get_thermo("vol")
            volume_per_CeO2 = volume / nCeO2
            volumes_per_CeO2.append(volume_per_CeO2)

            E_CeO2 = lmp.get_thermo("pe") * conv / nCeO2
            energies_per_CeO2.append(E_CeO2)
        except:
            print(f'Error with electric field {e} and volume {i:.2f}')
            continue
        # Store the data in the dictionary
        data_per_e_field[e] = (volumes_per_CeO2, energies_per_CeO2)

# Plotting all in one plot with different colors
for e, (volumes, energies) in data_per_e_field.items():
    plt.plot(volumes, energies, 'o-', label=f'E-field: {e} eV/Ang')

plt.xlabel('Volume per CeO2')
plt.ylabel('Energy per CeO2')
plt.title('Energy vs. Volume per CeO2 for different electric fields')
plt.legend()  # Display a legend to identify lines
plt.savefig('/mnt/c/Users/user/Documents/School/FFparamOpt/NewAttempt_supercell/energy_volume_with_efield/energy_vs_volume_all_efields.png')