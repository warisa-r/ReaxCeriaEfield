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

volumes_per_CeO2 = []
energies_per_CeO2 = []

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
        print(f'Error at volume {i:.2f}')
        continue

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(volumes_per_CeO2, energies_per_CeO2, '-o', color='blue')
plt.xlabel('Volume per CeO2 (Å³)', fontsize=14)
plt.ylabel('Energy eV', fontsize=14)
plt.title('Energy vs. Volume for CeO2', fontsize=16)
plt.grid(True)
plt.savefig(f'{plot_directory}energy_vs_volume.png')