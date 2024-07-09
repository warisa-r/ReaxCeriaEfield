import matplotlib.pyplot as plt
from lammps import lammps
import numpy as np

import os


# Define the path to the shell script
file_path = 'run-lammps.sh'

file_path2 = 'run-lammps_vacancy.sh'

conv = 0.0433641 # kcal/mol to eV
D_e_O2 = 4.876615876761253 # in eV Dissociation energy of O2 molecule fromm running oxygen_dis_energy.py

data_per_e_field = {}

e_intensities = np.arange(0, 150, 10)
E_slab_values = []
E_vacancy_values = []

for e_intensity in e_intensities:

    # change the magnitude of the electric field
    # Read the content of the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Modify the specific line
    with open(file_path, 'w') as file:
        for line in lines:
            if line.strip().startswith('variable fz'):
                file.write(f'variable fz atom q*{e_intensity:.1f}\n')
            elif line.strip().startswith('variable efz'):
                file.write(f'variable efz atom q*z*{e_intensity:.1f}\n')
            elif line.strip().startswith('dump'):
                file.write(f'dump            1 all custom 10000 steps-{e_intensity:.1f}.xyz id xu yu zu\n')
            else:
                file.write(line)
    
    lmp1 = lammps()
    lmp1.file("run-lammps.sh")
    E_slab = lmp1.get_thermo("pe") * conv
    E_slab_values.append(E_slab)

    with open(file_path2, 'r') as file:
        lines = file.readlines()

    # Modify the specific line
    with open(file_path2, 'w') as file:
        for line in lines:
            if line.strip().startswith('variable fz'):
                file.write(f'variable fz atom q*{e_intensity:.1f}\n')
            elif line.strip().startswith('variable efz'):
                file.write(f'variable efz atom q*z*{e_intensity:.1f}\n')
            elif line.strip().startswith('dump'):
                file.write(f'dump            1 all custom 10000 steps-vac-{e_intensity:.1f}.xyz id xu yu zu\n')
            else:
                file.write(line)

    lmp2 = lammps()
    lmp2.file("run-lammps_vacancy.sh")
    E_slab_vacancy = lmp2.get_thermo("pe") * conv

    E_vacancy = E_slab_vacancy - E_slab - D_e_O2 / 2
    E_vacancy_values.append(E_vacancy)


print(f'E_slab_values = {E_slab_values}')
print(f'E_vacancy_values = {E_vacancy_values}')
# Plotting
plt.plot(e_intensities, E_vacancy_values, 's-', label='E-vacancy')
plt.xlabel('Electric Field Intensity (eV/Ang)')
plt.ylabel('Energy (eV)')
plt.legend()
plt.savefig('energy_vacancy_all_efields.png')