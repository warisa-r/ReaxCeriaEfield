from lammps import lammps
import numpy as np
lmp = lammps()
lmp.file("run-lammps.sh")

#Function to scale data file

import os

def scale_lammps_data(input_file, scale_factor, output_file_dir, output_file_template):
    # Read the input file
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Find the line indices for box dimensions
    for i, line in enumerate(lines):
        if 'xlo xhi' in line:
            xlo_xhi_index = i
        elif 'ylo yhi' in line:
            ylo_yhi_index = i
        elif 'zlo zhi' in line:
            zlo_zhi_index = i

    # Parse box dimensions
    xlo, xhi = map(float, lines[xlo_xhi_index].split()[:2])
    ylo, yhi = map(float, lines[ylo_yhi_index].split()[:2])
    zlo, zhi = map(float, lines[zlo_zhi_index].split()[:2])

    # Calculate new box dimensions
    new_xhi = xhi * scale_factor
    new_yhi = yhi * scale_factor
    new_zhi = zhi * scale_factor

    # Replace box dimensions in the lines
    lines[xlo_xhi_index] = f"{xlo} {new_xhi} xlo xhi\n"
    lines[ylo_yhi_index] = f"{ylo} {new_yhi} ylo yhi\n"
    lines[zlo_zhi_index] = f"{zlo} {new_zhi} zlo zhi\n"

    # Find the start of the atoms section
    atoms_start = lines.index('Atoms\n') + 2

    # Scale the atom coordinates
    for i in range(atoms_start, len(lines)):
        parts = lines[i].split()
        if len(parts) == 6:
            id, type, charge = parts[:3]
            x, y, z = map(float, parts[3:])

            new_x = x * scale_factor
            new_y = y * scale_factor
            new_z = z * scale_factor

            lines[i] = f"{id} {type} {charge} {new_x} {new_y} {new_z}\n"

    # Determine the output file name
    output_file_name = f"{output_file_template}_{scale_factor:.2f}.CeO2"
    output_file = os.path.join(output_file_dir, output_file_name)

    # Create the output directory if it doesn't exist
    os.makedirs(output_file_dir, exist_ok=True)

    # Write the scaled data to the output file
    with open(output_file, 'w') as file:
        file.writelines(lines)

#This function is used to rewrite the run-modernreax-lammps for it to read data from another file.

import re

def replace_and_save(input_filename, string2, output_filename):
    # Read the content from the input file
    with open(input_filename, 'r') as file:
        content = file.read()
    
    # Define the regular expression pattern to find the part after "read_data"
    pattern = r"(read_data\s+)\S+"
    
    # Replace the part after "read_data" with string2
    updated_content = re.sub(pattern, r"\1" + string2, content)
    
    # Write the updated content to the output file
    with open(output_filename, 'w') as file:
        file.write(updated_content)

# Constants needed for calculations
conv = 0.0433641 # eV to kcal/mol
nCeO2 = 32

start_scale = 0.95
end_scale = 1.25
step = 0.01

size = len(np.arange(start_scale, end_scale, step))

input_file = "data.CeO2"
output_file_dir = "NewAttempt_supercell"
output_file_template = "data"
simulation_file = "run-lammps.sh"

volume_array = np.zeros(size)
energy_array = np.zeros(size)
index = 0

for scale_factor in np.arange(start_scale, end_scale, step):
    lmp.command("clear")
    scale_lammps_data("CeO2_flourite_with_charges.lmp", scale_factor, "/mnt/d/git-workfile/PA/FFparamOpt/NewAttempt_supercell" ,  output_file_template)
    new_Data_file = f'data_{scale_factor:.2f}.CeO2'
    new_simulation_file= f"run-modernreax-lammps{scale_factor:.2f}.sh"
    replace_and_save(simulation_file, new_Data_file, new_simulation_file)
    lmp.file(new_simulation_file)
    volume = lmp.get_thermo("vol")
    volume_per_CeO2 = volume/ nCeO2
    energy_per_CeO2 = lmp.get_thermo("pe")* conv / nCeO2

    volume_array[index] = volume_per_CeO2
    energy_array[index] = energy_per_CeO2
    index +=1

energy_array -= energy_array.min()
import matplotlib.pyplot as plt
plt.plot(volume_array, energy_array, 'o-')
plt.xlabel('Volume (Å³)')
plt.ylabel('Energy (eV/CeO2)')
plt.savefig('VEscan.png')