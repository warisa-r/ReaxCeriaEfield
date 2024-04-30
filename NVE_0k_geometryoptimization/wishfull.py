import subprocess
import os
import random
import numpy as np

# Function to generate random initial positions
def generate_data_file(num_atoms, box_size, filename):
    with open(filename, "w") as f:
        f.write("LAMMPS data file for a simple Lennard-Jones system\n")
        f.write("\n")
        f.write(f"{num_atoms} atoms\n")
        f.write("1 atom types\n")
        f.write("\n")
        f.write(f"0.0 {box_size} xlo xhi\n")
        f.write(f"0.0 {box_size} ylo yhi\n")
        f.write(f"0.0 {box_size} zlo zhi\n")
        f.write("\n")
        f.write("Masses\n")
        f.write("\n")
        f.write("1 1.0\n")
        f.write("\n")
        f.write("Atoms\n")
        f.write("\n")
        for i in range(num_atoms):
            x = random.uniform(0.0, box_size)
            y = random.uniform(0.0, box_size)
            z = random.uniform(0.0, box_size)
            f.write(f"{i+1} 1 {x:.6f} {y:.6f} {z:.6f}\n")

def extract_potential_energy(filename):
    potential_energy = None
    in_paragraph = False

    with open(filename, "r") as f:
        for line in f:
            if line.startswith("Minimization stats:"):
                in_paragraph = True
                continue  # Skip to the next line

            if in_paragraph:
                if line.strip():  # Check if the line is not empty
                    if line.startswith("  Energy initial, next-to-last, final ="):
                        # Extract potential energy from the next line
                        next_line = next(f)
                        potential_energy = float(next_line.split()[-1])
                        break  # We found the value, no need to continue
                else:
                    # End of the paragraph
                    in_paragraph = False

    return potential_energy

# Usage
potential_energy = extract_potential_energy("log.lammps")
if potential_energy is not None:
    print("Potential energy:", potential_energy)
else:
    print("Potential energy not found in the file.")


# Function to run LAMMPS simulation and return potential energy
def run_lammps(filename):
    command = f"lmp_serial -in simple_NVE.in"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    # Check if LAMMPS run was successful
    if process.returncode == 0:
        potential_energy = extract_potential_energy("log.lammps")
        return potential_energy
    else:
        print("Error: LAMMPS simulation failed.")
        return None


        

# Main script
if __name__ == "__main__":
    num_simulations = 100
    num_atoms = 13
    box_size = 10.0
    data_files = []
    min_energy = 0

    for i in range(num_simulations):
        # Generate data file with random initial positions
        data_file = f"data_13atoms.data"
        generate_data_file(num_atoms, box_size, data_file)
        data_files.append(data_file)

        # Run LAMMPS simulation
        potential_energy = run_lammps("data_13atoms.data")
        if potential_energy  < min_energy:
            min_energy = potential_energy
            new_file_name = f"data_13atoms_iteration{i+1}.data"

    # Read data from the original file and write it to the new file
        with open("data_13atoms.data", "r") as original_file:
            data = original_file.read()
            with open(new_file_name, "w") as new_file:
                new_file.write(data)

        print(f"Simulation {i+1}: Potential energy = {potential_energy}")

    # Visualize simulation with lowest potential energy

