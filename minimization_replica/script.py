import subprocess
import random

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

# Define the content for input.min.lammps (without create_atoms line)
content = """\
# 1) Initialization
units            lj
dimension 3
atom_style atomic
pair_style lj/cut 30
boundary f f f

# 2) System definition
region simulation_box block -5 5 -5 5 -5 5
create_box 1 simulation_box

# 3) Simulation settings
mass 1 1
pair_coeff 1 1 1.0 1.0

# 4) Visualization
thermo 100
thermo_style custom step temp pe ke etotal press
dump mydmp all atom 1 dump.lammpstrj

# 5) Run
minimize 1.0e-10 1.0e-12 10000 10000
"""

# Split content into lines
content_lines = content.split('\n')

# Find the index where to insert create_atoms line
insert_index = None
for i, line in enumerate(content_lines):
    if line.startswith("create_box 1 simulation_box"):
        insert_index = i + 1
        break

# Run simulation 1000 times with different seed numbers
min_pot = 0
min_seed = 0
for i in range(1000):
    
    # Generate a random seed number with 6 digits
    seed_number = random.randint(100000, 999999)

    # Insert create_atoms line with the random seed number at the appropriate index
    content_with_seed = content_lines[:insert_index] + [f"create_atoms 1 random 13 {seed_number} simulation_box"] + content_lines[insert_index:]

    # Write content to input.min.lammps file
    with open("input.min.lammps", "w") as file:
        file.write('\n'.join(content_with_seed))

    # Call lmp_serial with the input file
    subprocess.run(["lmp_serial", "-in", "input.min.lammps"])

    # Extract potential energy from log file
    potential_energy = extract_potential_energy("log.lammps")

    # Print the potential energy for each iteration
    if potential_energy is not None:
        print(f"Iteration {i+1}: Potential Energy = {potential_energy}")
        if potential_energy < min_pot:
            min_pot = potential_energy
            min_seed = seed_number
            
    else:
        print(f"Iteration {i+1}: Potential energy not found in the log file.")


print(f"Potential Energy {min_pot} with seed {min_seed}")