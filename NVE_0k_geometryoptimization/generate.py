import random

num_atoms = 13
box_size = 10.0  # Box size

with open("data_13atoms.data", "w") as f:
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
