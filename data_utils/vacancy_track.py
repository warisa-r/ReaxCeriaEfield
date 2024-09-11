def read_nearest_cerium_indices(vacant_oxygen_index):
    """
    Reads the indices of the 4 nearest cerium atoms from the file corresponding to the given vacant oxygen index
    and increments each index by 1 since lammps index start from 1 and cif's start from 0.

    Parameters:
    vacant_oxygen_index (int): The index of the vacant oxygen atom.

    Returns:
    list: A list of indices of the 4 nearest cerium atoms in lammps index.
    """
    filename = f'vacancy_indices_{vacant_oxygen_index}.txt'
    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("Indices of the 4 nearest cerium atoms:"):
                # Extract the list of indices from the line
                indices_str = line.split(":")[1].strip()
                # Convert the string representation of the list to an actual list of integers
                indices = eval(indices_str)
                # Increment each index by 1
                incremented_indices = [index + 1 for index in indices]
                return incremented_indices
            
def parse_lammpstrj_of_atoms(filename, atom_indices):
    """
    Parses the LAMMPS trajectory file to extract the coordinates of specified atoms for every time frame.

    Parameters:
    filename (str): The name of the LAMMPS trajectory file.
    atom_indices (list): A list of atom indices to extract coordinates for.

    Returns:
    list: A list of dictionaries containing the coordinates of the specified atoms for each time frame.
    """
    coordinates = []
    with open(filename, 'r') as file:
        while True:
            line = file.readline()
            if not line:
                break
            if line.startswith("ITEM: TIMESTEP"):
                timestep = int(file.readline().strip())
                file.readline()  # Skip "ITEM: NUMBER OF ATOMS"
                num_atoms = int(file.readline().strip())
                file.readline()  # Skip "ITEM: BOX BOUNDS"
                file.readline()  # Skip box bounds line 1
                file.readline()  # Skip box bounds line 2
                file.readline()  # Skip box bounds line 3
                file.readline()  # Skip "ITEM: ATOMS id type x y z vx vy vz q"
                frame_data = {"timestep": timestep, "atoms": {}}
                for _ in range(num_atoms):
                    atom_data = file.readline().strip().split()
                    atom_id = int(atom_data[0])
                    if atom_id in atom_indices:
                        frame_data["atoms"][atom_id] = {
                            "type": int(atom_data[1]),
                            "x": float(atom_data[2]),
                            "y": float(atom_data[3]),
                            "z": float(atom_data[4]),
                            "vx": float(atom_data[5]),
                            "vy": float(atom_data[6]),
                            "vz": float(atom_data[7]),
                            "q": float(atom_data[8])
                        }
                coordinates.append(frame_data)
    return coordinates

def highlight_vacancy_site(filename, atom_indices, output_filename):
    """
    Modifies the atom types of specified atoms in the LAMMPS trajectory file. This is to see the vacancy site better.
    Though, I (Warisa am not particularly sure how useful or practical this can be)

    Parameters:
    filename (str): The name of the LAMMPS trajectory file.
    atom_indices (list): A list of atom indices to modify.
    output_filename (str): The name of the output file to write the modified data to. We create new file so that the initial
    lammpstrj can remain intact.
    """
    # Assert that the input and output files are LAMMPS trajectory files
    assert filename.endswith('.lammpstrj'), "Input file must be a LAMMPS trajectory file with extension .lammpstrj"
    assert output_filename.endswith('.lammpstrj'), "Output file must be a LAMMPS trajectory file with extension .lammpstrj"
    
    with open(filename, 'r') as infile, open(output_filename, 'w') as outfile:
        while True:
            line = infile.readline()
            if not line:
                break
            outfile.write(line)
            if line.startswith("ITEM: ATOMS"):
                while True:
                    atom_line = infile.readline()
                    if not atom_line or atom_line.startswith("ITEM:"):
                        outfile.write(atom_line)
                        break
                    atom_data = atom_line.strip().split()
                    atom_id = int(atom_data[0])
                    atom_type = int(atom_data[1])
                    if atom_id in atom_indices:
                        if atom_type == 1:
                            atom_data[1] = '3'
                        else:
                            raise ValueError(f"Atom surrounding Oxygen {atom_id} is not of Cerium. Wrong index?")
                    outfile.write(" ".join(atom_data) + "\n")