def add_charges_to_atoms(input_file, output_file):
    """
    Adds predefined charges to atoms in a LAMMPS geometry file generated by atomsk from a CIF file, and writes the modified atom information to an output file.

    This function is specifically designed to work with LAMMPS files generated by the command 'atomsk <input>.cif lammps', where <input>
    is the name of the CIF file without the '.cif' extension. Before using this function, ensure that the '#atomic' comment line in the
    LAMMPS file generated by atomsk is removed, as this function does not process files with '#atomic' present and the scaling function
    may not work correctly with it.

    The function assigns charges based on atom type, with the current implementation supporting charges for Ce (4.0) and O (-2.0) atom types.
    It reads the input LAMMPS file, identifies the atom information section following the 'Atoms' marker, and adds the predefined charges
    to the appropriate atoms based on their type. The modified atom information is then written to the specified output file.

    Warning: This function assumes a specific format for the input LAMMPS file where atom information starts immediately after a line
    containing 'Atoms'. The atom information is expected to be in the format:
    atom_index atom_type x_coordinate y_coordinate z_coordinate
    Lines that do not conform to this format may cause the function to skip them or terminate unexpectedly.

    Parameters:
    - input_file (str): The path to the input LAMMPS file containing atom information.
    - output_file (str): The path to the output file where the modified atom information will be written.

    Returns:
    None
    """

    # Define the charges for each atom type
    atom_charges = {
        1: 4.0,  # Ce
        2: -2.0  # O
    }

    # Read the file and store its contents
    with open(input_file, "r") as file:
        lines = file.readlines()

    # Find the line number where atom information starts
    atom_info_line = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            atom_info_line = i + 1
            break

    # If the line number is found, add charges to each atom
    if atom_info_line is not None:
        # Iterate through atom lines and add charges
        for i in range(atom_info_line, len(lines)):
            atom_line = lines[i].split()
            if len(atom_line) < 2:
                continue  # Skip lines that do not have enough elements
            atom_index = atom_line[0]
            atom_type = int(atom_line[1])
            charge = atom_charges.get(atom_type)
            if charge is not None:
                # Modify the line to include charge
                lines[i] = f"         {atom_line[0]}         {atom_line[1]}   {charge}        {atom_line[2]}       {atom_line[3]}       {atom_line[4]}\n"

    # Write the modified contents back to the output file
    with open(output_file, "w") as file:
        file.writelines(lines)