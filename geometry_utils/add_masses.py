def add_masses_to_lmp(input_file, output_file):
    '''
        add masses to certain lmp file, if the mass information is lost when transforming
        from cif to lmp
        The masses is predefined as Ce and O
        Parameters:
            - input_file (str): The path to the input LAMMPS file containing atom information.
            - output_file (str): The path to the output file where the modified mass information will be written.
    '''

    # Define atom masses for Ce and O
    atom_masses = {
        1 :  [140.11600000, '# Ce'],
        2 :  [15.99900000, '# O'],
    }

    # Read the input file
    with open(input_file, "r") as file:
        lines = file.readlines()

    # Find the insert position
    insert_index = 0
    for i, line in enumerate(lines):
        if line.startswith('Atom'):
            insert_index = i
            break

    # Create the mass section
    mass_section = "Masses\n\n"
    for i in range(len(atom_masses)):
        mass_section += f"      {i+1}    {atom_masses.get(i+1,0)[0]}       {atom_masses.get(i+1,0)[1]}\n"
    mass_section += "\n"
    
    # Insert mass section into the file
    new_lines = lines[:insert_index] + [mass_section] + lines[insert_index:]

    # write the modified file
    with open(output_file, "w") as file:
        file.writelines(new_lines)

