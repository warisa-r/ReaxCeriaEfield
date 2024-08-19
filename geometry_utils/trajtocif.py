from ase.io import read, write

def trajtocif(traj_file, cif_file, frame=-1):
    '''
    convert one frame from trajectory into a cif file
    traj_file: input trajectory file name
    cif_file: output cif file name
    frame: to converted frame, default is the last frame
    '''

    # Read the entire trajectory
    trajectory = read(traj_file, index=':')

    # Extract the frame
    target_frame = trajectory[frame]

    # Manually adjust this part if needed
    atom_type_mapping = {
        1: 'Ce',
        2: 'O',
    }

    for atom in target_frame:
        old_type = atom.number  # or atom.symbol if atom types are initially symbols
        new_symbol = atom_type_mapping.get(old_type)
        if new_symbol:
            atom.symbol = new_symbol
        else:
            raise ValueError(f"Unmapped atom type: {old_type}")

    # Write the last frame to a CIF file
    write(cif_file, target_frame)