from ase.io import read, write
from pymatgen.io.ase import AseAtomsAdaptor

def trajtocif(traj_file, cif_file, frame=-1):
    '''
    convert one frame from trajectory into a cif file
    traj_file: input trajectory file name
    cif_file: output cif file name
    frame: desired frame, default is -1, last frame. ":" for all frame
    '''

    # Read the entire trajectory
    trajectory = read(traj_file, index=frame, format="lammps-dump-text")

    # Convert to pymatgen structure
    structure = AseAtomsAdaptor.get_structure(trajectory)

    # Manually adjust this part if needed
    atom_type_mapping = {
        1: 'Ce',
        2: 'O',
    }

    for site in structure:
        element = site.specie  # Get the specie at the site
        atomic_number = element.Z  # Get the atomic number

        # Map the atomic number to the desired symbol
        new_symbol = atom_type_mapping.get(atomic_number)
        if new_symbol:
            site.species = {new_symbol: 1.0}  # Set the new symbol with full occupancy
        else:
            raise ValueError(f"Unmapped atom type: {atomic_number}")

    # Write the last frame to a CIF file
    structure.to(fmt="cif", filename=cif_file)

