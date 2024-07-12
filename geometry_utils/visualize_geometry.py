from pymatgen.io.lammps.data import LammpsData
from pymatgen.io.cif import CifWriter


def convert_lammps_to_cif(lammps_geometry_file, atom_style, output_file):
    """
    Converts a LAMMPS geometry file to a CIF file. When you use the add_charges_to_atoms function, you must remove # charges or # atomic from the LAMMPS file. 
    In order to use this function, the opposite has to be done. You must correct add back # atom_style (# atomic, # charge, # full) to the LAMMPS file after the keyword
    'Atoms'.

    Parameters:
    - lammps_geometry_file (str): Path to the LAMMPS geometry file.
    - atom_style (str): The atom style used in the LAMMPS file. Common styles include 'atomic', 'charge', 'molecular'.
    - output_file (str): Path for the output CIF file. Ensures the file ends with '.cif'.

    Returns:
    - None

    Example:
    ```
    # From the file: visualize_geometry.py
    # Read the LAMMPS file
    lammps_file = 'data_0.95.CeO2'

    # Specify the atom style if different from 'full'
    # Common styles: 'atomic', 'charge', 'molecular'
    atom_style = 'charge'  # Replace with your actual atom style

    # Read the LAMMPS file with the specified atom style
    output_file = 'data_0.95'

    convert_lammps_to_cif(lammps_file, atom_style, output_file)
    ```
    """

    # Ensure the output file ends with '.cif'
    if not output_file.endswith('.cif'):
        output_file += '.cif'

    # Read the LAMMPS file with the specified atom style
    lammps_data = LammpsData.from_file(lammps_geometry_file, atom_style=atom_style)

    # Extract the structure
    structure = lammps_data.structure

    # Write to CIF file
    cif_writer = CifWriter(structure)
    cif_writer.write_file(output_file)
    print(f"Conversion complete. CIF file saved as: {output_file}")
