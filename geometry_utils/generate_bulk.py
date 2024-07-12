from ase import Atoms
from ase.build import make_supercell
from ase.io import write


def create_fluorite(supercell_dims, geometry_filename):
    """
    Creates a fluorite structure of CeO2 and writes it to a file.

    This function generates a CeO2 unit cell with a fluorite structure, applies the specified supercell dimensions to it, and then writes the resulting structure to a file in the specified format.

    Parameters:
    - supercell_dims (tuple of int): The dimensions (nx, ny, nz) to expand the unit cell into a supercell.
    - geometry_filename (str): The filename (including path and extension) where the structure should be written. The format is inferred from the file extension.

    Returns:
    None
    """

    # Define lattice parameters for fluorite CeO2
    a_fluorite = 5.46745  # Angstrom
    c_fluorite = 5.46745  # Angstrom
    
    # Create the CeO2 unit cell with 12 atoms
    fluorite_CeO2 = Atoms(
        symbols='Ce4O8',
        positions=[
            [0.0, 0.0, 0.0],  # Atom positions for Ce (original from the structure i got from MP)
            [0.0, 2.733725, 2.733725],
            [2.733725, 0.0, 2.733725],
            [2.733725, 2.733725, 0.0],
            [1.3668625, 1.3668625, 4.1005875],
            [1.3668625, 4.1005875, 4.1005875],
            [1.3668625, 4.1005875, 1.3668625],
            [1.3668625, 1.3668625, 1.3668625],
            [4.1005875, 1.3668625, 1.3668625],
            [4.1005875, 4.1005875, 1.3668625],
            [4.1005875, 4.1005875, 4.1005875],
            [4.1005875, 1.3668625, 4.1005875],
        ],
        cell=[
            [a_fluorite, 0, 0],
            [0,a_fluorite, 0],
            [0, 0, c_fluorite]
        ],
        pbc=True
    )


    # Create the supercell of size 2x2x2
    fluorite_CeO2_supercell = make_supercell(fluorite_CeO2, [supercell_dims[0],0,0],[0,supercell_dims[1],0],[0,0,supercell_dims[2]])

    # Define the CIF filename or specify your path here if u want this to be organized
    cif_filename = geometry_filename # Example: "CeO2_fluorite.cif"
    # Write the supercell to a CIF file
    write(cif_filename, fluorite_CeO2_supercell, format='cif')

def create_rutile(supercell_dims, geometry_filename):
    """
    Creates a rutile structure of CeO2 and writes it to a file.

    This function generates a CeO2 unit cell with a rutile structure, applies the specified supercell dimensions to it, and then writes the resulting structure to a file in the specified format.

    Parameters:
    - supercell_dims (tuple of int): The dimensions (nx, ny, nz) to expand the unit cell into a supercell.
    - geometry_filename (str): The filename (including path and extension) where the structure should be written. The format is inferred from the file extension.

    Returns:
    None
    """

    # Define lattice parameters for fluorite CeO2
    a_rutile = 5.078  # Angstrom
    c_rutile = 3.567  # Angstrom
    
    # Create the CeO2 unit cell with 12 atoms
    rutile_CeO2 = Atoms(
        symbols='Ce2O4',
        positions=[
            # Original 6 atoms from Shengming
            [0.0000, 0.0000, 0.0000],                      # Ce1
            [2.5673285, 2.5673285, 1.8182100],             # Ce2
            [0.9619266, 4.1727304, 1.8182100],             # O1
            [4.1727304, 0.9619266, 1.8182100],             # O2
            [3.5292551, 3.5292551, 0.0000],                # O3
            [1.6054019, 1.6054019, 0.0000],                # O4
        ],
        cell=[
            [a_rutile, 0, 0],
            [0, a_rutile, 0],
            [0, 0, c_rutile]
        ],
        pbc=True
    )


    # Create the supercell of size 2x2x2
    rutile_CeO2_supercell = make_supercell(rutile_CeO2, [supercell_dims[0],0,0],[0,supercell_dims[1],0],[0,0,supercell_dims[2]])

    # Define the CIF filename or specify your path here if u want this to be organized
    cif_filename = geometry_filename # Example: "CeO2_rutile.cif"
    # Write the supercell to a CIF file
    write(cif_filename, rutile_CeO2_supercell, format='cif')