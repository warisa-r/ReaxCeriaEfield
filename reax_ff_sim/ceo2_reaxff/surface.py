import pymatgen
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter

# load the structure cif file
structure = Structure.from_file("geometry_file\CeO2_mp-20194_conventional_standard.cif")

# create a slabGeneratir object
miller_index = (1, 1, 0)
min_slab_size = 15
min_vacuum_size = 15

slabgen = SlabGenerator(structure, miller_index, min_slab_size, min_vacuum_size)

slabs = slabgen.get_slab()

# export slab to a cif file
cif_writer = CifWriter(slabs)
cif_writer.write_file("110_slab.cif")