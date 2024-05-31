import pymatgen
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter
from ase.io import read, write
from ase.build import surface, make_supercell, sort

# load the structure cif file
#structure = Structure.from_file("geometry_file\CeO2_mp-20194_conventional_standard.cif")

# create a slabGeneratir object
#miller_index = (1, 1, 1)
#min_slab_size = 15
#min_vacuum_size = 15

#structure.make_supercell([4,4,1])
#slabgen = SlabGenerator(structure, miller_index, min_slab_size, min_vacuum_size,center_slab=True)

#slabs = slabgen.get_slab()
#print(len(slabs))
#slab = slabs[0]

# export slab to a cif file
#cif_writer = CifWriter(slabs)
#cif_writer.write_file("test_slab.cif")

# TODO: create supercell 4x4
#structure.make_supercell(scaling_matrix=[4,4,1])
#structure.to(fmt='cif', filename='test_slab.cif')

# with ase
structure = read("geometry_file\CeO2_mp-20194_conventional_standard.cif")

slab = surface(structure, (1,0,0), layers=8, vacuum=10.0)
slab.center(vacuum=15.0, axis=2)
#reconstruction

supercell = make_supercell(slab,[[4,0,0],[0,4,0],[0,0,1]])
write('100slab.cif', slab, format='cif')