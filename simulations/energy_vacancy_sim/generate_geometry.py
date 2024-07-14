# Prepare a slab geometry from the bulk geometry
from ase.io import read, write
from ase.build import surface, make_supercell

structure = read("bulk.cif") #TBD, download the bulk geometry from MP
slab = surface(structure, (1,1,1), layers=7, vacuum=40.0)
supercell = make_supercell(slab,[[1,0,0],[0,1,0],[0,0,1]])
write('111slab.cif', supercell, format='cif')