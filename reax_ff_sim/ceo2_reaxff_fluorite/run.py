# it's basically the same, but with old lammps
# I tried to get the exactly same structure as in the paper, which is in geometry_file/new110slab
# would like to get 110 correct first
# surface energy is not crazy huge anymore, it's around 30-60 compare to paper around 1
# assuming the mathformel is correct, and the unit conversion also correct
# there might be some trick we don't know while creating structure
from lammps import lammps
import numpy as np

# convertion factor
conv1 = 4184 # kcal/mol to J/mol
mol = 6.02214076e23 # J/mol to J
conv2 = 1e-20 # A^2 to m^2
lmp_bulk = lammps()
lmp_bulk.file("ceo2fluorite.lmp")
E_bulk = lmp_bulk.get_thermo("pe")


def get111():
    nBulk_111 = 224 #14 for 111slab
    A_111 = 828.4211067 #51.7762 for 111slab  A^2
    lmp111 = lammps()
    lmp111.file("ceo2slab111.lmp")
    E_slab = lmp111.get_thermo("pe")
    E_surf_111 = (E_slab - nBulk_111 * E_bulk) * conv1 / (2 * A_111 * conv2 * mol) 
    #print(E_bulk)
    #print(E_slab)
    print(E_surf_111)

def get110():
    nBulk_110 = 273 #21/4 for 'new110slab', 312 for '110supercell', 273 for 'new110supercell'
    A_110 = 618.2739 #79.2658 for 'new110slab', 676.4022963 for '110supercell',  618.2739 for 'new110supercell'
    lmp110 = lammps()
    lmp110.file("ceo2slab110.lmp")
    E_slab = lmp110.get_thermo("pe")
    E_surf_110 = (E_slab - nBulk_110 * E_bulk) * conv1 / (2 * A_110 * conv2 * mol)
    #print(E_bulk)
    #print(E_slab)
    print(E_surf_110)

if __name__=="__main__":
    #get111()
    get110()
    pass