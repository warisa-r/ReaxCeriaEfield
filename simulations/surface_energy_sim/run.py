# it's basically the same, but with old lammps
# I tried to get the exactly same structure as in the paper, which is in geometry_file/new110slab
# would like to get 110 correct first
# surface energy is not crazy huge anymore, it's around 30-60 compare to paper around 1
# assuming the mathformel is correct, and the unit conversion also correct
# there might be some trick we don't know while creating structure
from lammps import lammps
import numpy as np

# convertion factor
conv = 6.94769546e-21
conv2 = 1e-20 # A^2 to m^2

lmp_bulk = lammps()
lmp_bulk.file("ceo2bulk.lmp")
E_bulk = lmp_bulk.get_thermo("pe")
E_bulk /= 4

def get111():
    nBulk_111 = 224 
    A_111 = 207.1048714 # A^2
    lmp111 = lammps()
    lmp111.file("ceo2slab111.lmp")
    E_slab = lmp111.get_thermo("pe")
    E_surf_111 = (E_slab - nBulk_111 * E_bulk) * conv / (2 * A_111 * conv2) 
    #print(E_bulk)
    #print(E_slab)
    print(E_surf_111)

def get110():
    nBulk_110 = 640 
    A_110 = 338.2008389 # A^2
    lmp110 = lammps()
    lmp110.file("ceo2slab110.lmp")
    E_slab = lmp110.get_thermo("pe")
    E_surf_110 = (E_slab - nBulk_110 * E_bulk) * conv / (2 * A_110 * conv2)
    #print(E_bulk)
    #print(E_slab)
    print(E_surf_110)

if __name__=="__main__":
    get111()
    get110()
    pass