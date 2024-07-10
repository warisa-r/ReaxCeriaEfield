from lammps import lammps
import numpy as np

conv1 = 4184 # kcal/mol to J/mol
mol = 6.02214076e23
conv2 = 1e-20 # A^2 to m^2
lmp_bulk = lammps()
lmp_bulk.file("ceo2fluorite.lmp")
E_bulk = lmp_bulk.get_thermo("pe")


def get111():
    nBulk_111 = 224 #14
    A_111 = 828.4211067 #51.7762   A^2
    lmp111 = lammps()
    lmp111.file("ceo2slab111.lmp")
    E_slab = lmp111.get_thermo("pe")
    E_surf_111 = (E_slab - nBulk_111 * E_bulk) * conv1 / (2 * A_111 * conv2 * mol) 
    #E_surf_111 = (E_slab - nBulk_111 * E_bulk) / (2 * A_111) 
    print(E_bulk)
    print(E_slab)
    print(E_surf_111)

def get110():
    nBulk_110 = 312 #21/4 312
    A_110 = 676.4022963 #79.2658 676.4022963 A^2
    lmp110 = lammps()
    lmp110.file("ceo2slab110.lmp")
    E_slab = lmp110.get_thermo("pe")
    E_surf_110 = (E_slab - nBulk_110 * E_bulk) * conv1 / (2 * A_110 * conv2 * mol)
    print(E_bulk)
    print(E_slab)
    print(E_surf_110)

if __name__=="__main__":
    #get111()
    get110()
    #print("hello world")