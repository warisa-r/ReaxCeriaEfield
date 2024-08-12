# compare the surface energy with pbe+u

# ---------- Initialize Simulation --------------------- 
clear 
units real
boundary p p p  
atom_style charge

# ---------- create atoms ------------------------------------
read_data 111slab.lmp


# ---------- Define Interatomic Potential ---------------------
pair_style	reax 10.0 0 1 1.0e-6
pair_coeff	* * ffield.reax 4 3

neighbor	2 bin
neigh_modify	every 10 delay 0 check no

dump    1 all atom 1 111dump.lammpstrj 

minimize        1.0e-6 1.0e-8 1000 10000