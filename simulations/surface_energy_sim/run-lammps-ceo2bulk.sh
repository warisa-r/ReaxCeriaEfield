# compare the surface energy with pbe+u

# ---------- Initialize Simulation --------------------- 
clear 
units real
boundary p p p 
atom_style charge

# ---------- create atoms ------------------------------------
read_data bulk.lmp

# ---------- Define Interatomic Potential ---------------------
pair_style	reax/c NULL
pair_coeff	* * ffield.reax Ce O

fix             1 all qeq/reax 1 0.0 10.0 1e-6 reaxff

neighbor 2.0 bin 
neigh_modify every 10 delay 0 check no

minimize        1.0e-6 1.0e-8 1000 10000