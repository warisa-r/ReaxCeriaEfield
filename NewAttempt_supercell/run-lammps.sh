# REAX potential for CeO2 system
# .....
boundary       p p s
units		real

atom_style	charge
read_data data_1.05.CeO2

pair_style	reax 10.0 0 1 1.0e-6
pair_coeff	* * ffield.reax 4 3

neighbor	2 bin
neigh_modify	every 10 delay 0 check no

# Calculate force components due to the electric field
variable fx atom 0.0
variable fy atom 0.0
variable fz atom q*90.0

# Define the energy associated with the force (for minimization)
variable efz atom q*z*90.0

# Dump the charge (q) of each atom to a file
#dump 1 all custom 1 charges.dat id q

# Apply the force due to the electric field
fix 1 all addforce v_fx v_fy v_fz energy v_efz

#dump            1 all custom 10000 steps-ut.xyz id xu yu zu 
#dump            2 all custom 10000 steps-ut.CeO2 id xu yu zu q 

minimize        1.0e-6 1.0e-8 1000 10000