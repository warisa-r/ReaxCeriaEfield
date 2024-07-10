# REAX potential for CeO2 system
# .....
boundary       p p p
units		real

atom_style	charge
read_data geometry_files/O2_dis.lmp
pair_style	reax 10.0 0 1 1.0e-6
pair_coeff	* * ffield.reax 2

neighbor	2 bin
neigh_modify	every 10 delay 0 check no

# Dump the charge (q) of each atom to a file
#dump 1 all custom 1 charges.dat id q

# Apply the force due to the electric field


# Calculate force components due to the electric field
variable fx atom 0.0
variable fy atom 0.0
variable fz atom q*1000.0
 
# Define the energy associated with the force (for minimization)
variable efz atom q*z*1000.0
 
# Apply the force due to the electric field
fix 1 all addforce v_fx v_fy v_fz energy v_efz

#dump            1 all custom 10000 steps-ut.xyz id xu yu zu 
#dump            2 all custom 10000 steps-ut.CeO2 id xu yu zu q 

minimize        1.0e-6 1.0e-8 1000 10000