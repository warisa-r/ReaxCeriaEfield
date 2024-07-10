# REAX potential for CeO2 system
# .....
boundary       p p s
units		real

atom_style	charge
read_data geometry_files/111slab_with_charges.lmp

pair_style	reax 10.0 0 1 1.0e-6
pair_coeff	* * ffield.reax 2 1

neighbor	2 bin
neigh_modify	every 10 delay 0 check no

# Apply the electric field
fix 1 all efield 0.0 0.0 60.0


#dump            1 all custom 10000 steps-ut.xyz id xu yu zu 
#dump            2 all custom 10000 steps-ut.CeO2 id xu yu zu q 

minimize        1.0e-6 1.0e-8 1000 10000