# REAX potential for CeO2 system
# .....
boundary       p p p
units		real

atom_style	charge
read_data data_vacancy.CeO2_111slab

pair_style	reax 10.0 0 1 1.0e-6

pair_coeff * * ffield.reax 4 3

neighbor	2 bin
neigh_modify	every 10 delay 0 check no


dump            1 all custom 10000 steps-ut-vacancy.lammpstrj id xu yu zu  

minimize        1.0e-6 1.0e-8 1000 10000