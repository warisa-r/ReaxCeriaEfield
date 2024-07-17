# REAX potential for CeO2 system
# .....
boundary       p p p
units		real

atom_style	charge
read_data data_dis.O2

pair_style reax/c NULL checkqeq yes

pair_coeff * * ffield.reax O

neighbor	2 bin
neigh_modify	every 10 delay 0 check no

fix 1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

#dump            1 all custom 10000 steps-ut.xyz id xu yu zu 
#dump            2 all custom 10000 steps-ut.CeO2 id xu yu zu q 

minimize        1.0e-6 1.0e-8 1000 10000