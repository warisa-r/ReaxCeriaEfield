# REAX potential for CeO2 system
# .....
boundary       p p p
units		real

atom_style	charge
read_data data.CeO2

pair_style	reax 10.0 0 1 1.0e-6
pair_coeff	* * ffield.reax 4 3

neighbor	2 bin
neigh_modify	every 10 delay 0 check no

velocity        all create 300.0 12345 mom yes rot yes
# fix npt used Nose-Hoover thermostat and barostat
fix             1 all npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0 # T start = T end = 300 K, P start = P stop = 1 atm

timestep        1.0
thermo          100
thermo_style    custom step temp press lx ly lz

dump           1 all custom 1000 dump.lammpstrj id type x y z vx vy vz q
#run             10000