
#Script to create scaled geometry files. 


import numpy as np
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from geometry_utils.scale_geometry import scale_lammps_geometry


start_scale = 0.95
end_scale = 1.70
step = 0.01


#input_file = "data.CeO2"
input_file = "simulations/energy_volume_scan_sim/data.CeO2"
output_file_template = "data"

for scale_factor in np.arange(start_scale, end_scale, step):
    scale_lammps_geometry(input_file, scale_factor, "simulations/energy_volume_scan_sim" ,  output_file_template)