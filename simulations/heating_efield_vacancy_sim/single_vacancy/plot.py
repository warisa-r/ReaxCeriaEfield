import os
import sys
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..','..')))
from data_utils.extract_temperature import extract_temperature
from data_utils.plotting import single_plot, double_plot

x, y = extract_temperature("log.lammps")
#Reduce for better visualization
#x = x[:6000]
#y = y[:6000]
x_adjusted = x-x[0]

np.savez('single_vacancy.npz', x=x_adjusted, y=y)

single_plot(x_adjusted, y, 'Steps', 'Temperature', 'Single Vacancy', 'Singe_Vacancy.png')