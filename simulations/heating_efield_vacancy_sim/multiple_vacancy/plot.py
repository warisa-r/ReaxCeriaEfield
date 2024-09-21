import os
import sys
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..','..')))
from data_utils.extract_temperature import extract_temperature
from data_utils.plotting import single_plot, double_plot

x, y = extract_temperature("log.lammps")
#Reduce to avoid error due to inclomplete log file, reduce for better visulization
x = x[:15865]
y = y[:15865]
x_adjusted = x-x[0]

np.savez('multiple_vacancy.npz', x=x_adjusted, y=y)

single_plot(x_adjusted, y, 'Steps', 'Temperature', 'Multiple Vacancy', 'Multiple_Vacancy.png')