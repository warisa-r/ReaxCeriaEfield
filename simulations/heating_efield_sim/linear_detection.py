import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from data_utils.linear_detection import find_optimal_window_size, sliding_window_linear_detection, plot_sliding_window_linear_detection
from data_utils.extract_temperature import extract_temperature

x_complete, y_complete = extract_temperature("log.lammps")
#Reduced from around 7000 to 6000 for better visualization
x = x_complete[:6000]
y = y_complete[:6000]
window_sizes = range(50, 2000, 50)
optimal_window_size = find_optimal_window_size(x, y, window_sizes, r2_threshold=0.95, slope_tolerance=0.1)
segments = sliding_window_linear_detection(x, y, window_size=optimal_window_size, r2_threshold=0.95, slope_tolerance=0.1)
plot_sliding_window_linear_detection(x, y, segments)

