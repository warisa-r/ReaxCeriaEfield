import numpy as np
import matplotlib.pyplot as plt

# Function to extract temperature from log file
def extract_temperature(log_file):
    timesteps = []
    temperatures = []
    start_parsing = False
    with open(log_file, 'r') as file:
        for line in file:
            if "Step" in line and "Temp" in line:
                start_parsing = True
                continue
            if start_parsing:
                if line.strip() == "":
                    break
                data = line.split()
                try:
                    timesteps.append(int(data[0]))
                    temperatures.append(float(data[1]))
                except (IndexError, ValueError):
                    continue
    return np.array(timesteps), np.array(temperatures)

# Extract temperature data from log file
timesteps, temperatures = extract_temperature("log.lammps")

# Function to plot temperature vs. time with specified interval
def plot_temperature(timesteps, temperatures, start=0, end=None):
    if end is None:
        end = timesteps[-1]
    mask = (timesteps >= start) & (timesteps <= end)
    plt.figure()
    plt.plot(timesteps[mask], temperatures[mask], label='Temperature')
    plt.xlabel('Timestep')
    plt.ylabel('Temperature (K)')
    plt.title(f'Temperature vs. Timestep ({start} to {end})')
    plt.legend()
    plt.savefig(f'temperature_plot_{start}_to_{end}.png')
    plt.show()

plot_temperature(timesteps, temperatures, start=0, end=300000)