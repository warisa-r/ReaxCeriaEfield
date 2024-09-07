import numpy as np

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