from lammps import lammps
import matplotlib.pyplot as plt
import statistics

lmp = lammps()
lmp.file("run-lammps-npt.lmp")

# Total number of timesteps to run for NPT
total_steps = 10000
steps_per_run = 100

# Lists to store time (steps) and lattice constants
time_steps = []
lx_values = []
ly_values = []
lz_values = []
T_values = []

# Initialize old lattice constants with None
lx_old, ly_old, lz_old, T_old = None, None, None, None

# Run the simulation in chunks and collect lattice constants
for i in range(0, total_steps, steps_per_run):
    lmp.command(f"run {steps_per_run} pre no post no")  # Run in chunks of 500 steps
    
    # Get the box dimensions after each chunk of 500 steps
    lx_new = lmp.get_thermo("lx")
    ly_new = lmp.get_thermo("ly")
    lz_new = lmp.get_thermo("lz")
    T_new = lmp.get_thermo("temp")

    # Store the data
    time_steps.append(i + steps_per_run)
    lx_values.append(lx_new)
    ly_values.append(ly_new)
    lz_values.append(lz_new)
    T_values.append(T_new)

    # Update old lattice constants
    lx_old, ly_old, lz_old, T_old = lx_new, ly_new, lz_new, T_new


# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(time_steps, lx_values, label='lx', marker='o')
plt.plot(time_steps, ly_values, label='ly', marker='s')
plt.plot(time_steps, lz_values, label='lz', marker='^')
plt.xlabel('Time Steps')
plt.ylabel('Lattice Constant (Angstroms)')
plt.title('Lattice Constant as a Function of Time')
plt.legend()
plt.grid(True)
# Save the plot as a picture file
plt.savefig('plot4.png')
# Save the simulation results
#lmp.command("write_restart restart.NPT")  # Save the restart file

'''T_mean = statistics.fmean(T_values)
plt.figure()
plt.plot(time_steps, T_values)
plt.xlabel('Time Steps')
plt.ylabel('Temperatur')
plt.title('Temperature as a Funktion of Time')
plt.grid(True)
plt.savefig('plot_T.png')
print(T_mean)'''