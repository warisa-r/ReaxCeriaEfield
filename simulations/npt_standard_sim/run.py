from lammps import lammps
import matplotlib.pyplot as plt

# Initialize LAMMPS
lmp = lammps()

# Define your input script (load it from a file)
lmp.file("run-lammps.lmp")  # Ensure this script doesn't include a final "run" command

# Total number of timesteps to run
total_steps = 1000000
steps_per_run = 500

# Lists to store time (steps) and lattice constants
time_steps = []
lx_values = []
ly_values = []
lz_values = []

# Initialize old lattice constants with None
lx_old, ly_old, lz_old = None, None, None

# Run the simulation in chunks and collect lattice constants
for i in range(0, total_steps, steps_per_run):
    lmp.command(f"run {steps_per_run} pre no post no")  # Run in chunks of 500 steps
    
    # Get the box dimensions after each chunk of 500 steps
    lx_new = lmp.get_thermo("lx")
    ly_new = lmp.get_thermo("ly")
    lz_new = lmp.get_thermo("lz")

    # Store the data
    time_steps.append(i + steps_per_run)
    lx_values.append(lx_new)
    ly_values.append(ly_new)
    lz_values.append(lz_new)
    
    # Check if the change in lattice constants is less than 0.01
    if lx_old is not None and ly_old is not None and lz_old is not None:
        if abs(lx_new - lx_old) < 0.01 and abs(ly_new - ly_old) < 0.01 and abs(lz_new - lz_old) < 0.01:
            print("Convergence reached. Terminating the loop.")
            break
    
    # Update old lattice constants
    lx_old, ly_old, lz_old = lx_new, ly_new, lz_new

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
plt.savefig('plot_of_lattice_constant.png')