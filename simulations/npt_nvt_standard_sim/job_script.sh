#!/usr/bin/zsh 
#SBATCH --job-name=my_first_npt_job
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1

module load gcc
module load openmpi/4.1.4
module load CMake
module load mpi4py
module load Python/3.10

# Install matplotlib within the virtual environment

pip install mpi4py
pip install matplotlib

# Run the Python script
srun python run.py
