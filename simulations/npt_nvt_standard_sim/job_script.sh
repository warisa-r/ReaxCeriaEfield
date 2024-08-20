#!/usr/bin/zsh 
#SBATCH --job-name=my_first_npt_job
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

module load gcc
module load OpenMPI
module load CMake
module load Python/3.10

# Create a virtual environment in the job's temporary directory
python -m venv myenv

# Activate the virtual environment
source myenv/bin/activate

# Install matplotlib within the virtual environment
pip install matplotlib

# Run the Python script
srun python run.py

# Deactivate the virtual environment after the script has finished
deactivate