#!/bin/sh
#SBATCH --job-name=uq_abc18_1e4
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tstepien@ufl.edu # Where to send mail
#SBATCH --account=tstepien
#SBATCH --nodes=1                  # Use one node
#SBATCH --ntasks=1                 # Run a single task
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=48gb          # Memory per processor
#SBATCH --time=00:05:00             # Time limit hrs:min:sec
#SBATCH --output=logs/MATLAB_%j.txt    # Output and error log

module load matlab
uq/compile_abc18_1e4/uq_abc18_1e4

echo "job successfully submitted"
exit
