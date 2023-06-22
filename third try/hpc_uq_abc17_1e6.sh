#!/bin/sh
#SBATCH --job-name=uq_abc17_1e6
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tstepien@ufl.edu # Where to send mail
#SBATCH --account=tstepien
#SBATCH --qos=tstepien
#SBATCH --nodes=1                  # Use one node
#SBATCH --ntasks=1                 # Run a single task
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3gb          # Memory per processor
#SBATCH --time=24:00:00             # Time limit hrs:min:sec
#SBATCH --output=logs/MATLAB_%j.txt    # Output and error log
pwd; hostname; date

module load matlab
./uq_abc17_1e6

date

echo "job successfully submitted"
exit
