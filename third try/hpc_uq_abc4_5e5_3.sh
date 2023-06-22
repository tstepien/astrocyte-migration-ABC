#!/bin/sh
#SBATCH --job-name=uq_abc4_5e5_3
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tstepien@ufl.edu # Where to send mail
#SBATCH --account=tstepien
#SBATCH --qos=tstepien-b
#SBATCH --nodes=1                  # Use one node
#SBATCH --ntasks=1                 # Run a single task
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2gb          # Memory per processor
#SBATCH --time=24:00:00             # Time limit hrs:min:sec
#SBATCH --output=logs/MATLAB_%j.txt    # Output and error log
pwd; hostname; date

module load matlab/2022b
./uq_abc4_5e5_3

date

echo "job successfully submitted"
exit
