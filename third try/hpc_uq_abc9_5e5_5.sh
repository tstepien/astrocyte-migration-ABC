#!/bin/sh
#SBATCH --job-name=uq_abc9_5e5_5
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tstepien@ufl.edu # Where to send mail
#SBATCH --account=tstepien
#SBATCH --qos=tstepien-b
#SBATCH --nodes=1                  # Use one node
#SBATCH --ntasks=1                 # Run a single task
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2gb          # Memory per processor
#SBATCH --time=48:00:00             # Time limit hrs:min:sec
#SBATCH --output=logs/MATLAB_%j.txt    # Output and error log
pwd; hostname; date

module load matlab
./uq_abc9_5e5_5

date

echo "job successfully submitted"
exit
