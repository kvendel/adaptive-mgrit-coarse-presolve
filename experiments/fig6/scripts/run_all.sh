#!/bin/bash
#SBATCH --job-name=step_sizes
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:00:10
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null



srun ./src/lin-hermite_harmonic -eps0 1e-5 -eps 1e-7 -mc 8 -mi 30 -tol 1e-9
srun ./src/lin-germite_harmonic_cut -eps0 1e-5 -eps 1e-7 -mc 8 -mi 30 -tol 1e-9
