#!/bin/bash
#SBATCH --job-name=vdp-prop_IE
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:05:00
#SBATCH --array=0-200%20
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# eps = 10^(-11 + i/50)
eps=$(awk -v i="$SLURM_ARRAY_TASK_ID" \
   'BEGIN {printf "%.16e", 10^(-9.0+i/50.0)}')

# eps0 = 100 * eps
eps0=$(awk -v eps="$eps" \
   'BEGIN {printf "%.16e", 100.0*eps}')

srun ./src/vdp-hermite_harmonic_I -eps0 "$eps0" -eps "$eps" -mc 128 -mi 30 -tol 1e-9
srun -n 1 ./src/vdp-sequential_I -tol "$eps"
srun ./src/vdp-hermite_harmonic_exp -eps0 "$eps0" -eps "$eps" -mc 128 -mi 30 -tol 1e-9
srun -n 1 ./src/vdp-sequential_exp -tol "$eps"
srun ./src/vdp-hermite_harmonic -eps0 "$eps0" -eps "$eps" -mc 128 -mi 30 -tol 1e-9
srun -n 1 ./src/vdp-sequential -tol "$eps"
srun ./src/vdp-hermite_harmonic_PID -eps0 "$eps0" -eps "$eps" -mc 128 -mi 30 -tol 1e-9
srun -n 1 ./src/vdp-sequential_PID -tol "$eps"
