#!/bin/bash
#SBATCH --job-name=r_error
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:01:00
#SBATCH --array=0-26%5
#SBATCH --output=r_error.out
#SBATCH --error=r_error.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

# Fix tolerancia
eps="1e-9"

# Az array 27 r-erteke
r_values=(2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000)

# Az aktualis array-elemhez tartozo r
r="${r_values[$SLURM_ARRAY_TASK_ID]}"

# eps0 = eps * r^2
eps0=$(awk -v eps="$eps" -v r="$r" \
    'BEGIN {printf "%.16e", eps*r*r}')

srun -n 4 ../src/lin-hermite_harmonic -eps0 "$eps0" -eps "$eps" -mc 8 -mi 30 -tol 1e-9
srun -n 4 ../src/LV-hermite_harmonic -eps0 "$eps0" -eps "$eps" -mc 64 -mi 30 -tol 1e-9
srun -n 4 ../src/vdp-hermite_harmonic -eps0 "$eps0" -eps "$eps" -mc 64 -mi 30 -tol 1e-9
