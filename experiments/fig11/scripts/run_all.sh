#!/bin/bash
#SBATCH --job-name=step_sizes
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:05:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

srun ../src/lin-hermite_harmonic -eps0 1e-3 -eps 1e-5 -mc 8 -mi 30 -tol 1e-9
srun ../src/lin-linear -eps0 1e-3 -eps 1e-5 -mc 8 -mi 30 -tol 1e-9
srun -n 1 ../src/lin-sequential -tol 1e-5

srun ../src/lin-hermite_harmonic -eps0 1e-5 -eps 1e-7 -mc 8 -mi 30 -tol 1e-9
srun ../src/lin-linear -eps0 1e-5 -eps 1e-7 -mc 8 -mi 30 -tol 1e-9
srun -n 1 ../src/lin-sequential -tol 1e-7

srun ../src/LV-hermite_harmonic -eps0 1e-3 -eps 1e-5 -mc 64 -mi 30 -tol 1e-9
srun ../src/LV-linear -eps0 1e-3 -eps 1e-5 -mc 64 -mi 30 -tol 1e-9
srun -n 1 ../src/LV-sequential -tol 1e-5

srun ../src/LV-hermite_harmonic -eps0 1e-5 -eps 1e-7 -mc 64 -mi 30 -tol 1e-9
srun ../src/LV-linear -eps0 1e-5 -eps 1e-7 -mc 64 -mi 30 -tol 1e-9
srun -n 1 ../src/LV-sequential -tol 1e-7

srun ../src/vdp-hermite_harmonic -eps0 1e-3 -eps 1e-5 -mc 64 -mi 30 -tol 1e-9
srun ../src/vdp-linear -eps0 1e-3 -eps 1e-5 -mc 64 -mi 30 -tol 1e-9
srun -n 1 ../src/vdp-sequential -tol 1e-5

srun ../src/vdp-hermite_harmonic -eps0 1e-5 -eps 1e-7 -mc 64 -mi 30 -tol 1e-9
srun ../src/vdp-linear -eps0 1e-5 -eps 1e-7 -mc 64 -mi 30 -tol 1e-9
srun -n 1 ../src/vdp-sequential -tol 1e-7