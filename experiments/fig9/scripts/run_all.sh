#!/bin/bash
#SBATCH --job-name=step_sizes
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:05:00
#SBATCH --output=scripts/run_all.out
#SBATCH --error=scripts/run_all.err

srun ./src/vdp-hermite_harmonic_I -eps0 1e-3 -eps 1e-5 -mc 128 -mi 30 -tol 1e-9
srun ./src/vdp-hermite_harmonic_exp -eps0 1e-3 -eps 1e-5 -mc 128 -mi 30 -tol 1e-9
srun ./src/vdp-hermite_harmonic -eps0 1e-3 -eps 1e-5 -mc 128 -mi 30 -tol 1e-9
srun ./src/vdp-hermite_harmonic_PID -eps0 1e-3 -eps 1e-5 -mc 128 -mi 30 -tol 1e-9
srun -n 1 ./src/vdp-sequential_I -tol 1e-5
srun -n 1 ./src/vdp-sequential_exp -tol 1e-5
srun -n 1 ./src/vdp-sequential -tol 1e-5
srun -n 1 ./src/vdp-sequential_PID -tol 1e-5

srun ./src/vdp-hermite_harmonic_I -eps0 1e-5 -eps 1e-7 -mc 128 -mi 30 -tol 1e-9
srun ./src/vdp-hermite_harmonic_exp -eps0 1e-5 -eps 1e-7 -mc 128 -mi 30 -tol 1e-9
srun ./src/vdp-hermite_harmonic -eps0 1e-5 -eps 1e-7 -mc 128 -mi 30 -tol 1e-9
srun ./src/vdp-hermite_harmonic_PID -eps0 1e-5 -eps 1e-7 -mc 128 -mi 30 -tol 1e-9
srun -n 1 ./src/vdp-sequential_I -tol 1e-7
srun -n 1 ./src/vdp-sequential_exp -tol 1e-7
srun -n 1 ./src/vdp-sequential -tol 1e-7
srun -n 1 ./src/vdp-sequential_PID -tol 1e-7

srun ./src/vdp-hermite_harmonic_I -eps0 1e-7 -eps 1e-9 -mc 128 -mi 30 -tol 1e-9
srun ./src/vdp-hermite_harmonic_exp -eps0 1e-7 -eps 1e-9 -mc 128 -mi 30 -tol 1e-9
srun ./src/vdp-hermite_harmonic -eps0 1e-7 -eps 1e-9 -mc 128 -mi 30 -tol 1e-9
srun ./src/vdp-hermite_harmonic_PID -eps0 1e-7 -eps 1e-9 -mc 128 -mi 30 -tol 1e-9
srun -n 1 ./src/vdp-sequential_I -tol 1e-9
srun -n 1 ./src/vdp-sequential_exp -tol 1e-9
srun -n 1 ./src/vdp-sequential -tol 1e-9
srun -n 1 ./src/vdp-sequential_PID -tol 1e-9
