#!/bin/bash
#SBATCH --job-name=step_sizes
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:01:00
#SBATCH --output=scripts/run_all.out
#SBATCH --error=scripts/run_all.err

srun ./src/lin-original_IE -ntime 300 -eps 1e-7 -mc 16 -mi 30 -tol 1e-15
srun ./src/lin-original_BS -ntime 60 -eps 1e-9 -mc 16 -mi 30 -tol 1e-15
srun ./src/lin-article_IE -ntime 300 -eps 1e-7 -mc 16 -mi 30 -tol 1e-15
srun ./src/lin-article_BS -ntime 60 -eps 1e-9 -mc 16 -mi 30 -tol 1e-15
srun ./src/lin-IE -eps0 1e-5 -eps 1e-7 -mc 16 -mi 30 -tol 1e-15
srun ./src/lin-BS -eps0 1e-6 -eps 1e-9 -mc 16 -mi 30 -tol 1e-15
srun -n 1 ./src/lin-sequential_IE -tol 1e-7
srun -n 1 ./src/lin-sequential_BS -tol 1e-9

srun ./src/LV-original_IE -ntime 2000 -eps 1e-7 -mc 16 -mi 30 -tol 1e-15
srun ./src/LV-original_BS -ntime 300 -eps 1e-9 -mc 64 -mi 30 -tol 1e-15
srun ./src/LV-article_IE -ntime 2000 -eps 1e-7 -mc 16 -mi 30 -tol 1e-15
srun ./src/LV-article_BS -ntime 300 -eps 1e-9 -mc 64 -mi 30 -tol 1e-15
srun ./src/LV-IE -eps0 1e-5 -eps 1e-7 -mc 16 -mi 30 -tol 1e-15
srun ./src/LV-BS -eps0 1e-6 -eps 1e-9 -mc 64 -mi 30 -tol 1e-15
srun -n 1 ./src/LV-sequential_IE -tol 1e-7
srun -n 1 ./src/LV-sequential_BS -tol 1e-9

srun ./src/vdp-original_IE -ntime 5500 -eps 1e-7 -mc 128 -mi 30 -tol 1e-15
srun ./src/vdp-original_BS -ntime 800 -eps 1e-9 -mc 128 -mi 30 -tol 1e-15
srun ./src/vdp-article_IE -ntime 5500 -eps 1e-7 -mc 128 -mi 30 -tol 1e-15
srun ./src/vdp-article_BS -ntime 800 -eps 1e-9 -mc 128 -mi 30 -tol 1e-15
srun ./src/vdp-IE -eps0 1e-5 -eps 1e-7 -mc 128 -mi 30 -tol 1e-15
srun ./src/vdp-BS -eps0 1e-6 -eps 1e-9 -mc 128 -mi 30 -tol 1e-15
srun -n 1 ./src/vdp-sequential_IE -tol 1e-7
srun -n 1 ./src/vdp-sequential_BS -tol 1e-9
