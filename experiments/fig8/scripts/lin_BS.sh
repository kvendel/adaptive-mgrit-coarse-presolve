#!/bin/bash
#SBATCH --job-name=lin-prop_BS
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:01:00
#SBATCH --array=0-150%20
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# eps = 10^(-11 + i/50)
eps=$(awk -v i="$SLURM_ARRAY_TASK_ID" \
   'BEGIN {printf "%.16e", 10^(-11.0+i/50.0)}')

# eps0 = 100 * eps
eps0=$(awk -v eps="$eps" \
   'BEGIN {printf "%.16e", 1000.0*eps}')

nt=$(awk -v eps0="$eps0" \
    'BEGIN {
        value = 29.0 * sqrt(1.0e-5 / eps0)
        nt = int(value)
        if (nt < value) {
            nt++
        }
        printf "%d", nt
    }')
srun ../src/lin-original_BS  -ntime "$nt" -eps "$eps" -mc 8 -mi 30 -tol 1e-15
srun ../src/lin-article_BS  -ntime "$nt" -eps "$eps" -mc 8 -mi 30 -tol 1e-15
srun ../src/lin-BS -eps0 "$eps0" -eps "$eps" -mc 8 -mi 30 -tol 1e-9
srun -n 1 ../src/lin-sequential_BS -tol "$eps"
