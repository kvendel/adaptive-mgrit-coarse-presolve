#!/bin/bash
#SBATCH --job-name=LV-prop_IE
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:03:00
#SBATCH --array=0-200%20
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null




# eps = 10^(-11 + i/50)
eps=$(awk -v i="$SLURM_ARRAY_TASK_ID" \
   'BEGIN {printf "%.16e", 10^(-9.0+i/50.0)}')

# eps0 = 100 * eps
eps0=$(awk -v eps="$eps" \
   'BEGIN {printf "%.16e", 100.0*eps}')

nt=$(awk -v eps="$eps" \
    'BEGIN {
        value = 189.0 * sqrt(1.0e-5 / eps)
        nt = int(value)
        if (nt < value) {
            nt++
        }
        printf "%d", nt
    }')
srun ./src/LV-original_IE  -ntime "$nt" -eps "$eps" -mc 8 -mi 30 -tol 1e-15
srun ./src/LV-article_IE  -ntime "$nt" -eps "$eps" -mc 8 -mi 30 -tol 1e-15
srun ./src/LV-IE -eps0 "$eps0" -eps "$eps" -mc 8 -mi 30 -tol 1e-9
srun -n 1 ./src/LV-sequential_IE -tol "$eps"
