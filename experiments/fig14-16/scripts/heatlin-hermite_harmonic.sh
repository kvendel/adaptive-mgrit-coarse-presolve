#!/bin/bash
#SBATCH --job-name=heatlin-IE
#SBATCH --partition=cpu
#SBATCH --nodes=8
#SBATCH --ntasks=1024
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=scripts/heatlin-hermite_harmonic-%j.out
#SBATCH --error=scripts/heatlin-hermite_harmonic-%j.err

set -euo pipefail

EXECUTABLE="./src/heatlin-hermite_harmonic"
EPS="1e-9"
NX="1025"
MGRIT_TOL="1e-9"
MAX_ITER="20"
RESULT_DIR="results_heatlin-hermite_harmonic"

# Egy sikertelen srun ennyiszer indulhat ujra.
MAX_ATTEMPTS=3
# Varakozas ket probalkozas kozott masodpercben.
RETRY_DELAY=30

P_VALUES=(1 2 4 8 16 32 64 128 256 512 1024)
R_VALUES=(2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000)
REPETITIONS=(0 1 2)

mkdir -p "$RESULT_DIR"

if [[ ! -x "$EXECUTABLE" ]]; then
   echo "A futtathato program nem talalhato vagy nem futtathato: $EXECUTABLE" >&2
   exit 1
fi

TEMP_FILE=""

cleanup()
{
   if [[ -n "${TEMP_FILE:-}" ]]; then
      rm -f "$TEMP_FILE"
   fi
}

trap cleanup EXIT INT TERM

for TEST in "${REPETITIONS[@]}"; do
   for P in "${P_VALUES[@]}"; do
      # P < 8 eseten minden MPI-rank kulon node-ra kerul.
      # P >= 8 eseten mind a 8 lefoglalt node-ot hasznaljuk,
      # a rankokat egyenletesen elosztva.
      if (( P < 8 )); then
         RUN_NODES="$P"
         TASKS_PER_NODE=1
      else
         RUN_NODES=8
         TASKS_PER_NODE=$((P / 8))
      fi

      for R in "${R_VALUES[@]}"; do
         EPS0=$(awk -v eps="$EPS" -v r="$R" 'BEGIN {printf "%.17g", eps*r*r}')
         RESULT_FILE=$(printf "%s/r_%04d_p_%04d_test_%d.dat" \
            "$RESULT_DIR" "$R" "$P" "$TEST")

         # A mar elkeszult, nem ures eredmenyek kihagyasa.
         if [[ -s "$RESULT_FILE" ]]; then
            continue
         fi

         TEMP_FILE="${RESULT_FILE}.tmp.${SLURM_JOB_ID}"
         attempt=1
         success=0

         while (( attempt <= MAX_ATTEMPTS )); do
            rm -f "$TEMP_FILE"

            echo "Inditas: TEST=$TEST P=$P R=$R, probalkozas=$attempt/$MAX_ATTEMPTS" >&2

            if srun \
               --nodes="$RUN_NODES" \
               --ntasks="$P" \
               --ntasks-per-node="$TASKS_PER_NODE" \
               --cpus-per-task=1 \
               --distribution=block:block \
               --cpu-bind=cores \
               --exclusive \
               "$EXECUTABLE" \
               -eps0 "$EPS0" \
               -eps "$EPS" \
               -nx "$NX" \
               -tol "$MGRIT_TOL" \
               -mi "$MAX_ITER" \
               > "$TEMP_FILE"
            then
               success=1
               break
            fi

            echo "Sikertelen srun: TEST=$TEST P=$P R=$R, probalkozas=$attempt/$MAX_ATTEMPTS" >&2
            rm -f "$TEMP_FILE"

            if (( attempt < MAX_ATTEMPTS )); then
               echo "Ujraprobalkozas ${RETRY_DELAY} masodperc mulva." >&2
               sleep "$RETRY_DELAY"
            fi

            attempt=$((attempt + 1))
         done

         if (( success == 0 )); then
            echo "A futas $MAX_ATTEMPTS probalkozas utan is sikertelen: TEST=$TEST P=$P R=$R" >&2
            exit 1
         fi

         mv "$TEMP_FILE" "$RESULT_FILE"
         TEMP_FILE=""
      done
   done
done
