#!/bin/bash
#SBATCH --job-name=heatlin-art-orig
#SBATCH --partition=cpu
#SBATCH --nodes=8
#SBATCH --ntasks=1024
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=00:05:00
#SBATCH --output=heatlin-article-original-%j.out
#SBATCH --error=heatlin-article-original-%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

ARTICLE_EXECUTABLE="./heatlin-article"
ORIGINAL_EXECUTABLE="./heatlin-original"

NTIME="500"
EPS="1e-9"
NX="1025"
MGRIT_TOL="1e-9"
MAX_ITER="20"

ARTICLE_RESULT_DIR="results_heatlin-article"
ORIGINAL_RESULT_DIR="results_heatlin-original"

# Egy sikertelen srun ennyiszer indulhat ujra.
MAX_ATTEMPTS=3

P_VALUES=(1 2 4 8 16 32 64 128 256 512 1024)
REPETITIONS=(0 1 2)

mkdir -p "$ARTICLE_RESULT_DIR" "$ORIGINAL_RESULT_DIR"

for executable in "$ARTICLE_EXECUTABLE" "$ORIGINAL_EXECUTABLE"; do
   if [[ ! -x "$executable" ]]; then
      echo "A futtathato program nem talalhato vagy nem futtathato: $executable" >&2
      exit 1
   fi
done

TEMP_FILE=""

cleanup()
{
   if [[ -n "${TEMP_FILE:-}" ]]; then
      rm -f "$TEMP_FILE"
   fi
}

trap cleanup EXIT INT TERM

run_case()
{
   local executable="$1"
   local program_name="$2"
   local result_dir="$3"
   local test_id="$4"
   local process_count="$5"
   local run_nodes="$6"
   local tasks_per_node="$7"

   local result_file
   local attempt=1
   local success=0

   result_file=$(printf "%s/ntime_%04d_p_%04d_test_%d.dat" \
      "$result_dir" "$NTIME" "$process_count" "$test_id")

   # A mar elkeszult, nem ures eredmenyek kihagyasa.
   if [[ -s "$result_file" ]]; then
      echo "Kihagyva: PROGRAM=$program_name TEST=$test_id P=$process_count, az eredmeny mar letezik." >&2
      return 0
   fi

   TEMP_FILE="${result_file}.tmp.${SLURM_JOB_ID}"

   while (( attempt <= MAX_ATTEMPTS )); do
      rm -f "$TEMP_FILE"

      echo "Inditas: PROGRAM=$program_name TEST=$test_id P=$process_count NTIME=$NTIME, probalkozas=$attempt/$MAX_ATTEMPTS" >&2

      if srun \
         --nodes="$run_nodes" \
         --ntasks="$process_count" \
         --ntasks-per-node="$tasks_per_node" \
         --cpus-per-task=1 \
         --distribution=block:block \
         --cpu-bind=cores \
         --exclusive \
         "$executable" \
         -ntime "$NTIME" \
         -eps "$EPS" \
         -nx "$NX" \
         -tol "$MGRIT_TOL" \
         -mi "$MAX_ITER" \
         > "$TEMP_FILE"
      then
         success=1
         break
      fi

      echo "Sikertelen srun: PROGRAM=$program_name TEST=$test_id P=$process_count, probalkozas=$attempt/$MAX_ATTEMPTS" >&2
      rm -f "$TEMP_FILE"

      attempt=$((attempt + 1))
   done

   if (( success == 0 )); then
      echo "A futas $MAX_ATTEMPTS probalkozas utan is sikertelen: PROGRAM=$program_name TEST=$test_id P=$process_count" >&2
      return 1
   fi

   mv "$TEMP_FILE" "$result_file"
   TEMP_FILE=""
}

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

      # A ket program ugyanabban a SLURM-allokacioban, egymas utan fut.
      run_case "$ARTICLE_EXECUTABLE" \
         "heatlin-article" \
         "$ARTICLE_RESULT_DIR" \
         "$TEST" "$P" "$RUN_NODES" "$TASKS_PER_NODE"

      run_case "$ORIGINAL_EXECUTABLE" \
         "heatlin-original" \
         "$ORIGINAL_RESULT_DIR" \
         "$TEST" "$P" "$RUN_NODES" "$TASKS_PER_NODE"
   done
done
