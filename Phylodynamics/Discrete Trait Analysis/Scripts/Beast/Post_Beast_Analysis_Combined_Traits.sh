#!/bin/bash
#SBATCH --job-name=TreeCombiner
#SBATCH --partition=bahl_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G
#SBATCH --time=7-00:00:00
#SBATCH --output=combine.%j.out
#SBATCH --error=combine.%j.err
#SBATCH --mail-user=xy12345@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
module purge
ml Beast/1.10.4-GCC-11.3.0
cd "${SLURM_SUBMIT_DIR}"

BEAST_BIN="${EBROOTBEAST}/bin"
RUNS=(r1 r2 r3 r4 r5)
CLADES=(Habitat GeoCluster)
PREFIX="DTA11_rep1_equal_AL"
BURNIN=1000000
RESAMPLE=90000
TAG="equal"

for cl in "${CLADES[@]}"; do
    mkdir -p "combined_rates_log/${cl}" \
             "combined_history_log/${cl}" \
             "combined_history_trees/${cl}"
done
mkdir -p bssvs_rates

for run in "${RUNS[@]}"; do
    for cl in "${CLADES[@]}"; do
        cp "${run}/${PREFIX}.${cl}.rates.log" \
           "combined_rates_log/${cl}/${cl}_rates_${run}.log"
        cp "${run}/${cl}.history.log" \
           "combined_history_log/${cl}/${cl}_history_${run}.log"
        cp "${run}/${cl}_Subsample.history.trees" \
           "combined_history_trees/${cl}/${cl}_history_${run}.trees"
    done
    cp "${run}/${PREFIX}.trees" \
       "bssvs_rates/bssvs_${run}.trees"
done

for cl in "${CLADES[@]}"; do
    "${BEAST_BIN}/logcombiner" -burnin "${BURNIN}" \
        combined_rates_log/${cl}/${cl}_rates_*.log \
        combined_rates_log/${cl}/${cl}_rates_${TAG}_combined.log
    "${BEAST_BIN}/logcombiner" -burnin "${BURNIN}" \
        combined_history_log/${cl}/${cl}_history_*.log \
        combined_history_log/${cl}/${cl}_history_${TAG}_combined.log
done

for cl in "${CLADES[@]}"; do
    "${BEAST_BIN}/logcombiner" -trees -burnin "${BURNIN}" -resample "${RESAMPLE}" \
        combined_history_trees/${cl}/${cl}_history_*.trees \
        combined_history_trees/${cl}/${cl}_history_${TAG}_combined.trees
    "${BEAST_BIN}/treeannotator" -limit 0.95 -heights median \
        combined_history_trees/${cl}/${cl}_history_${TAG}_combined.trees \
        combined_history_trees/${cl}/${cl}_${TAG}_MCC.tree
done

echo "▶ Combining BSSVS trees"
"${BEAST_BIN}/logcombiner" -trees -burnin "${BURNIN}" -resample "${RESAMPLE}" \
    bssvs_rates/bssvs_*.trees \
    bssvs_rates/bssvs_rates_${TAG}_combined.trees

echo "▶ Annotating BSSVS MCC tree"
"${BEAST_BIN}/treeannotator" -limit 0 -heights median \
    bssvs_rates/bssvs_rates_${TAG}_combined.trees \
    bssvs_rates/bssvs_rates_${TAG}.tree