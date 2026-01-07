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
ml load Beast/1.10.4-GCC-12.3.0-CUDA-12.1.1
cd "${SLURM_SUBMIT_DIR}"

RUNS=(r1 r2 r3 r4 r5)
CLADES=(HG)
BURNIN=1000000
RESAMPLE=90000
TAG="Subsample1"

for cl in "${CLADES[@]}"; do
    mkdir -p "combined_rates_log/${cl}" \
             "combined_history_log/${cl}" \
             "combined_history_trees/${cl}"
done
mkdir -p bssvs_rates 

for run in "${RUNS[@]}"; do
    for cl in "${CLADES[@]}"; do
        cp "${run}/DTA22_combined_compact_subsampled_data1_AL.${cl}.rates.log" \
           "combined_rates_log/${cl}/${cl}_rates_${run}.log"

        cp "${run}/${cl}.history.log" \
           "combined_history_log/${cl}/${cl}_history_${run}.log"

        cp "${run}/${cl}.history.trees" \
           "combined_history_trees/${cl}/${cl}_history_${run}.trees"
    done

    cp "${run}/DTA22_combined_compact_subsampled_data1_AL.trees" \
       "bssvs_rates/bssvs_${run}.trees"
done

for cl in "${CLADES[@]}"; do
    logcombiner -burnin "${BURNIN}" \
        combined_rates_log/${cl}/${cl}_rates_r*.log \
        combined_rates_log/${cl}/${cl}_rates_${TAG}_combined.log

    logcombiner -burnin "${BURNIN}" \
        combined_history_log/${cl}/${cl}_history_r*.log \
        combined_history_log/${cl}/${cl}_history_${TAG}_combined.log
done

for cl in "${CLADES[@]}"; do
    logcombiner -trees -burnin "${BURNIN}" -resample "${RESAMPLE}" \
        combined_history_trees/${cl}/${cl}_history_r*.trees \
        combined_history_trees/${cl}/${cl}_history_${TAG}_combined.trees

    treeannotator -limit 0.95 -heights median \
        combined_history_trees/${cl}/${cl}_history_${TAG}_combined.trees \
        combined_history_trees/${cl}/${cl}_${TAG}_MCC.tree
done

echo "▶ Combining BSSVS trees"
logcombiner -trees -burnin "${BURNIN}" -resample "${RESAMPLE}" \
    bssvs_rates/bssvs_r*.trees \
    bssvs_rates/bssvs_rates_${TAG}_combined.trees

echo "▶ Annotating BSSVS MCC tree"
treeannotator -limit 0 -heights median \
    bssvs_rates/bssvs_rates_${TAG}_combined.trees \
    bssvs_rates/bssvs_rates_${TAG}.tree
