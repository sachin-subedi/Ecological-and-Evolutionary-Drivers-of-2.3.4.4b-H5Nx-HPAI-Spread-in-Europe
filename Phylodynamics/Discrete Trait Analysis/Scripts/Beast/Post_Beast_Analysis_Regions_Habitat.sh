#!/bin/bash
#SBATCH --job-name=GLM
#SBATCH --partition=bahl_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=7-00:00:00
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --mail-user=xy12345@uga.edu  
#SBATCH --mail-type=END,FAIL

set -euo pipefail
module purge
ml load  Beast/1.10.4-GCC-12.3.0-CUDA-12.1.1
cd "${SLURM_SUBMIT_DIR}"

RUNS=(r1 r2 r3 r4)
CLADES=(HG)
BURNIN=1000000
TAG="equal"

for cl in "${CLADES[@]}"; do
    mkdir -p "combined_rates_log/${cl}"
done

for run in "${RUNS[@]}"; do
    for cl in "${CLADES[@]}"; do
        cp "${run}/bssvs_rep1_equal_AL.${cl}.rates.log" \
           "combined_rates_log/${cl}/${cl}_rates_${run}.log"
    done
done

for cl in "${CLADES[@]}"; do
    echo "▶ Combining rates logs for ${cl}"
    logcombiner -burnin "${BURNIN}" \
        combined_rates_log/${cl}/${cl}_rates_r*.log \
        combined_rates_log/${cl}/${cl}_rates_${TAG}_combined.log
done

echo "✔ Done"