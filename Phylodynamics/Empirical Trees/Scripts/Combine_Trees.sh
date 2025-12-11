#!/bin/bash
#SBATCH --job-name=TreeCombiner
#SBATCH --partition=gpu_p
#SBATCH --gres=gpu:L4:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16  
#SBATCH --mem=10G
#SBATCH --time=7-00:00:00
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --mail-user=ss11645@uga.edu  
#SBATCH --mail-type=END,FAIL

cd $SLURM_SUBMIT_DIR

ml load Beast/1.10.4-GCC-12.3.0-CUDA-12.1.1
logcombiner -trees -burnin 10000000 -resample 50000 3.trees 4.trees 5.trees 6.trees 8.trees oldcombined_empirical_subsampled1.trees
