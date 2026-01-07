#!/bin/bash
#SBATCH --job-name=SkyGrid
#SBATCH --partition=gpu_p
#SBATCH --gres=gpu:L4:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16  
#SBATCH --mem=10G
#SBATCH --time=7-00:00:00
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --mail-user=xy12345@uga.edu  
#SBATCH --mail-type=END,FAIL

cd $SLURM_SUBMIT_DIR

ml load  Beast/1.10.4-GCC-12.3.0-CUDA-12.1.1
beast -threads 16 -beagle -beagle_GPU -overwrite emp_combined_compact_subsampled_data1_AL.xml