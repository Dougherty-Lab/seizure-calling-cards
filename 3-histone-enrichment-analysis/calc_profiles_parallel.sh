#!/bin/bash
#SBATCH --job-name=histone_profiles_py
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8  # Request 8 CPUs (adjust based on your cluster)
#SBATCH --ntasks=1          # Keep as a single task with multiple CPUs

# Activate conda environment if necessary
source ~/miniconda3/etc/profile.d/conda.sh
conda activate histone # or another environment with the required packages

# Run the parallelized Python script
python calculate_profiles_parallel.py