#!/usr/bin/bash
#SBATCH --job-name=chienlab-breseq             # Job name
#SBATCH --partition=cpu,cpu-preempt            # Partition (queue) name
#SBATCH -c 26                                  # Number of CPUs
#SBATCH --nodes=1                              # Number of nodes
#SBATCH --mem-per-cpu=4gb                      # Job memory request
#SBATCH --time=06:00:00                        # Time limit hrs:min:sec
#SBATCH --output=logs/chienlab-breseq_%j.log   # Standard output and error log

set -e

date;hostname;pwd


# Load modules
module load nextflow/24.04.3 conda/latest

# Run pipeline with slurm profile
nextflow run main.nf \
    --project "$proj_dir"  \
    -profile conda


