#!/bin/bash
#SBATCH --job-name=robu_app                             # Job name
#SBATCH --output=logs/app_%j.out                        # Standard output log
#SBATCH --error=logs/app_%j.err                         # Standard error log
#SBATCH --partition=short                               # Short partition
#SBATCH --time=1:30:00                                  # Time limit
#SBATCH --ntasks=1                                      # Number of tasks
#SBATCH --cpus-per-task=1                               # 1 core
#SBATCH --mem=8G                                        # 8GB RAM
#SBATCH --mail-type=END,FAIL                            # Send email when job ends or fails
#SBATCH --mail-user=anthony-alexander_christidis@hms.harvard.edu
  
echo "Starting Real Data Application Job at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"

# Ensure we're in the submission directory
cd $SLURM_SUBMIT_DIR

# Create necessary directories
mkdir -p logs
mkdir -p application/data
mkdir -p application/results

# Thread limits
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Load modules
module load gcc/14.2.0 R/4.4.2

# Execute the R script
echo "Starting R script at: $(date)"
Rscript application/run_application.R

echo "Job completed at: $(date)"