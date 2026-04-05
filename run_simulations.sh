#!/bin/bash
#SBATCH --job-name=robu_sims                            # Job name
#SBATCH --output=logs/sims_%j.out                       # Standard output log (%j = Job ID)
#SBATCH --error=logs/sims_%j.err                        # Standard error log
#SBATCH --partition=short                               # Short partition (max 12 hours on O2)
#SBATCH --time=12:00:00                                 # Time limit
#SBATCH --ntasks=1                                      # Number of tasks
#SBATCH --cpus-per-task=1                               # 1 core (script is sequential)
#SBATCH --mem=8G                                        # 8GB RAM is plenty for p=400
#SBATCH --mail-type=END,FAIL                            # Send email when job ends or fails
#SBATCH --mail-user=anthony-alexander_christidis@hms.harvard.edu
  
echo "Starting Main Simulations Job at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "Working directory: $SLURM_SUBMIT_DIR"

# Ensure we're in the submission directory (root of project)
cd $SLURM_SUBMIT_DIR

# Create necessary directories if they don't exist
mkdir -p logs
mkdir -p simulations/results

# Set thread limits for numerical libraries to prevent implicit threading clashes
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Load the gcc and R modules
module load gcc/14.2.0 R/4.4.2

# Execute the R script
echo "Starting R script at: $(date)"
Rscript simulations/run_simulations.R

echo "Job completed at: $(date)"