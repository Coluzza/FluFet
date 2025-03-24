#!/bin/bash
#SBATCH --job-name=vel_analysis         # Job name
#SBATCH --output=Standard_%j.out        # Std output
#SBATCH --error=Standard_%j.err         # Std error
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks (single job)
#SBATCH --cpus-per-task=10              # Use 10 threads
#SBATCH --nodelist=Cogitatore1          # Force job to run on Cogitatore1
#SBATCH --time=0                         # Unlimited time

# Load environment
source ~/.bashrc  # Load environment variables

# Change to the working directory
cd $SLURM_SUBMIT_DIR

# Define input velocity file and number of threads
VELOCITY_FILE="vel_srd_f1t"  # Change to your actual file
NUM_THREADS=10

echo "Running velocity profile analysis on $(hostname) with $NUM_THREADS threads..."
python3 ~/bin/Flow_Profile_Parallel.py $VELOCITY_FILE $NUM_THREADS

echo "✅ Job completed on $(hostname)."
