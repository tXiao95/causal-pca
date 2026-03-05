#!/bin/bash

# Define the experiments
EXPERIMENTS=("baseline" "additive_confounding" "weak_dim_signal")

# Define the SLURM script to execute (update this if your slurm file is named differently)
SLURM_SCRIPT="fit_simulation.slurm"

# Check if the SLURM script exists before running
if [ ! -f "$SLURM_SCRIPT" ]; then
    echo "Error: SLURM script '$SLURM_SCRIPT' not found!"
    exit 1
fi

echo "Starting SLURM job submissions..."
echo "=========================================="

# Loop through each experiment and launch an array of 500 jobs
for EXP in "${EXPERIMENTS[@]}"; do
    echo "Submitting experiment: $EXP"
    
    sbatch \
        --job-name="sim_${EXP}" \
        --array=1-100 \
        --output="/scratch/thsiao3/output/sim_${EXP}_%A_%a.out" \
        --error="/scratch/thsiao3/error/sim_${EXP}_%A_%a.err" \
        "$SLURM_SCRIPT" "$EXP"
        
    echo "------------------------------------------"
done

echo "All experiments have been submitted successfully!"
