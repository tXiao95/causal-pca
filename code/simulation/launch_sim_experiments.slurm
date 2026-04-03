#!/bin/bash

# Define the experiments
#NAME="final_results_nonsparse_SL_earth_nnet_100x50"   
NAME="final_results_nonsparse_nnet_128x64x32"   
EXPERIMENTS=("interaction")
#EXPERIMENTS=("interaction" "additive")
#EXPERIMENTS=("interaction" "additive" "weak_dim")
#EXPERIMENTS=("additive" "weak_dim")
#EXPERIMENTS=("weak_dim")

# Define the SLURM script to execute (update this if your slurm file is named differently)
SLURM_SCRIPT="fit_simulation.slurm"
COMPILE_SCRIPT="compile_results.slurm"

# Check if the SLURM script exists before running
if [ ! -f "$SLURM_SCRIPT" ]; then
    echo "Error: SLURM script '$SLURM_SCRIPT' not found!"
    exit 1
fi

# Check if the compile script exists before running
if [ ! -f "$COMPILE_SCRIPT" ]; then
    echo "Error: SLURM script '$COMPILE_SCRIPT' not found!"
    exit 1
fi

echo "Starting SLURM job submissions..."
echo "=========================================="

# Loop through each experiment and launch an array of 1000 jobs
for EXP in "${EXPERIMENTS[@]}"; do
    echo "Submitting experiment: $EXP"
    
    # 1. Submit the array job and capture its Job ID using --parsable
    ARRAY_JOB_ID=$(sbatch --parsable \
        --job-name="sim_${EXP}" \
        --array=1-1000 \
        --output="/scratch/thsiao3/output/sim_${EXP}_%A_%a.out" \
        --error="/scratch/thsiao3/error/sim_${EXP}_%A_%a.err" \
        "$SLURM_SCRIPT" "$NAME" "$EXP")
        
    echo "  -> Array Job Submitted (ID: $ARRAY_JOB_ID)"
    
    # 2. Submit the Compile job, dependent on the Array Job finishing
    COMP_JOB_ID=$(sbatch --parsable \
        --job-name="comp_${EXP}" \
        --dependency=afterany:${ARRAY_JOB_ID} \
        --output="/scratch/thsiao3/output/comp_${EXP}_%j.out" \
        --error="/scratch/thsiao3/error/comp_${EXP}_%j.err" \
        "$COMPILE_SCRIPT" "$NAME" "$EXP")
        
    echo "  -> Compile Job Submitted (ID: $COMP_JOB_ID) - Waiting on $ARRAY_JOB_ID"
    echo "------------------------------------------"
done

echo "All experiments and their compile jobs have been submitted successfully!"
