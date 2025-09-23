#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=48:00:00
#SBATCH --array=500
#SBATCH --partition=chang,day-long-cpu,week-long-cpu
#SBATCH --output=/scratch/thsiao3/output/fit_sim_%A_%a.out
#SBATCH --error=/scratch/thsiao3/error/fit_sim_%A_%a.err

# Load necessary modules or set any other environment variables if needed
module load singularity
# navigate to directory
cd /projects/hhchang/causal-pca

# Path to the Singularity image
SINGULARITY_IMAGE=~/singularity-images/geospatial_4.4.1.sif

# Path to the R script
R_SCRIPT=fit_sim_temp.R

# Launch Singularity Job
srun singularity exec --bind /projects:/projects --bind /scratch:/scratch $SINGULARITY_IMAGE Rscript $R_SCRIPT $1
