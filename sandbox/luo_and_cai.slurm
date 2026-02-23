#!/bin/bash
#SBATCH --job-name=luo_and_cai
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=thsiao3@emory.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --time=00:30:00
#SBATCH --partition=chang
#SBATCH --output=/scratch/%u/output/luo_and_cai_%j.log
#SBATCH --error=/scratch/%u/error/luo_and_cai_%j.log

# Load Singularity Module
module load singularity

# Execute R script using Singularity Image
singularity exec --bind /projects:/projects --bind /scratch:/scratch ~/singularity-images/geospatial_4.4.1.sif Rscript "2026-02-22_Reproducing_Luo_and_Cai_Table1.R" 
