#!/bin/bash
#SBATCH --job-name=RDA_BW
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=thsiao3@emory.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100GB
#SBATCH --time=168:00:00
#SBATCH --partition=chang,week-long-cpu
#SBATCH --output=/scratch/%u/output/RDA_BW_%j.log
#SBATCH --error=/scratch/%u/error/RDA_BW_%j.log

# Load Singularity Module
module load singularity

# Execute R script using Singularity Image
singularity exec --bind /projects:/projects --bind /scratch:/scratch ~/singularity-images/geospatial_4.4.1.sif Rscript "ATL_real_data_analysis.R" $1
