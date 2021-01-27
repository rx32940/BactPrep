#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=snakemake_pangenome
#SBATCH --ntasks=1                    	
#SBATCH --cpus-per-task=24             
#SBATCH --time=100:00:00
#SBATCH --mem=100G
#SBATCH --output=/scratch/rx32940/snakemake_prokka.%j.out       
#SBATCH --error=/scratch/rx32940/snakemake_prokka.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL


cd /home/rx32940/github/Bacterial_Genome_Analysis_Toolbox/

source activate /home/rx32940/miniconda3/envs/Bacterial_Genome_Analysis_Toolbox

snakemake --cores 24 --use-conda 

conda deactivate