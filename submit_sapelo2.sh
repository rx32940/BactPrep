#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=snakemake_intLineage4
#SBATCH --ntasks=1                    	
#SBATCH --cpus-per-task=10             
#SBATCH --time=150:00:00
#SBATCH --mem=100G
#SBATCH --output=/scratch/rx32940/snakemake_intLineage4.%j.out       
#SBATCH --error=/scratch/rx32940/snakemake_intLineage4.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL


cd /home/rx32940/github/Bacterial_Genome_Analysis_Toolbox/

source activate /home/rx32940/miniconda3/envs/Bacterial_Genome_Analysis_Toolbox 

snakemake --cores 10 --use-conda
snakemake --rulegraph| dot -Tpdf > dag.svg

conda deactivate