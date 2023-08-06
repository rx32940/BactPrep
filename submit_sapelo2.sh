#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=wgsRecomb
#SBATCH --ntasks=1                    	
#SBATCH --cpus-per-task=5      
#SBATCH --time=100:00:00
#SBATCH --mem=100G
#SBATCH --output=../%x.%j.out       
#SBATCH --error=../%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL


# ml Miniconda3/4.10.3
# conda env create -f workflow/env/install.yaml -n BactPrep python=3.7

source activate BactPrep


BACTPREP="/scratch/rx32940/whole_genome/analysis/Lint/BactPrep"
WORKPATH="/scratch/rx32940/whole_genome/analysis/Lint"
REF="/scratch/rx32940/whole_genome/ref"

cd $BACTPREP

time python start_analysis.py wgsRecomb \
-p Lint_recomb \
-o $WORKPATH -i $WORKPATH/assemblies \
-r $REF/GCF_000026665.1_ASM2666v1_genomic.fna \
-G " -f 30" -t 5 



conda deactivate
