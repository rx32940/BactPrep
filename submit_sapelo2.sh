#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=sub.sh
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

cd $SLURM_SUBMIT_DIR

WORKPATH="/scratch/rx32940/testBactPrep"

time python start_analysis.py ALL \
-p PMEN1.dated \
-o $WORKPATH -i $WORKPATH/assemblies \
-r $WORKPATH/GCF_000026665.1_ASM2666v1_genomic.fna \
-G " -f 30" -t 5 



conda deactivate
