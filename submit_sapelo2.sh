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

source activate snakemake

cd $SLURM_SUBMIT_DIR

WORKPATH="/scratch/rx32940/PMEN1"

time python start_analysis.py coreRecomb \
-p PMEN1.dated \
-o $WORKPATH -i $WORKPATH/assemblies \
-M \
-a $WORKPATH/PMEN1.dated.metadata.csv \
-m Year,Country \
-r $WORKPATH/GCF_000026665.1_ASM2666v1_genomic.fna \
-G " -f 30" -t 5 



conda deactivate
