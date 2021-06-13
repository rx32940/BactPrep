#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=sub.sh
#SBATCH --ntasks=1                    	
#SBATCH --cpus-per-task=10             
#SBATCH --time=150:00:00
#SBATCH --mem=30G
#SBATCH --output=../%x.%j.out       
#SBATCH --error=../%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

# ml Biopython/1.75-intel-2019b-Python-3.7.4
# ml snakemake/5.7.1-foss-2019b-Python-3.7.4
# ml Anaconda3/2020.02

# conda env create -f workflow/env/install.yaml -n snakemake

source activate snakemake

cd $SLURM_SUBMIT_DIR

GENPATH="/scratch/rx32940/pop.structure.lepto.int"
LINEAGEPATH="/scratch/rx32940/pop.structure.lepto.int/lineages/4.lineage"

python start_analysis.py ALL -p pop.struc.lineage.4 \
-t 10 \
-v $GENPATH/phage_region.bed \
-r $GENPATH/GCF_000092565.1_ASM9256v1_genomic.fna \
-i $LINEAGEPATH/assemblies 

conda deactivate
