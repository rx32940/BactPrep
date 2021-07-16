#!/bin/bash
#SBATCH --partition=bahl_salv_p
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

WORKPATH="/scratch/rx32940/PMEN1"

python start_analysis.py ALL \
-p PMEN1.dated \
-o $WORKPATH -g $WORKPATH/gff \
-r $WORKPATH/GCF_000026665.1_ASM2666v1_genomic.fna \
-R "-r -y -iv 1.5"
# -t 10 \
# -M \
# -a $WORKPATH/PMEN1.dated.metadata.csv \
# -m Year,Country


conda deactivate
