#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=sub.sh
#SBATCH --ntasks=1                    	
#SBATCH --cpus-per-task=10             
#SBATCH --time=150:00:00
#SBATCH --mem=10G
#SBATCH --output=../%x.%j.out       
#SBATCH --error=../%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL



cd $SLURM_SUBMIT_DIR

# ml Biopython/1.75-intel-2019b-Python-3.7.4
# ml snakemake/5.7.1-foss-2019b-Python-3.7.4
# ml Anaconda3/2020.02

# conda env create -f workflow/env/install.yaml -n snakemake

source activate snakemake

snakemake --cores 10 --use-conda
snakemake --rulegraph | dot -Tpdf > dag.pdf

conda deactivate
