#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=sub.sh
#SBATCH --ntasks=1                    	
#SBATCH --cpus-per-task=10             
#SBATCH --time=150:00:00
#SBATCH --mem=100G
#SBATCH --output=../%x.%j.out       
#SBATCH --error=../%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL



cd $SLURM_SUBMIT_DIR

ml snakemake/5.7.1-foss-2019b-Python-3.7.4
>>>>>>> 79c4a25a2ac90d80e107b90c0acb4bc4e9013152

snakemake --cores 10 --use-conda
snakemake --rulegraph| dot -Tpdf > dag.pdf
