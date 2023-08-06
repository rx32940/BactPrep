#BSUB -P BactPrep1
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=50GB]"
#BSUB -o BactPrep1.%J.out
#BSUB -n 4

source activate BactPrep

WORKPATH="/home/rxu28/random/BactPrep1"
REF="/home/rxu28/random/BactPrep1"

time start_analysis.py ALL \
-p BactPrep1 \
-o $WORKPATH -i $WORKPATH/assemblies \
-r $REF/GCF_000026665.1_ASM2666v1_genomic.fna \
-G " -f 30" -t 4 \
-M -a $WORKPATH/../PMEN1.dated.metadata.218.csv \
-s 1 \
-m Year,Country

conda deactivate
