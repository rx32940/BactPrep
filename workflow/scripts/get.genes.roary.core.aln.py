import os
from Bio import SeqIO
import sys

roary_aln_coreGenes = sys.argv[1]
roary_coreGenes_file = sys.argv[2]

with open(roary_coreGenes_file, "w") as w:
    print("Start writing core genes locus tag in Roary's core_gene_alignment.aln")


with open(roary_aln_coreGenes) as file:
    for line in file.readlines():
        if "locus_tag" in line:
                current_tag=line.split("=")[1]
                with open(roary_coreGenes_file, "a") as w:
                    w.write(current_tag)
