# Bacterial_Genome_Analysis_Toolbox

This pipeline is written specifically for annotating the **bacteria whole genome sequences (WGS)**. The pipeline handles multiple operations that are necessary for bacterial genome analysis. Including:

1) Annotating bacterial WGS 
2) Constructing a pangenome for bacterial WGS dataset
3) Identify core and accessory loci for bacterial WGS dataset (both intra- and inter- species)
4) Produce core gene concatenation alignment
4) Identify potential recombination regions (recent & ancestral)
5) Identify SNPs from conserved regions of the bacterial genomes

## Installation

1) Install conda in your local computer

2) make a working directory
    ```
    mkdir {Bacterial_Analysis}*
    ```
    _* this name can change base on your project_

3) clone the repository into local working directory 
    ```
    git clone git@github.com:rx32940/Bacterial_Genome_Assembly_Pipeline.git
    ```

4) create conda env for Snakemake and activate the env
    ```
    conda create -n {snakemake}*
    conda activate {snakemake}*
    ```
    _* this name can change base on your project_

5) install snakemake
    ```
    conda install -c bioconda snakemake
    ```

6) after running all your analysis, deactivate the env
    ```
    conda deactivate
    ```

## Instruction

1) change the params in the config files, including:
    1. output_dir


