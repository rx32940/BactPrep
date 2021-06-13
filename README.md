# Bacterial_Genome_Analysis_Toolbox

This pipeline is written specifically for annotating the **bacteria whole genome sequences (WGS)**. The pipeline handles multiple operations that are necessary for bacterial genome analysis. Including:

1) Annotating bacterial WGS 
2) Constructing a pangenome for bacterial WGS dataset
3) Identify core and accessory loci for bacterial WGS dataset (both intra- and inter- species)
4) Produce core gene concatenation alignment
4) Identify potential recombination regions (recent & ancestral)
5) Identify SNPs from conserved regions of the bacterial genomes

## Overall Workflow
<img width="723" alt="Screen Shot 2021-04-04 at 7 07 13 PM" src="https://user-images.githubusercontent.com/31255012/113523905-00af4500-9579-11eb-8264-cf4eab09491d.png">




## Installation

1) Install conda in your local computer

2) make a working directory
    ```
    mkdir {Bacterial_Analysis}*
    
    cd {Nacterial_Analysis}
    ```
    _* this name can change base on your project_

3) clone the repository into local working directory 
    ```
    git clone git@github.com:rx32940/Bacterial_Genome_Assembly_Pipeline.git
    ```

4) create conda env for Snakemake and activate the env
    ```
    conda env create -f workflow/env/install.yaml -n {snakemake}*
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
## Instruction #1 (Specify Options)

### Module Selection:
    
```Gubbins```: detect recombination from WGS alignment
    
```Roary```: construct bacteria pangenome
    
```fastGear```: will attempt to detect recombination for each gene in the  all genes in the  pangenome individually
                - predict recombinaiton among lineages detected by BAPs (can also provide your own lineage)
                - this module use gene loci detected by Roary, thus will also run module roary
                - please use fastGear_gene module for individual gene/alignment of interest

```fastGear_gene```: will detect recombination from a gene/alignment interested

```fastGear_core```: will dect recombinations only from the core genes detected by Roary
                - this is part of the ALL module
                - this module use gene loci detected by Roary, thus will also run module roary
                - will mask detected recombination region, and call SNPs from conserved region of core genome alignment
                - recombinations were detected for each gene individually
                - will also reconstruct phylogeny for the dataset based on the core clonal SNPs.


        
```ALL```: this module will attempt to run gubbins, roary, and fastGear_core module


```
python start_analysis.py -h
usage: start_analysis MODULE [options]

Please always specify the program to use in the first argument, or the whole
pipeline will attemp to run

positional arguments:
  {gubbins,roary,fastGear_core,fastGear,fastGear_gene,ALL}
                        Specify the module you would like to run

optional arguments:
  -h, --help            show this help message and exit

general arguments:
  -i , --input          path to input dir with assemlies
  -p , --name           provide name prefix for the output files
  -t , --thread         num of threads
  -o , --output         path to the output directory

arguments for if you would like to add metadata to output:
  -M, --addMetadata     must have the flag specify if want to allow annotation
  -a , --annotate       path to a csv file containing sample metadata
  -s , --sample         integer indicates which column the sample name is in
                        the metadata csv file
  -m , --metadata       metadata chosen to annotate ML tree/alignment after
                        the sample name

arguments for gubbins module:
  -r , --ref            reference (required for gubbins module)
  -v , --phage          phage region identified for masking (bed file)

arguments for roary module:
  -g , --gff            path to input dir with gff (this can replace input
                        assemblies dir in roary module Must be gff3 files)
  -c , --core           define core gene definition by percentage for roary
                        module (default=99)
  -k , --kingdom        specify the kingom of input assemlies for genome
                        annotation (default=Bacteria)
  -z , --genus          specify the genus of input assemlies for genome
                        annotation (default=Leptospira)

arguments for all three fastGear modules (fastGear_core, fastGear, fastGear_gene):
  --fg , --fastgear_param 
                        path to fastGear params
  --mcr_path            path to mcr runtime (need to install before use any of
                        the fastGear module
  --LD_LIBRARY_PATH     LD_LIBRARY_PATH, this can be find after mcr runtime
                        installation
  --fastgear-exe        path to the excutable of fastGear

arguments for fastGear_gene module:
  -n , --alignment      input alignment (either -n/-fl is required for
                        fastGear_gene module)
  -fl , --alnlist       input alignment list with path to gene alignments
                        (either -n/-fl is required for fastGear_gene module)

Enjoy the program! :)
```

**Run**:```python start_analysis.py ALL(roary/gubbins/fastGear)```

## Instruction #2 (Setting up config)

1) change the params in the ```config/config.yaml``` files, including:
    
    + Mandatory:
        1. **project_name** 
        2. **output_dir** (absolute path)
        3. **asm_dir** (not mandatory if prokka annotation is already available)
            + option 1:
                1) specify absolute path to the **prokka_dir** in the ```config/config.yaml``` file
            + option 2:
                1) create dir for genome annotation files (```.gff```, output for Prokka)
                    ```
                    cd {output_dir}
                    mkdir gff
                    ```
                2) copy ```.gff``` files you wish to analyze into the folder above
        
        4. **within_species** (TRUE/FALSE)
        5. **kingdom** and **genus**
    
    + Optional:
        1. if wish to detect recombination from core gene concatenation alignment
            - all resources below can be downloaded from [this link](https://users.ics.aalto.fi/~pemartti/fastGEAR/)
                - **fastGear_exe_path** (need to install [fastGear](https://mostowylab.com/news/fastgear))
                - **mcr_path** (need to install matlab complier runtime)
                - **LD_LIBRARY_PATH** (this will be provided by mcr installation after mcr is successfully installed, copy to config.yaml)
        2. if wish to detect recombination from whole genome sequences using [Gubbins](https://github.com/sanger-pathogens/gubbins) 
            - **reference** (absolute path)
            - **phage_region** (optional, absolute path to a bed file containing phage regions)
        3. threads_num (number of threads wish to use)

2) run the pipeline locally 
    
    1. make sure conda env is activated
    2. Run ```snakemake --cores {number of threads} --use-conda```

3) run the pipeline on cluster (example script for Sapelo2 cluster in UGA)

    1. edit ```submit_sapelo2.sh``` file
    2. submit the script ```sbatch submit_sapelo2.sh```

            



