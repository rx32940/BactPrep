# BactPrep

This pipeline is written specifically for annotating the **bacteria whole genome sequences (WGS)**. The pipeline handles multiple operations that are necessary for bacterial genome analysis. Including:

1) Annotating bacterial WGS 
2) Constructing a pangenome for bacterial WGS dataset
3) Identify core and accessory loci for bacterial WGS dataset 
4) Produce core gene concatenation alignment (with & without recombination detection)
4) Identify potential recombination regions (recent & ancestral) - (WGS wise & Per-gene)
5) Identify SNPs from conserved regions of the bacterial genomes
6) Reconstruct Phylogeny of input dataset (Maximum Liklihood)
7) Add Annotation to alignment and ML trees taxa (designed for BEAST Analysis)

## Overall Workflow

![pipeline_workflow](https://user-images.githubusercontent.com/31255012/126013330-8ffed1fb-af59-45f2-9393-2cfe6331b324.png)


## Installation

1) Install conda (Python3) in your local computer or on the computing cluster

2) make a working directory
    ```
    mkdir {BactPrep_dir}*
    
    cd {BactPrep_dir}
    ```
    _* this name can change base on your project_

3) clone the repository into local working directory or download the most recent released issue from github
    ```
    git clone git@github.com:rx32940/BactPrep.git
    ```
    
    **or**
    
    download BactPrep release from Github

4) If first time using the pipeline
    ```
    cd BactPrep
    
    conda conda create -n {BactPrep}* conda-forge::mamba

    conda activate {BactPrep}*
    
    mamba install --file workflow/env/install.yaml

    bash INSTALL.sh

    ```
    _* this name can change base on your project_

4.1) if used the pipeline before or has matlab runtime R2016b (MCR) **AND** fastGear executable installed on the machine, use flag ```--mcr_path``` and ```--fastgear_exe``` to specify the absolute path to MCR and fasrGear executable. IF these two software were installed during previous use of BactPrep. you can find them in the ```resources``` folder from the previous download (please see example #6 below for detail).

    ```
    conda activate {BactPrep}*

    python start_analysis.py panRecomb -p all_Lint_fastGear_pan \
    -o $OUTPATH \
    -t 10 \
    -i $INPATH/assemblies \
    --mcr_path {path_to_previous_BactPrep_folder}/resources/mcr \
    --fastgear_exe {path_to_previous_BactPrep_folder}/resources/fastGEARpackageLinux64bit 
    ```
    
5) If you have trouble installing fastGear with ```INSTALL.sh``` script. please follow the instruction below for installation.
    1. mcr has many versions, use the link to download the version compatible with **fastGear**:
    
        Download and install fastGear excutable:
            1. change directory to: ```{absolute_path_to_BactPrep}/resources/mcr```
            2. you can download mcr provided by fastGear developers: https://users.ics.aalto.fi/~pemartti/fastGEAR/ 
            ```wget --no-check-certificate https://users.ics.aalto.fi/~pemartti/fastGEAR/fastGEARpackageLinux64bit.tar.gz -P {absolute_path_to_BactPrep}/resources```)
            3. Unzip the downloaded file
            ```tar -zvxf fastGEARpackageLinux64bit.tar.gz```

        Download and install Matlab Runtime:
        1. Download MCR zip provided by fastGear developers:
        ```wget https://users.ics.aalto.fi/~pemartti/fastGEAR/MCRInstallerLinux64bit.zip -P {absolute_path_to_BactPrep}/resources --no-check-certificate```
        2. Unzip the downloaded file ```unzip MCRInstallerLinux64bit.zip```
            1. or download version **R2016a** from MATLAB: https://www.mathworks.com/products/compiler/matlab-runtime.html
        3. change directory after unzip the downloaded file ```cd MCRInstallerLinux64bit```
        4. install: ```./install -destinationFolder {absolute_path_to_BactPrep}/resources/mcr/ -mode silent -agreeToLicense yes```
            1. if you would like to install with a GUI interface, please allow **X11 display** at the terminial, do ```./install```, this will open the GUI installation, and will allow you to change the directory to install, please install to ```{absolute_path_to_BactPrep}/resources/mcr``` 

    2. if you already have mcr (R2016a) on your machine (or used this pipeline before), you do not need to reinstall mcr, 
    please specify the absolute path with ```--mcr_path``` flag, which leads to the absolute path of your installed mcr 
        
        1. ```--mcr_path```: ex. ```--mcr_path {absolute_path_to_BactPrep}resources/mcr/```
        

6) You are now good to go! 
    RUN: ```python start_analysis.py ALL(coreGen/wgsRecomb/panRecomb)```

5) after running all your analysis, deactivate the env
    ```
    conda deactivate
    ```
## Instruction #1 (Specify Options)

### Module Selection:

```ALL```: this module will attempt to run ```wgsRecomb```, ```coreGen```, and ```coreRecomb``` module
            - All options required for these three modules are also required for ```ALL``` module
    
```wgsRecomb```: detect recombination from WGS alignment
    
```coreGen```: construct bacteria pangenome
    
```panRecomb```: will attempt to detect recombination for each gene in the  all genes in the  pangenome individually
                - predict recombinaiton among lineages detected by BAPs (can also provide your own lineage)
                - this module use gene loci detected by Roary, thus will also run module coreGen
                - please use geneRecomb module for individual gene/alignment of interest

```geneRecomb```: will detect recombination from a gene/alignment interested

```coreRecomb```: will dect recombinations only from the core genes detected by coreGen module (Roary)
                - this is part of the ALL module
                - this module use gene loci detected by Roary, thus will also run module coreGen
                - will mask detected recombination region, and call SNPs from conserved region of core genome alignment
                - recombinations were detected for each gene individually
                - will also reconstruct phylogeny for the dataset based on the core clonal SNPs.



```
usage: start_analysis MODULE [options]

Please always specify the program to use in the first argument, or the whole pipeline will attemp to run

positional arguments:
  {ALL,wgsRecomb,coreGen,coreRecomb,panRecomb,geneRecomb}
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
  -s , --sample         integer indicates which column the sample name is in the metadata csv file
  -m , --metadata       metadata chosen to annotate ML tree/alignment after the sample name

arguments for wgsRecomb module:
  -r , --ref            reference (required for wgsRecomb module)
  -v , --phage          phage region identified for masking (bed file)
  -G , --gubbins        any additional Gubbins arguments (please refer to Gubbins manual)

arguments for coreGen module:
  -g , --gff            path to input dir with gff (this can replace input assemblies dir in coreGen module Must be gff3 files)
  -c , --core           define core gene definition by percentage for coreGen module (default=99)
  -k , --kingdom        specify the kingom of input assemlies for genome annotation (default=Bacteria)
  -R , --roary          any additional roary arguments (please refer to Roary manual)

arguments for all three fastGear modules (coreRecomb, panRecomb, geneRecomb):
  --mcr_path            path to mcr runtime (need to install before use any of the fastGear module
  --fastgear_exe        path to the excutable of fastGear
  --fg , --fastgear_param 
                        path to fastGear params

arguments for geneRecomb module:
  -n , --alignment      input alignment (either -n/-fl is required for geneRecomb module)
  -fl , --alnlist       input alignment list with path to gene alignments (either -n/-fl is required for geneRecomb module)

Enjoy the program! :)
```

**Run**: ```python start_analysis.py ALL(coreGen/wgsRecomb/panRecomb)```

## Output Files

### ALL module

**OUTPUT:**
```
    ├── fastgear_core (coreRecomb output)
        ├── core_loci_fastGear_out (FastGear recombinaition analysis for each core gene)
        ├── core_loci_list.txt (list of genes analyzed in coreRecomb module)
        ├── fastgear_iqtree (ML phylogeny reconstructed using recombinaition masked core SNPs)
        ├── fastGear_masked_coregene_aln.fasta (core gene concatenation aln w/recombinaition masked)
        ├── masked_coregene_aln (recombination masked gene alns - individual genes)
        ├── plot_coregenome (combined statistics and final plots for fastGear core)
        ├── PMEN1.dated_core_mask_snp.fasta (recomb-masked core SNPs alns used for fastGear_iqtree)
        ├── PMEN1.dated_core_mask_snp_meta.fasta (recomb-masked core SNPs alns with annotation added)
        └── roary_coreGeneAln_locustag.txt (same file as core_loci_list.txt)
    ├── gff (Prokka annotated gff)
    ├── gubbins (gubbins output - detailed explanantion for each file can be find in Gubbins manual)
        ├── iqtree (ML phylogeny reconstructed using recomb-Free SNPs)
        ├── PMEN1.dated.branch_base_reconstruction.embl
        ├── PMEN1.dated.filtered_polymorphic_sites.fasta
        ├── PMEN1.dated.filtered_polymorphic_sites.phylip
        ├── PMEN1.dated.final_tree.tre
        ├── PMEN1.dated_meta.recombFreeSnpsAtcg.fasta (taxa annotated alns for Gubbins detected recomb-free SNPs)
        ├── PMEN1.dated.node_labelled.final_tree.tre
        ├── PMEN1.dated_noref.filtered_polymorphic_sites.fasta
        ├── PMEN1.dated.per_branch_statistics.csv
        ├── PMEN1.dated.recombination_predictions.embl
        ├── PMEN1.dated.recombination_predictions.gff
        ├── PMEN1.dated.summary_of_snp_distribution.vcf
        └── snp-sites 
    ├── prokka (Prokka output - individual WGS assemblies's genome annotation)
    ├── roary (roary output)
        ├── accessory_binary_genes.fa
        ├── accessory_binary_genes.fa.newick
        ├── _accessory_clusters
        ├── _accessory_clusters.clstr
        ├── accessory_graph.dot
        ├── accessory.header.embl
        ├── accessory.tab
        ├── blast_identity_frequency.Rtab
        ├── _blast_results
        ├── _clustered
        ├── _clustered.clstr
        ├── clustered_proteins
        ├── _combined_files
        ├── _combined_files.groups
        ├── core_accessory_graph.dot
        ├── core_accessory.header.embl
        ├── core_accessory.tab
        ├── core_alignment_header.embl
        ├── core_gene_alignment.aln
        ├── gene_presence_absence.csv
        ├── gene_presence_absence.Rtab
        ├── _inflated_mcl_groups
        ├── _inflated_unsplit_mcl_groups
        ├── _labeled_mcl_groups
        ├── number_of_conserved_genes.Rtab
        ├── number_of_genes_in_pan_genome.Rtab
        ├── number_of_new_genes.Rtab
        ├── number_of_unique_genes.Rtab
        ├── pan_genome_reference.fa
        ├── pan_genome_sequences
        ├── PMEN1.dated_coreConcate_meta.fasta
        ├── roary_iqtree (ML phylogeny reconstructed with core gene concatenation alignment)
        ├── summary_statistics.txt
        └── _uninflated_mcl_groups
    └── snippy (snippy output - individual samples' snp calling with the reference genome)
```

## Examples

**Sample Dataset**: [218 Streptococcus pneumoniae  PMEN1 WGS assemblies collected from the year 1984 - 2008 from 22 unique countries globally]((https://zenodo.org/record/5603335#.YkxP9y9h1TY)
    - The sample dataset can be downloaded to your work directory by:
```
mkdir -p $INPATH/assemblies

cd $INPATH/assemblies

zenodo_get -d 10.5281/zenodo.5603335
```

**Reference genome** for _Streptococcus pneumoniae_ PMEN1 cab be downloaded from NCBI: [Streptococcus pneumoniae ATCC 700669 (firmicutes)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000026665.1/)
```
cd $INPATH/

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/665/GCF_000026665.1_ASM2666v1/GCF_000026665.1_ASM2666v1_genomic.fna.gz

gunzip GCF_000026665.1_ASM2666v1_genomic.fna.gz
```

**1) Get Start - How to run ALL Module**
if you would like to run "wgsRecomb", "coreGen", and "coreRecomb" modules all together, you can just use the "ALL" module. **Note: a reference genome (-r) is necessary to run "wgsRecomb" module**

EXAMPLE:
```
python start_analysis.py ALL -p PMEN1.dated \
-o $OUTPATH \
-i $INPATH/assemblies
-r $INPATH/GCF_000026665.1_ASM2666v1_genomic.fna
```

1.1) If you already have gff files obtained from previous analysis, **gff dir** can be used as input for "coreGen module". This will saves a lot time

```
python start_analysis.py ALL -p PMEN1.dated \
-o $OUTPATH \
-t 10 \
-g $INPATH/gff \ # gff dir as input
--mcr_path {path_to_previous_BactPrep_folder}/resources/mcr \ # absolute path to mcr R2016b
--fastgear_exe /home/user/SOFTWARE/fastGEARpackageLinux64bit # absolute path to fastGear excutable

```

**2) Obtain Annotated Outputs**
if you would like to obtain annotated phylogenies and alignments, please provide a CSV file with annotation of every isolates. Flag ```-M``` must be specified for annotation. ```-a``` is the path to the CSV metadata file. ```-s``` allows the you to specify the index of the column **matches with the input assemblies' file names**, default is 1. ```-m``` asks for the column names of the metadata you would like to add for annotations (comma separated). 

EXAMPLE CSV File:
ENA Accession | Strain | Year | Country
-- | -- | -- | --
ERS009226 | ARG 740 | 1995 | Argentina
ERS009778 | 3122 | 1994 | Canada
ERS009785 | 36148 | 2008 | Canada
ERS004773 | HK P1 | 2000 | China
ERS004775 | HK P38 | 2000 | China

EXAMPLE:
```
python start_analysis.py ALL -p PMEN1.dated \
-o $OUTPATH \
-i $INPATH/assemblies \
-r $INPATH/GCF_000026665.1_ASM2666v1_genomic.fna \
-M \
-a $INPATH/PMEN1.dated.metadata.csv
-s 1
-m Year,Country
```
**3) IF you would only like to run "wgsRecomb"**
please keep in mind, a reference genome must provided by user."wgsRecomb" module will call snps from the reference genome for each input WGS assemblies, and combine them into a multiple sequence alignment using **Snippy**, where genome regions shared by all isolates in the input dataset will be extracted. **Gubbins** will take the Snippy input to detect recombination regions from the multi-sequence alignment. At end of the pipeline, SNPs outside of the recombination regions will be used to reconstruct the input dataset's phylogeny with **IQTree**. Annotation will be added to phylogenies and SNPs alignments if a metadata file is provided by user (please see example 2 for details).

EXAMPLE:
```
python start_analysis.py wgsRecomb -p PMEN1.dated \
-o $OUTPATH \
-i $INPATH/assemblies \
-r $INPATH/GCF_000026665.1_ASM2666v1_genomic.fna 

```

**4) IF you would only like to run "coreGen"**
a reference file would be required for this module. All input WGS assemblies will be annotated by **Prokka**. Using prokka's gene annotations, **Roary** will 1) reconstruct the pangenome of the input dataset, and 2) identify core genes shared by 99% (this can be adjust by user by ```-c``` flag) of the isolates in the input dataset. Roary will also provide a core gene concatenation alignment, which will be used for phylogeny reconstruction using **IQTree** at end of the pipeline.  Annotation will be added to phylogenies and SNPs alignments if a metadata file is provided by user (please see example 2 for details).

EXAMPLE:
```
python start_analysis.py wgsRecomb -p PMEN1.dated \
-o $OUTPATH \
-i $INPATH/assemblies 

```

**5) IF you would only like to run "coreRecomb"**
pipeline implemented in the "coreGen" module will run first. "coreRecomb" module will identify homologous recombinaition from every core gene identified by **Roary**. The identified recombination regions will be masked in the gene alignments before all core genes' masked alignments are concatenated into a super-gene alignment. core SNPs outside of the recombination regions will be called, SNPs outside of the recombination regions will be used to reconstruct the input dataset's phylogeny with **IQTree**. Annotation will be added to phylogenies and SNPs alignments if a metadata file is provided by user (please see example 2 for details).

EXAMPLE:
```
python start_analysis.py coreRecomb \
-p PMEN1.dated \
-o $WORKPATH -i $WORKPATH/assemblies \
-r $WORKPATH/GCF_000026665.1_ASM2666v1_genomic.fna \
-t 10 \
-M \
-a $WORKPATH/PMEN1.dated.metadata.csv \
-m Year,Country 

```

**6) IF matlab runtime (MCR) version R2016a is installed or this is not the first time you are running this pipeline.**
if you have already installed MCR R2016a and fastGear executable before on your machine, or you have already installed these two dependencies the previous times you were using BactPrep. You can use flag ```--mcr_path``` and ```--fastgear_exe``` to avoid installing these two dependencies again. **you don't need to run ```INSTALL.sh``` script again if these two scripts is already installed, but a conda env still need to be created and activated to run BactPrep pipeline**

EXAMPLE:

```
conda env create -f workflow/env/install.yaml -n BactPrep

conda activate BactPrep

python start_analysis.py panRecomb -p PMEN1.dated_fastGear_pan \
-o $OUTPATH \
-t 10 \
-i $INPATH/assemblies \
--mcr_path {path_to_previous_BactPrep_folder}/matlab/v901 \
--fastgear_exe {path_to_previous_BactPrep_folder}/fastGEARpackageLinux64bit 
```

**7) IF you would like to inform wgsRecomb (gubbins) about already known phage region**
phage region can be provided to snippy before running gubbins. If you would like to provide known phage region while running gubbins, use ```-v``` or ```--phage``` to provide phage region in a BED file.

EXAMPLE:

```
python start_analysis.py wgsRecomb \
-p PMEN1.dated \
-o $WORKPATH -i $WORKPATH/assemblies \
-r $WORKPATH/GCF_000026665.1_ASM2666v1_genomic.fna \
-v $WORKPATH/phage_region.bed

```

**8) IF additional arguments need to be specificed for Roary and Gubbins when using "coreGEN", "wgsRecomb", or "ALL" module**
additional Roary and Gubbins arguments that is not specificed by BactPrep can be added by using the ```-R``` of ```-G``` flags, respectively. Dependencies used for these additional arguments need to be install by user.

**SPACE is necessary at the beginning of the string**

Example:
```
python start_analysis.py ALL \
-p PMEN1.dated \
-o $WORKPATH -g $WORKPATH/gff \
-r $WORKPATH/GCF_000026665.1_ASM2666v1_genomic.fna \
-R " -r -y -iv 1.5"

```
