import os
import argparse
import yaml
from pathlib import Path

modules = ['gubbins', 'roary', 'fastGear', 'fastGear_gene','ALL']


# declare an argparse variable 
general_parser = argparse.ArgumentParser(prog="start_analysis", usage='%(prog)s MODULE [options]', 
description='Please always specify the program to use in the first argument, \
\nor the whole pipeline will attemp to run',allow_abbrev=False,
epilog='Enjoy the program! :)')


# postional arguments required to provide by user
general_parser.add_argument('MODULE', action='store', type=str, 
help='Specify the module you would like to run', 
choices=['gubbins', 'roary', 'fastGear', 'fastGear_gene','ALL'])
curren_wd = Path(os.getcwd())

# General Arguments
general_arguments = general_parser.add_argument_group("general arguments")
general_arguments.add_argument("-i","--input", type=str, help= "path to input dir with assemlies", metavar='', default=os.path.join(curren_wd, "../assemblies"))
general_arguments.add_argument("-p", "--name", type=str, help= "provide name prefix for the output files",metavar='', default= "bacterial_analysis" )
general_arguments.add_argument('-t','--thread', type=int,  help = "num of threads", metavar='', default = 1)
general_arguments.add_argument('-o','--output', type=str,  help = "path to the output directory", metavar='', default=str(curren_wd.parent.absolute()))
general_arguments.add_argument("-a", "--annotate", type=str, help= "path to a csv file containing sample metadata", metavar='', default="")
general_arguments.add_argument("-s", "--sample", type=int, help= "integer indicates which column the sample name is in the metadata csv file", metavar='',default=1)
general_arguments.add_argument("-m", "--metadata", help= "metadata chosen to annotate ML tree/alignment after the sample name", nargs='*',metavar='',default="")



# gubbins arguments
gubbins_arguments = general_parser.add_argument_group("arguments for gubbins module")
gubbins_arguments.add_argument('-r','--ref', type=str, help = "reference (required for gubbins module)", metavar='',default="")
gubbins_arguments.add_argument("-v", "--phage", type=str,help= "phage region identified for masking (bed file)",metavar='')

# roary arguments
roary_arguments = general_parser.add_argument_group("arguments for roary module")
roary_arguments.add_argument("-g", "--gff", type=str, help= "path to input dir with gff (this can replace input assemblies dir in roary module Must be gff3 files)", metavar='')
roary_arguments.add_argument("-c", "--core", type=int,help= "define core gene definition by percentage for roary module (default=99)",metavar='',default=99)
roary_arguments.add_argument("-k", "--kingdom", type=str,help= "specify the kingom of input assemlies for genome annotation (default=Bacteria)",metavar='',default="Bacteria")
roary_arguments.add_argument("-z", "--genus", type=str,help= "specify the genus of input assemlies for genome annotation (default=Leptospira)",metavar='',default="Leptospira")

# fastgear_gene arguments
fastgear_gene_arguments = general_parser.add_argument_group("arguments for fastGear_gene module")
fastgear_gene_arguments.add_argument("-n", "--alignment", type=str,help= "input alignment\n(either -n/-fl is required for fastGear_gene module) ",metavar='',default="")
fastgear_gene_arguments.add_argument("-fl", "--alnlist", type=str,help= "input alignment list with path to gene alignments\n(either -n/-fl is required for fastGear_gene module) ",metavar='',default="")


args = general_parser.parse_args()

if args.MODULE == "fastGear_gene" and (args.alignment == "" and args.alnlist == ""):
    general_parser.error("gene alignment/alignment list (-n/-fl) must provided for fastGear_gene module") 
if args.alignment != "" and args.alnlist != "":
    general_parser.error("please do not specify gene file and list and same time")

MODULE = args.MODULE
NAME=args.name
INPUT=args.input
REF=args.ref
OUT=args.output
THREAD=args.thread
GFF=args.gff
ANOT=args.annotate
SAMPLE=args.sample
META=args.metadata
PHAGE=args.phage
CORE=args.core
KINGDOM=args.kingdom
GENUS=args.genus
ALN=args.alignment
FL=args.alnlist
SINGLE= True if ALN != "" else False



# open 

if FL != "":
    with open(FL) as f:
        FILELIST = f.readlines()
                
    FILELIST = [x.strip() for x in FILELIST] 
elif ALN !="":
    FILELIST = [ALN]
else:
    FILELIST = ""




def get_geneNames():
    gene_name=""
    if ALN != "":
        gene_name=os.path.basename(ALN).split(".")[0]
    elif FL != "": # list input
        with open(FL) as f:
            for line in f:
                gene_name=gene_name  + os.path.basename(line).split(".")[0]+ ","
        gene_name = gene_name[:-1]
        # you may also want to remove whitespace characters like `\n` at the end of each line
        
    else:
        gene_name=""
    
    return gene_name.split(",")

GENE_NAME=get_geneNames()


def get_annotated(module):
    anot_files=None
    if module == "roary":
        anot_files=str(os.path.join(OUT,"iqtree",str(NAME + "_meta.coreConcate.newick")))
    if module == "gubbins":
        anot_files=str(os.path.join(OUT, "gubbins",str(NAME + "_meta.recombFreeSnpsAtcg.fasta"))) +"," +str(os.path.join(OUT ,'gubbins' ,str(NAME + "_meta.GubbinsSNPs.newick")))
    if module == "fastgear":
        anot_files=str(os.path.join(OUT,"fastgear" , "fastgear_iqtree" ,  str(NAME + "_meta.coreSNPs.newick")))
    return anot_files

def get_output(module):
    output_files=""
    if module == "roary":
        output_files=str(os.path.join(OUT,"iqtree" , str(NAME +".treefile")))
        if ANOT:
            output_files=str(output_files) + ",",get_annotated(module)
    elif module == "gubbins":
        output_files=str(os.path.join(OUT , "gubbins", "iqtree" , str(NAME + ".recombFreeSnpsAtcg.treefile")))
        if ANOT:
            output_files=str(output_files) + ",",get_annotated(module)        
    elif module == "fastGear":
        output_files=str(os.path.join(OUT , "fastgear" , "fastgear_iqtree" + str(NAME + "_core_mask_snp.treefile")))
        if ANOT:
            output_files=str(output_files) + ",",get_annotated(module)  
    elif module == "fastGear_gene":
        if SINGLE:
            output_files=",".join([str(os.path.join(OUT,"fastgear_gene" ,GENE_NAME[0],str(GENE_NAME[0] + ".mat")))])
        else:
            output_files=",".join([str(os.path.join(OUT,"fastgear_gene" ,file,str(file + ".mat"))) for file in GENE_NAME])
    elif module == "ALL":
        output_files=str(os.path.join(OUT,"iqtree" , str(NAME +".treefile")))+ "," + str(os.path.join(OUT , "gubbins", "iqtree" , str(NAME + ".recombFreeSnpsAtcg.treefile")))+ "," + str(os.path.join(OUT , "fastgear" , "fastgear_iqtree" , str(NAME + "_core_mask_snp.treefile")))
        if ANOT:
            output_files=str(output_files) + ",",get_annotated(module)
    
    output={'output': output_files}      
    return output





# construct the config file 
config = {'project_name': NAME,
'asm_dir': INPUT,
'output_dir': OUT,
'reference': REF,
'threads_num': THREAD,
'sample_metadata': ANOT,
'metadata_include': META,
'biosample_column': SAMPLE,
'gff_dir': GFF,
'kingdom': KINGDOM,
'genus': GENUS,
'define_core': CORE,
'phage_region':PHAGE,
'LD_LIBRARY_PATH': '/home/rx32940/matlab/v901/runtime/glnxa64:/home/rx32940/matlab/v901/bin/glnxa64:/home/rx32940/matlab/v901/sys/os/glnxa64',
'fastGear_exe_path': '/home/rx32940/SOFTWARE/fastGEARpackageLinux64bit/',
'mcr_path': '/home/rx32940/matlab/v901/',
'fastGear_gene_alignment': ALN,
'fastGear_gene_alignment': FL,
'fastGear_batch':FL,
'fastgear_gene_file_list': FILELIST}


if MODULE == "roary":
    config.update(get_output(MODULE))
    with open("config/test.yaml", "w") as configfile:
        yaml.dump(config,configfile)
elif MODULE == "gubbins":
    config.update(get_output(MODULE))
    with open("config/test.yaml", "w") as configfile:
        yaml.dump(config,configfile)
elif MODULE == "fastgear":
    config.update(get_output(MODULE))   
    with open("config/test.yaml", "w") as configfile:
        yaml.dump(config,configfile)
else:
    config.update(get_output(MODULE))     
    with open("config/test.yaml", "w") as configfile:
        yaml.dump(config,configfile)


if __name__ == "__main__":
    os.system ("snakemake --cores 10 --use-conda")



