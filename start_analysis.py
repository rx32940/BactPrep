import os
import argparse
import yaml
from pathlib import Path

modules = ['wgsRecomb', 'coreGen', 'coreRecomb','panRecomb', 'geneRecomb','ALL']

CWD = os.getcwd()

# declare an argparse variable 
general_parser = argparse.ArgumentParser(prog="start_analysis", usage='%(prog)s MODULE [options]', 
description='Please always specify the program to use in the first argument, \
\nor the whole pipeline will attemp to run',allow_abbrev=False,
epilog='Enjoy the program! :)')


# postional arguments required to provide by user
general_parser.add_argument('MODULE', action='store', type=str, 
help='Specify the module you would like to run', 
choices=['ALL','wgsRecomb', 'coreGen', 'coreRecomb','panRecomb', 'geneRecomb'])
curren_wd = Path(os.getcwd())



# General Arguments
general_arguments = general_parser.add_argument_group("general arguments")
general_arguments.add_argument("-i","--input", type=str, help= "path to input dir with assemlies", metavar='', default=os.path.join(curren_wd, "../assemblies"))
general_arguments.add_argument("-p", "--name", type=str, help= "provide name prefix for the output files",metavar='', default= "bacterial_analysis" )
general_arguments.add_argument('-t','--thread', type=int,  help = "num of threads", metavar='', default = 1)
general_arguments.add_argument('-o','--output', type=str,  help = "path to the output directory", metavar='', default=os.path.join(CWD, 'results'))

# output annotation arguments
annotate_arguments = general_parser.add_argument_group("arguments for if you would like to add metadata to output")
annotate_arguments.add_argument("-M", "--addMetadata", default=False, action='store_true', help= "must have the flag specify if want to allow annotation")
annotate_arguments.add_argument("-a", "--annotate", type=str, help= "path to a csv file containing sample metadata", metavar='', default="")
annotate_arguments.add_argument("-s", "--sample", type=int, help= "integer indicates which column the sample name is in the metadata csv file", metavar='',default=1)
annotate_arguments.add_argument("-m", "--metadata", type=str,help= "metadata chosen to annotate ML tree/alignment after the sample name",metavar='',default="")


# wgsRecomb arguments
Recomb_arguments = general_parser.add_argument_group("arguments for wgsRecomb module")
Recomb_arguments.add_argument('-r','--ref', type=str, help = "reference (required for wgsRecomb module)", metavar='',default="")
Recomb_arguments.add_argument("-v", "--phage", type=str,help= "phage region identified for masking (bed file)",metavar='',default="")
Recomb_arguments.add_argument("-G", "--gubbins", type=str,help= "any additional Gubbins arguments (PLEASE start the string with space) ",metavar='',default="")

# coreGen arguments
roary_arguments = general_parser.add_argument_group("arguments for coreGen module")
roary_arguments.add_argument("-g", "--gff", type=str, help= "path to input dir with gff (this can replace input assemblies dir in coreGen module Must be gff3 files)", metavar='')
roary_arguments.add_argument("-c", "--core", type=int,help= "define core gene definition by percentage for coreGen module (default=99)",metavar='',default=99)
roary_arguments.add_argument("-k", "--kingdom", type=str,help= "specify the kingom of input assemlies for genome annotation (default=Bacteria)",metavar='',default="Bacteria")
roary_arguments.add_argument("-R", "--roary", type=str,help= "any additional roary arguments (PLEASE start the string with space)",metavar='',default="")

# fastGear modules alignments (for all three fastGear moudles)
fastgear_arguments = general_parser.add_argument_group("arguments for all three fastGear modules (coreRecomb, panRecomb, geneRecomb)")
fastgear_arguments.add_argument("--mcr_path", type=str, help="path to mcr runtime (need to install before use any of the fastGear module", metavar='', default=os.path.join(CWD, 'resources/mcr/v901/'))
fastgear_arguments.add_argument("--fastgear_exe", type=str, help="path to the excutable of fastGear", metavar='', default=str(os.path.join(CWD,'resources/fastGEARpackageLinux64bit/')))
fastgear_arguments.add_argument("--fg","--fastgear_param", type=str, help="path to fastGear params", metavar='', default="")

# geneRecomb arguments
fastgear_gene_arguments = general_parser.add_argument_group("arguments for geneRecomb module")
fastgear_gene_arguments.add_argument("-n", "--alignment", type=str,help= "input alignment\n(either -n/-fl is required for geneRecomb module) ",metavar='',default="")
fastgear_gene_arguments.add_argument("-fl", "--alnlist", type=str,help= "input alignment list with path to gene alignments\n(either -n/-fl is required for geneRecomb module) ",metavar='',default="")


args = general_parser.parse_args()

if args.MODULE == "geneRecomb" and (args.alignment == "" and args.alnlist == ""):
    general_parser.error("gene alignment/alignment list (-n/-fl) must provided for geneRecomb module") 
if args.alignment != "" and args.alnlist != "":
    general_parser.error("please do not specify gene file and list and same time")

MODULE = args.MODULE
NAME=args.name
INPUT=args.input
REF=args.ref
OUT=args.output
THREAD=args.thread
GFF=args.gff
ROARY=args.roary
ADDANOT=args.addMetadata
ANOT=args.annotate
SAMPLE=args.sample
META=args.metadata
PHAGE=args.phage
GUBBINS=args.gubbins
CORE=args.core
KINGDOM=args.kingdom
ALN=args.alignment
FL=args.alnlist
SINGLE= True if ALN != "" else False
FASTGEAR_EXE=args.fastgear_exe
FASTGEAR_PARAM=args.fg if args.fg != "" else str(os.path.join(FASTGEAR_EXE,'fG_input_specs.txt'))
MCR_PATH=args.mcr_path
LD_LIB_PATH=str(os.path.join(MCR_PATH,'runtime/glnxa64:')) + \
    str(os.path.join(MCR_PATH , 'bin/glnxa64:')) +\
        str(os.path.join(MCR_PATH , 'sys/os/glnxa64'))


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
    if module == "coreGen":
        anot_files=str(os.path.join(OUT,"roary","roary_iqtree",str(NAME + "_meta.coreConcate.newick"))) + "," + str(os.path.join(OUT,"roary",str(NAME + "_coreConcate_meta.fasta")))
    if module == "wgsRecomb":
        anot_files=str(os.path.join(OUT, "gubbins",str(NAME + "_meta.recombFreeSnpsAtcg.fasta"))) +"," +str(os.path.join(OUT ,'gubbins','iqtree' ,str(NAME + "_meta.GubbinsSNPs.newick")))
    if module == "coreRecomb":
        anot_files=str(os.path.join(OUT,"fastgear_core" , "fastgear_iqtree" ,  str(NAME + "_meta.coreSNPs.newick"))) + "," + str(os.path.join(OUT , "fastgear_core",str(NAME + "_core_mask_snp_meta.fasta")))
    if module == "ALL":
        anot_files=str(os.path.join(OUT,"roary","roary_iqtree",str(NAME + "_meta.coreConcate.newick"))) +","\
            + str(os.path.join(OUT,"roary",str(NAME + "_coreConcate_meta.fasta"))) + ","\
                + str(os.path.join(OUT, "gubbins",str(NAME + "_meta.recombFreeSnpsAtcg.fasta"))) +","\
                    +str(os.path.join(OUT ,'gubbins','iqtree' ,str(NAME + "_meta.GubbinsSNPs.newick"))) +","\
                        + str(os.path.join(OUT,"fastgear_core" , "fastgear_iqtree" ,  str(NAME + "_meta.coreSNPs.newick"))) + ","\
                            + str(os.path.join(OUT , "fastgear_core",str(NAME + "_core_mask_snp_meta.fasta")))
    if module == "geneRecomb":
        ADDANOT=Fasle
        anot_files=""
    if module == "panRecomb":
        ADDANOT=Fasle
        anot_files=""
    return anot_files

def get_output(module):

    if module == "coreGen":
        output_files=str(os.path.join(OUT,'roary',"roary_iqtree" , str(NAME +".treefile")))
        if ADDANOT:
            output_files=str(output_files) + ","+get_annotated(module)
    elif module == "wgsRecomb":
        output_files=str(os.path.join(OUT , "gubbins", "iqtree" , str(NAME + ".recombFreeSnpsAtcg.treefile")))
        if ADDANOT:
            output_files=str(output_files) + ","+get_annotated(module)        
    elif module == "panRecomb":
        output_files=str(os.path.join(OUT , "fastgear" , "plot_pangenome/pan_fastgear_plot_recombination_count.pdf")) 
        if ADDANOT:
            output_files=str(output_files) + ","+get_annotated(module) 
    elif module == "coreRecomb":
        output_files=str(os.path.join(OUT , "fastgear_core" , "plot_coregenome/core_fastgear_plot_recombination_count.pdf")) + ","\
            + str(os.path.join(OUT, "fastgear_core" , "fastgear_iqtree" , str(NAME + "_core_mask_snp.treefile")))  
    elif module == "geneRecomb":
        if SINGLE:
            output_files=",".join([str(os.path.join(OUT,"fastgear_gene" ,GENE_NAME[0],str(GENE_NAME[0] + ".mat")))])
        else:
            output_files=",".join([str(os.path.join(OUT,"fastgear_gene" ,file,str(file + ".mat"))) for file in GENE_NAME])
    elif module == "ALL":
        output_files=str(os.path.join(OUT,"roary","roary_iqtree" , str(NAME +".treefile")))+ ","\
            + str(os.path.join(OUT , "gubbins", "iqtree" , str(NAME + ".recombFreeSnpsAtcg.treefile")))+ ","\
                + str(os.path.join(OUT , "fastgear_core" , "plot_coregenome/core_fastgear_plot_recombination_count.pdf")) + ","\
                    + str(os.path.join(OUT, "fastgear_core" , "fastgear_iqtree" , str(NAME + "_core_mask_snp.treefile"))) + ","\
                        + str(os.path.join(OUT , "fastgear_core" , "plot_coregenome/core_fastgear_plot_recombination_count.pdf"))
        if ADDANOT:
            output_files=str(output_files) + "," + get_annotated(module)

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
'define_core': CORE,
'phage_region':PHAGE,
'LD_LIBRARY_PATH': LD_LIB_PATH,
'fastGear_exe_path': FASTGEAR_EXE,
'mcr_path': MCR_PATH,
'fastGear_params': FASTGEAR_PARAM,
'fastgear_gene_file_list': FILELIST,
'roary':ROARY,
'gubbins':GUBBINS}


config.update(get_output(MODULE))     
with open("config/config.yaml", "w") as configfile:
    yaml.dump(config,configfile)



if __name__ == "__main__":
    os.system ("snakemake --cores %d --use-conda"%THREAD)



