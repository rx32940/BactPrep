import os
import pandas as pd


# read recent recombinaitons
fastGear_recent = pd.read_csv(snakemake.input[0],sep="\s+",header=1)
fastGear_recent_out = fastGear_recent[["StrainName", "Start", "End"]] # format in bed file 


# every ancestral recombination involves multiple strains from two lineages (donor and recipient)
# we will list all the strains in these two lineages with the ancestral recombination region detected
# this is a rowwise function
def mask_strains_in_lineage(row, lineage):
    lineage1 = row["Lineage1"] # lineage 1
    lineage2 = row["Lineage2"] # lineage 2
    start = row["Start"] # start position
    end = row["End"] # end position
    lineages = lineage[lineage["Lineage"].isin([lineage1, lineage2])]["Name"] # strains in both lineage 1 and 2
    current_ances_row_bed=pd.DataFrame({"StrainName" : lineages.to_list(), "Start" : [int(start)] * len(lineages), "End" : [int(end)] * len(lineages)}) # format involved strains with their recombined regions
    return current_ances_row_bed

# read in ancstral recombinaition files
fastGear_ances = pd.read_csv(snakemake.input[1],sep="\s+",header=1)
lineage_file=pd.read_csv(snakemake.input[2],sep="\s+",header=0) # lineage files (include information about every strain)

fastGear_ances_out_df=pd.DataFrame(columns={"StrainName","Start","End"}) # declare an empty dataframe for ancestral recombinations
# iterate through each row to form the concatenation of dataframe in the BED format for all strains involved in ancestral recombination
fastGear_ances_out = fastGear_ances.apply(lambda row: mask_strains_in_lineage(row, lineage_file), axis=1)

if not fastGear_ances_out.empty:
    for item in fastGear_ances_out.iteritems():
        fastGear_ances_out_df = fastGear_ances_out_df.append(item[1])

fastGear_out=fastGear_ances_out_df.append(fastGear_recent_out,ignore_index=True,sort=False)[["StrainName", "Start", "End"]] # reorder columns

fastGear_out.to_csv(snakemake.output[0], header=False, index=False, sep="\t")






    