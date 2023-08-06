import pandas as pd
from Bio import SeqIO
import sys
import os

input_dir=sys.argv[1]
output_dir=sys.argv[2]

# read in roary's output file
roary_output=pd.read_csv(sys.argv[3])

# remove the first 14 columns
biosample_roaryID=roary_output.iloc[1,14:].to_frame()

# set the name for the columb with Roary ID 
biosample_roaryID.columns = ["RoaryID"]

# remove the cluster number after roary ID 
biosample_roaryID= biosample_roaryID["RoaryID"].str.split("_").str[0].to_frame()

# make rownames (BioSample Accesion) a column (instead of rownames)
biosample_roaryID.index.name = 'BioSampleID'

# reset the rwonames of the dataframe
biosample_roaryID.reset_index(inplace=True)

# turn dataframe into a dictionary between Roary ID and BioSample Accession 
biosample_roaryID_dict=biosample_roaryID.set_index('RoaryID').to_dict()['BioSampleID']

# declare a function to change the header from Roary output to Biosample ID
def change_header(gene_fasta, changed_fasta):
    with open(gene_fasta) as f, open(changed_fasta, "w") as w:
            records = SeqIO.parse(f,"fasta")
            for r in records:
                currentID=r.id
                roaryID = currentID.split("_")[0]
                BioSampleAcc=biosample_roaryID_dict[roaryID]
                r.id = BioSampleAcc
                r.description = ""
                SeqIO.write(r, w, "fasta")


# get directory
directory = os.fsencode(os.path.join(input_dir))

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    old_filename=os.path.join(input_dir,filename)
    new_filename=os.path.join(output_dir,filename)
    change_header(old_filename, new_filename)

