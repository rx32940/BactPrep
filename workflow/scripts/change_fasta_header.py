from Bio import SeqIO
import pandas as pd


metadata = pd.read_csv(snakemake.input[0])

dict_list=[]
for item in metadata_include:
    temp_dict=metadata[["BioSample.Accession",item]].set_index("BioSample.Accession").to_dict()[item]
    dict_list.append(temp_dict)

with open(output_file, "w") as w:
    print("start a new file, overwrriten the old file if exist")

with open(snakemake.input[1]) as f, open(snakemake.output[0], "w") as w:
    records = SeqIO.parse(f,"fasta")
    for r in records:
        old_id = r.id
        new_name=""
        for d in range(len(dict_list)):
            current_dict= dict_list[d]
            currnt_meta=current_dict[old_id]
            new_name = new_name + "|" + str(current_meta) 
        r.id = new_name
        r.description = ""
        SeqIO.write(r, w, "fasta")

