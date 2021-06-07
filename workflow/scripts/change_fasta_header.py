from Bio import SeqIO
import pandas as pd
import sys

print(sys.argv[1])
metadata = pd.read_csv(str(sys.argv[1]))

print(int(sys.argv[4]))
dict_key= metadata.columns[int(sys.argv[4])]

dict_list=[]
items_list=str(sys.argv[3]).split(",")
for item in items_list:
    temp_dict=metadata.loc[:,[dict_key,item]].set_index(dict_key).to_dict()[item]
    dict_list.append(temp_dict)

with open(str(sys.argv[5]), "w+") as w:
    print("start a new file, overwrriten the old file if exist")

with open(str(sys.argv[2])) as f, open(str(sys.argv[5]), "w") as w:
    records = SeqIO.parse(f,"fasta")
    for r in records:
        old_id = r.id
        new_name=r.id
        for d in range(len(dict_list)):
            current_dict= dict_list[d]
            current_meta=current_dict[old_id]
            new_name = new_name + "|" + str(current_meta) 
        r.id = new_name
        r.description = ""
        print(r.id)
        SeqIO.write(r, w, "fasta")

    



