import pandas as pd
import os
from Bio import Phylo
import sys

allDated_meta = pd.read_csv(str(sys.argv[1]))
dict_key= allDated_meta.columns[int(sys.argv[4])]

dict_list=[]
items_list=str(sys.argv[3]).split(",")
for item in items_list:
    temp_dict=allDated_meta.loc[:,[dict_key,item]].set_index(dict_key).to_dict()[item]
    dict_list.append(temp_dict)


tree = Phylo.read(str(sys.argv[2]),"newick")

for taxa in tree.get_terminals():
    newname=str(taxa)
    for d in range(len(dict_list)):
        newname=newname + "|"+ str(dict_list[d][str(taxa)]) 
    tree.find_any(taxa).name = newname

Phylo.write(tree,str(sys.argv[5]),"newick")




